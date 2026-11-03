use std::collections::HashMap;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::Interval;
use smallvec::SmallVec;

const BIN_SIZE: u64 = 16 * 1024;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) enum Strand {
    Forward,
    Reverse,
    Unknown,
}

impl Strand {
    pub(crate) fn parse(value: &str) -> Result<Self> {
        match value {
            "+" => Ok(Self::Forward),
            "-" => Ok(Self::Reverse),
            "." => Ok(Self::Unknown),
            _ => Err(RsomicsError::InvalidInput(format!(
                "invalid annotation strand {value:?}"
            ))),
        }
    }

    pub(crate) const fn as_str(self) -> &'static str {
        match self {
            Self::Forward => "+",
            Self::Reverse => "-",
            Self::Unknown => ".",
        }
    }
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum StrandMode {
    Unstranded,
    Sense,
    Antisense,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum CountLevel {
    MetaFeature,
    Feature,
}

#[derive(Clone, Debug)]
pub(crate) struct Block {
    pub(crate) chrom: Option<u32>,
    pub(crate) start: u64,
    pub(crate) end: u64,
}

pub(crate) type Blocks = SmallVec<[Block; 2]>;
pub(crate) type Hits = SmallVec<[Hit; 4]>;

impl Block {
    pub(crate) fn new(chrom: Option<u32>, start: u64, end: u64) -> Result<Self> {
        let interval = Interval::new((), start, end)
            .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
        if interval.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "alignment block has zero length".to_owned(),
            ));
        }
        Ok(Self {
            chrom,
            start: interval.start(),
            end: interval.end(),
        })
    }
}

#[derive(Clone, Debug)]
pub(crate) struct PartInput {
    pub(crate) id: String,
    pub(crate) chrom: String,
    pub(crate) start: u64,
    pub(crate) end: u64,
    pub(crate) strand: Strand,
}

#[derive(Clone, Debug)]
struct Part {
    interval: Interval<String>,
    chrom_id: u32,
    strand: Strand,
    row: usize,
}

#[derive(Clone, Debug)]
pub(crate) struct FeatureRow {
    id: String,
    parts: Vec<usize>,
    length: u64,
}

impl FeatureRow {
    #[must_use]
    pub(crate) fn id(&self) -> &str {
        &self.id
    }

    #[must_use]
    pub(crate) const fn length(&self) -> u64 {
        self.length
    }
}

#[derive(Clone, Debug)]
pub(crate) struct FeatureSet {
    rows: Vec<FeatureRow>,
    parts: Vec<Part>,
    chrom_ids: HashMap<String, u32>,
    index: HashMap<u32, HashMap<u64, Vec<usize>>>,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub(crate) struct Hit {
    pub(crate) row: usize,
    pub(crate) overlap_bases: u64,
    pub(crate) feature_overlap_bases: u64,
    pub(crate) feature_length: u64,
}

impl FeatureSet {
    pub(crate) fn build(inputs: Vec<PartInput>, level: CountLevel) -> Result<Self> {
        if inputs.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "annotation contains no selected features".to_owned(),
            ));
        }

        let mut rows = Vec::<FeatureRow>::new();
        let mut parts = Vec::<Part>::with_capacity(inputs.len());
        let mut groups = HashMap::<String, usize>::new();
        let mut chrom_ids = HashMap::<String, u32>::new();

        for input in inputs {
            let next_chrom_id = u32::try_from(chrom_ids.len()).map_err(|_| {
                RsomicsError::InvalidInput("too many annotation chromosomes".to_owned())
            })?;
            let chrom_id = *chrom_ids
                .entry(input.chrom.clone())
                .or_insert(next_chrom_id);
            let interval = Interval::new(input.chrom, input.start, input.end)
                .map_err(|error| RsomicsError::InvalidInput(error.to_string()))?;
            if interval.is_empty() {
                return Err(RsomicsError::InvalidInput(format!(
                    "feature {} has zero length",
                    input.id
                )));
            }
            let row = match level {
                CountLevel::Feature => {
                    let row = rows.len();
                    rows.push(FeatureRow {
                        id: input.id,
                        parts: Vec::new(),
                        length: 0,
                    });
                    row
                }
                CountLevel::MetaFeature => match groups.get(&input.id).copied() {
                    Some(row) => row,
                    None => {
                        let row = rows.len();
                        groups.insert(input.id.clone(), row);
                        rows.push(FeatureRow {
                            id: input.id,
                            parts: Vec::new(),
                            length: 0,
                        });
                        row
                    }
                },
            };
            let part = parts.len();
            parts.push(Part {
                interval,
                chrom_id,
                strand: input.strand,
                row,
            });
            rows[row].parts.push(part);
        }

        for row in &mut rows {
            row.length = union_length(row.parts.iter().map(|&part| &parts[part]));
        }

        let mut index: HashMap<u32, HashMap<u64, Vec<usize>>> = HashMap::new();
        for (part_index, part) in parts.iter().enumerate() {
            let first = part.interval.start() / BIN_SIZE;
            let last = (part.interval.end() - 1) / BIN_SIZE;
            let bins = index.entry(part.chrom_id).or_default();
            for bin in first..=last {
                bins.entry(bin).or_default().push(part_index);
            }
        }

        Ok(Self {
            rows,
            parts,
            chrom_ids,
            index,
        })
    }

    #[must_use]
    pub(crate) fn rows(&self) -> &[FeatureRow] {
        &self.rows
    }

    pub(crate) fn chrom_id(&self, name: &str) -> Option<u32> {
        self.chrom_ids.get(name).copied()
    }

    pub(crate) fn hits(
        &self,
        blocks: &[Block],
        strand_mode: StrandMode,
        reverse: bool,
        track_feature_fraction: bool,
    ) -> Hits {
        self.hits_grouped(&[(blocks, reverse, 0)], strand_mode, track_feature_fraction)
    }

    pub(crate) fn hits_grouped(
        &self,
        groups: &[(&[Block], bool, u32)],
        strand_mode: StrandMode,
        track_feature_fraction: bool,
    ) -> Hits {
        let mut candidates = SmallVec::<[usize; 16]>::new();
        for block in groups.iter().flat_map(|(blocks, _, _)| *blocks) {
            let Some(chrom) = block.chrom else {
                continue;
            };
            let Some(bins) = self.index.get(&chrom) else {
                continue;
            };
            let first = block.start / BIN_SIZE;
            let last = (block.end - 1) / BIN_SIZE;
            for bin in first..=last {
                if let Some(parts) = bins.get(&bin) {
                    candidates.extend_from_slice(parts);
                }
            }
        }
        candidates.sort_unstable();
        candidates.dedup();

        let mut intersections = SmallVec::<[(usize, u32, u64, u64); 16]>::new();
        let mut feature_overlaps = SmallVec::<[(usize, u64, u64); 16]>::new();
        for part_index in candidates {
            let part = &self.parts[part_index];
            let mut part_intersections = SmallVec::<[(u64, u64); 4]>::new();
            for &(blocks, reverse, partition) in groups {
                if !strand_matches(part.strand, strand_mode, reverse) {
                    continue;
                }
                for block in blocks
                    .iter()
                    .filter(|block| block.chrom == Some(part.chrom_id))
                {
                    let start = block.start.max(part.interval.start());
                    let end = block.end.min(part.interval.end());
                    if start < end {
                        intersections.push((part.row, partition, start, end));
                        if track_feature_fraction {
                            part_intersections.push((start, end));
                        }
                    }
                }
            }
            if track_feature_fraction {
                let overlap = union_interval_length(&mut part_intersections);
                if overlap > 0 {
                    feature_overlaps.push((part.row, overlap, part.interval.len()));
                }
            }
        }
        union_hits(intersections, feature_overlaps)
    }

    pub(crate) fn fields(&self, row: usize) -> (String, String, String, String) {
        let row = &self.rows[row];
        let mut chrom = Vec::with_capacity(row.parts.len());
        let mut start = Vec::with_capacity(row.parts.len());
        let mut end = Vec::with_capacity(row.parts.len());
        let mut strand = Vec::with_capacity(row.parts.len());
        for &part in &row.parts {
            let part = &self.parts[part];
            chrom.push(part.interval.chrom());
            start.push((part.interval.start() + 1).to_string());
            end.push(part.interval.end().to_string());
            strand.push(part.strand.as_str());
        }
        (
            chrom.join(";"),
            start.join(";"),
            end.join(";"),
            strand.join(";"),
        )
    }
}

fn strand_matches(feature: Strand, mode: StrandMode, reverse: bool) -> bool {
    match mode {
        StrandMode::Unstranded => true,
        StrandMode::Sense => matches!(
            (feature, reverse),
            (Strand::Forward, false) | (Strand::Reverse, true)
        ),
        StrandMode::Antisense => matches!(
            (feature, reverse),
            (Strand::Forward, true) | (Strand::Reverse, false)
        ),
    }
}

fn union_length<'a>(parts: impl Iterator<Item = &'a Part>) -> u64 {
    union_tuples(
        parts
            .map(|part| {
                (
                    part.interval.chrom().to_owned(),
                    part.interval.start(),
                    part.interval.end(),
                )
            })
            .collect(),
    )
}

fn union_tuples(mut intervals: Vec<(String, u64, u64)>) -> u64 {
    intervals.sort_unstable_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then(left.1.cmp(&right.1))
            .then(left.2.cmp(&right.2))
    });
    let mut total = 0;
    let mut current: Option<(String, u64, u64)> = None;
    for (chrom, start, end) in intervals {
        match current.as_mut() {
            Some((current_chrom, _, current_end))
                if *current_chrom == chrom && start <= *current_end =>
            {
                *current_end = (*current_end).max(end);
            }
            Some((_, current_start, current_end)) => {
                total += *current_end - *current_start;
                current = Some((chrom, start, end));
            }
            None => current = Some((chrom, start, end)),
        }
    }
    if let Some((_, start, end)) = current {
        total += end - start;
    }
    total
}

fn union_interval_length(intervals: &mut SmallVec<[(u64, u64); 4]>) -> u64 {
    intervals.sort_unstable();
    let mut total = 0;
    let mut current = None;
    for (start, end) in intervals.drain(..) {
        match current.as_mut() {
            Some((_, current_end)) if start <= *current_end => {
                *current_end = (*current_end).max(end);
            }
            Some((current_start, current_end)) => {
                total += *current_end - *current_start;
                current = Some((start, end));
            }
            None => current = Some((start, end)),
        }
    }
    if let Some((start, end)) = current {
        total += end - start;
    }
    total
}

fn union_hits(
    mut intervals: SmallVec<[(usize, u32, u64, u64); 16]>,
    feature_overlaps: SmallVec<[(usize, u64, u64); 16]>,
) -> Hits {
    let mut best = HashMap::<usize, (u64, u64)>::new();
    for (row, overlap, length) in feature_overlaps {
        let candidate = (overlap, length);
        let current = best.entry(row).or_insert(candidate);
        if u128::from(overlap) * u128::from(current.1) > u128::from(current.0) * u128::from(length)
        {
            *current = candidate;
        }
    }
    intervals.sort_unstable();
    let mut hits = Hits::new();
    let mut current: Option<(usize, u32, u64, u64)> = None;
    let mut row_total = 0;
    for (row, chrom, start, end) in intervals {
        match current.as_mut() {
            Some((current_row, current_chrom, _, current_end))
                if *current_row == row && *current_chrom == chrom && start <= *current_end =>
            {
                *current_end = (*current_end).max(end);
            }
            Some((current_row, _, current_start, current_end)) => {
                row_total += *current_end - *current_start;
                if *current_row != row {
                    let (feature_overlap_bases, feature_length) =
                        best.get(current_row).copied().unwrap_or((1, 1));
                    hits.push(Hit {
                        row: *current_row,
                        overlap_bases: row_total,
                        feature_overlap_bases,
                        feature_length,
                    });
                    row_total = 0;
                }
                current = Some((row, chrom, start, end));
            }
            None => current = Some((row, chrom, start, end)),
        }
    }
    if let Some((row, _, start, end)) = current {
        row_total += end - start;
        let (feature_overlap_bases, feature_length) = best.get(&row).copied().unwrap_or((1, 1));
        hits.push(Hit {
            row,
            overlap_bases: row_total,
            feature_overlap_bases,
            feature_length,
        });
    }
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    fn part(id: &str, start: u64, end: u64, strand: Strand) -> PartInput {
        PartInput {
            id: id.to_owned(),
            chrom: "chr1".to_owned(),
            start,
            end,
            strand,
        }
    }

    #[test]
    fn meta_feature_length_is_a_union() {
        let features = FeatureSet::build(
            vec![
                part("g", 300, 350, Strand::Forward),
                part("g", 340, 400, Strand::Forward),
            ],
            CountLevel::MetaFeature,
        )
        .unwrap();
        assert_eq!(features.rows()[0].length(), 100);
    }

    #[test]
    fn query_deduplicates_overlapping_parts_and_blocks() {
        let features = FeatureSet::build(
            vec![
                part("g", 100, 150, Strand::Forward),
                part("g", 140, 200, Strand::Forward),
            ],
            CountLevel::MetaFeature,
        )
        .unwrap();
        let hits = features.hits(
            &[Block::new(Some(0), 130, 170).unwrap()],
            StrandMode::Unstranded,
            false,
            false,
        );
        assert_eq!(hits[0].overlap_bases, 40);
    }

    #[test]
    fn strand_modes_filter_candidates() {
        let features = FeatureSet::build(
            vec![
                part("plus", 100, 200, Strand::Forward),
                part("minus", 100, 200, Strand::Reverse),
            ],
            CountLevel::MetaFeature,
        )
        .unwrap();
        let block = Block::new(Some(0), 120, 140).unwrap();
        assert_eq!(
            features.hits(
                std::slice::from_ref(&block),
                StrandMode::Sense,
                false,
                false,
            )[0]
            .row,
            0
        );
        assert_eq!(
            features.hits(&[block], StrandMode::Antisense, false, false)[0].row,
            1
        );
    }
}
