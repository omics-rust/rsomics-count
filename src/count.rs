use std::collections::{BTreeMap, HashMap};
use std::io;
use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};
use std::sync::mpsc;

use noodles::sam::header::record::value::map::header::{sort_order::QUERY_NAME, tag::SORT_ORDER};
use rsomics_common::{Result, RsomicsError};
use serde::Serialize;
use smallvec::SmallVec;

use crate::alignment::{NameMode, Reader, Segment};
use crate::annotation::{self, AnnotationOptions};
use crate::model::{Block, Blocks, FeatureSet, Hit, Hits, StrandMode};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum Contribution {
    Exclude,
    Full,
    Fractional,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SplitMode {
    Any,
    Only,
    Exclude,
}

#[derive(Clone, Debug)]
pub struct PairOptions {
    pub enabled: bool,
    pub require_both_mapped: bool,
    pub exclude_chimeric: bool,
    pub check_fragment_length: bool,
    pub min_fragment_length: u32,
    pub max_fragment_length: u32,
}

impl Default for PairOptions {
    fn default() -> Self {
        Self {
            enabled: false,
            require_both_mapped: false,
            exclude_chimeric: false,
            check_fragment_length: false,
            min_fragment_length: 50,
            max_fragment_length: 600,
        }
    }
}

#[derive(Clone, Debug)]
pub struct CountOptions {
    pub strands: Vec<StrandMode>,
    pub multi_overlap: Contribution,
    pub multi_map: Contribution,
    pub min_overlap: u64,
    pub frac_overlap: f64,
    pub frac_overlap_feature: f64,
    pub min_mapq: u8,
    pub primary_only: bool,
    pub ignore_duplicates: bool,
    pub split: SplitMode,
    pub pair: PairOptions,
    pub aliases: Option<PathBuf>,
    pub threads: Option<NonZeroUsize>,
}

impl Default for CountOptions {
    fn default() -> Self {
        Self {
            strands: vec![StrandMode::Unstranded],
            multi_overlap: Contribution::Exclude,
            multi_map: Contribution::Exclude,
            min_overlap: 1,
            frac_overlap: 0.0,
            frac_overlap_feature: 0.0,
            min_mapq: 0,
            primary_only: false,
            ignore_duplicates: false,
            split: SplitMode::Any,
            pair: PairOptions::default(),
            aliases: None,
            threads: None,
        }
    }
}

#[derive(Clone, Debug, Default, Eq, PartialEq, Serialize)]
pub struct CountSummary {
    pub assigned: u64,
    pub unassigned_unmapped: u64,
    pub unassigned_read_type: u64,
    pub unassigned_singleton: u64,
    pub unassigned_mapping_quality: u64,
    pub unassigned_chimera: u64,
    pub unassigned_fragment_length: u64,
    pub unassigned_duplicate: u64,
    pub unassigned_multi_mapping: u64,
    pub unassigned_secondary: u64,
    pub unassigned_non_split: u64,
    pub unassigned_split: u64,
    pub unassigned_no_features: u64,
    pub unassigned_overlapping_length: u64,
    pub unassigned_ambiguity: u64,
}

impl CountSummary {
    #[must_use]
    pub const fn total(&self) -> u64 {
        self.assigned
            + self.unassigned_unmapped
            + self.unassigned_read_type
            + self.unassigned_singleton
            + self.unassigned_mapping_quality
            + self.unassigned_chimera
            + self.unassigned_fragment_length
            + self.unassigned_duplicate
            + self.unassigned_multi_mapping
            + self.unassigned_secondary
            + self.unassigned_non_split
            + self.unassigned_split
            + self.unassigned_no_features
            + self.unassigned_overlapping_length
            + self.unassigned_ambiguity
    }

    fn merge(&mut self, other: &Self) {
        self.assigned += other.assigned;
        self.unassigned_unmapped += other.unassigned_unmapped;
        self.unassigned_read_type += other.unassigned_read_type;
        self.unassigned_singleton += other.unassigned_singleton;
        self.unassigned_mapping_quality += other.unassigned_mapping_quality;
        self.unassigned_chimera += other.unassigned_chimera;
        self.unassigned_fragment_length += other.unassigned_fragment_length;
        self.unassigned_duplicate += other.unassigned_duplicate;
        self.unassigned_multi_mapping += other.unassigned_multi_mapping;
        self.unassigned_secondary += other.unassigned_secondary;
        self.unassigned_non_split += other.unassigned_non_split;
        self.unassigned_split += other.unassigned_split;
        self.unassigned_no_features += other.unassigned_no_features;
        self.unassigned_overlapping_length += other.unassigned_overlapping_length;
        self.unassigned_ambiguity += other.unassigned_ambiguity;
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct LibraryReport {
    pub input: String,
    pub total: u64,
    pub assigned: u64,
}

#[derive(Clone, Debug, Serialize)]
pub struct CountReport {
    pub features: usize,
    pub libraries: Vec<LibraryReport>,
}

pub(crate) struct LibraryCounts {
    pub(crate) values: Vec<f64>,
    pub(crate) summary: CountSummary,
}

impl LibraryCounts {
    fn new(features: usize) -> Self {
        Self {
            values: vec![0.0; features],
            summary: CountSummary::default(),
        }
    }

    fn merge(&mut self, other: &Self) {
        for (value, addition) in self.values.iter_mut().zip(&other.values) {
            *value += addition;
        }
        self.summary.merge(&other.summary);
    }
}

pub fn count(
    inputs: &[PathBuf],
    annotation_path: &Path,
    counts_path: &Path,
    summary_path: &Path,
    annotation_options: &AnnotationOptions,
    options: &CountOptions,
) -> Result<CountReport> {
    validate(inputs, counts_path, summary_path, options)?;
    rsomics_common::reject_output_alias(
        counts_path,
        inputs
            .iter()
            .map(PathBuf::as_path)
            .chain(std::iter::once(annotation_path))
            .chain(options.aliases.as_deref())
            .chain(std::iter::once(summary_path)),
    )?;
    rsomics_common::reject_output_alias(
        summary_path,
        inputs
            .iter()
            .map(PathBuf::as_path)
            .chain(std::iter::once(annotation_path))
            .chain(options.aliases.as_deref())
            .chain(std::iter::once(counts_path)),
    )?;
    let features = annotation::load(annotation_path, annotation_options)?;
    let aliases = annotation::load_aliases(options.aliases.as_deref())?;
    let strands = expand_strands(&options.strands, inputs.len())?;
    let mut libraries = Vec::with_capacity(inputs.len());
    for (input, strand) in inputs.iter().zip(strands) {
        libraries.push(count_library(input, &features, &aliases, strand, options)?);
    }

    rsomics_common::write_atomic_pair(counts_path, summary_path, |counts, summary| {
        crate::output::write_counts(counts, &features, inputs, &libraries)?;
        crate::output::write_summary(summary, inputs, &libraries, options.split)
    })?;

    Ok(CountReport {
        features: features.rows().len(),
        libraries: inputs
            .iter()
            .zip(&libraries)
            .map(|(input, library)| LibraryReport {
                input: input.display().to_string(),
                total: library.summary.total(),
                assigned: library.summary.assigned,
            })
            .collect(),
    })
}

fn validate(
    inputs: &[PathBuf],
    counts_path: &Path,
    summary_path: &Path,
    options: &CountOptions,
) -> Result<()> {
    if inputs.is_empty() {
        return Err(RsomicsError::ConfigError(
            "at least one alignment input is required".to_owned(),
        ));
    }
    if inputs
        .iter()
        .filter(|path| path.as_path() == Path::new("-"))
        .count()
        > 1
    {
        return Err(RsomicsError::ConfigError(
            "standard input can name only one library".to_owned(),
        ));
    }
    if counts_path == Path::new("-") || summary_path == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "counts and summary outputs must be named files".to_owned(),
        ));
    }
    for (label, value) in [
        ("--frac-overlap", options.frac_overlap),
        ("--frac-overlap-feature", options.frac_overlap_feature),
    ] {
        if !value.is_finite() || !(0.0..=1.0).contains(&value) {
            return Err(RsomicsError::ConfigError(format!(
                "{label} must be finite and between 0 and 1"
            )));
        }
    }
    if options.min_overlap == 0 {
        return Err(RsomicsError::ConfigError(
            "--min-overlap must be at least 1".to_owned(),
        ));
    }
    if options.pair.min_fragment_length > options.pair.max_fragment_length {
        return Err(RsomicsError::ConfigError(
            "minimum fragment length exceeds maximum fragment length".to_owned(),
        ));
    }
    if options.pair.check_fragment_length && !options.pair.require_both_mapped {
        return Err(RsomicsError::ConfigError(
            "fragment-length checks require both mates to be mapped".to_owned(),
        ));
    }
    Ok(())
}

fn expand_strands(strands: &[StrandMode], libraries: usize) -> Result<Vec<StrandMode>> {
    match strands {
        [strand] => Ok(vec![*strand; libraries]),
        values if values.len() == libraries => Ok(values.to_vec()),
        values => Err(RsomicsError::ConfigError(format!(
            "--strand expects one value or one per library; found {} for {libraries} libraries",
            values.len()
        ))),
    }
}

fn count_library(
    input: &Path,
    features: &FeatureSet,
    aliases: &HashMap<String, String>,
    strand: StrandMode,
    options: &CountOptions,
) -> Result<LibraryCounts> {
    let threads = options
        .threads
        .unwrap_or_else(|| std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN));
    let parallel_count = threads.get() > 1;
    let mut reader = Reader::open(
        input,
        if parallel_count {
            Some(NonZeroUsize::MIN)
        } else {
            options.threads
        },
    )?;
    let header = reader.read_header()?;
    let mut result = LibraryCounts::new(features.rows().len());
    if options.pair.enabled {
        count_pairs(
            &mut reader,
            &header,
            CountContext {
                aliases,
                features,
                strand,
                options,
            },
            threads.get() - 1,
            &mut result,
        )?;
    } else if parallel_count {
        count_single_parallel(
            &mut reader,
            &header,
            CountContext {
                aliases,
                features,
                strand,
                options,
            },
            threads.get() - 1,
            &mut result,
        )?;
    } else {
        reader.visit(
            &header,
            aliases,
            &|name| features.chrom_id(name),
            NameMode::None,
            |segment| {
                process_unit(
                    Unit::single(segment),
                    features,
                    strand,
                    options,
                    &mut result,
                )
            },
        )?;
    }
    Ok(result)
}

fn count_pairs(
    reader: &mut Reader,
    header: &noodles::sam::Header,
    context: CountContext<'_>,
    workers: usize,
    result: &mut LibraryCounts,
) -> Result<()> {
    if is_queryname_sorted(header) && workers > 0 {
        count_queryname_pairs_parallel(reader, header, context, workers, result)
    } else if is_queryname_sorted(header) {
        let mut group_active = false;
        let mut buckets = PairBuckets::new();
        reader.visit(
            header,
            context.aliases,
            &|name| context.features.chrom_id(name),
            NameMode::Grouped,
            |segment| {
                if !segment.segmented || segment.first == segment.last {
                    result.summary.unassigned_read_type += 1;
                    return Ok(());
                }
                if !segment.name_present {
                    return Err(RsomicsError::InvalidInput(
                        "paired alignment has no read name".to_owned(),
                    ));
                }
                if segment.group_start {
                    if group_active {
                        process_pair_buckets(std::mem::take(&mut buckets), context, result)?;
                    }
                    group_active = true;
                }
                let hit_index = segment.hit_index.unwrap_or(1);
                insert_pair_segment(pair_bucket(&mut buckets, hit_index), segment, None)
            },
        )?;
        process_pair_buckets(buckets, context, result)
    } else {
        let mut pairs = BTreeMap::<(Vec<u8>, u32), PairBucket>::new();
        reader.visit(
            header,
            context.aliases,
            &|name| context.features.chrom_id(name),
            NameMode::Owned,
            |mut segment| {
                if !segment.segmented || segment.first == segment.last {
                    result.summary.unassigned_read_type += 1;
                    return Ok(());
                }
                let name = segment.name.take().ok_or_else(|| {
                    RsomicsError::InvalidInput("paired alignment has no read name".to_owned())
                })?;
                let hit_index = segment.hit_index.unwrap_or(1);
                insert_pair_segment(
                    pairs.entry((name.clone(), hit_index)).or_default(),
                    segment,
                    Some(&name),
                )
            },
        )?;
        for (_, pair) in pairs {
            process_pair(
                pair,
                context.features,
                context.strand,
                context.options,
                result,
            )?;
        }
        Ok(())
    }
}

fn count_queryname_pairs_parallel(
    reader: &mut Reader,
    header: &noodles::sam::Header,
    context: CountContext<'_>,
    workers: usize,
    result: &mut LibraryCounts,
) -> Result<()> {
    const BATCH_SIZE: usize = 4 * 1024;

    std::thread::scope(|scope| {
        let mut senders = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let (sender, receiver) = mpsc::sync_channel::<Vec<PairBucket>>(2);
            senders.push(sender);
            handles.push(scope.spawn(move || {
                let mut local = LibraryCounts::new(context.features.rows().len());
                for batch in receiver {
                    for pair in batch {
                        process_pair(
                            pair,
                            context.features,
                            context.strand,
                            context.options,
                            &mut local,
                        )?;
                    }
                }
                Ok::<_, RsomicsError>(local)
            }));
        }

        let mut group_active = false;
        let mut buckets = PairBuckets::new();
        let mut batch = Vec::with_capacity(BATCH_SIZE);
        let mut next_worker = 0;
        let read_result = reader.visit(
            header,
            context.aliases,
            &|name| context.features.chrom_id(name),
            NameMode::Grouped,
            |segment| {
                if !segment.segmented || segment.first == segment.last {
                    result.summary.unassigned_read_type += 1;
                    return Ok(());
                }
                if !segment.name_present {
                    return Err(RsomicsError::InvalidInput(
                        "paired alignment has no read name".to_owned(),
                    ));
                }
                if segment.group_start {
                    if group_active {
                        enqueue_pair_buckets(
                            std::mem::take(&mut buckets),
                            &mut batch,
                            &senders,
                            &mut next_worker,
                            BATCH_SIZE,
                        )?;
                    }
                    group_active = true;
                }
                let hit_index = segment.hit_index.unwrap_or(1);
                insert_pair_segment(pair_bucket(&mut buckets, hit_index), segment, None)
            },
        );
        if read_result.is_ok() {
            enqueue_pair_buckets(buckets, &mut batch, &senders, &mut next_worker, BATCH_SIZE)?;
            if !batch.is_empty() {
                send_batch(&senders, &mut next_worker, batch)?;
            }
        }
        drop(senders);

        let worker_result = merge_workers(handles, result);
        read_result?;
        worker_result
    })
}

fn enqueue_pair_buckets(
    mut buckets: PairBuckets,
    batch: &mut Vec<PairBucket>,
    senders: &[mpsc::SyncSender<Vec<PairBucket>>],
    next_worker: &mut usize,
    batch_size: usize,
) -> Result<()> {
    buckets.sort_unstable_by_key(|(hit_index, _)| *hit_index);
    for (_, pair) in buckets {
        batch.push(pair);
        if batch.len() == batch_size {
            let ready = std::mem::replace(batch, Vec::with_capacity(batch_size));
            send_batch(senders, next_worker, ready)?;
        }
    }
    Ok(())
}

fn send_batch<T>(
    senders: &[mpsc::SyncSender<Vec<T>>],
    next_worker: &mut usize,
    batch: Vec<T>,
) -> Result<()> {
    senders[*next_worker].send(batch).map_err(|_| {
        RsomicsError::Io(io::Error::new(
            io::ErrorKind::BrokenPipe,
            "counting worker stopped",
        ))
    })?;
    *next_worker = (*next_worker + 1) % senders.len();
    Ok(())
}

fn merge_workers(
    handles: Vec<std::thread::ScopedJoinHandle<'_, Result<LibraryCounts>>>,
    result: &mut LibraryCounts,
) -> Result<()> {
    let mut worker_error = None;
    for handle in handles {
        match handle.join() {
            Ok(Ok(local)) => result.merge(&local),
            Ok(Err(error)) if worker_error.is_none() => worker_error = Some(error),
            Err(_) if worker_error.is_none() => {
                worker_error = Some(RsomicsError::Io(io::Error::other(
                    "counting worker panicked",
                )));
            }
            _ => {}
        }
    }
    match worker_error {
        Some(error) => Err(error),
        None => Ok(()),
    }
}

fn is_queryname_sorted(header: &noodles::sam::Header) -> bool {
    header
        .header()
        .and_then(|value| value.other_fields().get(&SORT_ORDER))
        .is_some_and(|order| order == QUERY_NAME)
}

fn insert_pair_segment(
    bucket: &mut PairBucket,
    segment: Segment,
    name: Option<&[u8]>,
) -> Result<()> {
    let slot = if segment.first {
        &mut bucket.first
    } else {
        &mut bucket.last
    };
    if slot.replace(segment).is_some() {
        let identity = name.map_or_else(
            || "group".to_owned(),
            |name| format!("{:?}", String::from_utf8_lossy(name)),
        );
        return Err(RsomicsError::InvalidInput(format!(
            "paired alignment {identity} repeats the same segment and HI value"
        )));
    }
    Ok(())
}

fn process_pair_buckets(
    mut buckets: PairBuckets,
    context: CountContext<'_>,
    result: &mut LibraryCounts,
) -> Result<()> {
    buckets.sort_unstable_by_key(|(hit_index, _)| *hit_index);
    for (_, pair) in buckets {
        process_pair(
            pair,
            context.features,
            context.strand,
            context.options,
            result,
        )?;
    }
    Ok(())
}

type PairBuckets = SmallVec<[(u32, PairBucket); 2]>;

fn pair_bucket(buckets: &mut PairBuckets, hit_index: u32) -> &mut PairBucket {
    if let Some(index) = buckets
        .iter()
        .position(|(candidate, _)| *candidate == hit_index)
    {
        return &mut buckets[index].1;
    }
    buckets.push((hit_index, PairBucket::default()));
    &mut buckets.last_mut().unwrap().1
}

#[derive(Clone, Copy)]
struct CountContext<'a> {
    aliases: &'a HashMap<String, String>,
    features: &'a FeatureSet,
    strand: StrandMode,
    options: &'a CountOptions,
}

fn count_single_parallel(
    reader: &mut Reader,
    header: &noodles::sam::Header,
    context: CountContext<'_>,
    workers: usize,
    result: &mut LibraryCounts,
) -> Result<()> {
    const BATCH_SIZE: usize = 8 * 1024;

    std::thread::scope(|scope| {
        let mut senders = Vec::with_capacity(workers);
        let mut handles = Vec::with_capacity(workers);
        for _ in 0..workers {
            let (sender, receiver) = mpsc::sync_channel::<Vec<Segment>>(2);
            senders.push(sender);
            handles.push(scope.spawn(move || {
                let mut local = LibraryCounts::new(context.features.rows().len());
                for batch in receiver {
                    for segment in batch {
                        process_unit(
                            Unit::single(segment),
                            context.features,
                            context.strand,
                            context.options,
                            &mut local,
                        )?;
                    }
                }
                Ok::<_, RsomicsError>(local)
            }));
        }

        let mut batch = Vec::with_capacity(BATCH_SIZE);
        let mut next_worker = 0;
        let read_result = reader.visit(
            header,
            context.aliases,
            &|name| context.features.chrom_id(name),
            NameMode::None,
            |segment| {
                batch.push(segment);
                if batch.len() == BATCH_SIZE {
                    let ready = std::mem::replace(&mut batch, Vec::with_capacity(BATCH_SIZE));
                    send_batch(&senders, &mut next_worker, ready)?;
                }
                Ok(())
            },
        );
        if read_result.is_ok()
            && !batch.is_empty()
            && let Err(error) = send_batch(&senders, &mut next_worker, batch)
        {
            drop(senders);
            for handle in handles {
                let _ = handle.join();
            }
            return Err(error);
        }
        drop(senders);

        let worker_result = merge_workers(handles, result);
        read_result?;
        worker_result
    })
}

#[derive(Default)]
struct PairBucket {
    first: Option<Segment>,
    last: Option<Segment>,
}

fn process_pair(
    pair: PairBucket,
    features: &FeatureSet,
    strand: StrandMode,
    options: &CountOptions,
    result: &mut LibraryCounts,
) -> Result<()> {
    if options.pair.require_both_mapped
        && (pair.first.as_ref().is_none_or(|segment| segment.unmapped)
            || pair.last.as_ref().is_none_or(|segment| segment.unmapped))
    {
        result.summary.unassigned_singleton += 1;
        return Ok(());
    }
    if options.pair.exclude_chimeric
        && let (Some(first), Some(last)) = (pair.first.as_ref(), pair.last.as_ref())
        && !first.unmapped
        && !last.unmapped
        && (first.reference != last.reference || first.reverse == last.reverse)
    {
        result.summary.unassigned_chimera += 1;
        return Ok(());
    }
    let invalid_fragment_length = if options.pair.check_fragment_length {
        let length = pair
            .first
            .as_ref()
            .or(pair.last.as_ref())
            .map_or(0, |segment| segment.template_length.unsigned_abs());
        length < options.pair.min_fragment_length || length > options.pair.max_fragment_length
    } else {
        false
    };
    let mut unit = Unit::pair(pair);
    unit.invalid_fragment_length = invalid_fragment_length;
    process_unit(unit, features, strand, options, result)
}

struct Unit {
    block_groups: SmallVec<[BlockGroup; 2]>,
    unmapped: bool,
    secondary: bool,
    duplicate: bool,
    mapq: u8,
    nh: u32,
    split: bool,
    invalid_fragment_length: bool,
}

struct BlockGroup {
    blocks: Blocks,
    reverse: bool,
}

impl Unit {
    fn single(segment: Segment) -> Self {
        let mut block_groups = SmallVec::new();
        block_groups.push(BlockGroup {
            blocks: segment.blocks,
            reverse: segment.reverse,
        });
        Self {
            block_groups,
            unmapped: segment.unmapped,
            secondary: segment.secondary,
            duplicate: segment.duplicate,
            mapq: segment.mapq,
            nh: segment.nh,
            split: segment.split,
            invalid_fragment_length: false,
        }
    }

    fn pair(pair: PairBucket) -> Self {
        let segments = [pair.first, pair.last];
        let secondary = segments.iter().flatten().any(|segment| segment.secondary);
        let duplicate = segments.iter().flatten().any(|segment| segment.duplicate);
        let mapq = segments
            .iter()
            .flatten()
            .filter(|segment| !segment.unmapped)
            .map(|segment| segment.mapq)
            .min()
            .unwrap_or(0);
        let nh = segments
            .iter()
            .flatten()
            .map(|segment| segment.nh)
            .max()
            .unwrap_or(1);
        let split = segments.iter().flatten().any(|segment| segment.split);
        let mut block_groups = SmallVec::new();
        for segment in segments.into_iter().flatten() {
            if !segment.unmapped {
                block_groups.push(BlockGroup {
                    blocks: segment.blocks,
                    reverse: segment.reverse ^ segment.last,
                });
            }
        }
        Self {
            unmapped: block_groups.is_empty(),
            block_groups,
            secondary,
            duplicate,
            mapq,
            nh,
            split,
            invalid_fragment_length: false,
        }
    }
}

fn process_unit(
    unit: Unit,
    features: &FeatureSet,
    strand: StrandMode,
    options: &CountOptions,
    result: &mut LibraryCounts,
) -> Result<()> {
    if unit.unmapped {
        result.summary.unassigned_unmapped += 1;
        return Ok(());
    }
    if options.primary_only && unit.secondary {
        result.summary.unassigned_secondary += 1;
        return Ok(());
    }
    if unit.mapq < options.min_mapq {
        result.summary.unassigned_mapping_quality += 1;
        return Ok(());
    }
    if options.ignore_duplicates && unit.duplicate {
        result.summary.unassigned_duplicate += 1;
        return Ok(());
    }
    if unit.nh > 1 && options.multi_map == Contribution::Exclude {
        result.summary.unassigned_multi_mapping += 1;
        return Ok(());
    }
    match options.split {
        SplitMode::Only if !unit.split => {
            result.summary.unassigned_non_split += 1;
            return Ok(());
        }
        SplitMode::Exclude if unit.split => {
            result.summary.unassigned_split += 1;
            return Ok(());
        }
        _ => {}
    }

    let raw_hits = unit_hits(
        features,
        &unit.block_groups,
        strand,
        options.frac_overlap_feature > 0.0,
    );
    if raw_hits.is_empty() {
        result.summary.unassigned_no_features += 1;
        return Ok(());
    }
    let read_bases: u64 = unit
        .block_groups
        .iter()
        .flat_map(|group| &group.blocks)
        .map(Block::len)
        .sum();
    let overlap_hits: Hits = raw_hits
        .into_iter()
        .filter(|hit| {
            hit.overlap_bases >= options.min_overlap
                && ratio_at_least(hit.overlap_bases, read_bases, options.frac_overlap)
        })
        .collect();
    if overlap_hits.is_empty() {
        result.summary.unassigned_overlapping_length += 1;
        return Ok(());
    }
    let hits: Hits = overlap_hits
        .into_iter()
        .filter(|hit| {
            ratio_at_least(
                hit.feature_overlap_bases,
                hit.feature_length,
                options.frac_overlap_feature,
            )
        })
        .collect();
    if hits.is_empty() {
        result.summary.unassigned_no_features += 1;
        return Ok(());
    }
    if hits.len() > 1 && options.multi_overlap == Contribution::Exclude {
        result.summary.unassigned_ambiguity += 1;
        return Ok(());
    }
    if unit.invalid_fragment_length {
        result.summary.unassigned_fragment_length += 1;
        return Ok(());
    }

    let mapping_factor = if unit.nh > 1 && options.multi_map == Contribution::Fractional {
        1.0 / f64::from(unit.nh)
    } else {
        1.0
    };
    let overlap_factor = if hits.len() > 1 && options.multi_overlap == Contribution::Fractional {
        1.0 / hits.len() as f64
    } else {
        1.0
    };
    let contribution = mapping_factor * overlap_factor;
    for hit in hits {
        result.values[hit.row] += contribution;
    }
    result.summary.assigned += 1;
    Ok(())
}

fn ratio_at_least(numerator: u64, denominator: u64, threshold: f64) -> bool {
    threshold == 0.0 || numerator as f64 / denominator as f64 >= threshold
}

fn unit_hits(
    features: &FeatureSet,
    groups: &[BlockGroup],
    strand: StrandMode,
    track_feature_fraction: bool,
) -> Hits {
    let per_group: SmallVec<[Hits; 2]> = groups
        .iter()
        .map(|group| features.hits(&group.blocks, strand, group.reverse, track_feature_fraction))
        .collect();
    if per_group.len() <= 1 {
        return per_group.into_iter().next().unwrap_or_default();
    }

    let mut common: Option<std::collections::HashSet<usize>> = None;
    for hits in &per_group {
        let rows: std::collections::HashSet<_> = hits.iter().map(|hit| hit.row).collect();
        common = Some(match common {
            Some(current) => current.intersection(&rows).copied().collect(),
            None => rows,
        });
    }
    let common = common.unwrap_or_default();
    let restrict_to_common = !common.is_empty();
    let mut overlap = BTreeMap::<usize, Hit>::new();
    for hits in per_group {
        for hit in hits {
            if !restrict_to_common || common.contains(&hit.row) {
                overlap
                    .entry(hit.row)
                    .and_modify(|current| {
                        current.overlap_bases += hit.overlap_bases;
                        if u128::from(hit.feature_overlap_bases)
                            * u128::from(current.feature_length)
                            > u128::from(current.feature_overlap_bases)
                                * u128::from(hit.feature_length)
                        {
                            current.feature_overlap_bases = hit.feature_overlap_bases;
                            current.feature_length = hit.feature_length;
                        }
                    })
                    .or_insert(hit);
            }
        }
    }
    overlap.into_values().collect()
}

#[cfg(test)]
mod tests {
    use crate::model::{CountLevel, PartInput, Strand};

    use super::*;

    fn features() -> FeatureSet {
        FeatureSet::build(
            vec![
                PartInput {
                    id: "a".to_owned(),
                    chrom: "chr1".to_owned(),
                    start: 100,
                    end: 150,
                    strand: Strand::Forward,
                },
                PartInput {
                    id: "b".to_owned(),
                    chrom: "chr1".to_owned(),
                    start: 140,
                    end: 200,
                    strand: Strand::Reverse,
                },
            ],
            CountLevel::MetaFeature,
        )
        .unwrap()
    }

    fn unit() -> Unit {
        Unit {
            block_groups: smallvec::smallvec![BlockGroup {
                blocks: smallvec::smallvec![Block::new(Some(0), 145, 155).unwrap()],
                reverse: false,
            }],
            unmapped: false,
            secondary: false,
            duplicate: false,
            mapq: 60,
            nh: 1,
            split: false,
            invalid_fragment_length: false,
        }
    }

    #[test]
    fn ambiguity_is_one_accounted_unit() {
        let mut result = LibraryCounts {
            values: vec![0.0; 2],
            summary: CountSummary::default(),
        };
        process_unit(
            unit(),
            &features(),
            StrandMode::Unstranded,
            &CountOptions::default(),
            &mut result,
        )
        .unwrap();
        assert_eq!(result.summary.unassigned_ambiguity, 1);
        assert_eq!(result.summary.total(), 1);
    }

    #[test]
    fn fractional_multi_overlap_is_deterministic() {
        let mut result = LibraryCounts {
            values: vec![0.0; 2],
            summary: CountSummary::default(),
        };
        let options = CountOptions {
            multi_overlap: Contribution::Fractional,
            ..CountOptions::default()
        };
        process_unit(
            unit(),
            &features(),
            StrandMode::Unstranded,
            &options,
            &mut result,
        )
        .unwrap();
        assert_eq!(result.values, [0.5, 0.5]);
        assert_eq!(result.summary.assigned, 1);
    }
}
