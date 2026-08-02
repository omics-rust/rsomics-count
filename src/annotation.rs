use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use bstr::ByteSlice;
use flate2::read::MultiGzDecoder;
use noodles_gff as gff;
use noodles_gtf as gtf;
use rsomics_common::{Context, Result, RsomicsError};

use crate::model::{CountLevel, FeatureSet, PartInput, Strand};

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum AnnotationFormat {
    Gtf,
    Gff3,
    Saf,
}

#[derive(Clone, Debug)]
pub struct AnnotationOptions {
    pub format: AnnotationFormat,
    pub feature_types: Vec<String>,
    pub group_by: String,
    pub level: CountLevel,
}

pub(crate) fn load(path: &Path, options: &AnnotationOptions) -> Result<FeatureSet> {
    if options.feature_types.is_empty() {
        return Err(RsomicsError::ConfigError(
            "at least one annotation feature type is required".to_owned(),
        ));
    }
    if options.group_by.is_empty() {
        return Err(RsomicsError::ConfigError(
            "the grouping attribute cannot be empty".to_owned(),
        ));
    }
    let selected: HashSet<&str> = options.feature_types.iter().map(String::as_str).collect();
    let parts = match options.format {
        AnnotationFormat::Gtf => load_gtf(path, &selected, &options.group_by),
        AnnotationFormat::Gff3 => load_gff3(path, &selected, &options.group_by),
        AnnotationFormat::Saf => load_saf(path, &selected),
    }?;
    FeatureSet::build(parts, options.level)
        .rs_with_context(|| format!("loading annotation {}", path.display()))
}

fn load_gtf(path: &Path, selected: &HashSet<&str>, group_by: &str) -> Result<Vec<PartInput>> {
    let mut reader = gtf::io::Reader::new(open_text(path)?);
    let mut parts = Vec::new();
    for (record_index, result) in reader.record_bufs().enumerate() {
        let record = result.map_err(|error| annotation_error(path, record_index + 1, error))?;
        push_record(
            path,
            record_index + 1,
            &record,
            selected,
            group_by,
            &mut parts,
        )?;
    }
    Ok(parts)
}

fn load_gff3(path: &Path, selected: &HashSet<&str>, group_by: &str) -> Result<Vec<PartInput>> {
    let mut reader = gff::io::Reader::new(open_text(path)?);
    let mut parts = Vec::new();
    for (record_index, result) in reader.record_bufs().enumerate() {
        let record = result.map_err(|error| annotation_error(path, record_index + 1, error))?;
        push_record(
            path,
            record_index + 1,
            &record,
            selected,
            group_by,
            &mut parts,
        )?;
    }
    Ok(parts)
}

fn push_record(
    path: &Path,
    record_index: usize,
    record: &gff::feature::RecordBuf,
    selected: &HashSet<&str>,
    group_by: &str,
    parts: &mut Vec<PartInput>,
) -> Result<()> {
    let ty = record.ty().to_str().map_err(|error| {
        annotation_error(
            path,
            record_index,
            format!("feature type is not UTF-8: {error}"),
        )
    })?;
    if !selected.contains(ty) {
        return Ok(());
    }
    let raw_id = record
        .attributes()
        .get(group_by.as_bytes())
        .ok_or_else(|| {
            annotation_error(
                path,
                record_index,
                format!("missing grouping attribute {group_by:?}"),
            )
        })?;
    let id = raw_id.as_string().ok_or_else(|| {
        annotation_error(
            path,
            record_index,
            format!("grouping attribute {group_by:?} has multiple values"),
        )
    })?;
    let id = id.to_str().map_err(|error| {
        annotation_error(
            path,
            record_index,
            format!("grouping attribute {group_by:?} is not UTF-8: {error}"),
        )
    })?;
    if id.is_empty() {
        return Err(annotation_error(
            path,
            record_index,
            format!("grouping attribute {group_by:?} is empty"),
        ));
    }
    let chrom = record.reference_sequence_name().to_str().map_err(|error| {
        annotation_error(
            path,
            record_index,
            format!("reference sequence name is not UTF-8: {error}"),
        )
    })?;
    let strand = match record.strand() {
        gff::feature::record::Strand::Forward => Strand::Forward,
        gff::feature::record::Strand::Reverse => Strand::Reverse,
        gff::feature::record::Strand::None => Strand::Unknown,
        other => {
            return Err(annotation_error(
                path,
                record_index,
                format!("unsupported feature strand {other:?}"),
            ));
        }
    };
    parts.push(PartInput {
        id: id.to_owned(),
        chrom: chrom.to_owned(),
        start: usize::from(record.start()) as u64 - 1,
        end: usize::from(record.end()) as u64,
        strand,
    });
    Ok(())
}

fn load_saf(path: &Path, selected: &HashSet<&str>) -> Result<Vec<PartInput>> {
    if selected.len() != 1 || !selected.contains("exon") {
        return Err(RsomicsError::ConfigError(
            "--feature-type is not applicable to SAF annotations".to_owned(),
        ));
    }
    let reader = open_text(path)?;
    let mut lines = reader.lines().enumerate();
    let (header_line, header) = loop {
        let Some((index, result)) = lines.next() else {
            return Err(RsomicsError::InvalidInput(format!(
                "{}: empty SAF annotation",
                path.display()
            )));
        };
        let line = result.rs_with_context(|| format!("reading {}", path.display()))?;
        if !line.trim().is_empty() && !line.starts_with('#') {
            break (index + 1, line);
        }
    };
    if header.split('\t').collect::<Vec<_>>() != ["GeneID", "Chr", "Start", "End", "Strand"] {
        return Err(annotation_error(
            path,
            header_line,
            "expected SAF header GeneID, Chr, Start, End, Strand",
        ));
    }

    let mut parts = Vec::new();
    for (index, result) in lines {
        let line = result.rs_with_context(|| format!("reading {}", path.display()))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let line_number = index + 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() != 5 {
            return Err(annotation_error(
                path,
                line_number,
                format!("expected 5 SAF columns, found {}", fields.len()),
            ));
        }
        if fields[0].is_empty() || fields[1].is_empty() {
            return Err(annotation_error(
                path,
                line_number,
                "SAF GeneID and Chr must be nonempty",
            ));
        }
        let start: u64 = fields[2].parse().map_err(|_| {
            annotation_error(
                path,
                line_number,
                format!("invalid SAF start {:?}", fields[2]),
            )
        })?;
        let end: u64 = fields[3].parse().map_err(|_| {
            annotation_error(
                path,
                line_number,
                format!("invalid SAF end {:?}", fields[3]),
            )
        })?;
        if start == 0 || start > end {
            return Err(annotation_error(
                path,
                line_number,
                format!("invalid one-based SAF interval {start}..={end}"),
            ));
        }
        parts.push(PartInput {
            id: fields[0].to_owned(),
            chrom: fields[1].to_owned(),
            start: start - 1,
            end,
            strand: Strand::parse(fields[4])?,
        });
    }
    Ok(parts)
}

pub(crate) fn load_aliases(path: Option<&Path>) -> Result<HashMap<String, String>> {
    let Some(path) = path else {
        return Ok(HashMap::new());
    };
    let reader = open_text(path)?;
    let mut aliases = HashMap::new();
    for (index, result) in reader.lines().enumerate() {
        let line = result.rs_with_context(|| format!("reading {}", path.display()))?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<_> = line.split(',').map(str::trim).collect();
        if fields.len() != 2 || fields.iter().any(|field| field.is_empty()) {
            return Err(annotation_error(
                path,
                index + 1,
                "expected annotation,input chromosome alias pair",
            ));
        }
        match aliases.insert(fields[1].to_owned(), fields[0].to_owned()) {
            Some(previous) if previous != fields[0] => {
                return Err(annotation_error(
                    path,
                    index + 1,
                    format!("input chromosome {:?} has conflicting aliases", fields[1]),
                ));
            }
            _ => {}
        }
    }
    Ok(aliases)
}

fn open_text(path: &Path) -> Result<Box<dyn BufRead>> {
    let file =
        File::open(path).rs_with_context(|| format!("opening annotation {}", path.display()))?;
    let mut buffered = BufReader::new(file);
    let gzip = buffered
        .fill_buf()
        .rs_with_context(|| format!("reading annotation {}", path.display()))?
        .starts_with(&[0x1f, 0x8b]);
    if gzip {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(buffered))))
    } else {
        Ok(Box::new(buffered))
    }
}

fn annotation_error(path: &Path, record: usize, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{} record {record}: {error}", path.display()))
}

#[cfg(test)]
mod tests {
    use std::io::Write;

    use super::*;

    fn options(format: AnnotationFormat) -> AnnotationOptions {
        AnnotationOptions {
            format,
            feature_types: vec!["exon".to_owned()],
            group_by: "gene_id".to_owned(),
            level: CountLevel::MetaFeature,
        }
    }

    #[test]
    fn strict_saf_builds_union_lengths() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "GeneID\tChr\tStart\tEnd\tStrand").unwrap();
        writeln!(file, "g\tchr1\t101\t150\t+").unwrap();
        writeln!(file, "g\tchr1\t141\t200\t+").unwrap();
        let features = load(file.path(), &options(AnnotationFormat::Saf)).unwrap();
        assert_eq!(features.rows()[0].length(), 100);
    }

    #[test]
    fn malformed_saf_fails() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "GeneID\tChr\tStart\tEnd\tStrand").unwrap();
        writeln!(file, "g\tchr1\t0\t10\t+").unwrap();
        assert!(load(file.path(), &options(AnnotationFormat::Saf)).is_err());
    }

    #[test]
    fn aliases_map_input_names_to_annotation_names() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "chr1,1").unwrap();
        assert_eq!(load_aliases(Some(file.path())).unwrap()["1"], "chr1");
    }

    #[test]
    fn selected_records_require_the_grouping_attribute() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(
            file,
            "chr1\ttest\texon\t101\t120\t.\t+\t.\ttranscript_id \"t\";"
        )
        .unwrap();
        assert!(load(file.path(), &options(AnnotationFormat::Gtf)).is_err());
    }

    #[test]
    fn multiple_feature_types_share_the_declared_group() {
        let mut file = tempfile::NamedTempFile::new().unwrap();
        writeln!(file, "chr1\ttest\texon\t101\t120\t.\t+\t.\tgene_id \"g\";").unwrap();
        writeln!(file, "chr1\ttest\tCDS\t121\t140\t.\t+\t0\tgene_id \"g\";").unwrap();
        let features = load(
            file.path(),
            &AnnotationOptions {
                feature_types: vec!["exon".to_owned(), "CDS".to_owned()],
                ..options(AnnotationFormat::Gtf)
            },
        )
        .unwrap();
        assert_eq!(features.rows().len(), 1);
        assert_eq!(features.rows()[0].length(), 40);
    }
}
