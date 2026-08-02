use std::io::{BufWriter, Write};
use std::path::PathBuf;

use rsomics_common::{Result, RsomicsError};

use crate::count::{CountSummary, LibraryCounts, SplitMode};
use crate::model::FeatureSet;

type SummaryRow = (&'static str, fn(&CountSummary) -> u64);

pub(crate) fn write_counts(
    output: &mut std::fs::File,
    features: &FeatureSet,
    inputs: &[PathBuf],
    libraries: &[LibraryCounts],
) -> Result<()> {
    let mut output = BufWriter::with_capacity(256 * 1024, output);
    writeln!(
        output,
        "# Program:rsomics-count v{}",
        env!("CARGO_PKG_VERSION")
    )
    .map_err(RsomicsError::Io)?;
    write!(output, "Geneid\tChr\tStart\tEnd\tStrand\tLength").map_err(RsomicsError::Io)?;
    for input in inputs {
        write!(output, "\t{}", input.display()).map_err(RsomicsError::Io)?;
    }
    writeln!(output).map_err(RsomicsError::Io)?;

    for (row_index, row) in features.rows().iter().enumerate() {
        let (chrom, start, end, strand) = features.fields(row_index);
        write!(
            output,
            "{}\t{chrom}\t{start}\t{end}\t{strand}\t{}",
            row.id(),
            row.length()
        )
        .map_err(RsomicsError::Io)?;
        for library in libraries {
            let value = library.values[row_index];
            let rounded = value.round();
            if (value - rounded).abs() < 1e-9 {
                write!(output, "\t{rounded:.0}").map_err(RsomicsError::Io)?;
            } else {
                write!(output, "\t{value:.2}").map_err(RsomicsError::Io)?;
            }
        }
        writeln!(output).map_err(RsomicsError::Io)?;
    }
    output.flush().map_err(RsomicsError::Io)
}

pub(crate) fn write_summary(
    output: &mut std::fs::File,
    inputs: &[PathBuf],
    libraries: &[LibraryCounts],
    split: SplitMode,
) -> Result<()> {
    let mut output = BufWriter::new(output);
    write!(output, "Status").map_err(RsomicsError::Io)?;
    for input in inputs {
        write!(output, "\t{}", input.display()).map_err(RsomicsError::Io)?;
    }
    writeln!(output).map_err(RsomicsError::Io)?;

    let split_row: SummaryRow = match split {
        SplitMode::Exclude => ("Unassigned_Split", |summary| summary.unassigned_split),
        SplitMode::Any | SplitMode::Only => ("Unassigned_NonSplit", |summary| {
            summary.unassigned_non_split
        }),
    };
    let rows: [SummaryRow; 14] = [
        ("Assigned", |summary| summary.assigned),
        ("Unassigned_Unmapped", |summary| summary.unassigned_unmapped),
        ("Unassigned_Read_Type", |summary| {
            summary.unassigned_read_type
        }),
        ("Unassigned_Singleton", |summary| {
            summary.unassigned_singleton
        }),
        ("Unassigned_MappingQuality", |summary| {
            summary.unassigned_mapping_quality
        }),
        ("Unassigned_Chimera", |summary| summary.unassigned_chimera),
        ("Unassigned_FragmentLength", |summary| {
            summary.unassigned_fragment_length
        }),
        ("Unassigned_Duplicate", |summary| {
            summary.unassigned_duplicate
        }),
        ("Unassigned_MultiMapping", |summary| {
            summary.unassigned_multi_mapping
        }),
        ("Unassigned_Secondary", |summary| {
            summary.unassigned_secondary
        }),
        split_row,
        ("Unassigned_NoFeatures", |summary| {
            summary.unassigned_no_features
        }),
        ("Unassigned_Overlapping_Length", |summary| {
            summary.unassigned_overlapping_length
        }),
        ("Unassigned_Ambiguity", |summary| {
            summary.unassigned_ambiguity
        }),
    ];
    for (label, value) in rows {
        write!(output, "{label}").map_err(RsomicsError::Io)?;
        for library in libraries {
            write!(output, "\t{}", value(&library.summary)).map_err(RsomicsError::Io)?;
        }
        writeln!(output).map_err(RsomicsError::Io)?;
    }
    output.flush().map_err(RsomicsError::Io)
}
