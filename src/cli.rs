use std::ffi::OsString;
use std::num::NonZeroU64;
use std::path::{Path, PathBuf};

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{OutputArgs, Result, RsomicsError, ThreadArgs, ToolMeta};
use serde::Serialize;

use rsomics_count::{
    AnnotationFormat, AnnotationOptions, Contribution, CountLevel, CountOptions, CountReport,
    NormalizeMode, NormalizeOptions, NormalizeReport, PairOptions, SplitMode, StrandMode, count,
    normalize,
};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-count",
    version,
    about = "Assign alignments to genomic features and normalize count matrices",
    arg_required_else_help = true,
    subcommand_negates_reqs = true
)]
pub struct Cli {
    #[command(subcommand)]
    command: Option<Command>,

    #[command(flatten)]
    count: CountArgs,

    #[command(flatten)]
    threads: ThreadArgs,

    #[command(flatten)]
    pub output: OutputArgs,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Normalize a strict gene-by-sample count matrix.
    Normalize(NormalizeArgs),
}

#[derive(Debug, Args)]
struct CountArgs {
    /// SAM or BAM libraries; use - for a single standard-input library.
    #[arg(required = true, num_args = 1.., value_name = "ALIGNMENT")]
    inputs: Vec<PathBuf>,

    /// GTF, GFF3, or SAF annotation; gzip is detected by content.
    #[arg(short = 'a', long, required = true, value_name = "ANNOTATION")]
    annotation: Option<PathBuf>,

    /// Transactional featureCounts-compatible count matrix.
    #[arg(short, long, required = true, value_name = "COUNTS")]
    output: Option<PathBuf>,

    /// Transactional assignment summary; defaults to COUNTS.summary.
    #[arg(long, value_name = "SUMMARY")]
    summary: Option<PathBuf>,

    /// Annotation grammar.
    #[arg(short = 'F', long, value_enum, default_value_t = AnnotationFormatArg::Gtf)]
    annotation_format: AnnotationFormatArg,

    /// Selected annotation feature types.
    #[arg(long, value_delimiter = ',', default_value = "exon")]
    feature_type: Vec<String>,

    /// Attribute that groups features into meta-features.
    #[arg(short = 'g', long, default_value = "gene_id")]
    group_by: String,

    /// Count grouped meta-features or individual annotation records.
    #[arg(long, value_enum, default_value_t = LevelArg::Meta)]
    level: LevelArg,

    /// One strand mode for all libraries or a comma-separated mode per library.
    #[arg(
        short = 's',
        long,
        value_enum,
        value_delimiter = ',',
        default_value = "unstranded"
    )]
    strand: Vec<StrandArg>,

    /// Assignment of reads overlapping multiple targets.
    #[arg(long, value_enum, default_value_t = ContributionArg::Exclude)]
    multi_overlap: ContributionArg,

    /// Assignment of alignments with NH greater than one.
    #[arg(short = 'M', long, value_enum, default_value_t = ContributionArg::Exclude)]
    multi_map: ContributionArg,

    /// Minimum overlapping reference bases.
    #[arg(long, default_value_t = NonZeroU64::new(1).unwrap())]
    min_overlap: NonZeroU64,

    /// Minimum target overlap relative to original CIGAR read or fragment length.
    #[arg(long, default_value_t = 0.0, value_parser = ratio)]
    frac_overlap: f64,

    /// Minimum fraction of target bases overlapped.
    #[arg(long, default_value_t = 0.0, value_parser = ratio)]
    frac_overlap_feature: f64,

    /// Minimum mapping quality.
    #[arg(short = 'Q', long, default_value_t = 0)]
    min_mapq: u8,

    /// Exclude alignments carrying the secondary flag.
    #[arg(long)]
    primary_only: bool,

    /// Exclude PCR and optical duplicate alignments.
    #[arg(long)]
    ignore_duplicates: bool,

    /// Include any, only split, or only nonsplit alignments.
    #[arg(long, value_enum, default_value_t = SplitArg::Any)]
    split: SplitArg,

    /// Count read pairs as fragments.
    #[arg(short = 'p', long)]
    paired: bool,

    /// Require both mates to be mapped.
    #[arg(short = 'B', long, requires = "paired")]
    require_both_mapped: bool,

    /// Exclude fragments whose mapped mates use different references or orientations.
    #[arg(short = 'C', long, requires = "paired")]
    exclude_chimeric: bool,

    /// Enforce the configured fragment-length interval.
    #[arg(short = 'P', long, requires = "require_both_mapped")]
    check_fragment_length: bool,

    /// Minimum absolute template length.
    #[arg(
        short = 'd',
        long,
        default_value_t = 50,
        requires = "check_fragment_length"
    )]
    min_fragment_length: u32,

    /// Maximum absolute template length.
    #[arg(
        short = 'D',
        long,
        default_value_t = 600,
        requires = "check_fragment_length"
    )]
    max_fragment_length: u32,

    /// CSV pairs mapping annotation chromosome names to alignment names.
    #[arg(short = 'A', long, value_name = "ALIASES")]
    chromosome_aliases: Option<PathBuf>,
}

#[derive(Debug, Args)]
struct NormalizeArgs {
    /// Gene-by-sample count matrix.
    counts: PathBuf,

    /// Transactional normalized matrix.
    #[arg(short, long)]
    output: PathBuf,

    /// Separate gene and positive-length table when counts has no Length column.
    #[arg(short = 'l', long)]
    lengths: Option<PathBuf>,

    /// Normalization arithmetic.
    #[arg(short, long, value_enum, default_value_t = NormalizeModeArg::Tpm)]
    mode: NormalizeModeArg,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum AnnotationFormatArg {
    Gtf,
    Gff3,
    Saf,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum LevelArg {
    Meta,
    Feature,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum StrandArg {
    Unstranded,
    Sense,
    Antisense,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum ContributionArg {
    Exclude,
    Full,
    Fractional,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum SplitArg {
    Any,
    Only,
    Exclude,
}

#[derive(Clone, Copy, Debug, ValueEnum)]
enum NormalizeModeArg {
    Tpm,
    Fpkm,
    FpkmUq,
}

#[derive(Debug, Serialize)]
#[serde(tag = "operation", rename_all = "kebab-case")]
pub enum CommandReport {
    Count(CountReport),
    Normalize(NormalizeReport),
}

impl Cli {
    pub fn execute(self) -> Result<CommandReport> {
        match self.command {
            Some(Command::Normalize(args)) => execute_normalize(args).map(CommandReport::Normalize),
            None => self.execute_count().map(CommandReport::Count),
        }
    }

    fn execute_count(self) -> Result<CountReport> {
        let annotation = self.count.annotation.ok_or_else(|| {
            RsomicsError::ConfigError("--annotation is required for counting".to_owned())
        })?;
        let output = self.count.output.ok_or_else(|| {
            RsomicsError::ConfigError("--output is required for counting".to_owned())
        })?;
        let summary = self.count.summary.unwrap_or_else(|| summary_path(&output));
        let mut protected: Vec<&Path> = self.count.inputs.iter().map(PathBuf::as_path).collect();
        protected.push(annotation.as_path());
        if let Some(aliases) = self.count.chromosome_aliases.as_deref() {
            protected.push(aliases);
        }
        protected.push(summary.as_path());
        rsomics_common::reject_output_alias(&output, protected)?;
        rsomics_common::reject_output_alias(
            &summary,
            self.count
                .inputs
                .iter()
                .map(PathBuf::as_path)
                .chain(std::iter::once(annotation.as_path()))
                .chain(self.count.chromosome_aliases.as_deref()),
        )?;

        count(
            &self.count.inputs,
            &annotation,
            &output,
            &summary,
            &AnnotationOptions {
                format: self.count.annotation_format.into(),
                feature_types: self.count.feature_type,
                group_by: self.count.group_by,
                level: self.count.level.into(),
            },
            &CountOptions {
                strands: self.count.strand.into_iter().map(Into::into).collect(),
                multi_overlap: self.count.multi_overlap.into(),
                multi_map: self.count.multi_map.into(),
                min_overlap: self.count.min_overlap.get(),
                frac_overlap: self.count.frac_overlap,
                frac_overlap_feature: self.count.frac_overlap_feature,
                min_mapq: self.count.min_mapq,
                primary_only: self.count.primary_only,
                ignore_duplicates: self.count.ignore_duplicates,
                split: self.count.split.into(),
                pair: PairOptions {
                    enabled: self.count.paired,
                    require_both_mapped: self.count.require_both_mapped,
                    exclude_chimeric: self.count.exclude_chimeric,
                    check_fragment_length: self.count.check_fragment_length,
                    min_fragment_length: self.count.min_fragment_length,
                    max_fragment_length: self.count.max_fragment_length,
                },
                aliases: self.count.chromosome_aliases,
                threads: self.threads.requested(),
            },
        )
    }
}

fn execute_normalize(args: NormalizeArgs) -> Result<NormalizeReport> {
    rsomics_common::reject_output_alias(
        &args.output,
        std::iter::once(args.counts.as_path()).chain(args.lengths.as_deref()),
    )?;
    normalize(
        &args.counts,
        &args.output,
        &NormalizeOptions {
            mode: args.mode.into(),
            lengths: args.lengths,
        },
    )
}

fn summary_path(output: &Path) -> PathBuf {
    let mut path = OsString::from(output.as_os_str());
    path.push(".summary");
    PathBuf::from(path)
}

fn ratio(value: &str) -> std::result::Result<f64, String> {
    let ratio: f64 = value
        .parse()
        .map_err(|_| format!("expected a ratio from 0 through 1, found {value:?}"))?;
    if ratio.is_finite() && (0.0..=1.0).contains(&ratio) {
        Ok(ratio)
    } else {
        Err(format!(
            "expected a ratio from 0 through 1, found {value:?}"
        ))
    }
}

impl From<AnnotationFormatArg> for AnnotationFormat {
    fn from(value: AnnotationFormatArg) -> Self {
        match value {
            AnnotationFormatArg::Gtf => Self::Gtf,
            AnnotationFormatArg::Gff3 => Self::Gff3,
            AnnotationFormatArg::Saf => Self::Saf,
        }
    }
}

impl From<LevelArg> for CountLevel {
    fn from(value: LevelArg) -> Self {
        match value {
            LevelArg::Meta => Self::MetaFeature,
            LevelArg::Feature => Self::Feature,
        }
    }
}

impl From<StrandArg> for StrandMode {
    fn from(value: StrandArg) -> Self {
        match value {
            StrandArg::Unstranded => Self::Unstranded,
            StrandArg::Sense => Self::Sense,
            StrandArg::Antisense => Self::Antisense,
        }
    }
}

impl From<ContributionArg> for Contribution {
    fn from(value: ContributionArg) -> Self {
        match value {
            ContributionArg::Exclude => Self::Exclude,
            ContributionArg::Full => Self::Full,
            ContributionArg::Fractional => Self::Fractional,
        }
    }
}

impl From<SplitArg> for SplitMode {
    fn from(value: SplitArg) -> Self {
        match value {
            SplitArg::Any => Self::Any,
            SplitArg::Only => Self::Only,
            SplitArg::Exclude => Self::Exclude,
        }
    }
}

impl From<NormalizeModeArg> for NormalizeMode {
    fn from(value: NormalizeModeArg) -> Self {
        match value {
            NormalizeModeArg::Tpm => Self::Tpm,
            NormalizeModeArg::Fpkm => Self::Fpkm,
            NormalizeModeArg::FpkmUq => Self::FpkmUq,
        }
    }
}

#[cfg(test)]
mod tests {
    use clap::CommandFactory;

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        rsomics_help::command::<Cli>().debug_assert();
        Cli::command().debug_assert();
    }
}
