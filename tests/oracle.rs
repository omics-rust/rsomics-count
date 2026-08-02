use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-count"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn run(input: &str, extra: &[&str]) -> (String, String) {
    let directory = tempfile::tempdir().unwrap();
    let counts = directory.path().join("counts.tsv");
    let summary = directory.path().join("summary.tsv");
    let mut command = Command::new(binary());
    command
        .arg(fixture(input))
        .args(["-a", fixture("model.gtf").to_str().unwrap()])
        .args(["-o", counts.to_str().unwrap()])
        .args(["--summary", summary.to_str().unwrap()])
        .args(["--threads", "1"])
        .args(extra);
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    (
        fs::read_to_string(counts).unwrap(),
        fs::read_to_string(summary).unwrap(),
    )
}

fn count_rows(output: &str) -> Vec<Vec<&str>> {
    output
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').collect())
        .collect()
}

fn counts(output: &str) -> Vec<&str> {
    count_rows(output)
        .into_iter()
        .map(|row| *row.last().unwrap())
        .collect()
}

fn summary_values(output: &str) -> Vec<u64> {
    output
        .lines()
        .skip(1)
        .map(|line| line.split('\t').nth(1).unwrap().parse().unwrap())
        .collect()
}

#[test]
fn default_matches_featurecounts_2_1_1() {
    let (matrix, summary) = run("single.sam", &[]);
    assert_eq!(counts(&matrix), ["3", "1", "1"]);
    assert_eq!(count_rows(&matrix)[2][5], "100");
    assert_eq!(
        summary_values(&summary),
        [5, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 2]
    );
}

#[test]
fn overlap_fraction_strand_and_filters_match_featurecounts_2_1_1() {
    assert_eq!(
        counts(&run("single.sam", &["--multi-overlap", "full"]).0),
        ["5", "3", "1"]
    );
    assert_eq!(
        counts(&run("single.sam", &["--multi-map", "fractional"]).0),
        ["3.50", "1", "1"]
    );
    assert_eq!(
        counts(&run("single.sam", &["--multi-overlap", "fractional"],).0),
        ["4", "2", "1"]
    );
    assert_eq!(
        counts(&run("single.sam", &["--strand", "sense"]).0),
        ["5", "1", "1"]
    );
    assert_eq!(
        counts(&run("single.sam", &["--strand", "antisense"]).0),
        ["0", "2", "0"]
    );
    assert_eq!(
        counts(&run("single.sam", &["--min-mapq", "10", "--ignore-duplicates"]).0),
        ["1", "1", "1"]
    );
}

#[test]
fn feature_and_paired_modes_match_featurecounts_2_1_1() {
    assert_eq!(
        counts(&run("single.sam", &["--level", "feature"]).0),
        ["3", "0", "1", "0", "0"]
    );
    let (paired, summary) = run("paired.sam", &["--paired"]);
    assert_eq!(counts(&paired), ["2", "1", "0"]);
    assert_eq!(
        summary_values(&summary),
        [3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    );
    let (coordinate, coordinate_summary) = run("paired-coordinate.sam", &["--paired"]);
    assert_eq!(counts(&coordinate), counts(&paired));
    assert_eq!(
        summary_values(&coordinate_summary),
        summary_values(&summary)
    );
    let (checked, summary) = run(
        "paired.sam",
        &["--paired", "--require-both-mapped", "--exclude-chimeric"],
    );
    assert_eq!(counts(&checked), ["1", "1", "0"]);
    assert_eq!(
        summary_values(&summary),
        [2, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
    );
}

#[test]
fn overlap_threshold_and_split_modes_match_featurecounts_2_1_1() {
    let (matrix, summary) = run(
        "single.sam",
        &["--multi-overlap", "fractional", "--multi-map", "fractional"],
    );
    assert_eq!(counts(&matrix), ["4.50", "2", "1"]);
    assert_eq!(
        summary_values(&summary),
        [8, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0]
    );

    let (matrix, summary) = run("single.sam", &["--min-overlap", "15"]);
    assert_eq!(counts(&matrix), ["3", "1", "1"]);
    assert_eq!(
        summary_values(&summary),
        [5, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1]
    );

    let (matrix, summary) = run("single.sam", &["--frac-overlap", "0.75"]);
    assert_eq!(counts(&matrix), ["3", "2", "1"]);
    assert_eq!(
        summary_values(&summary),
        [6, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 1]
    );

    let (matrix, summary) = run("single.sam", &["--frac-overlap-feature", "0.25"]);
    assert_eq!(counts(&matrix), ["3", "1", "1"]);
    assert_eq!(
        summary_values(&summary),
        [5, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 2, 0, 1]
    );

    let (_, summary) = run("single.sam", &["--split", "only"]);
    assert_eq!(
        summary_values(&summary),
        [0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 7, 0, 0, 1]
    );
    let (matrix, summary) = run("single.sam", &["--split", "exclude"]);
    assert_eq!(counts(&matrix), ["3", "1", "1"]);
    assert!(summary.contains("Unassigned_Split\t1"));
    assert_eq!(
        summary_values(&summary),
        [5, 1, 0, 0, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1]
    );
}

#[test]
fn paired_filter_and_strand_modes_match_featurecounts_2_1_1() {
    let (matrix, summary) = run(
        "paired.sam",
        &[
            "--paired",
            "--require-both-mapped",
            "--check-fragment-length",
            "-d",
            "100",
            "-D",
            "200",
        ],
    );
    assert_eq!(counts(&matrix), ["1", "0", "0"]);
    assert_eq!(
        summary_values(&summary),
        [1, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1]
    );

    let (matrix, summary) = run("paired.sam", &["--paired", "--strand", "sense"]);
    assert_eq!(counts(&matrix), ["3", "1", "0"]);
    assert_eq!(
        summary_values(&summary),
        [4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );

    let (matrix, summary) = run("paired-multimap.sam", &["--paired", "--multi-map", "full"]);
    assert_eq!(counts(&matrix), ["1", "0", "1"]);
    assert_eq!(
        summary_values(&summary),
        [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );

    let (matrix, summary) = run(
        "paired-multimap.sam",
        &["--paired", "--multi-map", "fractional"],
    );
    assert_eq!(counts(&matrix), ["0.50", "0", "0.50"]);
    assert_eq!(
        summary_values(&summary),
        [2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    );
}

#[test]
fn cigar_and_primary_modes_match_featurecounts_2_1_1() {
    let (matrix, summary) = run("edge.sam", &[]);
    assert_eq!(counts(&matrix), ["10", "0", "0"]);
    assert_eq!(
        summary_values(&summary),
        [10, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0]
    );

    let (matrix, summary) = run("edge.sam", &["--primary-only"]);
    assert_eq!(counts(&matrix), ["9", "0", "0"]);
    assert_eq!(
        summary_values(&summary),
        [9, 1, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0]
    );
}
