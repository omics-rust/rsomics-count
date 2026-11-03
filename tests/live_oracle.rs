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

#[test]
fn featurecounts_2_1_1_contract() {
    let Some(oracle) = std::env::var_os("FEATURECOUNTS") else {
        return;
    };
    let version = Command::new(&oracle).arg("-v").output().unwrap();
    let version_text = format!(
        "{}{}",
        String::from_utf8_lossy(&version.stdout),
        String::from_utf8_lossy(&version.stderr)
    );
    assert!(version_text.contains("v2.1.1"), "{version_text}");

    let cases: &[(&str, &[&str], &[&str])] = &[
        ("default", &[], &[]),
        ("overlap", &["-O"], &["--multi-overlap", "full"]),
        (
            "overlap-fraction",
            &["-O", "--fraction"],
            &["--multi-overlap", "fractional"],
        ),
        (
            "multimap-fraction",
            &["-M", "--fraction"],
            &["--multi-map", "fractional"],
        ),
        (
            "combined-fraction",
            &["-O", "-M", "--fraction"],
            &["--multi-overlap", "fractional", "--multi-map", "fractional"],
        ),
        (
            "minimum-overlap",
            &["--minOverlap", "15"],
            &["--min-overlap", "15"],
        ),
        (
            "read-overlap-fraction",
            &["--fracOverlap", "0.75"],
            &["--frac-overlap", "0.75"],
        ),
        (
            "feature-overlap-fraction",
            &["--fracOverlapFeature", "0.25"],
            &["--frac-overlap-feature", "0.25"],
        ),
        ("split-only", &["--splitOnly"], &["--split", "only"]),
        (
            "nonsplit-only",
            &["--nonSplitOnly"],
            &["--split", "exclude"],
        ),
        ("sense", &["-s", "1"], &["--strand", "sense"]),
        ("antisense", &["-s", "2"], &["--strand", "antisense"]),
        (
            "filters",
            &["-Q", "10", "--ignoreDup"],
            &["--min-mapq", "10", "--ignore-duplicates"],
        ),
        ("feature", &["-f"], &["--level", "feature"]),
    ];
    for (name, oracle_args, ours_args) in cases {
        compare_case(&oracle, name, "single.sam", oracle_args, ours_args);
    }
    compare_case(
        &oracle,
        "paired",
        "paired.sam",
        &["-p", "--countReadPairs"],
        &["--paired"],
    );
    compare_case(
        &oracle,
        "paired-coordinate",
        "paired-coordinate.sam",
        &["-p", "--countReadPairs"],
        &["--paired"],
    );
    compare_case(
        &oracle,
        "paired-checks",
        "paired.sam",
        &["-p", "--countReadPairs", "-B", "-C"],
        &["--paired", "--require-both-mapped", "--exclude-chimeric"],
    );
    compare_case(
        &oracle,
        "paired-fragment-length",
        "paired.sam",
        &[
            "-p",
            "--countReadPairs",
            "-B",
            "-P",
            "-d",
            "100",
            "-D",
            "200",
        ],
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
    compare_case(
        &oracle,
        "paired-sense",
        "paired.sam",
        &["-p", "--countReadPairs", "-s", "1"],
        &["--paired", "--strand", "sense"],
    );
    compare_case(
        &oracle,
        "paired-multimap-full",
        "paired-multimap.sam",
        &["-p", "--countReadPairs", "-M"],
        &["--paired", "--multi-map", "full"],
    );
    compare_case(
        &oracle,
        "paired-multimap-fraction",
        "paired-multimap.sam",
        &["-p", "--countReadPairs", "-M", "--fraction"],
        &["--paired", "--multi-map", "fractional"],
    );
    compare_case(&oracle, "edge-cigar", "edge.sam", &[], &[]);
    compare_case(
        &oracle,
        "deletion-gap",
        "deletion-gap.sam",
        &["-O"],
        &["--multi-overlap", "full"],
    );
    compare_case(
        &oracle,
        "edge-primary",
        "edge.sam",
        &["--primary"],
        &["--primary-only"],
    );
}

fn compare_case(
    oracle: &std::ffi::OsStr,
    name: &str,
    input: &str,
    oracle_args: &[&str],
    ours_args: &[&str],
) {
    let directory = tempfile::tempdir().unwrap();
    let oracle_counts = directory.path().join(format!("{name}.oracle.tsv"));
    let ours_counts = directory.path().join(format!("{name}.rsomics.tsv"));
    let ours_summary = directory.path().join(format!("{name}.rsomics.summary"));
    let input = fixture(input);

    let result = Command::new(oracle)
        .args(oracle_args)
        .args(["-a", fixture("model.gtf").to_str().unwrap()])
        .args(["-o", oracle_counts.to_str().unwrap()])
        .arg(&input)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "oracle {name}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let result = Command::new(binary())
        .arg(&input)
        .args(["-a", fixture("model.gtf").to_str().unwrap()])
        .args(["-o", ours_counts.to_str().unwrap()])
        .args(["--summary", ours_summary.to_str().unwrap()])
        .args(["--threads", "1"])
        .args(ours_args)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "rsomics {name}: {}",
        String::from_utf8_lossy(&result.stderr)
    );

    assert_eq!(
        normalized_counts(&fs::read_to_string(&ours_counts).unwrap()),
        normalized_counts(&fs::read_to_string(&oracle_counts).unwrap()),
        "count mismatch for {name}"
    );
    assert_eq!(
        normalized_summary(&fs::read_to_string(&ours_summary).unwrap()),
        normalized_summary(
            &fs::read_to_string(format!("{}.summary", oracle_counts.display())).unwrap()
        ),
        "summary mismatch for {name}"
    );
}

fn normalized_counts(input: &str) -> Vec<Vec<String>> {
    input
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

fn normalized_summary(input: &str) -> Vec<(String, String)> {
    input
        .lines()
        .skip(1)
        .map(|line| {
            let fields: Vec<_> = line.split('\t').collect();
            (fields[0].to_owned(), fields[1].to_owned())
        })
        .collect()
}
