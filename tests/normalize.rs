use std::fs;
use std::path::PathBuf;
use std::process::Command;

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-count"))
}

#[test]
fn embedded_lengths_normalize_to_tpm() {
    let directory = tempfile::tempdir().unwrap();
    let counts = directory.path().join("counts.tsv");
    let output = directory.path().join("tpm.tsv");
    fs::write(
        &counts,
        "# provenance\nGeneid\tChr\tStart\tEnd\tStrand\tLength\ts1\ng1\tchr1\t1\t100\t+\t100\t3\ng2\tchr1\t1\t80\t+\t80\t1\ng3\tchr1\t1\t100\t+\t100\t1\n",
    )
    .unwrap();
    let result = Command::new(binary())
        .args(["normalize", counts.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let sum: f64 = fs::read_to_string(output)
        .unwrap()
        .lines()
        .skip(1)
        .map(|line| line.split('\t').nth(1).unwrap().parse::<f64>().unwrap())
        .sum();
    assert!((sum - 1_000_000.0).abs() < 1e-4);
}

#[test]
fn malformed_matrix_preserves_existing_output() {
    let directory = tempfile::tempdir().unwrap();
    let counts = directory.path().join("counts.tsv");
    let lengths = directory.path().join("lengths.tsv");
    let output = directory.path().join("normalized.tsv");
    fs::write(&counts, "gene\ts1\ng1\tbad\n").unwrap();
    fs::write(&lengths, "g1\t100\n").unwrap();
    fs::write(&output, "existing\n").unwrap();
    let result = Command::new(binary())
        .args(["normalize", counts.to_str().unwrap()])
        .args(["-l", lengths.to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read_to_string(output).unwrap(), "existing\n");
}

#[test]
fn fpkm_and_upper_quartile_match_frozen_references() {
    let directory = tempfile::tempdir().unwrap();
    let counts = directory.path().join("counts.tsv");
    let lengths = directory.path().join("lengths.tsv");
    fs::write(&counts, "gene\ts1\ts2\ng1\t10\t0\ng2\t20\t30\ng3\t0\t10\n").unwrap();
    fs::write(&lengths, "gene\tLength\ng1\t1000\ng2\t2000\ng3\t500\n").unwrap();

    for (mode, expected) in [
        (
            "fpkm",
            [
                [333_333.333_333, 0.0],
                [333_333.333_333, 375_000.0],
                [0.0, 500_000.0],
            ],
        ),
        (
            "fpkm-uq",
            [
                [571_428.571_429, 0.0],
                [571_428.571_429, 600_000.0],
                [0.0, 800_000.0],
            ],
        ),
    ] {
        let output = directory.path().join(format!("{mode}.tsv"));
        let result = Command::new(binary())
            .args(["normalize", counts.to_str().unwrap()])
            .arg("-l")
            .arg(&lengths)
            .arg("-o")
            .arg(&output)
            .args(["--mode", mode])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        let actual: Vec<Vec<f64>> = fs::read_to_string(output)
            .unwrap()
            .lines()
            .skip(1)
            .map(|line| {
                line.split('\t')
                    .skip(1)
                    .map(|value| value.parse().unwrap())
                    .collect()
            })
            .collect();
        for (actual_row, expected_row) in actual.iter().zip(expected) {
            for (actual, expected) in actual_row.iter().zip(expected_row) {
                assert!((actual - expected).abs() < 1e-5, "{actual} != {expected}");
            }
        }
    }
}
