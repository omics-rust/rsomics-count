use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use flate2::{Compression, write::GzEncoder};
use noodles::sam::alignment::io::Write as _;
use noodles::{bam, sam};

fn binary() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-count"))
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/fixtures")
        .join(name)
}

fn output_counts(path: &Path) -> Vec<Vec<String>> {
    fs::read_to_string(path)
        .unwrap()
        .lines()
        .filter(|line| !line.starts_with('#'))
        .skip(1)
        .map(|line| line.split('\t').map(str::to_owned).collect())
        .collect()
}

fn sam_to_bam(source: &Path, destination: &Path) {
    let mut reader = sam::io::Reader::new(BufReader::new(File::open(source).unwrap()));
    let header = reader.read_header().unwrap();
    let records = reader
        .record_bufs(&header)
        .collect::<std::io::Result<Vec<_>>>()
        .unwrap();
    let mut writer = bam::io::Writer::new(File::create(destination).unwrap());
    writer.write_header(&header).unwrap();
    for record in records {
        writer.write_alignment_record(&header, &record).unwrap();
    }
    writer.try_finish().unwrap();
}

#[test]
fn invalid_projection_options_preserve_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("counts.tsv");
    fs::write(&output, "keep").unwrap();
    for args in [["--read2pos", "4"], ["--read-extension5", "2147483648"]] {
        let result = Command::new(binary())
            .arg(fixture("projection.sam"))
            .arg("-a")
            .arg(fixture("projection.saf"))
            .args(["-F", "saf"])
            .arg("-o")
            .arg(&output)
            .args(args)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert_eq!(fs::read_to_string(&output).unwrap(), "keep");
    }
}

#[test]
fn projection_preserves_section_endpoints_and_original_fragment_length() {
    let directory = tempfile::tempdir().unwrap();
    let single = fixture("projection.sam");
    let bam = directory.path().join("projection.bam");
    sam_to_bam(&single, &bam);
    let paired = fixture("projection-paired.sam");
    for (name, input, annotation, extra, expected) in [
        (
            "extensions",
            &single,
            "projection.saf",
            vec![],
            vec!["2", "2", "2", "1", "1"],
        ),
        (
            "five",
            &bam,
            "projection.saf",
            vec!["--read2pos", "5"],
            vec!["2", "0", "2", "0", "1"],
        ),
        (
            "three",
            &single,
            "projection.saf",
            vec!["--read2pos", "3"],
            vec!["1", "2", "0", "1", "1"],
        ),
        (
            "paired-low",
            &paired,
            "projection-paired.saf",
            vec![
                "--paired",
                "--strand",
                "sense",
                "--read2pos",
                "5",
                "--frac-overlap",
                ".05",
            ],
            vec!["1", "1"],
        ),
        (
            "paired-high",
            &paired,
            "projection-paired.saf",
            vec![
                "--paired",
                "--strand",
                "sense",
                "--read2pos",
                "5",
                "--frac-overlap",
                ".1",
            ],
            vec!["0", "0"],
        ),
    ] {
        let output = directory.path().join(format!("{name}.tsv"));
        let result = Command::new(binary())
            .arg(input)
            .arg("-a")
            .arg(fixture(annotation))
            .arg("-o")
            .arg(&output)
            .args([
                "-F",
                "saf",
                "--threads",
                "2",
                "--multi-overlap",
                "full",
                "--read-extension5",
                "10",
                "--read-extension3",
                "7",
            ])
            .args(extra)
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            output_counts(&output)
                .iter()
                .map(|row| row.last().unwrap().as_str())
                .collect::<Vec<_>>(),
            expected,
            "{name}"
        );
    }
}

#[test]
fn read_fraction_uses_original_cigar_and_fragment_unions() {
    let directory = tempfile::tempdir().unwrap();
    let single = fixture("read-fraction-single.sam");
    let paired = fixture("read-fraction-paired.sam");
    let bam = directory.path().join("paired.bam");
    sam_to_bam(&paired, &bam);
    for (name, input, fraction, is_paired, expected) in [
        (
            "single",
            &single,
            ".45",
            false,
            ["1", "0", "0", "0", "0", "0"],
        ),
        (
            "paired-low",
            &bam,
            ".3",
            true,
            ["0", "1", "1", "1", "1", "1"],
        ),
        (
            "paired-mid",
            &paired,
            ".45",
            true,
            ["0", "0", "1", "0", "1", "1"],
        ),
        (
            "paired-high",
            &bam,
            ".7",
            true,
            ["0", "0", "0", "0", "0", "0"],
        ),
    ] {
        let output = directory.path().join(format!("{name}.tsv"));
        let mut command = Command::new(binary());
        command
            .arg(input)
            .arg("-a")
            .arg(fixture("read-fraction.saf"))
            .arg("-o")
            .arg(&output)
            .args(["-F", "saf", "--threads", "2", "--frac-overlap", fraction]);
        if is_paired {
            command.arg("--paired");
        }
        let result = command.output().unwrap();
        assert!(
            result.status.success(),
            "{name}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            output_counts(&output)
                .iter()
                .map(|row| row.last().unwrap().as_str())
                .collect::<Vec<_>>(),
            expected,
            "{name}"
        );
    }
}

#[test]
fn deletion_gaps_do_not_overlap_features_in_sam_or_bam() {
    let directory = tempfile::tempdir().unwrap();
    let sam = fixture("deletion-gap.sam");
    let bam = directory.path().join("deletion-gap.bam");
    sam_to_bam(&sam, &bam);
    for input in [&sam, &bam] {
        let output = directory.path().join(format!(
            "{}.tsv",
            input.extension().unwrap().to_string_lossy()
        ));
        let result = Command::new(binary())
            .arg(input)
            .arg("-a")
            .arg(fixture("model.gtf"))
            .arg("-o")
            .arg(&output)
            .args(["--threads", "2", "--multi-overlap", "full"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            output_counts(&output)
                .iter()
                .map(|row| row.last().unwrap().as_str())
                .collect::<Vec<_>>(),
            ["0", "0", "0"],
            "{}",
            input.display()
        );
    }
}

#[test]
fn bam_uses_parallel_reader_and_matches_sam() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("single.bam");
    let output = directory.path().join("counts.tsv");
    sam_to_bam(&fixture("single.sam"), &bam);
    let result = Command::new(binary())
        .arg(&bam)
        .args(["-a", fixture("model.gtf").to_str().unwrap()])
        .args(["-o", output.to_str().unwrap()])
        .args(["--threads", "2"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let rows = output_counts(&output);
    assert_eq!(
        rows.iter()
            .map(|row| row.last().unwrap().as_str())
            .collect::<Vec<_>>(),
        ["3", "1", "1"]
    );
}

#[test]
fn parallel_fractional_counting_is_byte_stable() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("single.bam");
    sam_to_bam(&fixture("single.sam"), &bam);
    let mut outputs = Vec::new();
    for threads in [1, 4] {
        let counts = directory.path().join(format!("counts-{threads}.tsv"));
        let summary = directory.path().join(format!("summary-{threads}.tsv"));
        let result = Command::new(binary())
            .arg(&bam)
            .arg("-a")
            .arg(fixture("model.gtf"))
            .arg("-o")
            .arg(&counts)
            .arg("--summary")
            .arg(&summary)
            .args(["--threads", &threads.to_string()])
            .args(["--multi-overlap", "fractional"])
            .args(["--multi-map", "fractional"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        outputs.push((fs::read(counts).unwrap(), fs::read(summary).unwrap()));
    }
    assert_eq!(outputs[0], outputs[1]);
}

#[test]
fn raw_bam_cigar_and_flag_projection_matches_sam() {
    let directory = tempfile::tempdir().unwrap();
    let sam = fixture("edge.sam");
    let bam = directory.path().join("edge.bam");
    sam_to_bam(&sam, &bam);
    let mut outputs = Vec::new();
    for input in [&sam, &bam] {
        let counts = directory.path().join(format!(
            "{}.tsv",
            input.extension().unwrap().to_string_lossy()
        ));
        let summary = directory.path().join(format!(
            "{}.summary",
            input.extension().unwrap().to_string_lossy()
        ));
        let result = Command::new(binary())
            .arg(input)
            .arg("-a")
            .arg(fixture("model.gtf"))
            .arg("-o")
            .arg(&counts)
            .arg("--summary")
            .arg(&summary)
            .args(["--threads", "2"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        outputs.push((
            output_counts(&counts),
            fs::read_to_string(summary)
                .unwrap()
                .lines()
                .skip(1)
                .map(str::to_owned)
                .collect::<Vec<_>>(),
        ));
    }
    assert_eq!(outputs[0], outputs[1]);
}

#[test]
fn coordinate_sorted_bam_pairs_by_name() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("paired-coordinate.bam");
    sam_to_bam(&fixture("paired-coordinate.sam"), &bam);
    let counts = directory.path().join("counts.tsv");
    let result = Command::new(binary())
        .arg(&bam)
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-o")
        .arg(&counts)
        .args(["--threads", "2", "--paired"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        output_counts(&counts)
            .iter()
            .map(|row| row.last().unwrap().as_str())
            .collect::<Vec<_>>(),
        ["2", "1", "0"]
    );
}

#[test]
fn sam_standard_input_uses_the_same_assignment_path() {
    let directory = tempfile::tempdir().unwrap();
    let counts = directory.path().join("counts.tsv");
    let summary = directory.path().join("summary.tsv");
    let mut child = Command::new(binary())
        .arg("-")
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-o")
        .arg(&counts)
        .arg("--summary")
        .arg(&summary)
        .args(["--threads", "4"])
        .stdin(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture("single.sam")).unwrap())
        .unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        output_counts(&counts)
            .iter()
            .map(|row| row.last().unwrap().as_str())
            .collect::<Vec<_>>(),
        ["3", "1", "1"]
    );
}

#[test]
fn multiple_libraries_share_one_matrix_and_summary() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("counts.tsv");
    let summary = directory.path().join("summary.tsv");
    let result = Command::new(binary())
        .args([fixture("single.sam"), fixture("single.sam")])
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-o")
        .arg(&output)
        .arg("--summary")
        .arg(&summary)
        .args(["--threads", "1"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let rows = output_counts(&output);
    assert_eq!(&rows[0][6..], ["3", "3"]);
    assert_eq!(
        fs::read_to_string(summary)
            .unwrap()
            .lines()
            .next()
            .unwrap()
            .split('\t')
            .count(),
        3
    );
}

#[test]
fn saf_gff3_gzip_and_chromosome_aliases_are_real_paths() {
    let directory = tempfile::tempdir().unwrap();
    for (annotation, format) in [("model.saf", "saf"), ("model.gff3", "gff3")] {
        let output = directory.path().join(format!("{format}.tsv"));
        let result = Command::new(binary())
            .arg(fixture("single.sam"))
            .args(["-a", fixture(annotation).to_str().unwrap()])
            .args(["-F", format])
            .args(["-o", output.to_str().unwrap()])
            .args(["--threads", "1"])
            .output()
            .unwrap();
        assert!(
            result.status.success(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(output_counts(&output)[0].last().unwrap(), "3");
    }

    let compressed = directory.path().join("model.gtf.gz");
    let mut encoder = GzEncoder::new(File::create(&compressed).unwrap(), Compression::default());
    encoder
        .write_all(&fs::read(fixture("model.gtf")).unwrap())
        .unwrap();
    encoder.finish().unwrap();
    let gzip_output = directory.path().join("gzip.tsv");
    let result = Command::new(binary())
        .arg(fixture("single.sam"))
        .args(["-a", compressed.to_str().unwrap()])
        .args(["-o", gzip_output.to_str().unwrap()])
        .args(["--threads", "1"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(output_counts(&gzip_output)[0].last().unwrap(), "3");

    let alias_output = directory.path().join("alias.tsv");
    let result = Command::new(binary())
        .arg(fixture("aliases.sam"))
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-A")
        .arg(fixture("aliases.csv"))
        .arg("-o")
        .arg(&alias_output)
        .args(["--threads", "1"])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(output_counts(&alias_output)[0].last().unwrap(), "1");
}

#[test]
fn missing_second_output_parent_preserves_existing_counts() {
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("counts.tsv");
    fs::write(&output, "existing\n").unwrap();
    let summary = directory.path().join("missing/summary.tsv");
    let result = Command::new(binary())
        .arg(fixture("single.sam"))
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-o")
        .arg(&output)
        .arg("--summary")
        .arg(summary)
        .args(["--threads", "1"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read_to_string(output).unwrap(), "existing\n");
}

#[test]
fn truncated_bam_fails_without_replacing_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let bam = directory.path().join("truncated.bam");
    sam_to_bam(&fixture("single.sam"), &bam);
    let file = OpenOptions::new().write(true).open(&bam).unwrap();
    file.set_len(file.metadata().unwrap().len() / 2).unwrap();

    let counts = directory.path().join("counts.tsv");
    let summary = directory.path().join("summary.tsv");
    fs::write(&counts, "existing counts\n").unwrap();
    fs::write(&summary, "existing summary\n").unwrap();
    let result = Command::new(binary())
        .arg(&bam)
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-o")
        .arg(&counts)
        .arg("--summary")
        .arg(&summary)
        .args(["--threads", "2"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read_to_string(counts).unwrap(), "existing counts\n");
    assert_eq!(fs::read_to_string(summary).unwrap(), "existing summary\n");
}

#[test]
fn malformed_cigar_fails_without_replacing_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("invalid.sam");
    fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\nbad\t0\tchr1\t105\t60\t0M\t*\t0\t0\t*\t*\n",
    )
    .unwrap();
    let counts = directory.path().join("counts.tsv");
    let summary = directory.path().join("summary.tsv");
    fs::write(&counts, "existing counts\n").unwrap();
    fs::write(&summary, "existing summary\n").unwrap();
    let result = Command::new(binary())
        .arg(&sam)
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-o")
        .arg(&counts)
        .arg("--summary")
        .arg(&summary)
        .args(["--threads", "1"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read_to_string(counts).unwrap(), "existing counts\n");
    assert_eq!(fs::read_to_string(summary).unwrap(), "existing summary\n");
}

#[test]
fn paired_bam_requires_read_names() {
    let directory = tempfile::tempdir().unwrap();
    let sam = directory.path().join("missing-name.sam");
    fs::write(
        &sam,
        "@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n*\t99\tchr1\t105\t60\t20M\t=\t205\t120\tAAAAAAAAAAAAAAAAAAAA\t*\n*\t147\tchr1\t205\t60\t20M\t=\t105\t-120\tAAAAAAAAAAAAAAAAAAAA\t*\n",
    )
    .unwrap();
    let bam = directory.path().join("missing-name.bam");
    sam_to_bam(&sam, &bam);
    let result = Command::new(binary())
        .arg(&bam)
        .arg("-a")
        .arg(fixture("model.gtf"))
        .arg("-o")
        .arg(directory.path().join("counts.tsv"))
        .args(["--threads", "2", "--paired"])
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(
        String::from_utf8_lossy(&result.stderr).contains("no read name"),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
}

#[test]
fn invalid_bam_assignment_tags_fail() {
    let directory = tempfile::tempdir().unwrap();
    let cases = [
        (
            "negative-nh",
            "r\t0\tchr1\t105\t60\t20M\t*\t0\t0\tAAAAAAAAAAAAAAAAAAAA\t*\tNH:i:-1\n",
            &[][..],
        ),
        (
            "zero-hi",
            "p\t99\tchr1\t105\t60\t20M\t=\t205\t120\tAAAAAAAAAAAAAAAAAAAA\t*\tHI:i:0\np\t147\tchr1\t205\t60\t20M\t=\t105\t-120\tAAAAAAAAAAAAAAAAAAAA\t*\tHI:i:0\n",
            &["--paired"][..],
        ),
    ];
    for (name, records, options) in cases {
        let sam = directory.path().join(format!("{name}.sam"));
        fs::write(
            &sam,
            format!("@HD\tVN:1.6\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n{records}"),
        )
        .unwrap();
        let bam = directory.path().join(format!("{name}.bam"));
        sam_to_bam(&sam, &bam);
        let result = Command::new(binary())
            .arg(&bam)
            .arg("-a")
            .arg(fixture("model.gtf"))
            .arg("-o")
            .arg(directory.path().join(format!("{name}.tsv")))
            .args(["--threads", "2"])
            .args(options)
            .output()
            .unwrap();
        assert!(!result.status.success(), "{name}");
    }
}
