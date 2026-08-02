use std::env;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

const UNIQUE_GENES: usize = 450;
const OVERLAP_PAIRS: usize = 25;
const REFERENCE_LENGTH: usize = 8_000_000;

fn main() -> io::Result<()> {
    let mut args = env::args_os().skip(1);
    let annotation = args
        .next()
        .expect("usage: generate_fixture <annotation.gtf> <reads>");
    let reads: usize = args
        .next()
        .expect("usage: generate_fixture <annotation.gtf> <records> [single|paired]")
        .to_string_lossy()
        .parse()
        .expect("record count must be an integer");
    let mode = args.next().unwrap_or_else(|| "single".into());
    assert!(args.next().is_none(), "too many arguments");

    write_annotation(Path::new(&annotation))?;
    match mode.to_str() {
        Some("single") => write_single_alignments(reads),
        Some("paired") => write_paired_alignments(reads),
        _ => panic!("mode must be single or paired"),
    }
}

fn write_annotation(path: &Path) -> io::Result<()> {
    let mut output = BufWriter::new(File::create(path)?);
    for gene in 0..UNIQUE_GENES {
        let base = 1_000 + gene * 15_000;
        write_gene(&mut output, gene, base, gene % 2 == 0)?;
    }
    for pair in 0..OVERLAP_PAIRS {
        let base = 7_000_000 + pair * 20_000;
        write_gene(&mut output, UNIQUE_GENES + pair * 2, base, true)?;
        write_gene(&mut output, UNIQUE_GENES + pair * 2 + 1, base + 100, false)?;
    }
    output.flush()
}

fn write_gene(output: &mut impl Write, gene: usize, base: usize, forward: bool) -> io::Result<()> {
    let strand = if forward { '+' } else { '-' };
    for offset in [0, 500, 1_000] {
        writeln!(
            output,
            "chr1\tfixture\texon\t{}\t{}\t.\t{strand}\t.\tgene_id \"gene{gene:04}\"; transcript_id \"tx{gene:04}\";",
            base + offset + 1,
            base + offset + 200
        )?;
    }
    Ok(())
}

fn write_single_alignments(reads: usize) -> io::Result<()> {
    let mut output = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    writeln!(output, "@HD\tVN:1.6\tSO:unsorted")?;
    writeln!(output, "@SQ\tSN:chr1\tLN:{REFERENCE_LENGTH}")?;
    for read in 0..reads {
        let category = read % 20;
        let sequence = sequence(read as u64);
        let quality = "I".repeat(100);
        if category == 0 {
            writeln!(
                output,
                "r{read:08}\t4\t*\t0\t0\t*\t*\t0\t0\t{sequence}\t{quality}"
            )?;
            continue;
        }
        let (position, cigar) = match category {
            2 => {
                let pair = (read / 20) % OVERLAP_PAIRS;
                (7_000_000 + pair * 20_000 + 151, "100M")
            }
            3 => (7_900_000 + (read / 20) % 10_000, "100M"),
            4..=7 => {
                let gene = (read / 20) % UNIQUE_GENES;
                (1_000 + gene * 15_000 + 151, "50M300N50M")
            }
            _ => {
                let gene = (read / 20) % UNIQUE_GENES;
                (1_000 + gene * 15_000 + 51, "100M")
            }
        };
        let nh = if category == 1 { "\tNH:i:2" } else { "" };
        writeln!(
            output,
            "r{read:08}\t0\tchr1\t{position}\t60\t{cigar}\t*\t0\t0\t{sequence}\t{quality}{nh}"
        )?;
    }
    output.flush()
}

fn write_paired_alignments(fragments: usize) -> io::Result<()> {
    let mut output = BufWriter::with_capacity(1024 * 1024, io::stdout().lock());
    writeln!(output, "@HD\tVN:1.6\tSO:queryname")?;
    writeln!(output, "@SQ\tSN:chr1\tLN:{REFERENCE_LENGTH}")?;
    let quality = "I".repeat(100);
    for fragment in 0..fragments {
        let category = fragment % 20;
        let name = format!("p{fragment:08}");
        let first_sequence = sequence((fragment * 2) as u64);
        let last_sequence = sequence((fragment * 2 + 1) as u64);
        match category {
            0 => {
                writeln!(
                    output,
                    "{name}\t77\t*\t0\t0\t*\t*\t0\t0\t{first_sequence}\t{quality}"
                )?;
                writeln!(
                    output,
                    "{name}\t141\t*\t0\t0\t*\t*\t0\t0\t{last_sequence}\t{quality}"
                )?;
            }
            1 => {
                let gene = (fragment / 20) % UNIQUE_GENES;
                let first = 1_000 + gene * 15_000 + 51;
                writeln!(
                    output,
                    "{name}\t73\tchr1\t{first}\t60\t100M\t*\t0\t0\t{first_sequence}\t{quality}"
                )?;
                writeln!(
                    output,
                    "{name}\t133\t*\t0\t0\t*\tchr1\t{first}\t0\t{last_sequence}\t{quality}"
                )?;
            }
            _ => {
                let (first, last, first_cigar) = match category {
                    2 => {
                        let pair = (fragment / 20) % OVERLAP_PAIRS;
                        let first = 7_000_000 + pair * 20_000 + 151;
                        (first, first + 10, "100M")
                    }
                    3 => {
                        let first = 7_900_000 + (fragment / 20) % 10_000;
                        (first, first + 150, "100M")
                    }
                    4..=7 => {
                        let gene = (fragment / 20) % UNIQUE_GENES;
                        let base = 1_000 + gene * 15_000;
                        (base + 151, base + 1_051, "50M300N50M")
                    }
                    _ => {
                        let gene = (fragment / 20) % UNIQUE_GENES;
                        let base = 1_000 + gene * 15_000;
                        (base + 51, base + 551, "100M")
                    }
                };
                let template_length = last + 99 - first;
                let nh = if category == 8 { "\tNH:i:2" } else { "" };
                writeln!(
                    output,
                    "{name}\t99\tchr1\t{first}\t60\t{first_cigar}\t=\t{last}\t{template_length}\t{first_sequence}\t{quality}{nh}"
                )?;
                writeln!(
                    output,
                    "{name}\t147\tchr1\t{last}\t60\t100M\t=\t{first}\t-{template_length}\t{last_sequence}\t{quality}{nh}"
                )?;
            }
        }
    }
    output.flush()
}

fn sequence(mut state: u64) -> String {
    let mut sequence = String::with_capacity(100);
    for _ in 0..100 {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        sequence.push(match state & 3 {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            _ => 'T',
        });
    }
    sequence
}
