use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use rsomics_common::{Context, Result, RsomicsError};
use serde::Serialize;

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum NormalizeMode {
    Tpm,
    Fpkm,
    FpkmUq,
}

#[derive(Clone, Debug)]
pub struct NormalizeOptions {
    pub mode: NormalizeMode,
    pub lengths: Option<PathBuf>,
}

#[derive(Clone, Debug, Serialize)]
pub struct NormalizeReport {
    pub genes: usize,
    pub samples: usize,
    pub mode: &'static str,
}

pub fn normalize(
    counts_path: &Path,
    output_path: &Path,
    options: &NormalizeOptions,
) -> Result<NormalizeReport> {
    rsomics_common::reject_output_alias(
        output_path,
        std::iter::once(counts_path).chain(options.lengths.as_deref()),
    )?;
    let matrix = Matrix::read(counts_path, options.lengths.as_deref())?;
    let values = matrix.normalize(options.mode)?;
    rsomics_common::write_atomic(output_path, |output| matrix.write(output, &values))?;
    Ok(NormalizeReport {
        genes: matrix.genes.len(),
        samples: matrix.samples.len(),
        mode: match options.mode {
            NormalizeMode::Tpm => "tpm",
            NormalizeMode::Fpkm => "fpkm",
            NormalizeMode::FpkmUq => "fpkm-uq",
        },
    })
}

struct Matrix {
    genes: Vec<String>,
    samples: Vec<String>,
    lengths: Vec<f64>,
    counts: Vec<Vec<f64>>,
}

impl Matrix {
    fn read(path: &Path, lengths_path: Option<&Path>) -> Result<Self> {
        let file = File::open(path)
            .rs_with_context(|| format!("opening count matrix {}", path.display()))?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines().enumerate();
        let (header_line, header) = loop {
            let Some((index, result)) = lines.next() else {
                return Err(RsomicsError::InvalidInput(format!(
                    "{}: empty count matrix",
                    path.display()
                )));
            };
            let line = result.rs_with_context(|| format!("reading {}", path.display()))?;
            if !line.starts_with('#') && !line.trim().is_empty() {
                break (index + 1, line);
            }
        };
        let header: Vec<&str> = header.split('\t').collect();
        let embedded = header.len() >= 7
            && header[..6] == ["Geneid", "Chr", "Start", "End", "Strand", "Length"];
        let sample_start = if embedded { 6 } else { 1 };
        if header.len() <= sample_start {
            return Err(matrix_error(
                path,
                header_line,
                "matrix has no sample columns",
            ));
        }
        if !embedded && header[0] != "Geneid" && header[0] != "gene" {
            return Err(matrix_error(
                path,
                header_line,
                "first column must be Geneid or gene",
            ));
        }
        let samples: Vec<String> = header[sample_start..]
            .iter()
            .map(|sample| (*sample).to_owned())
            .collect();
        if samples.iter().any(String::is_empty) {
            return Err(matrix_error(
                path,
                header_line,
                "sample names must be nonempty",
            ));
        }
        let unique_samples: HashSet<_> = samples.iter().collect();
        if unique_samples.len() != samples.len() {
            return Err(matrix_error(path, header_line, "duplicate sample name"));
        }

        let supplied_lengths = lengths_path.map(read_lengths).transpose()?;
        if !embedded && supplied_lengths.is_none() {
            return Err(RsomicsError::ConfigError(
                "--lengths is required when the matrix has no Length column".to_owned(),
            ));
        }

        let mut genes = Vec::new();
        let mut lengths = Vec::new();
        let mut counts = Vec::new();
        let mut seen = HashSet::new();
        for (index, result) in lines {
            let line = result.rs_with_context(|| format!("reading {}", path.display()))?;
            if line.trim().is_empty() {
                return Err(matrix_error(path, index + 1, "blank matrix row"));
            }
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() != header.len() {
                return Err(matrix_error(
                    path,
                    index + 1,
                    format!("expected {} columns, found {}", header.len(), fields.len()),
                ));
            }
            let gene = fields[0];
            if gene.is_empty() || !seen.insert(gene.to_owned()) {
                return Err(matrix_error(
                    path,
                    index + 1,
                    format!("empty or duplicate gene identifier {gene:?}"),
                ));
            }
            let length = match supplied_lengths.as_ref() {
                Some(values) => values.get(gene).copied().ok_or_else(|| {
                    matrix_error(path, index + 1, format!("missing length for gene {gene:?}"))
                })?,
                None => parse_positive(fields[5], path, index + 1, "length")?,
            };
            let row = fields[sample_start..]
                .iter()
                .enumerate()
                .map(|(sample, value)| parse_count(value, path, index + 1, &samples[sample], gene))
                .collect::<Result<Vec<_>>>()?;
            genes.push(gene.to_owned());
            lengths.push(length);
            counts.push(row);
        }
        if genes.is_empty() {
            return Err(RsomicsError::InvalidInput(format!(
                "{}: count matrix has no data rows",
                path.display()
            )));
        }
        Ok(Self {
            genes,
            samples,
            lengths,
            counts,
        })
    }

    fn normalize(&self, mode: NormalizeMode) -> Result<Vec<Vec<f64>>> {
        let samples = self.samples.len();
        let mut denominators = vec![0.0; samples];
        match mode {
            NormalizeMode::Tpm => {
                for (row, length) in self.counts.iter().zip(&self.lengths) {
                    for (sample, count) in row.iter().enumerate() {
                        denominators[sample] += count / (length / 1_000.0);
                    }
                }
            }
            NormalizeMode::Fpkm => {
                for row in &self.counts {
                    for (sample, count) in row.iter().enumerate() {
                        denominators[sample] += count;
                    }
                }
            }
            NormalizeMode::FpkmUq => {
                for (sample, denominator) in denominators.iter_mut().enumerate() {
                    let mut expressed: Vec<f64> = self
                        .counts
                        .iter()
                        .map(|row| row[sample])
                        .filter(|count| *count > 0.0)
                        .collect();
                    *denominator = upper_quartile(&mut expressed);
                }
            }
        }
        for (sample, denominator) in denominators.iter().enumerate() {
            if *denominator == 0.0 {
                return Err(RsomicsError::InvalidInput(format!(
                    "sample {:?} has a zero {:?} denominator",
                    self.samples[sample], mode
                )));
            }
        }

        Ok(self
            .counts
            .iter()
            .enumerate()
            .map(|(gene, row)| {
                row.iter()
                    .enumerate()
                    .map(|(sample, count)| match mode {
                        NormalizeMode::Tpm => {
                            count / (self.lengths[gene] / 1_000.0) / denominators[sample]
                                * 1_000_000.0
                        }
                        NormalizeMode::Fpkm | NormalizeMode::FpkmUq => {
                            count * 1_000_000_000.0 / (self.lengths[gene] * denominators[sample])
                        }
                    })
                    .collect()
            })
            .collect())
    }

    fn write(&self, output: &mut File, values: &[Vec<f64>]) -> Result<()> {
        let mut output = BufWriter::new(output);
        write!(output, "Geneid").map_err(RsomicsError::Io)?;
        for sample in &self.samples {
            write!(output, "\t{sample}").map_err(RsomicsError::Io)?;
        }
        writeln!(output).map_err(RsomicsError::Io)?;
        for (gene, row) in self.genes.iter().zip(values) {
            write!(output, "{gene}").map_err(RsomicsError::Io)?;
            for value in row {
                write!(output, "\t{value:.6}").map_err(RsomicsError::Io)?;
            }
            writeln!(output).map_err(RsomicsError::Io)?;
        }
        output.flush().map_err(RsomicsError::Io)
    }
}

fn read_lengths(path: &Path) -> Result<HashMap<String, f64>> {
    let file =
        File::open(path).rs_with_context(|| format!("opening lengths {}", path.display()))?;
    let reader = BufReader::new(file);
    let mut lengths = HashMap::new();
    for (index, result) in reader.lines().enumerate() {
        let line = result.rs_with_context(|| format!("reading {}", path.display()))?;
        if line.starts_with('#') || line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if index == 0
            && fields.len() == 2
            && matches!(fields[0], "Geneid" | "gene")
            && matches!(fields[1], "Length" | "length")
        {
            continue;
        }
        if fields.len() != 2 || fields[0].is_empty() {
            return Err(matrix_error(
                path,
                index + 1,
                "expected gene and positive length columns",
            ));
        }
        let value = parse_positive(fields[1], path, index + 1, "length")?;
        if lengths.insert(fields[0].to_owned(), value).is_some() {
            return Err(matrix_error(
                path,
                index + 1,
                format!("duplicate length for gene {:?}", fields[0]),
            ));
        }
    }
    if lengths.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "{}: no gene lengths",
            path.display()
        )));
    }
    Ok(lengths)
}

fn parse_positive(value: &str, path: &Path, line: usize, label: &str) -> Result<f64> {
    let value: f64 = value
        .parse()
        .map_err(|_| matrix_error(path, line, format!("invalid {label} value {value:?}")))?;
    if !value.is_finite() || value <= 0.0 {
        return Err(matrix_error(
            path,
            line,
            format!("{label} must be finite and positive"),
        ));
    }
    Ok(value)
}

fn parse_count(value: &str, path: &Path, line: usize, sample: &str, gene: &str) -> Result<f64> {
    let count: f64 = value.parse().map_err(|_| {
        matrix_error(
            path,
            line,
            format!("invalid count {value:?} for gene {gene:?}, sample {sample:?}"),
        )
    })?;
    if !count.is_finite() || count < 0.0 {
        return Err(matrix_error(
            path,
            line,
            format!("count for gene {gene:?}, sample {sample:?} must be finite and nonnegative"),
        ));
    }
    Ok(count)
}

fn upper_quartile(values: &mut [f64]) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    values.sort_by(f64::total_cmp);
    let position = (values.len() - 1) as f64 * 0.75;
    let lower = position.floor() as usize;
    let fraction = position.fract();
    values[lower]
        + values
            .get(lower + 1)
            .map_or(0.0, |upper| fraction * (*upper - values[lower]))
}

fn matrix_error(path: &Path, line: usize, error: impl std::fmt::Display) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{} line {line}: {error}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quartile_matches_linear_interpolation() {
        let mut values = [1.0, 2.0, 3.0, 4.0];
        assert_eq!(upper_quartile(&mut values), 3.25);
    }

    #[test]
    fn zero_library_fails_loud() {
        let matrix = Matrix {
            genes: vec!["g".to_owned()],
            samples: vec!["s".to_owned()],
            lengths: vec![100.0],
            counts: vec![vec![0.0]],
        };
        assert!(matrix.normalize(NormalizeMode::Tpm).is_err());
    }
}
