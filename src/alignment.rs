use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufReader, Read};
use std::num::NonZeroUsize;
use std::path::Path;

use bstr::ByteSlice;
use noodles::{bgzf, sam};
use noodles_util::alignment;
use rsomics_bamio::raw::{RecordReader as RawRecordReader, RecordRef};
use rsomics_common::{Context, Result, RsomicsError};

use crate::model::{Block, Blocks};

enum Inner {
    Bam(rsomics_bamio::ParallelBamReader),
    General(alignment::io::Reader<Box<dyn Read>>),
}

pub(crate) struct Reader {
    inner: Inner,
    input: String,
}

#[derive(Clone, Debug)]
pub(crate) struct Segment {
    pub(crate) name: Option<Vec<u8>>,
    pub(crate) name_present: bool,
    pub(crate) group_start: bool,
    pub(crate) blocks: Blocks,
    pub(crate) reverse: bool,
    pub(crate) segmented: bool,
    pub(crate) first: bool,
    pub(crate) last: bool,
    pub(crate) unmapped: bool,
    pub(crate) secondary: bool,
    pub(crate) duplicate: bool,
    pub(crate) mapq: u8,
    pub(crate) nh: u32,
    pub(crate) hit_index: Option<u32>,
    pub(crate) split: bool,
    pub(crate) template_length: i32,
    pub(crate) reference: Option<usize>,
}

#[derive(Clone, Copy)]
pub(crate) enum NameMode {
    None,
    Owned,
    Grouped,
}

struct NameTracker {
    mode: NameMode,
    previous: Vec<u8>,
}

impl NameTracker {
    fn new(mode: NameMode) -> Self {
        Self {
            mode,
            previous: Vec::new(),
        }
    }

    fn read(&mut self, raw: Option<&[u8]>) -> (Option<Vec<u8>>, bool, bool) {
        let Some(raw) = raw.filter(|name| !name.is_empty() && *name != b"*") else {
            return (None, false, false);
        };
        match self.mode {
            NameMode::None => (None, true, false),
            NameMode::Owned => (Some(raw.to_vec()), true, false),
            NameMode::Grouped => {
                let group_start = self.previous != raw;
                if group_start {
                    self.previous.clear();
                    self.previous.extend_from_slice(raw);
                }
                (None, true, group_start)
            }
        }
    }
}

impl Reader {
    pub(crate) fn open(input: &Path, workers: Option<NonZeroUsize>) -> Result<Self> {
        let display = input.display().to_string();
        if input == Path::new("-") {
            if workers.is_some_and(|workers| workers.get() != 1) {
                return Err(RsomicsError::ConfigError(
                    "standard-input alignment decoding requires --threads 1".to_owned(),
                ));
            }
            let inner = alignment::io::reader::Builder::default()
                .build_from_reader(Box::new(io::stdin()) as Box<dyn Read>)
                .rs_context("opening alignment input -")?;
            return Ok(Self {
                inner: Inner::General(inner),
                input: display,
            });
        }

        match detect(input)? {
            Format::Cram => Err(RsomicsError::InvalidInput(format!(
                "CRAM is not supported in this release: {}",
                input.display()
            ))),
            Format::Bam => Ok(Self {
                inner: Inner::Bam(rsomics_bamio::open_with_workers(
                    input,
                    workers.unwrap_or_else(|| {
                        std::thread::available_parallelism().unwrap_or(NonZeroUsize::MIN)
                    }),
                )?),
                input: display,
            }),
            Format::Sam => {
                if workers.is_some_and(|workers| workers.get() != 1) {
                    return Err(RsomicsError::ConfigError(format!(
                        "SAM input requires --threads 1: {}",
                        input.display()
                    )));
                }
                let file = File::open(input)
                    .rs_with_context(|| format!("opening alignment input {}", input.display()))?;
                let inner = alignment::io::reader::Builder::default()
                    .build_from_reader(Box::new(file) as Box<dyn Read>)
                    .rs_with_context(|| format!("opening alignment input {}", input.display()))?;
                Ok(Self {
                    inner: Inner::General(inner),
                    input: display,
                })
            }
        }
    }

    pub(crate) fn read_header(&mut self) -> Result<sam::Header> {
        match &mut self.inner {
            Inner::Bam(reader) => reader.read_header(),
            Inner::General(reader) => reader.read_header(),
        }
        .map_err(|error| {
            RsomicsError::InvalidInput(format!(
                "reading alignment header from {}: {error}",
                self.input
            ))
        })
    }

    pub(crate) fn visit(
        &mut self,
        header: &sam::Header,
        aliases: &HashMap<String, String>,
        chrom_id: &dyn Fn(&str) -> Option<u32>,
        name_mode: NameMode,
        mut visit: impl FnMut(Segment) -> Result<()>,
    ) -> Result<()> {
        let mut names = NameTracker::new(name_mode);
        match &mut self.inner {
            Inner::Bam(reader) => {
                let chrom_by_ref = reference_ids(header, aliases, chrom_id)?;
                let mut records = RawRecordReader::new(reader.get_mut());
                while let Some(record) = records
                    .next()
                    .rs_with_context(|| format!("reading alignment record from {}", self.input))?
                {
                    visit(parse_raw_record(&record, &chrom_by_ref, &mut names)?)?;
                }
            }
            Inner::General(reader) => {
                for result in reader.records(header) {
                    let record = result.map_err(|error| record_error(&self.input, error))?;
                    visit(parse_record(
                        header,
                        record.as_ref(),
                        aliases,
                        chrom_id,
                        &mut names,
                    )?)?;
                }
            }
        }
        Ok(())
    }
}

fn parse_record(
    header: &sam::Header,
    record: &dyn sam::alignment::Record,
    aliases: &HashMap<String, String>,
    chrom_id: &dyn Fn(&str) -> Option<u32>,
    names: &mut NameTracker,
) -> Result<Segment> {
    use sam::alignment::record::data::field::Tag;

    let flags = record
        .flags()
        .map_err(|error| RsomicsError::InvalidInput(format!("invalid alignment flags: {error}")))?;
    let (name, name_present, group_start) = names.read(record.name().map(|name| name.as_ref()));
    let (reference, chrom) = reference_info(
        record.reference_sequence_id(header),
        header,
        aliases,
        chrom_id,
    )?;
    let (blocks, split) = if flags.is_unmapped() {
        (Blocks::new(), false)
    } else {
        reference.ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "mapped record {} has no reference sequence",
                display_name(name.as_deref())
            ))
        })?;
        cigar_blocks(record, chrom)?
    };
    let data = record.data();
    let nh = integer_tag(data.get(&Tag::ALIGNMENT_HIT_COUNT), "NH")?.unwrap_or(1);
    if nh == 0 {
        return Err(RsomicsError::InvalidInput(format!(
            "record {} has NH=0",
            display_name(name.as_deref())
        )));
    }
    let hit_index = integer_tag(data.get(&Tag::HIT_INDEX), "HI")?;
    if hit_index == Some(0) {
        return Err(RsomicsError::InvalidInput(format!(
            "record {} has HI=0",
            display_name(name.as_deref())
        )));
    }
    let mapq = record
        .mapping_quality()
        .transpose()
        .map_err(|error| RsomicsError::InvalidInput(format!("invalid mapping quality: {error}")))?
        .map_or(0, |quality| quality.get());

    Ok(Segment {
        name,
        name_present,
        group_start,
        blocks,
        reverse: flags.is_reverse_complemented(),
        segmented: flags.is_segmented(),
        first: flags.is_first_segment(),
        last: flags.is_last_segment(),
        unmapped: flags.is_unmapped(),
        secondary: flags.is_secondary(),
        duplicate: flags.is_duplicate(),
        mapq,
        nh,
        hit_index,
        split,
        template_length: record.template_length().map_err(|error| {
            RsomicsError::InvalidInput(format!("invalid template length: {error}"))
        })?,
        reference,
    })
}

fn cigar_blocks(record: &dyn sam::alignment::Record, chrom: Option<u32>) -> Result<(Blocks, bool)> {
    use sam::alignment::record::cigar::op::Kind;

    let start = record
        .alignment_start()
        .transpose()
        .map_err(|error| RsomicsError::InvalidInput(format!("invalid alignment start: {error}")))?
        .ok_or_else(|| {
            RsomicsError::InvalidInput("mapped record has no alignment start".to_owned())
        })?;
    let mut cursor = u64::try_from(usize::from(start) - 1)
        .map_err(|_| RsomicsError::InvalidInput("alignment coordinate overflow".to_owned()))?;
    let mut block_start = cursor;
    let mut blocks = Blocks::new();
    let mut split = false;
    for result in record.cigar().iter() {
        let op = result
            .map_err(|error| RsomicsError::InvalidInput(format!("invalid CIGAR: {error}")))?;
        if op.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "CIGAR operation has zero length".to_owned(),
            ));
        }
        match op.kind() {
            Kind::Skip | Kind::Deletion => {
                split |= op.kind() == Kind::Skip;
                if block_start < cursor {
                    blocks.push(Block::new(chrom, block_start, cursor)?);
                }
                cursor = cursor.checked_add(op.len() as u64).ok_or_else(|| {
                    RsomicsError::InvalidInput("CIGAR coordinate overflow".to_owned())
                })?;
                block_start = cursor;
            }
            kind if kind.consumes_reference() => {
                cursor = cursor.checked_add(op.len() as u64).ok_or_else(|| {
                    RsomicsError::InvalidInput("CIGAR coordinate overflow".to_owned())
                })?;
            }
            _ => {}
        }
    }
    if block_start < cursor {
        blocks.push(Block::new(chrom, block_start, cursor)?);
    }
    if blocks.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "mapped record CIGAR consumes no reference bases".to_owned(),
        ));
    }
    Ok((blocks, split))
}

fn parse_raw_record(
    record: &RecordRef<'_>,
    chrom_by_ref: &[Option<u32>],
    names: &mut NameTracker,
) -> Result<Segment> {
    let flags = record.flags();
    let (name, name_present, group_start) = names.read(Some(record.name()));
    let reference = raw_reference_id(record.reference_sequence_id(), chrom_by_ref.len())?;
    let chrom = reference.and_then(|id| chrom_by_ref[id]);
    let unmapped = flags & 0x4 != 0;
    let (blocks, split) = if unmapped {
        (Blocks::new(), false)
    } else {
        reference.ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "mapped record {} has no reference sequence",
                display_name(name.as_deref())
            ))
        })?;
        raw_cigar_blocks(record, chrom)?
    };
    let nh = raw_integer_tag(record, *b"NH")?.unwrap_or(1);
    if nh == 0 {
        return Err(RsomicsError::InvalidInput(format!(
            "record {} has NH=0",
            display_name(name.as_deref())
        )));
    }
    let hit_index = raw_integer_tag(record, *b"HI")?;
    if hit_index == Some(0) {
        return Err(RsomicsError::InvalidInput(format!(
            "record {} has HI=0",
            display_name(name.as_deref())
        )));
    }

    Ok(Segment {
        name,
        name_present,
        group_start,
        blocks,
        reverse: flags & 0x10 != 0,
        segmented: flags & 0x1 != 0,
        first: flags & 0x40 != 0,
        last: flags & 0x80 != 0,
        unmapped,
        secondary: flags & 0x100 != 0,
        duplicate: flags & 0x400 != 0,
        mapq: record.mapping_quality(),
        nh,
        hit_index,
        split,
        template_length: record.template_length(),
        reference,
    })
}

fn raw_reference_id(value: i32, count: usize) -> Result<Option<usize>> {
    if value == -1 {
        return Ok(None);
    }
    let id = usize::try_from(value).map_err(|_| {
        RsomicsError::InvalidInput(format!("invalid reference sequence identifier {value}"))
    })?;
    if id >= count {
        return Err(RsomicsError::InvalidInput(format!(
            "reference sequence identifier {id} is absent from header"
        )));
    }
    Ok(Some(id))
}

fn raw_cigar_blocks(record: &RecordRef<'_>, chrom: Option<u32>) -> Result<(Blocks, bool)> {
    let start = record.alignment_start();
    let mut cursor = u64::try_from(start).map_err(|_| {
        RsomicsError::InvalidInput("mapped record has no alignment start".to_owned())
    })?;
    let mut block_start = cursor;
    let mut blocks = Blocks::new();
    let mut split = false;
    for (kind, len) in record.cigar_ops() {
        if len == 0 {
            return Err(RsomicsError::InvalidInput(
                "CIGAR operation has zero length".to_owned(),
            ));
        }
        match kind {
            2 | 3 => {
                split |= kind == 3;
                if block_start < cursor {
                    blocks.push(Block::new(chrom, block_start, cursor)?);
                }
                cursor = cursor.checked_add(u64::from(len)).ok_or_else(|| {
                    RsomicsError::InvalidInput("CIGAR coordinate overflow".to_owned())
                })?;
                block_start = cursor;
            }
            0 | 7 | 8 => {
                cursor = cursor.checked_add(u64::from(len)).ok_or_else(|| {
                    RsomicsError::InvalidInput("CIGAR coordinate overflow".to_owned())
                })?;
            }
            1 | 4 | 5 | 6 => {}
            _ => {
                return Err(RsomicsError::InvalidInput(format!(
                    "invalid BAM CIGAR operation code {kind}"
                )));
            }
        }
    }
    if block_start < cursor {
        blocks.push(Block::new(chrom, block_start, cursor)?);
    }
    if blocks.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "mapped record CIGAR consumes no reference bases".to_owned(),
        ));
    }
    Ok((blocks, split))
}

fn raw_integer_tag(record: &RecordRef<'_>, tag: [u8; 2]) -> Result<Option<u32>> {
    let Some(kind) = record.aux_type(tag) else {
        return Ok(None);
    };
    let value = record.aux_value(tag).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "invalid {}{} alignment tag",
            char::from(tag[0]),
            char::from(tag[1])
        ))
    })?;
    let signed = match (kind, value) {
        (b'c', [value]) => i64::from(i8::from_le_bytes([*value])),
        (b'C', [value]) => i64::from(*value),
        (b's', [a, b]) => i64::from(i16::from_le_bytes([*a, *b])),
        (b'S', [a, b]) => i64::from(u16::from_le_bytes([*a, *b])),
        (b'i', [a, b, c, d]) => i64::from(i32::from_le_bytes([*a, *b, *c, *d])),
        (b'I', [a, b, c, d]) => i64::from(u32::from_le_bytes([*a, *b, *c, *d])),
        _ => {
            return Err(RsomicsError::InvalidInput(format!(
                "{}{} alignment tag is not an integer",
                char::from(tag[0]),
                char::from(tag[1])
            )));
        }
    };
    u32::try_from(signed).map(Some).map_err(|_| {
        RsomicsError::InvalidInput(format!(
            "{}{} alignment tag is out of range",
            char::from(tag[0]),
            char::from(tag[1])
        ))
    })
}

fn display_name(name: Option<&[u8]>) -> Cow<'_, str> {
    name.map(String::from_utf8_lossy)
        .unwrap_or(Cow::Borrowed("<unnamed>"))
}

fn integer_tag(
    value: Option<io::Result<sam::alignment::record::data::field::Value<'_>>>,
    tag: &str,
) -> Result<Option<u32>> {
    let Some(value) = value else {
        return Ok(None);
    };
    let value = value.map_err(|error| {
        RsomicsError::InvalidInput(format!("invalid {tag} alignment tag: {error}"))
    })?;
    let value = value.as_int().ok_or_else(|| {
        RsomicsError::InvalidInput(format!("{tag} alignment tag is not an integer"))
    })?;
    u32::try_from(value)
        .map(Some)
        .map_err(|_| RsomicsError::InvalidInput(format!("{tag} alignment tag is out of range")))
}

fn reference_info(
    id: Option<io::Result<usize>>,
    header: &sam::Header,
    aliases: &HashMap<String, String>,
    chrom_id: &dyn Fn(&str) -> Option<u32>,
) -> Result<(Option<usize>, Option<u32>)> {
    let Some(id) = id.transpose().map_err(|error| {
        RsomicsError::InvalidInput(format!("invalid reference sequence identifier: {error}"))
    })?
    else {
        return Ok((None, None));
    };
    let (name, _) = header.reference_sequences().get_index(id).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "reference sequence identifier {id} is absent from header"
        ))
    })?;
    let name = name.to_str().map_err(|error| {
        RsomicsError::InvalidInput(format!("reference sequence name is not UTF-8: {error}"))
    })?;
    let canonical = aliases.get(name).map_or(name, String::as_str);
    Ok((Some(id), chrom_id(canonical)))
}

fn reference_ids(
    header: &sam::Header,
    aliases: &HashMap<String, String>,
    chrom_id: &dyn Fn(&str) -> Option<u32>,
) -> Result<Vec<Option<u32>>> {
    header
        .reference_sequences()
        .keys()
        .map(|name| {
            let name = name.to_str().map_err(|error| {
                RsomicsError::InvalidInput(format!("reference sequence name is not UTF-8: {error}"))
            })?;
            let canonical = aliases.get(name).map_or(name, String::as_str);
            Ok(chrom_id(canonical))
        })
        .collect()
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
enum Format {
    Sam,
    Bam,
    Cram,
}

fn detect(path: &Path) -> Result<Format> {
    let mut source = BufReader::new(
        File::open(path)
            .rs_with_context(|| format!("opening alignment input {}", path.display()))?,
    );
    let mut magic = [0; 4];
    source.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting alignment format for {}: {error}",
            path.display()
        ))
    })?;
    if magic == *b"CRAM" {
        return Ok(Format::Cram);
    }
    if magic == *b"BAM\x01" {
        return Ok(Format::Bam);
    }
    if magic[..2] != [0x1f, 0x8b] {
        return Ok(Format::Sam);
    }
    let mut reader = bgzf::io::Reader::new(
        File::open(path)
            .rs_with_context(|| format!("opening alignment input {}", path.display()))?,
    );
    reader.read_exact(&mut magic).map_err(|error| {
        RsomicsError::InvalidInput(format!(
            "detecting BGZF alignment format for {}: {error}",
            path.display()
        ))
    })?;
    Ok(if magic == *b"BAM\x01" {
        Format::Bam
    } else {
        Format::Sam
    })
}

fn record_error(input: &str, error: io::Error) -> RsomicsError {
    RsomicsError::InvalidInput(format!("reading alignment record from {input}: {error}"))
}

#[cfg(test)]
mod tests {
    use noodles::core::Position;
    use noodles::sam::alignment::{
        RecordBuf,
        record::cigar::{Op, op::Kind},
    };

    use super::*;

    #[test]
    fn cigar_skip_creates_reference_blocks() {
        let record = RecordBuf::builder()
            .set_alignment_start(Position::try_from(101).unwrap())
            .set_cigar(
                [
                    Op::new(Kind::Match, 20),
                    Op::new(Kind::Skip, 50),
                    Op::new(Kind::SequenceMatch, 10),
                ]
                .into_iter()
                .collect(),
            )
            .build();
        let (blocks, split) = cigar_blocks(&record, Some(0)).unwrap();
        assert_eq!((blocks[0].start, blocks[0].end), (100, 120));
        assert_eq!((blocks[1].start, blocks[1].end), (170, 180));
        assert!(split);
    }
}
