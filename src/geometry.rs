use rsomics_common::{Result, RsomicsError};
use smallvec::SmallVec;

use crate::model::{Block, Blocks};

pub(crate) struct Geometry {
    pub(crate) blocks: Blocks,
    pub(crate) split: bool,
    pub(crate) extent: Extent,
}

#[derive(Clone, Debug, Default)]
pub(crate) struct Extent {
    intervals: SmallVec<[Span; 2]>,
}

#[derive(Clone, Copy, Debug)]
struct Span {
    start: u64,
    end: u64,
    insertion: Option<(u64, u64)>,
}

impl Span {
    fn length(self) -> u64 {
        self.end - self.start + self.insertion.map_or(0, |(_, length)| length)
    }

    fn prefix(&mut self, end: u64) -> u64 {
        let mut length = end - self.start;
        if self.insertion.is_some_and(|(anchor, _)| anchor < end) {
            length += self.insertion.take().map_or(0, |(_, length)| length);
        }
        self.start = end;
        length
    }
}

impl Extent {
    pub(crate) fn length(&self) -> u64 {
        self.intervals.iter().copied().map(Span::length).sum()
    }

    pub(crate) fn union_length(&self, other: &Self) -> u64 {
        let mut left = self.intervals.iter().copied();
        let mut right = other.intervals.iter().copied();
        let mut a = left.next();
        let mut b = right.next();
        let mut length = 0;
        while let (Some(mut x), Some(mut y)) = (a, b) {
            if (x.start <= y.start && x.end > y.start) || (y.start <= x.start && y.end > x.start) {
                let start = x.start.max(y.start);
                let end = x.end.min(y.end);
                length += x.prefix(start) + y.prefix(start) + end - start;
                if let (Some(ix), Some(iy)) = (x.insertion, y.insertion) {
                    if ix == iy {
                        length += ix.1;
                    }
                    if ix.0 <= iy.0 {
                        x.insertion = None;
                    }
                    if iy.0 <= ix.0 {
                        y.insertion = None;
                    }
                }
                a = if x.end > end {
                    x.start = end;
                    Some(x)
                } else {
                    left.next()
                };
                b = if y.end > end {
                    y.start = end;
                    Some(y)
                } else {
                    right.next()
                };
            } else if x.start > y.start {
                length += y.length();
                b = right.next();
            } else {
                length += x.length();
                a = left.next();
            }
        }
        length
            + a.map_or(0, Span::length)
            + b.map_or(0, Span::length)
            + left.map(Span::length).sum::<u64>()
            + right.map(Span::length).sum::<u64>()
    }
}

fn advance(position: u64, length: u64) -> Result<u64> {
    position
        .checked_add(length)
        .ok_or_else(|| RsomicsError::InvalidInput("CIGAR coordinate overflow".to_owned()))
}

pub(crate) fn parse(
    start: u64,
    chrom: Option<u32>,
    operations: impl IntoIterator<Item = Result<(u8, u64)>>,
    track_extent: bool,
) -> Result<Geometry> {
    let mut blocks = Blocks::new();
    let mut extent = Extent::default();
    let mut cursor = start;
    let mut block_start = start;
    let mut extent_start = if track_extent { advance(start, 1)? } else { 0 };
    let mut read_consumed = false;
    let mut trailing_clip = 0;
    let mut split = false;
    for operation in operations {
        let (kind, length) = operation?;
        if length == 0 {
            return Err(RsomicsError::InvalidInput(
                "CIGAR operation has zero length".to_owned(),
            ));
        }
        match kind {
            0 | 7 | 8 => {
                cursor = advance(cursor, length)?;
                read_consumed = true;
            }
            1..=3 => {
                if cursor > block_start {
                    blocks.push(Block::new(chrom, block_start, cursor)?);
                }
                if track_extent {
                    let end = advance(cursor, 1)?;
                    extent.intervals.push(Span {
                        start: extent_start,
                        end,
                        insertion: (kind == 1).then_some((end, length)),
                    });
                }
                if kind == 1 {
                    read_consumed = true;
                } else {
                    cursor = advance(cursor, length)?;
                    split |= kind == 3;
                }
                block_start = cursor;
                if track_extent {
                    extent_start = advance(cursor, 1)?;
                }
            }
            4 => {
                if read_consumed {
                    trailing_clip = length;
                } else {
                    read_consumed = true;
                    if track_extent {
                        extent_start = extent_start.saturating_sub(length);
                    }
                }
            }
            5 | 6 => {}
            _ => {
                return Err(RsomicsError::InvalidInput(format!(
                    "invalid BAM CIGAR operation code {kind}"
                )));
            }
        }
    }
    if cursor > block_start {
        blocks.push(Block::new(chrom, block_start, cursor)?);
        if track_extent {
            extent.intervals.push(Span {
                start: extent_start,
                end: advance(advance(cursor, 1)?, trailing_clip)?,
                insertion: None,
            });
        }
    }
    if blocks.is_empty() {
        return Err(RsomicsError::InvalidInput(
            "mapped record CIGAR consumes no reference bases".to_owned(),
        ));
    }
    Ok(Geometry {
        blocks,
        split,
        extent,
    })
}
