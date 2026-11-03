use rsomics_common::{Result, RsomicsError};

use crate::model::{Block, Blocks};

/// The biological read end retained as a single reference base.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum ReadPosition {
    FivePrime,
    ThreePrime,
}

/// Strand-aware reference extensions, applied before optional end reduction.
#[derive(Clone, Copy, Debug, Default)]
pub struct ProjectionOptions {
    pub extension5: u32,
    pub extension3: u32,
    pub read_position: Option<ReadPosition>,
}

#[derive(Clone, Copy, Debug, Default)]
pub(crate) struct Ends {
    pub(crate) first_start: u64,
    pub(crate) first_end: u64,
    pub(crate) last_start: u64,
    pub(crate) last_end: u64,
}

impl ProjectionOptions {
    pub(crate) fn validate(&self) -> Result<()> {
        if self.extension5 > i32::MAX as u32 || self.extension3 > i32::MAX as u32 {
            return Err(RsomicsError::ConfigError(
                "read extensions must not exceed 2147483647 bases".to_owned(),
            ));
        }
        Ok(())
    }

    pub(crate) fn is_active(&self) -> bool {
        self.extension5 != 0 || self.extension3 != 0 || self.read_position.is_some()
    }

    pub(crate) fn apply(
        &self,
        blocks: &mut Blocks,
        ends: Ends,
        reference_length: u64,
        reverse: bool,
    ) -> Result<()> {
        if !self.is_active() {
            return Ok(());
        }
        let invalid_extent = || {
            RsomicsError::InvalidInput(
                "mapped alignment extent is outside its reference sequence".to_owned(),
            )
        };
        let chrom = blocks.first().ok_or_else(invalid_extent)?.chrom;
        if ends.last_end > reference_length {
            return Err(invalid_extent());
        }
        let left_extension = u64::from(if reverse {
            self.extension3
        } else {
            self.extension5
        });
        let right_extension = u64::from(if reverse {
            self.extension5
        } else {
            self.extension3
        });
        let start = ends.first_start.saturating_sub(left_extension);
        let end = ends.last_end + right_extension.min(reference_length - ends.last_end);

        if let Some(position) = self.read_position {
            let retain_left = matches!(position, ReadPosition::FivePrime) != reverse;
            let base = if retain_left {
                start
            } else {
                end.checked_sub(1).ok_or_else(invalid_extent)?
            };
            if base >= reference_length {
                return Err(invalid_extent());
            }
            let block = Block::new(chrom, base, base + 1)?;
            blocks.clear();
            blocks.push(block);
            return Ok(());
        }

        if start < ends.first_start {
            if ends.first_start < ends.first_end {
                if let Some(first) = blocks.first_mut() {
                    first.start = start;
                }
            } else {
                blocks.insert(0, Block::new(chrom, start, ends.first_end)?);
            }
        }
        if end > ends.last_end {
            if ends.last_start < ends.last_end {
                if let Some(last) = blocks.last_mut() {
                    last.end = end;
                }
            } else {
                blocks.push(Block::new(chrom, ends.last_start, end)?);
            }
        }
        Ok(())
    }
}
