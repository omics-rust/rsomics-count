//! Genomic feature counting and count-matrix normalization.

#![forbid(unsafe_code)]

mod alignment;
mod annotation;
mod count;
mod model;
mod normalize;
mod output;

pub use annotation::{AnnotationFormat, AnnotationOptions};
pub use count::{
    Contribution, CountOptions, CountReport, CountSummary, LibraryReport, PairOptions, SplitMode,
    count,
};
pub use model::{CountLevel, StrandMode};
pub use normalize::{NormalizeMode, NormalizeOptions, NormalizeReport, normalize};
