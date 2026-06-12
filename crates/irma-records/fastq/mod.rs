//! Operations for editing FASTQ reads before analysis.
//!
//! This includes hard clipping, masking, adapter or barcode removal, primer
//! handling, poly-G cleanup, canonical base recoding, and read-quality
//! summaries. These operations are exposed through the [`ReadTransforms`]
//! trait.

mod transforms;

pub use transforms::*;

#[cfg(test)]
mod test;
