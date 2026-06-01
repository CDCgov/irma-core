mod aligned;
mod cigar;
mod merge_pairs;
mod paired_stats;

pub use merge_pairs::*;
pub use paired_stats::*;

pub(crate) use aligned::*;
pub(crate) use cigar::*;

#[cfg(test)]
mod test;
