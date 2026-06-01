pub mod zip_reads;

mod deinterleave;
mod header_error;
mod id_side;

pub use deinterleave::{DeinterleaveError, DeinterleavedPairedReads, DeinterleavedPairedReadsExt};
pub use header_error::PairedHeaderError;
pub use id_side::{ReadSide, check_paired_headers, get_molecular_id_side};
pub use zip_reads::{ZipPairedReadsError, ZipPairedReadsExt, ZipReadsError};

#[cfg(test)]
mod test;
