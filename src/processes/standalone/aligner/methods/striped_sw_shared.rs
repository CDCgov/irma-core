//! The striped Smith Waterman algorithm

use crate::aligner::methods::AlignmentMethod;
use std::{borrow::Borrow, fmt::Display};
use zoe::{
    alignment::{Alignment, MaybeAligned, SharedProfiles},
    data::{matrices::WeightMatrix, records::SequenceReadable},
    prelude::ProfileSets,
};

/// An [`AlignmentMethod`] for performing Striped Smith Waterman with
/// [`SharedProfiles`] (profiles that are not shared between threads).
///
/// A register width of 256 is assumed.
///
/// ## Parameters
///
/// * `W` - The type for the weight matrix (either owned or a reference)
/// * `S` - The alphabet size
pub struct StripedSmithWatermanShared<W, const S: usize>
where
    W: Borrow<WeightMatrix<'static, i8, S>>, {
    weight_matrix: W,
    gap_open:      i8,
    gap_extend:    i8,
}

impl<W, const S: usize> StripedSmithWatermanShared<W, S>
where
    W: Borrow<WeightMatrix<'static, i8, S>>,
{
    /// Creates a new [`StripedSmithWatermanShared`] alignment method with the
    /// provided scoring parameters.
    pub fn new(weight_matrix: W, gap_open: i8, gap_extend: i8) -> Self {
        Self {
            weight_matrix,
            gap_open,
            gap_extend,
        }
    }
}

impl<W, const S: usize> AlignmentMethod for StripedSmithWatermanShared<W, S>
where
    W: Borrow<WeightMatrix<'static, i8, S>> + Sync + Send + 'static,
{
    type Profile<'a>
        = SharedProfiles<'a, 32, 16, 8, S>
    where
        Self: 'a;

    type Score = u32;

    fn build_profile<'a, R>(&'a self, profile_seq: &'a R) -> std::io::Result<Self::Profile<'a>>
    where
        R: SequenceReadable + Display, {
        SharedProfiles::new_with_w256(
            Box::from(profile_seq.sequence_bytes()),
            self.weight_matrix.borrow(),
            self.gap_open,
            self.gap_extend,
        )
        .map_err(|e| {
            std::io::Error::other(format!(
                "The following sequence failed to be preprocessed:\n\n{profile_seq}\n\nThis was caused by:\n{e}"
            ))
        })
    }

    fn align_one(&self, profile: &Self::Profile<'_>, regular_seq: &[u8]) -> std::io::Result<Option<Alignment<Self::Score>>> {
        // TODO: Offer starting integer option?
        match profile.smith_waterman_alignment_from_i8(regular_seq) {
            MaybeAligned::Some(alignment) => Ok(Some(alignment)),
            MaybeAligned::Overflowed => Err(std::io::Error::other("The score exceeded the capacity of i32!")),
            MaybeAligned::Unmapped => Ok(None),
        }
    }
}
