//! The striped Smith Waterman algorithm

use crate::aligner::methods::AlignmentMethod;
use std::borrow::Borrow;
use zoe::{
    alignment::{Alignment, LocalProfiles, MaybeAligned},
    data::{
        err::{ErrorWithContext, ResultWithErrorContext},
        matrices::WeightMatrix,
        records::{HeaderReadable, SequenceReadable},
    },
    prelude::{ProfileSets, SeqSrc},
};

/// An [`AlignmentMethod`] for performing Striped Smith Waterman with
/// [`LocalProfiles`] (profiles that are not shared between threads).
///
/// A register width of 256 is assumed.
///
/// ## Parameters
///
/// * `M` - The type for the weight matrix (either owned or a reference)
/// * `S` - The alphabet size
pub struct StripedSmithWatermanLocal<M, const S: usize>
where
    M: Borrow<WeightMatrix<'static, i8, S>>, {
    weight_matrix: M,
    gap_open:      i8,
    gap_extend:    i8,
}

impl<M, const S: usize> StripedSmithWatermanLocal<M, S>
where
    M: Borrow<WeightMatrix<'static, i8, S>>,
{
    /// Creates a new [`StripedSmithWatermanLocal`] alignment method with the
    /// provided scoring parameters.
    pub fn new(weight_matrix: M, gap_open: i8, gap_extend: i8) -> Self {
        Self {
            weight_matrix,
            gap_open,
            gap_extend,
        }
    }
}

impl<M, const S: usize> AlignmentMethod for StripedSmithWatermanLocal<M, S>
where
    M: Borrow<WeightMatrix<'static, i8, S>> + Sync + Send + 'static,
{
    type Profile<'a>
        = LocalProfiles<'a, 32, 16, 8, S>
    where
        Self: 'a;

    type Score = u32;

    fn build_profile<'a, R>(&'a self, profile_seq: &'a R) -> Result<Self::Profile<'a>, ErrorWithContext>
    where
        R: SequenceReadable + HeaderReadable, {
        LocalProfiles::new_with_w256(
            profile_seq.sequence_bytes(),
            self.weight_matrix.borrow(),
            self.gap_open,
            self.gap_extend,
        )
        .with_context(format!(
            "The sequence with the following header failed to be preprocessed: {}",
            profile_seq.header()
        ))
    }

    #[inline]
    fn align_one(
        &self, profile: &Self::Profile<'_>, regular_seq: SeqSrc<&[u8]>,
    ) -> std::io::Result<Option<Alignment<Self::Score>>> {
        // TODO: Offer starting integer option?
        match profile.sw_align_from_i8(regular_seq) {
            MaybeAligned::Some(alignment) => Ok(Some(alignment)),
            MaybeAligned::Overflowed => Err(std::io::Error::other("The score exceeded the capacity of i32!")),
            MaybeAligned::Unmapped => Ok(None),
        }
    }
}
