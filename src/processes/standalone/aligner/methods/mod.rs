//! Traits and structs for the different methods of pairwise alignment used by
//! aligner

use std::fmt::Display;
use zoe::{
    alignment::Alignment,
    data::{fasta::FastaSeq, records::SequenceReadable},
    math::AnyInt,
    prelude::NucleotidesView,
};

#[cfg(not(feature = "dev_no_rayon"))]
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

mod striped_sw_local;
mod striped_sw_shared;

pub use striped_sw_local::*;
pub use striped_sw_shared::*;

use crate::aligner::Strand;

/// A trait unifying the different strategies for pairwise sequence alignment
/// used in aligner.
pub trait AlignmentMethod: Sync + Send {
    /// The type of the profile (built by preprocessing either the query or the
    /// reference)
    type Profile<'a>: Send
    where
        Self: 'a;

    /// The type of the output alignment
    type Score: AnyInt;

    /// Builds a profile from a sequence
    fn build_profile<'a, R>(&'a self, profile_seq: &'a R) -> std::io::Result<Self::Profile<'a>>
    where
        R: SequenceReadable + Display;

    /// Aligns a preprocessed sequence against a regular sequence.
    ///
    /// The alignment returned assumes that the profile corresponds to the
    /// query. If it corresponds to the reference, then the alignment needs to
    /// be inverted.
    ///
    /// If the alignment is unmapped, [`None`] should be returned.
    fn align_one(&self, profile: &Self::Profile<'_>, regular_seq: &[u8]) -> std::io::Result<Option<Alignment<Self::Score>>>;

    /// Aligns a preprocessed sequence against a regular sequence and the
    /// reverse complement of that sequence (if passed).
    ///
    /// The alignment returned assumes that the profile corresponds to the
    /// query. If it corresponds to the reference, then the alignment needs to
    /// be inverted.
    ///
    /// If the alignment is unmapped, [`None`] is returned.
    fn align(
        &self, profile: &Self::Profile<'_>, regular_seq: &[u8], rc_regular_seq: Option<&Vec<u8>>,
    ) -> std::io::Result<Option<(Alignment<Self::Score>, Strand)>> {
        let forward_alignment = self.align_one(profile, regular_seq)?;

        let revcomp_alignment = match rc_regular_seq {
            Some(rc_regular_seq) => self.align_one(profile, rc_regular_seq)?,
            None => {
                if let Some(forward_alignment) = forward_alignment {
                    return Ok(Some((forward_alignment, Strand::Forward)));
                } else {
                    return Ok(None);
                }
            }
        };

        let out = match (forward_alignment, revcomp_alignment) {
            (None, None) => None,
            (None, Some(revcomp)) => Some((revcomp, Strand::Reverse)),
            (Some(forward), None) => Some((forward, Strand::Forward)),
            (Some(forward), Some(revcomp)) => {
                if revcomp.score > forward.score {
                    Some((revcomp, Strand::Reverse))
                } else {
                    Some((forward, Strand::Forward))
                }
            }
        };

        Ok(out)
    }

    /// Builds all the profiles for the sequences, returning a vector where each
    /// element is a tuple of the original record and the corresponding profile.
    fn zip_with_profiles<'a>(
        &'a self, profile_seqs: &'a [FastaSeq],
    ) -> std::io::Result<Vec<(&'a FastaSeq, Self::Profile<'a>)>> {
        profile_seqs
            .maybe_par_iter()
            .map(|profile_seq| -> std::io::Result<_> { Ok((profile_seq, self.build_profile(profile_seq)?)) })
            .collect()
    }

    /// If `rev_comp` is true, gets the reverse complement of all the
    /// `regular_seqs`, returning a vector where each element is a tuple of the
    /// original sequence and the reverse complement. If `rev_comp` is false,
    /// use `None` for the second tuple value.
    fn maybe_zip_with_revcomp<'a>(
        &'a self, regular_seqs: &'a [FastaSeq], rev_comp: bool,
    ) -> Vec<(&'a FastaSeq, Option<Vec<u8>>)> {
        if rev_comp {
            regular_seqs
                .maybe_par_iter()
                .map(|seq| {
                    (
                        seq,
                        Some(
                            NucleotidesView::from(seq.sequence.as_slice())
                                .to_reverse_complement()
                                .into_vec(),
                        ),
                    )
                })
                .collect()
        } else {
            regular_seqs.iter().map(|seq| (seq, None)).collect()
        }
    }
}

/// An extension trait providing [`maybe_par_iter`], which calls [`par_iter`]
/// when the `one-thread` feature is not set and [`into_iter`] otherwise.
///
/// [`par_iter`]: rayon::iter::IntoParallelRefIterator::par_iter
/// [`into_iter`]: IntoIterator::into_iter
#[cfg(not(feature = "dev_no_rayon"))]
trait MaybeParIter<'a>: IntoParallelRefIterator<'a> {
    /// Calls [`par_iter`] if the `one-thread` feature is not set and
    /// [`into_iter`] otherwise.
    ///
    /// [`par_iter`]: rayon::iter::IntoParallelRefIterator::par_iter
    /// [`into_iter`]: IntoIterator::into_iter
    #[inline]
    fn maybe_par_iter(&'a self) -> Self::Iter {
        self.par_iter()
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl<'a, T: ?Sized + IntoParallelRefIterator<'a>> MaybeParIter<'a> for T {}

/// An extension trait providing [`maybe_par_iter`], which calls [`par_iter`]
/// when the `one-thread` feature is not set and [`into_iter`] otherwise.
///
/// [`maybe_par_iter`]: MaybeParIter::maybe_par_iter
/// [`par_iter`]: rayon::iter::IntoParallelRefIterator::par_iter
/// [`into_iter`]: IntoIterator::into_iter
#[cfg(feature = "dev_no_rayon")]
trait MaybeParIter<'a>
where
    &'a Self: IntoIterator,
    Self: 'a, {
    /// Calls [`par_iter`] if the `one-thread` feature is not set and
    /// [`into_iter`] otherwise.
    ///
    /// [`par_iter`]: rayon::iter::IntoParallelRefIterator::par_iter
    /// [`into_iter`]: IntoIterator::into_iter
    #[inline]
    fn maybe_par_iter(&'a self) -> <&'a Self as IntoIterator>::IntoIter {
        self.into_iter()
    }
}

#[cfg(feature = "dev_no_rayon")]
impl<'a, T: 'a + ?Sized> MaybeParIter<'a> for T where &'a T: IntoIterator {}
