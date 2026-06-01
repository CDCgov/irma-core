use crate::sam::merge_pairs::SamExpandableAlignment;
use std::ops::Range;
use zoe::data::sam::SamData;

/// Struct holding the alignment information of a SAM query sequence to some
/// reference.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub(crate) struct SamAligned {
    /// Query bases aligned to the reference (insertions removed, deletions
    /// added).
    pub(crate) aligned:           Vec<u8>,
    /// Query quality scores aligned to the reference (insertions removed,
    /// deletions added).
    pub(crate) qaligned:          Vec<u8>,
    /// Start of range of reference indices. Start is inclusive.
    pub(crate) ref_start:         usize,
    /// End of range of reference indices. The `ref_end` is *exclusive*, thus
    /// the final aligned query base will be at `ref_end - 1`.
    pub(crate) ref_end:           usize,
    /// Vector of query insertion relative to reference. Ordered by reference
    /// indices.
    pub(crate) insertions:        Vec<SamInsertion>,
    /// Number of bases hard/soft clipped at the start
    pub(crate) num_clipped_start: usize,
    /// Number of bases hard/soft clipped at the end
    pub(crate) num_clipped_end:   usize,
}

impl SamAligned {
    /// Initializes a new [`SamAligned`] object from the fields.
    #[inline]
    pub(crate) fn new(
        aligned: Vec<u8>, qaligned: Vec<u8>, ref_start: usize, ref_end: usize, insertions: Vec<SamInsertion>,
        num_clipped_start: usize, num_clipped_end: usize,
    ) -> Self {
        SamAligned {
            aligned,
            qaligned,
            ref_start,
            ref_end,
            insertions,
            num_clipped_start,
            num_clipped_end,
        }
    }

    /// At the given reference index, provides the reference-aligned query's
    /// nucleotide and encoded (ASCII) quality score as an optional tuple of
    /// `(base, qs)`, otherwise, `None` is returned.
    #[must_use]
    pub(crate) fn get_base_and_quality(&self, at_index: usize) -> Option<(u8, u8)> {
        // Safety: reference range is the same length as `aligned` and
        // `qaligned`; thus, subtracting the start of the range will render
        // appropriate 0-based index. If the reference index it outside the
        // range, then `None` is returned.
        if self.ref_range().contains(&at_index) {
            Some((
                self.aligned[at_index - self.ref_start],
                self.qaligned[at_index - self.ref_start],
            ))
        } else {
            None
        }
    }

    /// At the given reference index, provides the reference-aligned query's
    /// nucleotide as an `Option`, otherwise, `None` is returned.
    #[must_use]
    pub(crate) fn get_base(&self, at_index: usize) -> Option<u8> {
        // Safety: reference range is the same length as `aligned` and
        // `qaligned`; thus, subtracting the start of the range will render
        // appropriate 0-based index. If the reference index it outside the
        // range, then `None` is returned.
        if self.ref_range().contains(&at_index) {
            Some(self.aligned[at_index - self.ref_start])
        } else {
            None
        }
    }

    /// Return a half-open reference range as `Range<uize>` for the aligned query.
    #[must_use]
    pub(crate) fn ref_range(&self) -> Range<usize> {
        self.ref_start..self.ref_end
    }

    /// Merges two aligned queries' reference ranges such that the combined
    /// range spans both aligned regions.
    #[must_use]
    pub(crate) fn merge_ref_range(&self, other: &SamAligned) -> Range<usize> {
        std::cmp::min(self.ref_start, other.ref_start)..std::cmp::max(self.ref_end, other.ref_end)
    }

    /// Checks if the `SamAligned` contains an insertion after the 0-based
    /// reference index and returns the query's unaligned range.
    #[must_use]
    pub(crate) fn get_insert_after(&self, reference_index: usize) -> Option<SamInsertion> {
        self.insertions
            .iter()
            .find(|insert| insert.ref_index == reference_index)
            .cloned()
    }
}

impl From<&SamData> for SamAligned {
    /// Builds a [`SamAligned`] from [`SamData`].
    ///
    /// ## Panics
    ///
    /// Currently panics on invalid cigar states. To be removed in the future
    /// for either a `Result` or type-state validated CIGAR strings.
    ///
    /// This also has a chance of panicking in debug mode if an insertion
    /// appears at the start of the alignment.
    #[inline]
    fn from(row: &SamData) -> Self {
        row.get_aligned()
    }
}

/// Represents a single insertion in a [`SamAligned`] record.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub(crate) struct SamInsertion {
    /// The reference index (0-based) before which the insertion occurs.
    pub(crate) ref_index:   usize,
    /// The starting index (0-based) of the insertion in the query sequence.
    pub(crate) query_start: usize,
    /// The ending index (0-based, non-inclusive) of the insertion in the query
    ///   sequence.
    pub(crate) query_end:   usize,
}

impl SamInsertion {
    #[inline]
    #[must_use]
    pub(crate) fn query_range(&self) -> Range<usize> {
        self.query_start..self.query_end
    }
}
