use crate::aligner::{
    AlignmentAndSeqs, AlignmentMethod, References,
    arg_parsing::{AlignerConfig, NumPasses, WhichSequence},
};
use std::sync::atomic::{AtomicU64, Ordering};
use zoe::{alignment::sw::max_score_for_int_type, data::WeightMatrix};

/// A collection of non-blocking thread-safe tallies for the queries
/// encountered.
#[derive(Debug, Default)]
pub struct QueryTallies {
    /// The number of queries.
    num_queries:         AtomicU64,
    /// The number of queries of length at most 300.
    queries_at_most_300: AtomicU64,
}

impl QueryTallies {
    pub fn tally(&self, query: &[u8]) {
        self.num_queries.fetch_add(1, Ordering::Relaxed);
        if query.len() <= 300 {
            self.queries_at_most_300.fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn num_queries(&self) -> u64 {
        self.num_queries.load(Ordering::Relaxed)
    }

    pub fn queries_at_most_300(&self) -> u64 {
        self.queries_at_most_300.load(Ordering::Relaxed)
    }
}

/// A collection information regarding the reference sequences.
#[derive(Debug)]
pub struct RefTallies {
    /// The total number of references being aligned against.
    num_refs:      usize,
    /// The length of the first reference sequence.
    first_ref_len: usize,
}

impl RefTallies {
    pub fn new<const S: usize>(references: &References<'_, S>) -> Self {
        Self {
            num_refs:      references.0.len(),
            first_ref_len: references.0.first().map_or(0, |seq| seq.forward.sequence.len()),
        }
    }
}

/// A collection of non-blocking thread-safe tallies for information regarding
/// the outcomes of the alignments performed.
#[derive(Debug, Default)]
pub struct AlignmentTallies {
    /// The number of alignments performed. The forward and reverse complement
    /// are counted as one alignment together.
    num_alignments:    AtomicU64,
    /// The number of scores fitting an `i8` (at most 254).
    ///
    /// This only includes the counts for the alignments against the proper
    /// strand. When `--rev-comp` is enabled, the alignment for the worse strand
    /// is not counted.
    scores_fitting_i8: AtomicU64,
}

impl AlignmentTallies {
    pub fn tally<const S: usize>(&self, alignment: &AlignmentAndSeqs<'_, '_>, matrix: &WeightMatrix<i8, S>) {
        self.num_alignments.fetch_add(1, Ordering::Relaxed);
        if let Some(alignment) = &alignment.mapping {
            if alignment.inner.score <= max_score_for_int_type::<i8, i8, S>(matrix) {
                self.scores_fitting_i8.fetch_add(1, Ordering::Relaxed);
            }
        } else {
            // Unmapped also fits
            self.scores_fitting_i8.fetch_add(1, Ordering::Relaxed);
        }
    }

    pub fn num_alignments(&self) -> u64 {
        self.num_alignments.load(Ordering::Relaxed)
    }

    pub fn scores_fitting_i8(&self) -> u64 {
        self.scores_fitting_i8.load(Ordering::Relaxed)
    }
}

pub fn pick_alignment_method(
    _query_tallies: &QueryTallies, ref_tallies: &RefTallies, alignment_tallies: &AlignmentTallies, config: &AlignerConfig,
) -> AlignmentMethod {
    let mut num_alignments = alignment_tallies.num_alignments();
    let mut scores_fitting_i8 = alignment_tallies.scores_fitting_i8();

    if config.rev_comp {
        // Assume that all incorrect-strand alignments fit an i8. This may not
        // be strictly true for long-read strand chimera or ITRs
        scores_fitting_i8 += num_alignments;
        num_alignments *= 2;
    }

    let num_passes = config.method.unwrap_or_else(|| {
        // Validity: this is a best effort snapshot and may be slightly skewed
        // based on relaxed ordering. This affects only performance and not
        // correctness.
        if num_alignments >= 50 && scores_fitting_i8 as f64 / num_alignments as f64 >= 0.66 {
            NumPasses::OnePass
        } else {
            NumPasses::ThreePass
        }
    });

    let profile_from = config.profile_from.unwrap_or(match num_passes {
        NumPasses::OnePass => {
            if ref_tallies.first_ref_len < 600 {
                WhichSequence::Query
            } else {
                WhichSequence::Reference
            }
        }
        NumPasses::ThreePass => WhichSequence::Query,
    });

    match (num_passes, profile_from) {
        (NumPasses::OnePass, WhichSequence::Query) => AlignmentMethod::OnePassQueryProfile,
        (NumPasses::OnePass, WhichSequence::Reference) => AlignmentMethod::OnePassRefProfile,
        (NumPasses::ThreePass, WhichSequence::Query) => AlignmentMethod::ThreePassQueryProfile,
        (NumPasses::ThreePass, WhichSequence::Reference) => AlignmentMethod::ThreePassRefProfile,
    }
}

pub struct AllTallies {
    /// The number of queries.
    pub num_queries:           u64,
    /// The number of queries of length at most 300.
    pub queries_at_most_300:   u64,
    /// The total number of references being aligned against.
    pub num_refs:              usize,
    /// The length of the first reference sequence.
    pub first_ref_len:         usize,
    /// The number of alignments performed (including both forward and reverse
    /// complement).
    pub num_alignments:        u64,
    /// The estimated number of scores fitting `i8` (at most 254). This is
    /// counted exactly for the highest scoring strand, it is assumed that the
    /// other strand (if `--rev-comp` is enabled) will always fit an `i8`.
    pub est_scores_fitting_i8: u64,
}

impl AllTallies {
    pub fn new(
        query_tallies: &QueryTallies, ref_tallies: &RefTallies, alignment_tallies: &AlignmentTallies, config: &AlignerConfig,
    ) -> Self {
        let mut num_alignments = alignment_tallies.num_alignments();
        let mut est_scores_fitting_i8 = alignment_tallies.scores_fitting_i8();

        if config.rev_comp {
            // Assume that all incorrect-strand alignments fit an i8. This may
            // not be strictly true for long-read strand chimera or ITRs
            est_scores_fitting_i8 += num_alignments;
            num_alignments *= 2;
        }

        Self {
            num_queries: query_tallies.num_queries(),
            queries_at_most_300: query_tallies.queries_at_most_300(),
            num_refs: ref_tallies.num_refs,
            first_ref_len: ref_tallies.first_ref_len,
            num_alignments,
            est_scores_fitting_i8,
        }
    }
}
