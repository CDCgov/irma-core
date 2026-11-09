use crate::aligner::{
    arg_parsing::{AlignerConfig, Alphabet, EitherParams, NumPasses, ParsedAlignerArgs, WhichSequence, parse_aligner_args},
    writers::{write_alignment, write_header},
};
use clap::{Args, builder::RangedI64ValueParser};
use irma_records::io::{FastX, FastXReader, Finish, OutputOptions, ReadFileZipInThread, ValidatePaths, WithContext};
use std::{cell::RefCell, cmp::Ordering, path::PathBuf};
use zoe::{
    alignment::{Alignment, LocalProfiles, MaybeAligned, ScoreParams, SharedProfiles, TracebackBuffer},
    data::{
        err::{ResultWithErrorContext, ResultWithSubitem},
        fasta::Fasta,
    },
    prelude::{NucleotidesView, ProfileSets, SeqSrc},
};

#[cfg(not(feature = "dev_no_rayon"))]
use irma_records::io::WriterThreaded;

#[cfg(not(feature = "dev_no_rayon"))]
use rayon::iter::{ParallelBridge, ParallelIterator};

/// The capacity to use for [`WriterThreaded::with_batch_capacity`].
#[cfg(not(feature = "dev_no_rayon"))]
const WRITE_BATCH_CAPACITY: usize = 64 * 1024;

mod arg_parsing;
mod writers;

/// A type alias for the query reader used by `aligner`.
type QueryReader = WithContext<FastXReader<ReadFileZipInThread>>;

thread_local! {
    /// A reusable traceback buffer for 1-pass alignment to reduce impact of
    /// repeated large allocations.
    static TRACEBACK_BUFFER: RefCell<TracebackBuffer> = RefCell::new(TracebackBuffer::new());
}

/// A type alias for the writer being used for the SAM file, which depends on
/// whether `dev_no_rayon` is set.
#[cfg(not(feature = "dev_no_rayon"))]
type SamWriter = WriterThreaded;

/// A type alias for the writer being used for the SAM file, which depends on
/// whether `dev_no_rayon` is set.
#[cfg(feature = "dev_no_rayon")]
type SamWriter = irma_records::io::WriteFileStdout;

/// The command line arguments for `aligner`
#[derive(Args, Debug)]
pub struct AlignerArgs {
    /// Path to the FASTA file containing the reference sequence(s)
    ref_file: PathBuf,

    /// Path to the FASTQ or FASTA file containing the query sequence(s)
    query_file: PathBuf,

    #[arg(long, alias = "out")]
    /// Output filepath for alignments. Alignments print to STDOUT if not
    /// provided
    output: Option<PathBuf>,

    #[arg(short = 'r', long)]
    /// Also align against the reverse complement, reporting the alignment with
    /// the better score
    rev_comp: bool,

    #[arg(short = 'm', long, value_parser = RangedI64ValueParser::<u8>::new().range(0..=127))]
    /// The score for a match, in [0, 127]
    matching: Option<u8>,

    #[arg(short = 'x', long, value_parser = RangedI64ValueParser::<u8>::new().range(0..=127))]
    /// The penalty for a mismatch, expressed as a nonnegative value in [0, 127]
    mismatch: Option<u8>,

    #[arg(short = 'o', long, default_value_t = 10, value_parser = RangedI64ValueParser::<u8>::new().range(0..=127))]
    /// The penalty for opening a gap, expressed as a nonnegative value in [0,
    /// 127]
    gap_open: u8,

    #[arg(short = 'e', long, default_value_t = 1, value_parser = RangedI64ValueParser::<u8>::new().range(0..=127))]
    /// The penalty for extending a gap, expressed as a nonnegative value in [0,
    /// 127]
    gap_extend: u8,

    #[arg(long, conflicts_with_all = ["matching", "mismatch", "ignore_n"])]
    /// The protein substitution matrix to use, specified by name. This defaults
    /// to `blosum62` when `alphabet` is `AA`
    matrix: Option<String>,

    #[arg(long)]
    /// If specified, any base involving N has a score of 0. This is only
    /// allowed when alphabet is DNA
    ignore_n: bool,

    #[arg(long)]
    /// The alphabet to use. [defaults: DNA, if --matrix, then AA]
    alphabet: Option<Alphabet>,

    #[arg(long, conflicts_with = "profile_from_query")]
    /// Builds the profile from the reference sequences. This is default when
    /// --method 3pass is used or left unspecified
    profile_from_ref: bool,

    #[arg(long, conflicts_with = "profile_from_ref")]
    /// Builds the profile from the query sequences. This is default when
    /// --method 1pass is used
    profile_from_query: bool,

    #[arg(long, default_value_t = NumPasses::ThreePass)]
    /// The method to use for alignment
    method: NumPasses,

    #[arg(long)]
    /// Excludes the unmapped alignments from the final alignment
    exclude_unmapped: bool,

    #[arg(long)]
    /// Only output the best scoring alignment for each query
    best_match: bool,

    #[arg(long)]
    /// Set the code to use only one thread for performing alignments
    single_thread: bool,

    #[arg(long)]
    /// Include the SAM header line
    header: bool,
}

impl ValidatePaths for AlignerArgs {
    fn inputs(&self) -> impl IntoIterator<Item = &PathBuf> {
        [&self.ref_file, &self.query_file]
    }

    fn outputs(&self) -> impl IntoIterator<Item = &PathBuf> {
        self.output.iter()
    }
}

/// Sub-program for performing sequence alignment
pub fn aligner_process(args: AlignerArgs) -> std::io::Result<()> {
    args.validate_paths()?;

    let ParsedAlignerArgs {
        query_reader,
        references,
        params,
        header,
        config,
    } = parse_aligner_args(args)?;

    #[cfg(not(feature = "dev_no_rayon"))]
    if config.single_thread {
        rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    }

    let mut writer = OutputOptions::new_from_opt_path(config.output.as_ref())
        .use_file_or_stdout()
        .open()?;

    if header {
        write_header(&mut writer, &references)?;
    }

    #[cfg(not(feature = "dev_no_rayon"))]
    let writer = WriterThreaded::with_batch_capacity(writer, WRITE_BATCH_CAPACITY);

    dispatch_alphabet(query_reader, references, writer, params, &config)?;

    Ok(())
}

/// Dispatches the aligner based on the alphabet and weight matrix.
///
/// ## Errors
///
/// Errors while reading the queries, performing the alignment, and writing the
/// alignment are propagated. Context containing the header(s) is added for
/// failed alignment.
fn dispatch_alphabet(
    query_reader: QueryReader, references: Vec<Fasta>, writer: SamWriter, params: EitherParams, config: &AlignerConfig,
) -> std::io::Result<()> {
    match params {
        EitherParams::Dna(params) => dispatch_method(query_reader, references, writer, params, config),
        EitherParams::Aa(params) => dispatch_method(query_reader, references, writer, params, config),
    }
}

/// Dispatches the aligner based on whether the profile is built from the query
/// or reference, whether all alignments or just the best match alignment are
/// being reported, and whether the three-pass algorithm is used.
///
/// ## Errors
///
/// Errors while reading the queries, performing the alignment, and writing the
/// alignment are propagated. Context containing the header(s) is added for
/// failed alignment.
fn dispatch_method<const S: usize>(
    query_reader: QueryReader, references: Vec<Fasta>, writer: SamWriter, params: ScoreParams<'static, i8, S>,
    config: &AlignerConfig,
) -> std::io::Result<()> {
    let references = References::new(&references, &params, config.rev_comp);

    if config.best_match {
        align_queries(query_reader, writer, |writer, query| {
            align_query_best_match(query, &references, writer, &params, config)
        })
    } else {
        align_queries(query_reader, writer, |writer, query| {
            align_query_all(query, &references, writer, &params, config)
        })
    }
}

/// Aligns a query against all references, writing the outputs to `writer`.
///
/// ## Errors
///
/// Errors while performing and writing the alignment are propagated. Context
/// containing the header is added for failed alignment.
fn align_query_all<'r, const S: usize>(
    query: FastX, references: &References<'r, S>, writer: &mut SamWriter, params: &ScoreParams<i8, S>,
    config: &AlignerConfig,
) -> std::io::Result<()> {
    match (config.method, config.profile_from) {
        (NumPasses::OnePass, WhichSequence::Query) => {
            let query = QueryWithProfile::new(&query, params);

            for reference in references {
                let alignment = query.sw_1pass_query_profile(reference)?;
                write_alignment(writer, alignment, config)?;
            }
        }
        (NumPasses::OnePass, WhichSequence::Reference) => {
            let query = QueryWithRc::new(&query, config.rev_comp);

            for reference in references {
                let alignment = reference.sw_1pass_ref_profile(&query)?;
                write_alignment(writer, alignment, config)?;
            }
        }
        (NumPasses::ThreePass, WhichSequence::Query) => {
            let query = QueryWithProfile::new(&query, params);

            for reference in references {
                let alignment = query.sw_3pass_query_profile(reference)?;
                write_alignment(writer, alignment, config)?;
            }
        }
        (NumPasses::ThreePass, WhichSequence::Reference) => {
            let query = QueryWithRc::new(&query, config.rev_comp);

            for reference in references {
                let alignment = reference.sw_3pass_ref_profile(&query)?;
                write_alignment(writer, alignment, config)?;
            }
        }
    }

    Ok(())
}

/// Aligns a query against all references, writing the best match to `writer`.
///
/// ## Errors
///
/// Errors while performing and writing the alignment are propagated. Context
/// containing the header is added for failed alignment.
fn align_query_best_match<'r, const S: usize>(
    query: FastX, references: &References<'r, S>, writer: &mut SamWriter, params: &ScoreParams<i8, S>,
    config: &AlignerConfig,
) -> std::io::Result<()> {
    let best_alignment = match (config.method, config.profile_from) {
        (NumPasses::OnePass, WhichSequence::Query) => {
            let query = QueryWithProfile::new(&query, params);
            align_best_ref(references, |reference| query.sw_1pass_query_profile(reference))?
        }
        (NumPasses::OnePass, WhichSequence::Reference) => {
            let query = QueryWithRc::new(&query, config.rev_comp);
            align_best_ref(references, |reference| reference.sw_1pass_ref_profile(&query))?
        }
        (NumPasses::ThreePass, WhichSequence::Query) => {
            let query = QueryWithProfile::new(&query, params);
            align_best_ref(references, |reference| query.sw_3pass_query_profile(reference))?
        }
        (NumPasses::ThreePass, WhichSequence::Reference) => {
            let query = QueryWithRc::new(&query, config.rev_comp);
            align_best_ref(references, |reference| reference.sw_3pass_ref_profile(&query))?
        }
    };

    write_alignment(writer, best_alignment, config)
}

/// Performs all alignments as indicated by closure `f`, using either a parallel
/// iterator (`par_bridge`) or a serial iterator depending on the `dev_no_rayon`
/// feature.
///
/// This implementation is for the serial case.
///
/// ## Errors
///
/// Any IO errors occurring within the [`QueryReader`], while writing to
/// `writer`, while finalizing `writer`, or while calling `f` are propagated.
#[inline]
#[cfg(feature = "dev_no_rayon")]
fn align_queries<F>(query_reader: QueryReader, mut writer: SamWriter, f: F) -> std::io::Result<()>
where
    F: Fn(&mut SamWriter, FastX) -> std::io::Result<()> + Sync + Send, {
    let mut query_reader = query_reader;
    query_reader.try_for_each(|record| f(&mut writer, record?))?;
    writer.finish()
}

/// Performs all alignments as indicated by closure `f`, using either a parallel
/// iterator (`par_bridge`) or a serial iterator depending on the `dev_no_rayon`
/// feature.
///
/// This implementation is for the parallel case.
///
/// ## Errors
///
/// Any IO errors occurring within the [`QueryReader`], within the writer
/// thread, while finalizing `writer`, or while calling `f` are propagated.
#[inline]
#[cfg(not(feature = "dev_no_rayon"))]
fn align_queries<F>(query_reader: QueryReader, writer: WriterThreaded, f: F) -> std::io::Result<()>
where
    F: Fn(&mut WriterThreaded, FastX) -> std::io::Result<()> + Sync + Send, {
    query_reader
        .par_bridge()
        .try_for_each_with(writer.clone(), |w, record| f(w, record?))?;
    writer.finish()
}

/// An enum to track which strand the alignment mapped to
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, Default)]
pub enum Strand {
    /// The alignment was with the forward sequence
    #[default]
    Forward,
    /// The alignment was with the reverse complement
    Reverse,
}

/// The reverse complement of a sequence, if `--rev-comp` was passed.
#[derive(Clone, Eq, PartialEq, Debug)]
pub struct MaybeRevComp(Option<Vec<u8>>);

impl MaybeRevComp {
    /// Computes the reverse complement of a sequence, if `rev_comp` is true.
    pub fn new(seq: &[u8], rev_comp: bool) -> Self {
        MaybeRevComp(rev_comp.then(|| NucleotidesView::from(seq).to_reverse_complement().into_vec()))
    }
}

/// A query record together with a set of profiles for alignment.
pub struct QueryWithProfile<'q, const S: usize> {
    /// The record for the query.
    forward: &'q FastX,
    /// The profile set for the query sequence.
    profile: LocalProfiles<'q, 32, 16, 8, S>,
}

impl<'q, const S: usize> QueryWithProfile<'q, S> {
    /// Bundles a query record together with a corresponding profile set for use
    /// in alignment.
    ///
    /// ## Errors
    ///
    /// An error containing the query header as context is returned if a profile
    /// fails to be made from the sequence.
    pub fn new(query: &'q FastX, params: &'q ScoreParams<i8, S>) -> Self {
        let forward = query;

        let forward_seq = forward.sequence.as_slice();
        let profile = LocalProfiles::new(forward_seq, params);

        Self { forward, profile }
    }

    /// Aligns the query profile against the provided reference using the 1-pass
    /// algorithm.
    ///
    /// This is a higher-level abstraction around [`AlignerMethods::sw_1pass`]
    /// that handles the reverse complement alignment (if `--rev-comp` is used)
    /// and returns an [`AlignmentAndSeqs`].
    ///
    /// ## Errors
    ///
    /// If the alignment fails (due to overflow), context with the query and
    /// reference header is added. If it was the reverse complement alignment
    /// that failed, context is also added mentioning this.
    pub fn sw_1pass_query_profile<'r>(&self, reference: &Reference<'r, S>) -> std::io::Result<AlignmentAndSeqs<'q, 'r>> {
        let mapping = align_maybe_rc(SeqSrc::Reference(&reference.forward.sequence), &reference.reverse, |seq| {
            self.profile.sw_1pass(seq)
        })
        .with_context("Failed to align the sequences with the following headers:")
        .with_subitem(format!(
            "Query: {q_header}\nReference: {r_header}",
            q_header = self.forward.header,
            r_header = reference.forward.header
        ))?;

        Ok(AlignmentAndSeqs {
            mapping,
            query: self.forward,
            reference: reference.forward,
        })
    }

    /// Aligns the query profile against the provided reference using the 3-pass
    /// algorithm.
    ///
    /// This is a higher-level abstraction around [`AlignerMethods::sw_3pass`]
    /// that handles the reverse complement alignment (if `--rev-comp` is used)
    /// and returns an [`AlignmentAndSeqs`].
    ///
    /// ## Errors
    ///
    /// If the alignment fails (due to overflow), context with the query and
    /// reference header is added. If it was the reverse complement alignment
    /// that failed, context is also added mentioning this.
    pub fn sw_3pass_query_profile<'r>(&self, reference: &Reference<'r, S>) -> std::io::Result<AlignmentAndSeqs<'q, 'r>> {
        let mapping = align_maybe_rc(SeqSrc::Reference(&reference.forward.sequence), &reference.reverse, |seq| {
            self.profile.sw_3pass(seq)
        })
        .with_context("Failed to align the sequences with the following headers:")
        .with_subitem(format!(
            "Query: {q_header}\nReference: {r_header}",
            q_header = self.forward.header,
            r_header = reference.forward.header
        ))?;

        Ok(AlignmentAndSeqs {
            mapping,
            query: self.forward,
            reference: reference.forward,
        })
    }
}

/// The query record together with its reverse complement (if `--rev-comp` was
/// passed).
pub struct QueryWithRc<'q, const S: usize> {
    /// The record for the query.
    forward: &'q FastX,
    /// The reverse complement of the query, if `--rev-comp` was passed.
    reverse: MaybeRevComp,
}

impl<'q, const S: usize> QueryWithRc<'q, S> {
    /// Bundles a query record together with its reverse complement, if
    /// `rev_comp` is true.
    pub fn new(query: &'q FastX, rev_comp: bool) -> Self {
        let forward = query;

        let forward_seq = forward.sequence.as_slice();
        let reverse = MaybeRevComp::new(forward_seq, rev_comp);

        Self { forward, reverse }
    }

    /// Returns whether `--rev-comp` is being used.
    pub fn rev_comp(&self) -> bool {
        self.reverse.0.is_some()
    }
}

/// A reference record, together with its reverse complement and a profile set
/// for alignment.
pub struct Reference<'r, const S: usize> {
    /// The record for the reference.
    forward: &'r Fasta,
    /// The reverse complement of the reference, if `--rev-comp` was passed.
    reverse: MaybeRevComp,
    /// The profile set for the forward reference sequence.
    profile: SharedProfiles<'r, 32, 16, 8, S>,
}

impl<'r, const S: usize> Reference<'r, S> {
    /// Bundles a reference together with its reverse complement (if `rev_comp`
    /// is true) and a profile set for use in alignment.
    pub fn new(reference: &'r Fasta, params: &'r ScoreParams<'static, i8, S>, rev_comp: bool) -> Self {
        let forward = reference;

        let forward_seq = forward.sequence.as_slice();
        let reverse = MaybeRevComp::new(forward_seq, rev_comp);
        let profile = SharedProfiles::new(forward_seq, params);

        Self {
            forward,
            reverse,
            profile,
        }
    }

    /// Returns whether `--rev-comp` is being used.
    pub fn rev_comp(&self) -> bool {
        self.reverse.0.is_some()
    }

    /// Aligns the reference profile against the provided query using the 1-pass
    /// algorithm.
    ///
    /// This is a higher-level abstraction around [`AlignerMethods::sw_1pass`]
    /// that handles the reverse complement alignment (if `--rev-comp` is used)
    /// and returns an [`AlignmentAndSeqs`].
    ///
    /// ## Errors
    ///
    /// If the alignment fails (due to overflow), context with the query and
    /// reference header is added. If it was the reverse complement alignment
    /// that failed, context is also added mentioning this.
    pub fn sw_1pass_ref_profile<'q>(&self, query: &QueryWithRc<'q, S>) -> std::io::Result<AlignmentAndSeqs<'q, 'r>> {
        let mapping = align_maybe_rc(SeqSrc::Query(&query.forward.sequence), &query.reverse, |seq| {
            self.profile.sw_1pass(seq)
        })
        .with_context("Failed to align the sequences with the following headers:")
        .with_subitem(format!(
            "Query: {q_header}\nReference: {r_header}",
            q_header = query.forward.header,
            r_header = self.forward.header
        ))?;

        Ok(AlignmentAndSeqs {
            mapping,
            query: query.forward,
            reference: self.forward,
        })
    }

    /// Aligns the reference profile against the provided query using the 3-pass
    /// algorithm.
    ///
    /// This is a higher-level abstraction around [`AlignerMethods::sw_3pass`]
    /// that handles the reverse complement alignment (if `--rev-comp` is used)
    /// and returns an [`AlignmentAndSeqs`].
    ///
    /// ## Errors
    ///
    /// If the alignment fails (due to overflow), context with the query and
    /// reference header is added. If it was the reverse complement alignment
    /// that failed, context is also added mentioning this.
    pub fn sw_3pass_ref_profile<'q>(&self, query: &QueryWithRc<'q, S>) -> std::io::Result<AlignmentAndSeqs<'q, 'r>> {
        let mapping = align_maybe_rc(SeqSrc::Query(&query.forward.sequence), &query.reverse, |seq| {
            self.profile.sw_3pass(seq)
        })
        .with_context("Failed to align the sequences with the following headers:")
        .with_subitem(format!(
            "Query: {q_header}\nReference: {r_header}",
            q_header = query.forward.header,
            r_header = self.forward.header
        ))?;

        Ok(AlignmentAndSeqs {
            mapping,
            query: query.forward,
            reference: self.forward,
        })
    }
}

/// A collection of references to align against.
pub struct References<'r, const S: usize>(Vec<Reference<'r, S>>);

impl<'r, const S: usize> References<'r, S> {
    /// Bundles the `references` with reverse complement and profile information
    /// for alignment.
    pub fn new(references: &'r [Fasta], params: &'r ScoreParams<'static, i8, S>, rev_comp: bool) -> Self {
        References(
            references
                .iter()
                .map(|reference| Reference::new(reference, params, rev_comp))
                .collect(),
        )
    }

    /// Returns an iterator over the references.
    pub fn iter(&self) -> std::slice::Iter<'_, Reference<'r, S>> {
        self.0.iter()
    }
}

impl<'r, 'c, const S: usize> IntoIterator for &'c References<'r, S> {
    type Item = &'c Reference<'r, S>;
    type IntoIter = std::slice::Iter<'c, Reference<'r, S>>;

    #[inline]
    fn into_iter(self) -> Self::IntoIter {
        self.0.iter()
    }
}

/// A trait extending [`ProfileSets`] with methods that add context to errors.
pub trait AlignerMethods<'a, const S: usize>: ProfileSets<'a, 32, 16, 8, S> {
    /// Performs the 1-pass Striped Smith Waterman alignment algorithm starting
    /// from an `i8`, using [`ProfileSets::sw_align_from_i8_with_buf`].
    ///
    /// ## Errors
    ///
    /// If the alignment overflows, then this is returned as an error.
    fn sw_1pass(&self, seq: SeqSrc<&[u8]>) -> std::io::Result<Option<Alignment<u32>>> {
        let seq = seq.map(AsRef::as_ref);
        // Within a single thread, this is only borrowed one at a time, so this
        // will not panic
        let alignment = TRACEBACK_BUFFER.with_borrow_mut(|buffer| self.sw_align_from_i8_with_buf(seq, buffer));
        maybe_aligned_to_option(alignment)
    }

    /// Performs the 1-pass Striped Smith Waterman alignment algorithm starting
    /// from an `i8`, similar to [`ProfileSets::sw_align_from_i8`].
    ///
    /// ## Errors
    ///
    /// If the alignment overflows, then this is returned as an error.
    fn sw_3pass(&self, seq: SeqSrc<&[u8]>) -> std::io::Result<Option<Alignment<u32>>> {
        let seq = seq.map(AsRef::as_ref);
        let alignment = self.sw_align_from_i8_3pass(seq);
        maybe_aligned_to_option(alignment)
    }
}

/// A helper function for converting [`MaybeAligned`] to a result of an option,
/// assuming that a profile set was used to form the alignment.
///
/// [`MaybeAligned::Unmapped`] is converted to `None`.
///
/// ## Errors
///
/// [`MaybeAligned::Overflowed`] is converted to an error.
fn maybe_aligned_to_option(alignment: MaybeAligned<Alignment<u32>>) -> std::io::Result<Option<Alignment<u32>>> {
    match alignment {
        MaybeAligned::Some(alignment) => Ok(Some(alignment)),
        MaybeAligned::Overflowed => Err(std::io::Error::other("The score exceeded the capacity of i32!")),
        MaybeAligned::Unmapped => Ok(None),
    }
}

impl<'a, const S: usize> AlignerMethods<'a, S> for LocalProfiles<'a, 32, 16, 8, S> {}
impl<'a, const S: usize> AlignerMethods<'a, S> for SharedProfiles<'a, 32, 16, 8, S> {}

/// Performs an alignment involving `seq`, as well as its reverse complement if
/// it was precomputed in `seq_rc`.
///
/// The alignment to perform is given by `f`, which is a closure accepting the
/// sequence to align as an argument (either the forward or reverse complement
/// of `seq`).
///
/// In the case of a tie, the forward strand is preferred.
///
/// ## Errors
///
/// Any errors from the forward alignment are propagated without additional
/// context. For the reverse alignment, errors are added with context specifying
/// that the reverse complement alignment failed.
pub fn align_maybe_rc<T, F>(seq: SeqSrc<&T>, seq_rc: &MaybeRevComp, f: F) -> std::io::Result<Option<AlignmentAndStrand>>
where
    T: AsRef<[u8]> + ?Sized,
    F: Fn(SeqSrc<&[u8]>) -> std::io::Result<Option<Alignment<u32>>>, {
    let alignment_forward = f(seq.map(AsRef::as_ref))?;

    let (alignment, strand) = if let Some(seq_rc) = &seq_rc.0 {
        let seq_rc = seq.map(|_| seq_rc.as_slice());
        let mut alignment_rc = f(seq_rc).with_context("Failed to perform the reverse complement alignment")?;

        if alignment_rc > alignment_forward {
            if seq.is_reference()
                && let Some(alignment) = alignment_rc.as_mut()
            {
                alignment.make_reverse();
            }

            (alignment_rc, Strand::Reverse)
        } else {
            (alignment_forward, Strand::Forward)
        }
    } else {
        (alignment_forward, Strand::Forward)
    };

    Ok(alignment.map(|alignment| AlignmentAndStrand {
        inner: alignment,
        strand,
    }))
}

/// Performs all the alignments against the provided `reference`, returning the
/// one with the best score.
///
/// The alignment to perform is given by `f`, which is a closure accepting the
/// reference to align against as an argument.
///
/// In the case of a tie, the last reference is preferred.
///
/// ## Errors
///
/// Any errors from the alignments are propagated without context.
///
/// ## Panics
///
/// The `references` provided must be non-empty.
pub fn align_best_ref<'q, 'r, F, const S: usize>(
    references: &References<'r, S>, f: F,
) -> std::io::Result<AlignmentAndSeqs<'q, 'r>>
where
    F: Fn(&Reference<'r, S>) -> std::io::Result<AlignmentAndSeqs<'q, 'r>>, {
    let mut references = references.iter();

    let first_reference = references.next().expect("The references field should be non-empty");
    let mut best_alignment = f(first_reference)?;

    for reference in references {
        let alignment = f(reference)?;
        match alignment.partial_cmp(&best_alignment) {
            Some(Ordering::Greater) | None => best_alignment = alignment,
            _ => {}
        }
    }

    Ok(best_alignment)
}

/// An [`Alignment`] together with the [`Strand`] of the alignment.
#[derive(PartialEq, Debug)]
pub struct AlignmentAndStrand {
    /// The inner alignment, using a score of `u32`.
    pub inner:  Alignment<u32>,
    /// The strand corresponding to the alignment.
    pub strand: Strand,
}

/// An [`Alignment`] together with the [`Strand`] of the alignment, a reference
/// to the query record, and a reference to the reference record. Unmapped
/// alignments are also represented.
///
/// ## Parameters
///
/// - `'q`: The lifetime of the query record
/// - `'r`: The lifetime of the reference record
#[derive(PartialEq)]
pub struct AlignmentAndSeqs<'q, 'r> {
    /// The [`Alignment`] and [`Strand`], if the alignment is mapped. If
    /// unmapped, then this is `None`.
    pub mapping:   Option<AlignmentAndStrand>,
    /// A reference to the query record.
    pub query:     &'q FastX,
    /// A reference to the reference record.
    pub reference: &'r Fasta,
}

impl PartialOrd for AlignmentAndStrand {
    /// This method returns an ordering between `self` and `other` values if one
    /// exists.
    ///
    /// The ordering is based on the `score` field of the alignment, and
    /// equivalent scores either produce `None` or [`Ordering::Equal`] (if all
    /// other fields are identical).
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        match self.inner.partial_cmp(&other.inner) {
            Some(core::cmp::Ordering::Equal) => {
                if self.strand == other.strand {
                    Some(Ordering::Equal)
                } else {
                    None
                }
            }
            ord => ord,
        }
    }
}

impl PartialOrd for AlignmentAndSeqs<'_, '_> {
    /// This method returns an ordering between `self` and `other` values if one
    /// exists.
    ///
    /// The ordering is based on the `score` field of the alignment, and
    /// equivalent scores either produce `None` or [`Ordering::Equal`] (if all
    /// other fields are identical).
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        match self.mapping.partial_cmp(&other.mapping) {
            Some(core::cmp::Ordering::Equal) => {
                if self == other {
                    Some(Ordering::Equal)
                } else {
                    None
                }
            }
            ord => ord,
        }
    }
}

#[cfg(test)]
mod traceback_tests {
    use super::*;

    #[test]
    fn reuse_across_queries_matches_allocating_api() {
        let params = ScoreParams::new_dna_params(2, -5, -10, -1, None).unwrap();
        let long = b"ACGT".repeat(100);
        let mismatch = b"T".repeat(45);
        // Growth, smaller prefixes, unmapped results, and i8 overflow retries.
        for query in [b"ACGT".as_slice(), long.as_slice(), mismatch.as_slice(), b"ACGT"] {
            let local = LocalProfiles::new_with_w256(query, &params);
            let shared = SharedProfiles::new_with_w256(query, &params);
            for target in [b"ACGT".as_slice(), long.as_slice(), mismatch.as_slice(), b""] {
                for seq in [SeqSrc::Reference(target), SeqSrc::Query(target)] {
                    let expected = maybe_aligned_to_option(local.sw_align_from_i8(seq)).unwrap();
                    assert_eq!(local.sw_1pass(seq).unwrap(), expected);
                    assert_eq!(shared.sw_1pass(seq).unwrap(), expected);
                }
            }
        }
    }
}
