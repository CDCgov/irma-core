use crate::{
    aligner::{
        arg_parsing::{AlignerConfig, Alphabet, AnyMatrix, NumPasses, ParsedAlignerArgs, WhichSequence, parse_aligner_args},
        writers::{AlignmentWriter, write_header},
    },
    io::{FastX, FastXReader, IterWithContext, OutputOptions, ReadFileZipPipe},
};
use clap::{Args, builder::RangedI64ValueParser};
use std::{cmp::Ordering, io::Write, path::PathBuf, sync::atomic::AtomicU64};
use zoe::{
    alignment::{Alignment, LocalProfiles, MaybeAligned, SharedProfiles, sw::max_score_for_int_type},
    data::{err::ResultWithErrorContext, fasta::FastaSeq, matrices::WeightMatrix},
    prelude::{NucleotidesView, ProfileSets, SeqSrc},
};

#[cfg(not(feature = "dev_no_rayon"))]
use crate::aligner::writers::{AlignmentWriterThreaded, ThreadedWriteError};
#[cfg(not(feature = "dev_no_rayon"))]
use rayon::iter::{ParallelBridge, ParallelIterator};

#[cfg(feature = "dev_no_rayon")]
use std::io::Write;

mod arg_parsing;
mod writers;

/// A type alias for the query reader used by `aligner`.
type QueryReader = IterWithContext<FastXReader<ReadFileZipPipe>>;

/// A type alias for the writer being used for the SAM file, which depends on
/// whether `dev_no_rayon` is set.
#[cfg(not(feature = "dev_no_rayon"))]
type SamWriter = AlignmentWriterThreaded;

/// A type alias for the writer being used for the SAM file, which depends on
/// whether `dev_no_rayon` is set.
#[cfg(feature = "dev_no_rayon")]
type SamWriter = crate::io::WriteFileZipStdout;

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
    /// Builds the profile from the reference sequences instead of the queries
    profile_from_ref: bool,

    #[arg(long, conflicts_with = "profile_from_ref")]
    /// Builds the profile from the query sequences (currently the default)
    profile_from_query: bool,

    #[arg(long)]
    /// The method to use for alignment. If not specified, the 1pass algorithm
    /// is used
    method: Option<NumPasses>,

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

    #[arg(long)]
    /// The file to print tally diagnostics to
    tally_diagnostics: Option<PathBuf>,
}

/// Sub-program for performing sequence alignment
pub fn aligner_process(args: AlignerArgs) -> std::io::Result<()> {
    let ParsedAlignerArgs {
        query_reader,
        references,
        weight_matrix,
        header,
        tally_diagnostics,
        config,
    } = parse_aligner_args(args)?;

    #[cfg(not(feature = "dev_no_rayon"))]
    if config.single_thread {
        rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    }

    let mut writer = OutputOptions::new_from_opt_path(config.output.as_ref())
        .use_file_zip_or_stdout()
        .open()?;

    if header {
        write_header(&mut writer, &references)?;
    }

    #[cfg(not(feature = "dev_no_rayon"))]
    let writer = AlignmentWriterThreaded::from_writer(writer);

    // Validity: No context is added to the result
    let tallies = dispatch_alphabet(query_reader, references, writer, weight_matrix, &config)?;

    if let Some(path) = tally_diagnostics {
        let mut tally_diagnostics = OutputOptions::new_from_path(&path).use_file().open()?;

        let UnpackedTallies {
            satisfying: scores_fitting_i8,
            failing: scores_exceeding_i8,
        } = tallies.scores_fitting_i8.unpack();

        let UnpackedTallies {
            satisfying: queries_at_most_300,
            failing: queries_over_300,
        } = tallies.queries_at_most_300.unpack();

        let first_ref_len = tallies.first_ref_len;

        writeln!(tally_diagnostics, "Scores fitting i8: {scores_fitting_i8}")?;
        writeln!(tally_diagnostics, "Scores exceeding i8: {scores_exceeding_i8}")?;
        writeln!(tally_diagnostics, "Queries at most length 300: {queries_at_most_300}")?;
        writeln!(tally_diagnostics, "Queries over length 300: {queries_over_300}")?;
        writeln!(tally_diagnostics, "First reference length: {first_ref_len}")?;
    }

    Ok(())
}

/// Dispatches the aligner based on the alphabet and weight matrix.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. Context containing the
/// header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn dispatch_alphabet(
    query_reader: QueryReader, references: Vec<FastaSeq>, writer: SamWriter, weight_matrix: AnyMatrix<'static, i8>,
    config: &AlignerConfig,
) -> std::io::Result<AlignmentTallies> {
    // Validity: No context is added to the results
    match weight_matrix {
        AnyMatrix::Dna(weight_matrix) => dispatch_method(query_reader, references, writer, &weight_matrix, config),
        AnyMatrix::AaNamed(weight_matrix) => dispatch_method(query_reader, references, writer, weight_matrix, config),
        AnyMatrix::AaSimple(weight_matrix) => dispatch_method(query_reader, references, writer, &weight_matrix, config),
    }
}

/// Dispatches the aligner based on whether the profile is built from the query
/// or reference, whether all alignments or just the best match alignment are
/// being reported, and whether the three-pass algorithm is used.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. Context containing the
/// header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn dispatch_method<const S: usize>(
    query_reader: QueryReader, references: Vec<FastaSeq>, writer: SamWriter, weight_matrix: &WeightMatrix<'static, i8, S>,
    config: &AlignerConfig,
) -> std::io::Result<AlignmentTallies> {
    let references = References::new(
        &references,
        weight_matrix,
        config.gap_open,
        config.gap_extend,
        config.rev_comp,
    )?;

    if config.best_match {
        align_best_match(query_reader, references, writer, weight_matrix, config)
    } else {
        align_all(query_reader, references, writer, weight_matrix, config)
    }
}

/// Aligns all the queries in `query_reader` to the `references`, writing the
/// outputs to `writer`. The method used is specified by the first argument.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. Context containing the
/// header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn align_all<'r, const S: usize>(
    query_reader: QueryReader, references: References<'r, S>, writer: SamWriter,
    weight_matrix: &WeightMatrix<'static, i8, S>, config: &AlignerConfig,
) -> std::io::Result<AlignmentTallies> {
    let tallies = AlignmentTallies::new(&references);

    align_queries(query_reader, writer, |writer, query| {
        let query = query?;
        let method = tallies.pick_alignment_method(config);

        match method {
            AlignmentMethod::OnePassQueryProfile => {
                let query = QueryWithProfile::new(&query, weight_matrix, config.gap_open, config.gap_extend)?;

                for reference in &references {
                    let alignment = query.sw_1pass_query_profile(reference)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    writer.write_alignment(alignment, config)?;
                }
            }
            AlignmentMethod::OnePassRefProfile => {
                let query = QueryWithRc::new(&query, config.rev_comp);

                for reference in &references.0 {
                    let alignment = reference.sw_1pass_ref_profile(&query)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    writer.write_alignment(alignment, config)?;
                }
            }
            AlignmentMethod::ThreePassQueryProfile => {
                let query = QueryWithProfile::new(&query, weight_matrix, config.gap_open, config.gap_extend)?;

                for reference in references.0.iter() {
                    let alignment = query.sw_3pass_query_profile(reference)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    writer.write_alignment(alignment, config)?;
                }
            }
            AlignmentMethod::ThreePassRefProfile => {
                let query = QueryWithRc::new(&query, config.rev_comp);

                for reference in references.0.iter() {
                    let alignment = reference.sw_3pass_ref_profile(&query)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    writer.write_alignment(alignment, config)?;
                }
            }
        }

        Ok(())
    })?;

    Ok(tallies)
}

/// Aligns all the queries in `query_reader` to the `references`, picking the
/// best reference for each and writing that alignment to `writer`. The method
/// used is specified by the first argument.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. Context containing the
/// header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn align_best_match<'r, const S: usize>(
    query_reader: QueryReader, references: References<'r, S>, writer: SamWriter,
    weight_matrix: &WeightMatrix<'static, i8, S>, config: &AlignerConfig,
) -> std::io::Result<AlignmentTallies> {
    let tallies = AlignmentTallies::new(&references);

    align_queries(query_reader, writer, |writer, query| {
        let query = query?;
        let method = tallies.pick_alignment_method(config);

        // Each match statement ends with a write, which appears redundant.
        // However, this is needed since the lifetime of the query is limited to
        // the match statement scope, and hence the alignment will not live long
        // enough to move this after

        match method {
            AlignmentMethod::OnePassQueryProfile => {
                let query = QueryWithProfile::new(&query, weight_matrix, config.gap_open, config.gap_extend)?;

                let best_alignment = align_best_ref(&references, |reference| {
                    let alignment = query.sw_1pass_query_profile(reference)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    Ok(alignment)
                })?;

                writer.write_alignment(best_alignment, config)?;
            }
            AlignmentMethod::OnePassRefProfile => {
                let query = QueryWithRc::new(&query, config.rev_comp);

                let best_alignment = align_best_ref(&references, |reference| {
                    let alignment = reference.sw_1pass_ref_profile(&query)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    Ok(alignment)
                })?;

                writer.write_alignment(best_alignment, config)?;
            }
            AlignmentMethod::ThreePassQueryProfile => {
                let query = QueryWithProfile::new(&query, weight_matrix, config.gap_open, config.gap_extend)?;

                let best_alignment = align_best_ref(&references, |reference| {
                    let alignment = query.sw_3pass_query_profile(reference)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    Ok(alignment)
                })?;

                writer.write_alignment(best_alignment, config)?;
            }
            AlignmentMethod::ThreePassRefProfile => {
                let query = QueryWithRc::new(&query, config.rev_comp);

                let best_alignment = align_best_ref(&references, |reference| {
                    let alignment = reference.sw_3pass_ref_profile(&query)?;
                    tallies.tally_alignment(&alignment, query.forward, weight_matrix);
                    Ok(alignment)
                })?;

                writer.write_alignment(best_alignment, config)?;
            }
        }

        Ok(())
    })?;

    Ok(tallies)
}

/// Performs all alignments as indicated by closure `f`, using either a parallel
/// iterator (`par_bridge`) or a serial iterator depending on the `dev_no_rayon`
/// feature.
///
/// This implementation is for the serial case.
///
/// ## Errors
///
/// Any IO errors occurring while calling `f` or flusing `writer` are
/// propagated.
#[inline]
#[cfg(feature = "dev_no_rayon")]
fn align_queries<F>(query_reader: QueryReader, mut writer: SamWriter, f: F) -> std::io::Result<()>
where
    F: Fn(&mut SamWriter, std::io::Result<FastX>) -> std::io::Result<()> + Sync + Send, {
    let mut query_reader = query_reader;
    query_reader.try_for_each(|query| f(&mut writer, query))?;
    writer.flush()
}

/// Performs all alignments as indicated by closure `f`, using either a parallel
/// iterator (`par_bridge`) or a serial iterator depending on the `dev_no_rayon`
/// feature.
///
/// This implementation is for the parallel case.
///
/// ## Errors
///
/// Any IO errors occurring within the writer thread, while flushing `writer`,
/// or while calling `f` are propagated. If a
/// [`ThreadedWriteError::ReceiverDeallocated`] occurs, and there is no IO error
/// found which could've caused this, then an error with a custom message is
/// thrown.
#[inline]
#[cfg(not(feature = "dev_no_rayon"))]
fn align_queries<F>(query_reader: QueryReader, writer: AlignmentWriterThreaded, f: F) -> std::io::Result<()>
where
    F: Fn(&mut AlignmentWriterThreaded, std::io::Result<FastX>) -> Result<(), ThreadedWriteError> + Sync + Send, {
    let res = query_reader
        .par_bridge()
        .try_for_each_with(writer.clone(), |w, record| f(w, record));

    match res {
        Ok(()) => writer.flush(),
        Err(ThreadedWriteError::IoError(e)) => Err(e),
        Err(ThreadedWriteError::ReceiverDeallocated) => Err(writer.flush().err().unwrap_or(std::io::Error::other(
            "The receiver in the writing thread unexpectedly closed",
        ))),
    }
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
    pub fn new(
        query: &'q FastX, matrix: &'q WeightMatrix<'q, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> std::io::Result<Self> {
        let forward = query;

        let forward_seq = forward.sequence.as_slice();
        let header = &forward.header;
        let profile = LocalProfiles::make_profile(forward_seq, header, matrix, gap_open, gap_extend)?;

        Ok(Self { forward, profile })
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
    pub fn sw_1pass_query_profile<'r>(&'q self, reference: &Reference<'r, S>) -> std::io::Result<AlignmentAndSeqs<'q, 'r>> {
        let mapping = align_maybe_rc(SeqSrc::Reference(&reference.forward.sequence), &reference.reverse, |seq| {
            self.profile.sw_1pass(seq)
        })
        .with_context(format!(
            "Failed to align the sequences with the following headers:\n    | Query: {q_header}\n    | Reference: {r_header}",
            q_header=self.forward.header, r_header=reference.forward.name
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
    pub fn sw_3pass_query_profile<'r>(&'q self, reference: &Reference<'r, S>) -> std::io::Result<AlignmentAndSeqs<'q, 'r>> {
        let mapping = align_maybe_rc(SeqSrc::Reference(&reference.forward.sequence), &reference.reverse, |seq| {
            self.profile.sw_3pass(seq)
        }).with_context(format!(
            "Failed to align the sequences with the following headers:\n    | Query: {q_header}\n    | Reference: {r_header}",
            q_header=self.forward.header, r_header=reference.forward.name
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
}

/// A reference record, together with its reverse complement and a profile set
/// for alignment.
#[derive(Clone, Debug)]
pub struct Reference<'r, const S: usize> {
    /// The record for the reference.
    forward: &'r FastaSeq,
    /// The reverse complement of the reference, if `--rev-comp` was passed.
    reverse: MaybeRevComp,
    /// The profile set for the forward reference sequence.
    profile: SharedProfiles<'r, 32, 16, 8, S>,
}

impl<'r, const S: usize> Reference<'r, S> {
    /// Bundles a reference together with its reverse complement (if `rev_comp`
    /// is true) and a profile set for use in alignment.
    ///
    /// ## Errors
    ///
    /// An error containing the reference header as context is returned if a
    /// profile fails to be made from the sequence.
    pub fn new(
        reference: &'r FastaSeq, matrix: &'r WeightMatrix<'r, i8, S>, gap_open: i8, gap_extend: i8, rev_comp: bool,
    ) -> std::io::Result<Self> {
        let forward = reference;

        let forward_seq = forward.sequence.as_slice();
        let reverse = MaybeRevComp::new(forward_seq, rev_comp);
        let header = &forward.name;
        let profile = SharedProfiles::make_profile(forward_seq, header, matrix, gap_open, gap_extend)?;

        Ok(Self {
            forward,
            reverse,
            profile,
        })
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
        }).with_context(format!(
            "Failed to align the sequences with the following headers:\n    | Query: {q_header}\n    | Reference: {r_header}",
            q_header=query.forward.header, r_header=self.forward.name
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
        }).with_context(format!(
            "Failed to align the sequences with the following headers:\n    | Query: {q_header}\n    | Reference: {r_header}",
            q_header=query.forward.header, r_header=self.forward.name
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
    ///
    /// ## Errors
    ///
    /// If any of the profiles fail to build, an error is returned with context
    /// including the header.
    pub fn new(
        references: &'r [FastaSeq], matrix: &'r WeightMatrix<'r, i8, S>, gap_open: i8, gap_extend: i8, rev_comp: bool,
    ) -> std::io::Result<Self> {
        references
            .iter()
            .map(|reference| Reference::new(reference, matrix, gap_open, gap_extend, rev_comp))
            .collect::<Result<_, _>>()
            .map(References)
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

/// An enum containing the different alignment methods which can be used in
/// [`align_all`] or [`align_best_match`].
#[allow(clippy::enum_variant_names)]
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
enum AlignmentMethod {
    OnePassQueryProfile,
    OnePassRefProfile,
    ThreePassQueryProfile,
    ThreePassRefProfile,
}

/// A trait extending [`ProfileSets`] with methods that add context to errors.
pub trait AlignerMethods<'a, const S: usize>: ProfileSets<'a, 32, 16, 8, S> {
    /// Constructs a new profile, similar to [`LocalProfiles::new`] or
    /// [`SharedProfiles::new`].
    ///
    /// ## Errors
    ///
    /// Context is added that includes the header of the sequence if it fails to
    /// be preprocessed.
    fn make_profile(
        seq: &'a [u8], header: &str, matrix: &'a WeightMatrix<'a, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> std::io::Result<Self>;

    /// Performs the 1-pass Striped Smith Waterman alignment algorithm starting
    /// from an `i8`, similar to [`ProfileSets::sw_align_from_i8`].
    ///
    /// ## Errors
    ///
    /// If the alignment overflows, then this is returned as an error.
    fn sw_1pass(&self, seq: SeqSrc<&[u8]>) -> std::io::Result<Option<Alignment<u32>>> {
        let seq = seq.map(AsRef::as_ref);
        let alignment = self.sw_align_from_i8(seq);
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

impl<'a, const S: usize> AlignerMethods<'a, S> for LocalProfiles<'a, 32, 16, 8, S> {
    fn make_profile(
        seq: &'a [u8], header: &str, matrix: &'a WeightMatrix<'a, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> std::io::Result<Self> {
        Ok(Self::new(seq, matrix, gap_open, gap_extend).with_context(format!(
            "The sequence with the following header failed to be preprocessed: {header}"
        ))?)
    }
}

impl<'a, const S: usize> AlignerMethods<'a, S> for SharedProfiles<'a, 32, 16, 8, S> {
    fn make_profile(
        seq: &'a [u8], header: &str, matrix: &'a WeightMatrix<'a, i8, S>, gap_open: i8, gap_extend: i8,
    ) -> std::io::Result<Self> {
        Ok(Self::new(seq, matrix, gap_open, gap_extend).with_context(format!(
            "The sequence with the following header failed to be preprocessed: {header}"
        ))?)
    }
}

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

    let alignment = if let Some(seq_rc) = &seq_rc.0 {
        let seq_rc = seq.map(|_| seq_rc.as_slice());
        let alignment_rc = f(seq_rc).with_context("Failed to perform the reverse complement alignment")?;

        if alignment_rc > alignment_forward {
            alignment_rc.map(|mut alignment| {
                if seq.is_reference() {
                    alignment.make_reverse();
                }
                AlignmentAndStrand {
                    inner:  alignment,
                    strand: Strand::Reverse,
                }
            })
        } else {
            alignment_forward.map(|alignment| AlignmentAndStrand {
                inner:  alignment,
                strand: Strand::Forward,
            })
        }
    } else {
        alignment_forward.map(|alignment| AlignmentAndStrand {
            inner:  alignment,
            strand: Strand::Forward,
        })
    };

    Ok(alignment)
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
#[derive(PartialEq)]
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
    pub reference: &'r FastaSeq,
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

/// The unpacked tallies for the number of alignments satisfying and failing a
/// particular condition.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
struct UnpackedTallies {
    satisfying: u64,
    failing:    u64,
}

/// A tally for the number of alignments satisfying a condition and failing a
/// condition.
///
/// To prevent skew in the update of the tallies, the two are packed into a
/// single atomic (which is guaranteed to have an ordering which is consistent
/// with the program's order). The lower 32 bits are the tally for the number
/// satisfying the condition, and the upper 32 bits are for the tally failing
/// the condition.
#[repr(transparent)]
#[derive(Debug, Default)]
struct AlignmentTally(AtomicU64);

impl AlignmentTally {
    /// Unpacks the two counters contained in the [`AtomicU64`] for easier
    /// handling and computation.
    fn unpack(&self) -> UnpackedTallies {
        let bits = self.0.load(std::sync::atomic::Ordering::Relaxed);
        let (satisfying, failing) = (bits & 0xFFFF_FFFF, bits >> 32);
        UnpackedTallies { satisfying, failing }
    }

    /// Increments one of the tallies, depending on whether `condition` is
    /// `true` (satisfying) or `false` (failing).
    fn tally(&self, condition: bool) {
        let increment = if condition { 1 } else { 1u64 << 32 };
        self.0.fetch_add(increment, std::sync::atomic::Ordering::Relaxed);
    }

    /// Returns whether the specified percent (in `0..=100`) of alignments
    /// satisfy the given condition. If fewer than `min_total` alignments are
    /// tallied, then `false` is returned by default.
    fn percent_meets_condition(&self, percent: u64, min_total: u64) -> bool {
        assert!(percent <= 100);

        let UnpackedTallies { satisfying, failing } = self.unpack();
        let total = satisfying + failing;

        total >= min_total && satisfying * 100 >= total * percent
    }
}

/// A collection of non-blocking thread-safe tallies for use in adaptively
/// determining the proper alignment algorithm to use.
///
/// The tallies use [`AlignmentTally`] ([`AtomicU64`] with two tallies packed
/// into each atomic) to prevent possible skew in the update order.
#[derive(Debug)]
struct AlignmentTallies {
    /// Tallies for the number of alignments fitting within the capacity of an
    /// `i8`.
    scores_fitting_i8:   AlignmentTally,
    /// Tallies for the number of queries of length at most 300.
    queries_at_most_300: AlignmentTally,
    /// The length of the first reference sequence.
    first_ref_len:       usize,
}

impl AlignmentTallies {
    /// Initializes new tallies with everything starting at 0 (and stores the
    /// length of the first reference).
    fn new<const S: usize>(references: &References<'_, S>) -> Self {
        Self {
            scores_fitting_i8:   AlignmentTally::default(),
            queries_at_most_300: AlignmentTally::default(),
            first_ref_len:       references.0.first().map_or(0, |reference| reference.forward.sequence.len()),
        }
    }

    /// Tallies the characteristics about an alignment.
    fn tally_alignment<const S: usize>(&self, alignment: &AlignmentAndSeqs, query: &FastX, matrix: &WeightMatrix<i8, S>) {
        if let Some(mapping) = &alignment.mapping {
            self.scores_fitting_i8
                .tally(mapping.inner.score <= max_score_for_int_type::<i8, i8, S>(matrix));
        }

        self.queries_at_most_300.tally(query.sequence.len() <= 300);
    }

    /// Picks an alignment method to use based on the tallies.
    fn pick_alignment_method(&self, config: &AlignerConfig) -> AlignmentMethod {
        let num_passes = config.method.unwrap_or_else(|| {
            if self.scores_fitting_i8.percent_meets_condition(66, 50) {
                NumPasses::OnePass
            } else {
                NumPasses::ThreePass
            }
        });

        let profile_from = config.profile_from.unwrap_or(match num_passes {
            NumPasses::OnePass => {
                if self.first_ref_len < 600 {
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
}
