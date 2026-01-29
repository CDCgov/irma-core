use crate::{
    aligner::{
        arg_parsing::{AlignerConfig, Alphabet, AnyMatrix, ParsedAlignerArgs, parse_aligner_args},
        methods::{AlignmentMethod, StripedSmithWatermanLocal, StripedSmithWatermanShared},
        writers::{AlignmentWriter, write_header},
    },
    io::{FastX, FastXReader, IterWithContext, OutputOptions, ReadFileZipPipe},
};
use clap::{Args, builder::RangedI64ValueParser};
use std::{borrow::Borrow, path::PathBuf};
use zoe::{
    data::{err::ResultWithErrorContext, fasta::FastaSeq, matrices::WeightMatrix},
    math::AnyInt,
    prelude::{NucleotidesView, SeqSrc},
};

#[cfg(not(feature = "dev_no_rayon"))]
use crate::aligner::writers::{AlignmentWriterThreaded, ThreadedWriteError};
#[cfg(not(feature = "dev_no_rayon"))]
use rayon::iter::{ParallelBridge, ParallelIterator};

#[cfg(feature = "dev_no_rayon")]
use std::io::Write;

mod arg_parsing;
mod methods;
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

    #[arg(long)]
    /// Builds the profile from the reference sequences instead of the queries
    profile_from_ref: bool,

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

/// Sub-program for performing sequence alignment
pub fn aligner_process(args: AlignerArgs) -> std::io::Result<()> {
    let ParsedAlignerArgs {
        query_reader,
        references,
        weight_matrix,
        header,
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
        let res = write_header(&mut writer, &references);

        if let Some(output) = &config.output {
            res.with_file_context("Failed to write SAM output to file", output)?;
        } else {
            res.with_context("Failed to write SAM output")?;
        }
    }

    #[cfg(not(feature = "dev_no_rayon"))]
    let writer = AlignmentWriterThreaded::from_writer(writer);

    // Validity: No context is added to the result
    dispatch_alphabet(query_reader, references, writer, weight_matrix, &config)
}

/// Dispatches the aligner based on the alphabet and weight matrix.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. See the documentation
/// for the particular [`AlignmentMethod`] for more details. Context containing
/// the header(s) is added for failed profile building or alignment.
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
) -> std::io::Result<()> {
    // Validity: No context is added to the results
    match weight_matrix {
        AnyMatrix::Dna(weight_matrix) => dispatch_runner(query_reader, references, writer, weight_matrix, config),
        AnyMatrix::AaNamed(weight_matrix) => dispatch_runner(query_reader, references, writer, weight_matrix, config),
        AnyMatrix::AaSimple(weight_matrix) => dispatch_runner(query_reader, references, writer, weight_matrix, config),
    }
}

/// Dispatches the aligner based on whether the profile is built from the query
/// or reference, and whether all alignments or just the best match alignment
/// are being reported.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. See the documentation
/// for the particular [`AlignmentMethod`] for more details. Context containing
/// the header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn dispatch_runner<M, const S: usize>(
    query_reader: QueryReader, references: Vec<FastaSeq>, writer: SamWriter, weight_matrix: M, config: &AlignerConfig,
) -> std::io::Result<()>
where
    M: Borrow<WeightMatrix<'static, i8, S>> + Sync + Send + 'static, {
    match (config.profile_from_ref, config.best_match) {
        (false, false) => {
            let method = StripedSmithWatermanLocal::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            // Validity: No context is added to this result
            align_all_profile_from_query(method, query_reader, references, writer, config)
        }
        (true, false) => {
            let method = StripedSmithWatermanShared::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            // Validity: No context is added to this result
            align_all_profile_from_ref(method, query_reader, references, writer, config)
        }
        (false, true) => {
            let method = StripedSmithWatermanLocal::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            // Validity: No context is added to this result
            align_best_match_profile_from_query(method, query_reader, references, writer, config)
        }
        (true, true) => {
            let method = StripedSmithWatermanShared::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            // Validity: No context is added to this result
            align_best_match_profile_from_ref(method, query_reader, references, writer, config)
        }
    }
}

/// Performs all alignments using the provided `method`, with profiles built
/// from the queries.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. See the documentation
/// for the particular [`AlignmentMethod`] for more details. Context containing
/// the header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn align_all_profile_from_query<A>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: SamWriter, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: AlignmentMethod, {
    // Compute all the reverse complements of the references ahead of time, if
    // rev_comp is true
    let references = method.maybe_zip_with_revcomp(&references, config.rev_comp);

    align_all(query_reader, writer, |writer, query| {
        let query = query?;
        let profile = method.build_profile(&query)?;

        for (reference, rc_reference) in &references {
            // Validity: This error message formatting requires that no callers
            // add additional context
            let mut alignment = method
                .align(&profile, SeqSrc::Reference(&reference.sequence), rc_reference.as_ref())
                .with_context(format!(
                    "Failed to align the sequences with the following headers:\n    | Query: {}\n    | Reference: {}",
                    query.header, reference.name
                ))?;
            if let Some((alignment, Strand::Reverse)) = alignment.as_mut() {
                // We have the reverse complement of the reference, but we want
                // the alignment to correspond to the reverse complement of the
                // query
                alignment.make_reverse();
            }

            writer.write_alignment(alignment, &query, reference, config)?;
        }

        Ok(())
    })
}

/// Performs all alignments using the provided `method`, with profiles built
/// from the references.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. See the documentation
/// for the particular [`AlignmentMethod`] for more details. Context containing
/// the header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn align_all_profile_from_ref<A>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: SamWriter, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: for<'a> AlignmentMethod<Profile<'a>: Sync>, {
    // Build profiles first, which requires SharedProfile
    let ref_profiles = method.zip_with_profiles(&references)?;

    align_all(query_reader, writer, |writer, query| {
        let query = query?;

        // Compute the reverse complement of the query, if rev_comp is true
        let rc_query = config.rev_comp.then(|| {
            NucleotidesView::from(query.sequence.as_slice())
                .to_reverse_complement()
                .into_vec()
        });

        for (reference, ref_profile) in &ref_profiles {
            // Validity: This error message formatting requires that no callers
            // add additional context
            let alignment = method
                .align(ref_profile, SeqSrc::Query(&query.sequence), rc_query.as_ref())
                .with_context(format!(
                    "Failed to align the sequences with the following headers:\n    | Query: {}\n    | Reference: {}",
                    query.header, reference.name
                ))?;

            writer.write_alignment(alignment, &query, reference, config)?;
        }

        Ok(())
    })
}

/// Performs all alignments between the streamed queries and the slurped
/// references using the provided alignment method, but only keeping one
/// alignment with the best score for each query.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. See the documentation
/// for the particular [`AlignmentMethod`] for more details. Context containing
/// the header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn align_best_match_profile_from_query<A>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: SamWriter, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: AlignmentMethod, {
    // Compute all the reverse complements of the references ahead of time, if
    // rev_comp is true
    let references = method.maybe_zip_with_revcomp(&references, config.rev_comp);

    align_all(query_reader, writer, |writer, query| {
        let query = query?;
        let profile = method.build_profile(&query)?;

        let (best_reference, mut best_alignment) = references
            .iter()
            .map(|(reference, rc_reference)| {
                // Validity: This error message formatting requires that no
                // callers add additional context
                let alignment = method
                    .align(&profile, SeqSrc::Reference(&reference.sequence), rc_reference.as_ref())
                    .with_context(format!(
                        "Failed to align the sequences with the following headers:\n    | Query: {}\n    | Reference: {}",
                        query.header, reference.name
                    ))?;
                std::io::Result::Ok((reference, alignment))
            })
            .max_by_key(|res| match res {
                Ok((_reference, alignment)) => alignment
                    .as_ref()
                    .map_or(A::Score::ZERO, |(alignment, _strand)| alignment.score),
                // Map errors to maximum value, so they are guaranteed to propagate
                Err(_) => A::Score::MAX,
            })
            .expect("The references field should be non-empty")?;

        if let Some((best_alignment, Strand::Reverse)) = best_alignment.as_mut() {
            // We have the reverse complement of the reference, but we want
            // the alignment to correspond to the reverse complement of the
            // query
            best_alignment.make_reverse();
        }

        writer.write_alignment(best_alignment, &query, best_reference, config)
    })
}

/// Performs all alignments between the streamed queries and the slurped
/// references using the provided alignment method, but only keeping one
/// alignment with the best score for each query.
///
/// ## Errors
///
/// Errors while reading the queries, building the profiles, performing the
/// alignment, and writing the alignment are propagated. See the documentation
/// for the particular [`AlignmentMethod`] for more details. Context containing
/// the header(s) is added for failed profile building or alignment.
///
/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn align_best_match_profile_from_ref<A>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: SamWriter, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: for<'a> AlignmentMethod<Profile<'a>: Sync>, {
    // Build profiles first, which requires `SharedProfile`
    let ref_profiles = method.zip_with_profiles(&references)?;

    align_all(query_reader, writer, |writer, query| {
        let query = query?;

        // Compute the reverse complement of the query, if rev_comp is true
        let rc_query = config.rev_comp.then(|| {
            NucleotidesView::from(query.sequence.as_slice())
                .to_reverse_complement()
                .into_vec()
        });

        let (best_reference, best_alignment) = ref_profiles
            .iter()
            .map(|(reference, ref_profile)| {
                // Validity: This error message formatting requires that no
                // callers add additional context
                let alignment = method
                    .align(ref_profile, SeqSrc::Query(&query.sequence), rc_query.as_ref())
                    .with_context(format!(
                        "Failed to align the sequences with the following headers:\n    | Query: {}\n    | Reference: {}",
                        query.header, reference.name
                    ))?;
                std::io::Result::Ok((reference, alignment))
            })
            .max_by_key(|res| match res {
                Ok((_reference, alignment)) => alignment
                    .as_ref()
                    .map_or(A::Score::ZERO, |(alignment, _strand)| alignment.score),
                // Map errors to maximum value, so they are guaranteed to propagate
                Err(_) => A::Score::MAX,
            })
            .expect("The references field should be non-empty")?;

        writer.write_alignment(best_alignment, &query, best_reference, config)
    })
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
fn align_all<F>(query_reader: QueryReader, mut writer: SamWriter, f: F) -> std::io::Result<()>
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
fn align_all<F>(query_reader: QueryReader, writer: AlignmentWriterThreaded, f: F) -> std::io::Result<()>
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
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum Strand {
    /// The alignment was with the forward sequence
    Forward,
    /// The alignment was with the reverse complement
    Reverse,
}

impl Default for Strand {
    #[inline]
    fn default() -> Self {
        Self::Forward
    }
}
