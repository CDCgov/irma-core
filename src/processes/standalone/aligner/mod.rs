#[cfg(not(feature = "dev_no_rayon"))]
use crate::aligner::writers::AlignmentWriterThreaded;
use crate::{
    aligner::{
        arg_parsing::{AlignerConfig, Alphabet, AnyMatrix, ParsedAlignerArgs, parse_aligner_args},
        methods::{AlignmentMethod, StripedSmithWatermanLocal, StripedSmithWatermanShared},
        writers::{AdditionalBounds, AlignmentWriter, write_header},
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
use rayon::iter::{ParallelBridge, ParallelIterator};
#[cfg(feature = "dev_no_rayon")]
use std::io::Write;

mod arg_parsing;
mod methods;
mod writers;

/// A type alias for the query reader used by `aligner`.
type QueryReader = IterWithContext<FastXReader<ReadFileZipPipe>>;

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
    let mut writer = AlignmentWriterThreaded::from_writer(writer);

    dispatch_alphabet(query_reader, references, &mut writer, weight_matrix, &config)?;

    // Must be called to either flush the bufwriter, or properly terminate the
    // thread
    let res = writer.flush();

    if let Some(output) = config.output {
        res.with_file_context("Failed to write SAM output to file", output)?;
    } else {
        res.with_context("Failed to write SAM output")?;
    }

    Ok(())
}

/// Dispatches the aligner based on the alphabet and weight matrix
fn dispatch_alphabet<W>(
    query_reader: QueryReader, references: Vec<FastaSeq>, writer: &mut W, weight_matrix: AnyMatrix<'static, i8>,
    config: &AlignerConfig,
) -> std::io::Result<()>
where
    W: AlignmentWriter + AdditionalBounds, {
    match weight_matrix {
        AnyMatrix::Dna(weight_matrix) => dispatch_runner(query_reader, references, writer, weight_matrix, config),
        AnyMatrix::AaNamed(weight_matrix) => dispatch_runner(query_reader, references, writer, weight_matrix, config),
        AnyMatrix::AaSimple(weight_matrix) => dispatch_runner(query_reader, references, writer, weight_matrix, config),
    }
}

/// Dispatches the aligner based on profile_from and best_match
fn dispatch_runner<W, M, const S: usize>(
    query_reader: QueryReader, references: Vec<FastaSeq>, writer: &mut W, weight_matrix: M, config: &AlignerConfig,
) -> std::io::Result<()>
where
    W: AlignmentWriter + AdditionalBounds,
    M: Borrow<WeightMatrix<'static, i8, S>> + Sync + Send + 'static, {
    match (config.profile_from_ref, config.best_match) {
        (false, false) => {
            let method = StripedSmithWatermanLocal::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_all_profile_from_query(method, query_reader, references, writer, config)
        }
        (true, false) => {
            let method = StripedSmithWatermanShared::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_all_profile_from_ref(method, query_reader, references, writer, config)
        }
        (false, true) => {
            let method = StripedSmithWatermanLocal::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_best_match_profile_from_query(method, query_reader, references, writer, config)
        }
        (true, true) => {
            let method = StripedSmithWatermanShared::<M, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_best_match_profile_from_ref(method, query_reader, references, writer, config)
        }
    }
}

/// Performs all alignments using the provided `method`, with profiles built
/// from the queries.
fn align_all_profile_from_query<A, W>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: &mut W, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: AlignmentMethod,
    W: AlignmentWriter + AdditionalBounds, {
    // Compute all the reverse complements of the references ahead of time, if
    // rev_comp is true
    let references = method.maybe_zip_with_revcomp(&references, config.rev_comp);

    align_all(query_reader, writer, |writer, query| -> std::io::Result<()> {
        let query = query?;
        let profile = method.build_profile(&query)?;

        for (reference, rc_reference) in &references {
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

            let res = writer.write_alignment(alignment, &query, reference, config);

            if let Some(output) = &config.output {
                res.with_file_context("Failed to write SAM output to file", output)?;
            } else {
                res.with_context("Failed to write SAM output")?;
            }
        }

        Ok(())
    })
}

/// Performs all alignments using the provided `method`, with profiles built
/// from the references.
fn align_all_profile_from_ref<A, W>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: &mut W, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: for<'a> AlignmentMethod<Profile<'a>: Sync>,
    W: AlignmentWriter + AdditionalBounds, {
    // Build profiles first, which requires SharedProfile
    let ref_profiles = method.zip_with_profiles(&references)?;

    align_all(query_reader, writer, |writer, query| -> std::io::Result<()> {
        let query = query?;

        // Compute the reverse complement of the query, if rev_comp is true
        let rc_query = config.rev_comp.then(|| {
            NucleotidesView::from(query.sequence.as_slice())
                .to_reverse_complement()
                .into_vec()
        });

        for (reference, ref_profile) in &ref_profiles {
            let alignment = method
                .align(ref_profile, SeqSrc::Query(&query.sequence), rc_query.as_ref())
                .with_context(format!(
                    "Failed to align the sequences with the following headers:\n    | Query: {}\n    | Reference: {}",
                    query.header, reference.name
                ))?;

            let res = writer.write_alignment(alignment, &query, reference, config);

            if let Some(output) = &config.output {
                res.with_file_context("Failed to write SAM output to file", output)?;
            } else {
                res.with_context("Failed to write SAM output")?;
            }
        }

        Ok(())
    })
}

/// Performs all alignments between the streamed queries and the slurped
/// references using the provided alignment method, but only keeping one
/// alignment with the best score for each query.
fn align_best_match_profile_from_query<A, W>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: &mut W, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: AlignmentMethod,
    W: AlignmentWriter + AdditionalBounds, {
    // Compute all the reverse complements of the references ahead of time, if
    // rev_comp is true
    let references = method.maybe_zip_with_revcomp(&references, config.rev_comp);

    align_all(query_reader, writer, |writer, query| -> std::io::Result<()> {
        let query = query?;
        let profile = method.build_profile(&query)?;

        let (best_reference, mut best_alignment) = references
            .iter()
            .map(|(reference, rc_reference)| {
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

        let res = writer.write_alignment(best_alignment, &query, best_reference, config);

        if let Some(output) = &config.output {
            Ok(res.with_file_context("Failed to write SAM output to file", output)?)
        } else {
            Ok(res.with_context("Failed to write SAM output")?)
        }
    })
}

/// Performs all alignments between the streamed queries and the slurped
/// references using the provided alignment method, but only keeping one
/// alignment with the best score for each query.
fn align_best_match_profile_from_ref<A, W>(
    method: A, query_reader: QueryReader, references: Vec<FastaSeq>, writer: &mut W, config: &AlignerConfig,
) -> std::io::Result<()>
where
    A: for<'a> AlignmentMethod<Profile<'a>: Sync>,
    W: AlignmentWriter + AdditionalBounds, {
    // Build profiles first, which requires `SharedProfile`
    let ref_profiles = method.zip_with_profiles(&references)?;

    align_all(query_reader, writer, |writer, query| -> std::io::Result<()> {
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

        let res = writer.write_alignment(best_alignment, &query, best_reference, config);

        if let Some(output) = &config.output {
            Ok(res.with_file_context("Failed to write SAM output to file", output)?)
        } else {
            Ok(res.with_context("Failed to write SAM output")?)
        }
    })
}

/// Performs all alignments as indicated by closure `f`, using either a parallel
/// iterator (`par_bridge`) or a serial iterator depending on the `one-thread`
/// feature.
#[inline]
fn align_all<W, F>(query_reader: QueryReader, writer: &mut W, f: F) -> std::io::Result<()>
where
    W: AlignmentWriter + AdditionalBounds,
    F: Fn(&mut W, std::io::Result<FastX>) -> std::io::Result<()> + Sync + Send, {
    #[cfg(not(feature = "dev_no_rayon"))]
    query_reader.par_bridge().try_for_each_with(writer.clone(), f)?;

    #[cfg(feature = "dev_no_rayon")]
    {
        let mut query_reader = query_reader;
        query_reader.try_for_each(|query| f(writer, query))?;
    }

    Ok(())
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
