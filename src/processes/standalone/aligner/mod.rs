use crate::{
    aligner::{
        arg_parsing::{AlignerConfig, Alphabet, AnyMatrix, ParsedAlignerArgs, parse_aligner_args},
        methods::{AlignmentMethod, StripedSmithWatermanLocal, StripedSmithWatermanShared},
        writers::send_alignment,
    },
    io::FastXReader,
};
use clap::{Args, builder::RangedI64ValueParser};
use rayon::iter::{ParallelBridge, ParallelIterator};
use std::{
    borrow::Borrow,
    io::{Read, Write},
    path::PathBuf,
    sync::mpsc,
    thread,
};
use zoe::{
    data::{fasta::FastaSeq, matrices::WeightMatrix},
    math::AnyInt,
    prelude::NucleotidesView,
};

mod arg_parsing;
mod methods;
mod writers;

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
}

/// Sub-program for performing sequence alignment
pub fn aligner_process(args: AlignerArgs) -> std::io::Result<()> {
    let ParsedAlignerArgs {
        query_reader,
        references,
        mut output,
        weight_matrix,
        config,
    } = parse_aligner_args(args)?;

    if config.single_thread {
        rayon::ThreadPoolBuilder::new().num_threads(1).build_global().unwrap();
    }

    let (sender, receiver) = mpsc::channel();
    let writer_thread = thread::spawn(move || -> std::io::Result<()> {
        while let Ok(string) = receiver.recv() {
            writeln!(output, "{string}")?;
        }
        Ok(())
    });

    dispatch_alphabet(query_reader, references, sender, weight_matrix, config)?;

    writer_thread.join().unwrap()
}

/// Dispatches the aligner based on the alphabet and weight matrix
fn dispatch_alphabet<R: Read + Send>(
    query_reader: FastXReader<R>, references: Vec<FastaSeq>, sender: mpsc::Sender<String>,
    weight_matrix: AnyMatrix<'static, i8>, config: AlignerConfig,
) -> std::io::Result<()> {
    match weight_matrix {
        AnyMatrix::Dna(weight_matrix) => dispatch_runner(query_reader, references, sender, weight_matrix, config),
        AnyMatrix::AaNamed(weight_matrix) => dispatch_runner(query_reader, references, sender, weight_matrix, config),
        AnyMatrix::AaSimple(weight_matrix) => dispatch_runner(query_reader, references, sender, weight_matrix, config),
    }
}

/// Dispatches the aligner based on profile_from and best_match
fn dispatch_runner<R, W, const S: usize>(
    query_reader: FastXReader<R>, references: Vec<FastaSeq>, sender: mpsc::Sender<String>, weight_matrix: W,
    config: AlignerConfig,
) -> std::io::Result<()>
where
    R: Read + Send,
    W: Borrow<WeightMatrix<'static, i8, S>> + Sync + Send + 'static, {
    match (config.profile_from_ref, config.best_match) {
        (false, false) => {
            let method = StripedSmithWatermanLocal::<W, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_all_profile_from_query(method, query_reader, references, sender, config)
        }
        (true, false) => {
            let method = StripedSmithWatermanShared::<W, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_all_profile_from_ref(method, query_reader, references, sender, config)
        }
        (false, true) => {
            let method = StripedSmithWatermanLocal::<W, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_best_match_profile_from_query(method, query_reader, references, sender, config)
        }
        (true, true) => {
            let method = StripedSmithWatermanShared::<W, S>::new(weight_matrix, config.gap_open, config.gap_extend);
            align_best_match_profile_from_ref(method, query_reader, references, sender, config)
        }
    }
}

/// Performs all alignments using the provided `method`, with profiles built
/// from the queries.
///
/// The results are sent to `sender`.
fn align_all_profile_from_query<R, A>(
    method: A, query_reader: FastXReader<R>, references: Vec<FastaSeq>, sender: mpsc::Sender<String>, config: AlignerConfig,
) -> std::io::Result<()>
where
    R: Read + Send,
    A: AlignmentMethod, {
    // Compute all the reverse complements of the references ahead of time, if
    // rev_comp is true
    let references = method.maybe_zip_with_revcomp(&references, config.rev_comp);

    query_reader
        .par_bridge()
        .try_for_each_with(sender, |sender, query| -> std::io::Result<()> {
            let query = query?;
            let profile = method.build_profile(&query)?;

            for (reference, rc_reference) in &references {
                let mut alignment = method.align(&profile, &reference.sequence, rc_reference.as_ref())?;
                if let Some((alignment, Strand::Reverse)) = alignment.as_mut() {
                    // We have the reverse complement of the reference, but we
                    // want the alignment to correspond to the reverse
                    // complement of the query
                    alignment.make_reverse();
                }
                send_alignment(sender, alignment, &query, reference, &config);
            }

            Ok(())
        })
}

/// Performs all alignments using the provided `method`, with profiles built
/// from the references.
///
/// The results are sent to `sender`.
fn align_all_profile_from_ref<R, A>(
    method: A, query_reader: FastXReader<R>, references: Vec<FastaSeq>, sender: mpsc::Sender<String>, config: AlignerConfig,
) -> std::io::Result<()>
where
    R: Read + Send,
    A: for<'a> AlignmentMethod<Profile<'a>: Sync>, {
    // Build profiles first, which requires SharedProfile
    let ref_profiles = method.zip_with_profiles(&references)?;

    query_reader
        .par_bridge()
        .try_for_each_with(sender, |sender, query| -> std::io::Result<()> {
            let query = query?;

            // Compute the reverse complement of the query, if rev_comp is true
            let rc_query = config.rev_comp.then(|| {
                NucleotidesView::from(query.sequence.as_slice())
                    .to_reverse_complement()
                    .into_vec()
            });

            for (reference, ref_profile) in &ref_profiles {
                // Invert the alignment, since we built the profile from the
                // reference
                let alignment = method
                    .align(ref_profile, &query.sequence, rc_query.as_ref())?
                    .map(|(alignment, strand)| (alignment.invert(), strand));

                send_alignment(sender, alignment, &query, reference, &config);
            }

            Ok(())
        })
}

/// Performs all alignments between the streamed queries and the slurped
/// references using the provided alignment method, but only keeping one
/// alignment with the best score for each query.
///
/// `profile_seqs` corresponds to the queries, and `regular_seqs` corresponds to
/// the references. The results are sent to `sender`.
fn align_best_match_profile_from_query<R, A>(
    method: A, query_reader: FastXReader<R>, references: Vec<FastaSeq>, sender: mpsc::Sender<String>, config: AlignerConfig,
) -> std::io::Result<()>
where
    R: Read + Send,
    A: AlignmentMethod, {
    // Compute all the reverse complements of the references ahead of time, if
    // rev_comp is true
    let references = method.maybe_zip_with_revcomp(&references, config.rev_comp);

    query_reader
        .par_bridge()
        .try_for_each_with(sender, |sender, query| -> std::io::Result<()> {
            let query = query?;
            let profile = method.build_profile(&query)?;

            let (best_reference, best_alignment) = references
                .iter()
                .map(|(reference, rc_reference)| {
                    let alignment = method.align(&profile, &reference.sequence, rc_reference.as_ref())?;
                    std::io::Result::Ok((reference, alignment))
                })
                .max_by_key(|res| match res {
                    Ok((_reference, alignment)) => alignment
                        .as_ref()
                        .map_or(A::Score::ZERO, |(alignment, _strand)| alignment.score),
                    // Map errors to maximum value, so they are guaranteed to propagate
                    Err(_) => A::Score::MAX,
                })
                .ok_or(std::io::Error::other("No references were specified!"))??;

            send_alignment(sender, best_alignment, &query, best_reference, &config);
            Ok(())
        })
}

/// Performs all alignments between the streamed queries and the slurped
/// references using the provided alignment method, but only keeping one
/// alignment with the best score for each query.
///
/// `profile_seqs` corresponds to the references, and `regular_seqs` corresponds
/// to the queries. The results are sent to `sender`.
fn align_best_match_profile_from_ref<R, A>(
    method: A, query_reader: FastXReader<R>, references: Vec<FastaSeq>, sender: mpsc::Sender<String>, config: AlignerConfig,
) -> std::io::Result<()>
where
    R: Read + Send,
    A: for<'a> AlignmentMethod<Profile<'a>: Sync>, {
    // Build profiles first, which requires `SharedProfile`
    let ref_profiles = method.zip_with_profiles(&references)?;

    query_reader
        .par_bridge()
        .try_for_each_with(sender, |sender, query| -> std::io::Result<()> {
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
                    let alignment = method.align(ref_profile, &query.sequence, rc_query.as_ref())?;
                    std::io::Result::Ok((reference, alignment))
                })
                .max_by_key(|res| match res {
                    Ok((_reference, alignment)) => alignment
                        .as_ref()
                        .map_or(A::Score::ZERO, |(alignment, _strand)| alignment.score),
                    // Map errors to maximum value, so they are guaranteed to propagate
                    Err(_) => A::Score::MAX,
                })
                .ok_or(std::io::Error::other("No references were specified!"))??;

            // Invert the alignment, since we built profile from reference
            let best_alignment = best_alignment.map(|(alignment, strand)| (alignment.invert(), strand));

            send_alignment(sender, best_alignment, &query, best_reference, &config);
            Ok(())
        })
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
