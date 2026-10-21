use crate::{
    aligner::{AlignerArgs, QueryReader},
    args::abort_clap,
};
use clap::{ValueEnum, error::ErrorKind};
use irma_records::io::InputOptions;
use std::{fmt::Display, path::PathBuf};
use zoe::{
    alignment::{ProfileError, ScoreParams},
    data::{
        AA_ALL_AMBIG_PROFILE_MAP_WITH_STOP, WeightMatrix,
        fasta::Fasta,
        matrices::{BLOSUM_62, aa_mat_from_name},
    },
};

/// The parsed and validated command line arguments for `aligner`
pub struct ParsedAlignerArgs {
    /// The streamed query sequences
    pub query_reader: QueryReader,
    /// The slurped reference sequences
    ///
    /// ## Validity
    ///
    /// This field must be non-empty.
    pub references:   Vec<Fasta>,
    /// The parameters for scoring the alignment
    pub params:       EitherParams,
    /// Whether to write the SAM header
    pub header:       bool,
    /// Any additional configuration
    pub config:       AlignerConfig,
}

/// The parsed and validated configuration options for `aligner`
pub struct AlignerConfig {
    /// Whether to also align the reverse complements
    pub rev_comp:         bool,
    /// Any override for which sequence to build the profile from
    pub profile_from:     WhichSequence,
    /// Any override for the number of passes to use
    pub method:           NumPasses,
    /// Whether to exclude unmapped alignments from the final output
    pub exclude_unmapped: bool,
    /// Whether to perform best match alignment
    pub best_match:       bool,
    /// The output path for the alignments (included in the config so that error
    /// context can be added)
    pub output:           Option<PathBuf>,
    /// Whether to set the Rayon number of threads to one
    #[cfg(not(feature = "dev_no_rayon"))]
    pub single_thread:    bool,
}

/// Parses and validates the arguments for `aligner` from the clap struct.
///
/// [`abort_clap`] will be called if:
///
/// - [`EitherParams::parse_from_clap`] fails (see the docs)
/// - The alphabet is [`Aa`] and `rev_comp` is true
///
/// ## Errors
///
/// Any IO errors from opening the queries, references, and output file are
/// propagated, with context containing the file path. If there is an invalid
/// record in the reference file, an error with the file path as context is
/// returned.
///
/// Any invalid records in the query file do not immediately produce errors
/// (since the reader is lazy), but any errors later produced will contain the
/// file path as context.
///
/// [`Aa`]: Alphabet::Aa
#[allow(unused_mut)]
pub fn parse_aligner_args(mut args: AlignerArgs) -> std::io::Result<ParsedAlignerArgs> {
    let params = EitherParams::parse_from_clap(
        args.alphabet,
        args.matrix,
        args.matching,
        args.mismatch,
        args.gap_open,
        args.gap_extend,
        args.ignore_n,
    );

    if matches!(params, EitherParams::Aa(_)) && args.rev_comp {
        abort_clap(
            ErrorKind::ArgumentConflict,
            "`--rev-comp` or `-r` cannot be specified with an amino acid alphabet",
            Some("aligner"),
        );
    }

    let query_reader = InputOptions::new_from_path(&args.query_file)
        .use_file_or_zip()
        .decode_in_thread()
        .parse_fastx()
        .open()?;

    let references = InputOptions::new_from_path(&args.ref_file)
        .use_file_or_zip()
        .parse_fasta()
        .open()?
        .collect::<Result<Vec<_>, _>>()?;

    // Validity: references field is required to be non-empty
    if references.is_empty() {
        return Err(std::io::Error::other(format!(
            "Empty reference file: {}",
            args.ref_file.display()
        )));
    }

    let mut profile_from = if args.profile_from_query {
        WhichSequence::Query
    } else if args.profile_from_ref {
        WhichSequence::Reference
    } else {
        match args.method {
            NumPasses::OnePass => WhichSequence::Reference,
            NumPasses::ThreePass => WhichSequence::Query,
        }
    };

    Ok(ParsedAlignerArgs {
        query_reader,
        references,
        params,
        header: args.header,
        config: AlignerConfig {
            rev_comp: args.rev_comp,
            profile_from,
            method: args.method,
            exclude_unmapped: args.exclude_unmapped,
            best_match: args.best_match,
            output: args.output,
            #[cfg(not(feature = "dev_no_rayon"))]
            single_thread: args.single_thread,
        },
    })
}

/// A clap enum for specifying the alphabet.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, ValueEnum)]
pub enum Alphabet {
    Dna,
    Aa,
}

/// A clap enum for specifying the method to use when aligning (one-pass or
/// three-pass).
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug, ValueEnum)]
pub enum NumPasses {
    #[value(name = "1pass", alias = "1-pass")]
    OnePass,
    #[value(name = "3pass", alias = "3-pass")]
    ThreePass,
}

impl Display for NumPasses {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            NumPasses::OnePass => write!(f, "1pass"),
            NumPasses::ThreePass => write!(f, "3pass"),
        }
    }
}

/// A clap enum for specifying which sequence to build the profile from.
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum WhichSequence {
    Query,
    Reference,
}

/// An enum abstracting over the different weight matrices supported by
/// `aligner`.
#[allow(clippy::large_enum_variant)]
pub enum EitherParams {
    /// Scoring parameters for a DNA alphabet
    Dna(ScoreParams<'static, i8, 5>),
    /// Scoring parameters for a protein alphabet
    Aa(ScoreParams<'static, i8, 25>),
}

impl EitherParams {
    /// Parses a [`EitherParams`] from clap arguments, determining the
    /// appropriate enum variant.
    ///
    /// The alphabet will default to [`Dna`], unless `matrix` is `Some`. When
    /// the alphabet is DNA, the scores default to `match = 2` and `mismatch =
    /// -5`. If `alphabet` is passed as [`Aa`] and no matrix name or scores are
    /// provided, [`BLOSUM_62`] is used.
    ///
    /// [`abort_clap`] will be called if:
    ///
    /// - `gap_extend` is a larger penalty than `gap_open`
    /// - If a matrix is loaded by name, the name must be valid using
    ///   [`aa_mat_from_name`]
    /// - `matrix` cannot be specified for a [`Dna`] alphabet
    /// - `ignore_n` cannot be specified when the alphabet is inferred to be
    ///   [`Aa`]
    ///
    /// ## Panics
    ///
    /// It is expected that clap will ensure the following:
    ///
    /// - `matching` and `mismatch` must be at most 127 if specified
    /// - Both `matrix` and `matching` cannot be specified
    ///
    /// [`Dna`]: Alphabet::Dna
    /// [`Aa`]: Alphabet::Aa
    fn parse_from_clap(
        alphabet: Option<Alphabet>, matrix: Option<String>, matching: Option<u8>, mismatch: Option<u8>, gap_open: u8,
        gap_extend: u8, ignore_n: bool,
    ) -> Self {
        // If either matching or mismatch is specified, then use the defaults
        // for the other
        let scores = match (matching, mismatch) {
            (None, None) => None,
            (Some(matching), None) => Some((matching, 5)),
            (None, Some(mismatch)) => Some((2, mismatch)),
            (Some(matching), Some(mismatch)) => Some((matching, mismatch)),
        };

        // Convert scores/penalties to signed weights
        let scores = scores.map(|(matching, mismatch)| {
            let matching = i8::try_from(matching).expect("Validated by clap");
            let mismatch = 0i8.checked_sub_unsigned(mismatch).expect("Validated by clap");
            (matching, mismatch)
        });
        let gap_open = 0i8.checked_sub_unsigned(gap_open).expect("Validated by clap");
        let gap_extend = 0i8.checked_sub_unsigned(gap_extend).expect("Validated by clap");

        let ignoring = if ignore_n { Some(b'N') } else { None };

        let result = match (alphabet, matrix, scores) {
            // Simple DNA weight matrix with default weights
            (None | Some(Alphabet::Dna), None, None) => {
                ScoreParams::new_dna_params(2, -5, gap_open, gap_extend, ignoring).map(EitherParams::Dna)
            }
            // Simple DNA weight matrix with user-specified weights
            (None | Some(Alphabet::Dna), None, Some((matching, mismatch))) => {
                ScoreParams::new_dna_params(matching, mismatch, gap_open, gap_extend, ignoring).map(EitherParams::Dna)
            }
            // BLOSUM62 default named protein weight matrix
            (Some(Alphabet::Aa), None, None) => ScoreParams::new(&BLOSUM_62, gap_open, gap_extend).map(EitherParams::Aa),
            // Named protein weight matrix
            (None | Some(Alphabet::Aa), Some(name), None) => {
                let Some(matrix) = aa_mat_from_name(&name) else {
                    abort_clap(
                        ErrorKind::InvalidValue,
                        format!("{name} is not a valid value for --matrix"),
                        Some("aligner"),
                    );
                };
                ScoreParams::new(matrix, gap_open, gap_extend).map(EitherParams::Aa)
            }
            // Simple protein weight matrix with user-specified weights
            (Some(Alphabet::Aa), None, Some((matching, mismatch))) => {
                let matrix = WeightMatrix::new(&AA_ALL_AMBIG_PROFILE_MAP_WITH_STOP, matching, mismatch, None);
                ScoreParams::new(matrix, gap_open, gap_extend).map(EitherParams::Aa)
            }

            // Invalid combinations
            (Some(Alphabet::Dna), Some(_), _) => abort_clap(
                ErrorKind::ArgumentConflict,
                "--matrix cannot be specified with a DNA alphabet",
                Some("aligner"),
            ),
            (_, Some(_), Some(_)) => {
                unreachable!("Validated by clap")
            }
        };

        let params = match result {
            Ok(params) => params,
            Err(ProfileError::BadGapWeights { gap_open, gap_extend }) => abort_clap(
                ErrorKind::InvalidValue,
                format!(
                    "The gap open penalty must be greater than or equal to the gap extend penalty, but {gap_open} (gap open) and {gap_extend} (gap extend) were provided",
                    gap_open = gap_open.abs(),
                    gap_extend = gap_extend.abs(),
                ),
                Some("aligner"),
            ),
            Err(_) => unreachable!("Validated by clap"),
        };

        if ignore_n && matches!(params, EitherParams::Aa(_)) {
            abort_clap(
                ErrorKind::ArgumentConflict,
                "--ignore-n cannot be specified with an amino acid alphabet",
                Some("aligner"),
            )
        }

        params
    }
}
