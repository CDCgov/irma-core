use crate::{
    aligner::{AlignerArgs, QueryReader},
    args::abort_clap,
    io::{FastXReader, IterFromFilename, ReadFileZip, ReadFileZipPipe},
};
use clap::{ValueEnum, builder::PossibleValue, error::ErrorKind};
use std::{fmt::Display, path::PathBuf};
use zoe::{
    data::{
        AA_ALL_AMBIG_PROFILE_MAP_WITH_STOP, WeightMatrix,
        fasta::FastaSeq,
        matrices::{BLOSUM_62, aa_mat_from_name},
    },
    math::AnyInt,
    prelude::FastaReader,
};

/// The parsed and validated command line arguments for `aligner`
pub struct ParsedAlignerArgs {
    /// The streamed query sequences
    pub query_reader:  QueryReader,
    /// The slurped reference sequences
    ///
    /// ## Validity
    ///
    /// This field must be non-empty.
    pub references:    Vec<FastaSeq>,
    /// The output path for the alignments
    pub output:        Option<PathBuf>,
    /// The weight matrix to use for the alignment
    pub weight_matrix: AnyMatrix<'static, i8>,
    /// Whether to write the SAM header
    pub header:        bool,
    /// Any additional configuration
    pub config:        AlignerConfig,
}

/// The parsed and validated configuration options for `aligner`
pub struct AlignerConfig {
    /// The gap open weight (should be non-positive)
    pub gap_open:         i8,
    /// The gap extend weight (should be non-positive)
    pub gap_extend:       i8,
    /// Whether to also align the reverse complements
    pub rev_comp:         bool,
    /// Whether the profiles are built from the references
    pub profile_from_ref: bool,
    /// Whether to exclude unmapped alignments from the final output
    pub exclude_unmapped: bool,
    /// Whether to perform best match alignment
    pub best_match:       bool,
    /// Whether to set the Rayon number of threads to one
    #[cfg(not(feature = "dev_no_rayon"))]
    pub single_thread:    bool,
}

/// Parses and validates the arguments for `aligner` from the clap struct.
///
/// [`abort_clap`] will be called if:
///
/// - [`AnyMatrix::parse_from_clap`] fails (see the docs)
/// - The alphabet is [`Aa`] and `rev_comp` is true
/// - The gap open penalty is smaller than the gap extend penalty
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
pub fn parse_aligner_args(args: AlignerArgs) -> std::io::Result<ParsedAlignerArgs> {
    let weight_matrix = AnyMatrix::parse_from_clap(args.alphabet, args.matrix, args.matching, args.mismatch, args.ignore_n);

    if weight_matrix.alphabet() == Alphabet::Aa && args.rev_comp {
        abort_clap(
            ErrorKind::ArgumentConflict,
            "`--rev-comp` or `-r` cannot be specified with an amino acid alphabet",
            Some("aligner"),
        );
    }

    let gap_open = -(args.gap_open as i8);
    let gap_extend = -(args.gap_extend as i8);

    if gap_extend < gap_open {
        abort_clap(
            ErrorKind::InvalidValue,
            format!(
                "The gap open penalty must be greater than or equal to the gap extend penalty, but {gap_open} (gap open) and {gap_extend} (gap extend) were provided",
                gap_open = args.gap_open,
                gap_extend = args.gap_extend
            ),
            Some("aligner"),
        )
    }

    let query_reader = FastXReader::<ReadFileZipPipe>::from_filename(&args.query_file)?;

    let references = FastaReader::<ReadFileZip>::from_filename(&args.ref_file)?.collect::<Result<Vec<_>, _>>()?;

    // Validity: references field is required to be non-empty
    if references.is_empty() {
        return Err(std::io::Error::other(format!(
            "Empty reference file: {}",
            args.ref_file.display()
        )));
    }

    Ok(ParsedAlignerArgs {
        query_reader,
        references,
        output: args.output,
        weight_matrix,
        header: args.header,
        config: AlignerConfig {
            gap_open,
            gap_extend,
            rev_comp: args.rev_comp,
            profile_from_ref: args.profile_from_ref,
            exclude_unmapped: args.exclude_unmapped,
            best_match: args.best_match,
            #[cfg(not(feature = "dev_no_rayon"))]
            single_thread: args.single_thread,
        },
    })
}

/// A clap enum for specifying the alphabet
#[derive(Debug, Clone, Copy, Eq, PartialEq, Hash)]
pub enum Alphabet {
    Dna,
    Aa,
}

impl Display for Alphabet {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Alphabet::Dna => write!(f, "DNA"),
            Alphabet::Aa => write!(f, "AA"),
        }
    }
}

impl ValueEnum for Alphabet {
    #[inline]
    fn value_variants<'a>() -> &'a [Self] {
        &[Self::Dna, Self::Aa]
    }

    #[inline]
    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            Self::Dna => Some(PossibleValue::new("DNA").alias("dna")),
            Self::Aa => Some(PossibleValue::new("AA").alias("aa")),
        }
    }
}

/// An enum abstracting over the different weight matrices supported by
/// `aligner`.
///
/// This is required because different alphabets require different const
/// generics.
pub enum AnyMatrix<'a, T: AnyInt + 'static> {
    /// A weight matrix for a DNA alphabet
    Dna(WeightMatrix<'a, T, 5>),
    /// A named weight matrix for a protein alphabet, obtained from Zoe
    AaNamed(&'static WeightMatrix<'static, T, 25>),
    /// A simple weight matrix for a protein alphabet
    AaSimple(WeightMatrix<'a, T, 25>),
}

impl<'a, T: AnyInt + 'static> AnyMatrix<'a, T> {
    /// Gets the alphabet associated with the weight matrix variant
    #[inline]
    #[must_use]
    fn alphabet(&self) -> Alphabet {
        match self {
            AnyMatrix::Dna(_) => Alphabet::Dna,
            AnyMatrix::AaNamed(_) => Alphabet::Aa,
            AnyMatrix::AaSimple(_) => Alphabet::Aa,
        }
    }
}

impl<'a> AnyMatrix<'a, i8> {
    /// Parses a [`AnyMatrix`] from clap arguments, determining the appropriate
    /// enum variant.
    ///
    /// The alphabet will default to [`Dna`], unless `matrix` is `Some`. When
    /// the alphabet is DNA, the scores default to `match = 2` and `mismatch =
    /// -5`. If `alphabet` is passed as [`Aa`] and no matrix name or scores are
    /// provided, [`BLOSUM_62`] is used.
    ///
    /// [`abort_clap`] will be called if:
    ///
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
        alphabet: Option<Alphabet>, matrix: Option<String>, matching: Option<u8>, mismatch: Option<u8>, ignore_n: bool,
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

        let matrix: AnyMatrix<'_, i8> = match (alphabet, matrix, scores) {
            // Simple DNA weight matrix with default weights
            (None | Some(Alphabet::Dna), None, None) => {
                WeightMatrix::new_dna_matrix(2, -5, if ignore_n { Some(b'N') } else { None }).into()
            }
            // Simple DNA weight matrix with user-specified weights
            (None | Some(Alphabet::Dna), None, Some((matching, mismatch))) => {
                WeightMatrix::new_dna_matrix(matching, mismatch, if ignore_n { Some(b'N') } else { None }).into()
            }
            // BLOSUM62 default named protein weight matrix
            (Some(Alphabet::Aa), None, None) => (&BLOSUM_62).into(),
            // Named protein weight matrix
            (None | Some(Alphabet::Aa), Some(name), None) => {
                let Some(matrix) = aa_mat_from_name(&name) else {
                    abort_clap(
                        ErrorKind::InvalidValue,
                        format!("{name} is not a valid value for --matrix"),
                        Some("aligner"),
                    );
                };
                matrix.into()
            }
            // Simple protein weight matrix with user-specified weights
            (Some(Alphabet::Aa), None, Some((matching, mismatch))) => {
                WeightMatrix::new(&AA_ALL_AMBIG_PROFILE_MAP_WITH_STOP, matching, mismatch, None).into()
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

        if ignore_n && matrix.alphabet() == Alphabet::Aa {
            abort_clap(
                ErrorKind::ArgumentConflict,
                "--ignore-n cannot be specified with an amino acid alphabet",
                Some("aligner"),
            )
        }

        matrix
    }
}

impl<'a, T: AnyInt> From<WeightMatrix<'a, T, 5>> for AnyMatrix<'a, T> {
    #[inline]
    fn from(value: WeightMatrix<'a, T, 5>) -> Self {
        Self::Dna(value)
    }
}

impl<'a, T: AnyInt + 'static> From<&'static WeightMatrix<'a, T, 25>> for AnyMatrix<'a, T> {
    #[inline]
    fn from(value: &'static WeightMatrix<'a, T, 25>) -> Self {
        Self::AaNamed(value)
    }
}

impl<'a, T: AnyInt> From<WeightMatrix<'a, T, 25>> for AnyMatrix<'a, T> {
    #[inline]
    fn from(value: WeightMatrix<'a, T, 25>) -> Self {
        Self::AaSimple(value)
    }
}
