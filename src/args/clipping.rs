use crate::utils::get_hasher;
use clap::{Args, ValueEnum, builder::PossibleValue};
use foldhash::fast::SeedableRandomState;
use std::{fmt::Debug, num::NonZeroUsize, path::PathBuf};
use zoe::{
    kmer::ThreeBitKmerSet,
    prelude::{CheckNucleotides, FastaReader, IsValidDNA, Nucleotides},
};

pub const MAX_KMER_LENGTH: usize = 21;

#[derive(Args, Debug)]
pub struct ClippingArgs {
    #[arg(short = 's', long, alias = "preserve_fastq")]
    /// Preserve original formatting: disables encoding sequences into uppercase
    /// canonical bases (A, C, T, G, N)
    pub preserve_bases: bool,

    #[arg(short = 'G', long)]
    /// Trims multiple consecutive G's (a possible artifact of Illumina
    /// sequencing) found at the ends of a sequence
    pub polyg_trim: Option<NonZeroUsize>,

    #[arg(long, default_value = "b", requires = "polyg_trim")]
    /// Specifies the end of the sequence for poly-g trimming : 'l' (left), 'r'
    /// (right), or 'b' (both)
    pub g_polyg_end: TrimEnd,

    #[arg(long, requires = "polyg_trim")]
    /// Overrides the threshold for consecutive G's on the left end of the
    /// sequence for poly-g trimming
    pub g_polyg_left: Option<NonZeroUsize>,

    #[arg(long, requires = "polyg_trim")]
    /// Overrides the threshold for consecutive G's on the right end of the
    /// sequence for poly-g trimming
    pub g_polyg_right: Option<NonZeroUsize>,

    #[arg(short = 'B', long, value_parser = validate_acgtn, group = "adapter_vs_barcode")]
    /// Trim barcodes and their reverse complements from sequence using string
    /// matching. Requires literal barcode as argument
    pub barcode_trim: Option<Nucleotides>,

    #[arg(long, value_enum, default_value = "b", requires = "barcode_trim")]
    /// Specifies the end of the sequence for barcode trimming : 'l' (left), 'r'
    /// (right), or 'b' (both)
    pub b_end: TrimEnd,

    #[arg(long, requires = "barcode_trim")]
    /// Restriction window size for barcode trimming on both ends of the
    /// sequence. If no restriction is provided, full scan is performed
    pub b_restrict: Option<NonZeroUsize>,

    #[arg(long, requires = "barcode_trim")]
    /// Restriction window for trimming barcodes on the left end of the
    /// sequence. Overrides --b-restrict
    pub b_restrict_left: Option<NonZeroUsize>,

    #[arg(long, requires = "barcode_trim")]
    /// Restriction window for trimming barcodes on the right end of the
    /// sequence. Overrides --b_restrict
    pub b_restrict_right: Option<NonZeroUsize>,

    #[arg(long, value_parser = validate_b_hdist, default_value = "0", requires = "barcode_trim")]
    /// Accepted Hamming distance for fuzzy barcode matching and trimming,
    /// between 0 and 3
    pub b_hdist: usize,

    #[arg(short = 'A', long, value_parser = validate_acgtn, group = "adapter_vs_barcode")]
    /// Trim adapters and their reverse complements from sequence. Requires
    /// literal adapter as argument
    pub adapter_trim: Option<Nucleotides>,

    #[arg(long, requires = "adapter_trim")]
    /// Allow up to one mismatch during adapter matching and trimming
    pub a_fuzzy: bool,

    #[arg(short = 'P', long, requires = "p_kmer_length")]
    /// Trim primers from sequence using k-mer matching. Requires path to primer
    /// fasta file and a kmer length
    pub primer_trim: Option<PathBuf>,

    #[arg(long, requires = "primer_trim")]
    /// Enables fuzzy matching (one mismatch) for k-mer searching of primers
    pub p_fuzzy: bool,

    #[arg(long, value_parser = validate_kmer_length, requires = "primer_trim")]
    /// Length of k-mer used for matching primers.
    pub p_kmer_length: Option<usize>,

    #[arg(long, default_value = "b", requires = "primer_trim")]
    /// Specifies the end of the sequence for primer trimming : 'l' (left), 'r'
    /// (right), or 'b' (both)
    pub p_end: TrimEnd,

    #[arg(long, default_value = "30", requires = "primer_trim")]
    /// Restriction window size for primer trimming on both ends of the sequence
    pub p_restrict: NonZeroUsize,

    #[arg(long, requires = "primer_trim")]
    /// Restriction window for trimming primer on the left end of the sequence
    /// Overrides --p_restrict
    pub p_restrict_left: Option<NonZeroUsize>,

    #[arg(long, requires = "primer_trim")]
    /// Restriction window for trimming barcodes on the right end of the
    /// sequence. Overrides --p_restrict
    pub p_restrict_right: Option<NonZeroUsize>,

    #[arg(short = 'H', long)]
    /// Hard trim from each end the specified number of bases
    pub hard_trim: Option<usize>,

    #[arg(long)]
    /// Hard trim range for only the left end of the sequence. Overrides
    /// hard-trim
    pub h_left: Option<usize>,

    #[arg(long)]
    /// Hard trim range for only the right end of the sequence. Overrides
    /// hard-trim
    pub h_right: Option<usize>,
}

/// Enum for trimming end options
#[derive(Debug, Clone, Copy)]
pub enum TrimEnd {
    L, // Left
    R, // Right
    B, // Both
}

// Allows case insensitivity for trim ends
impl ValueEnum for TrimEnd {
    #[inline]
    fn value_variants<'a>() -> &'a [Self] {
        &[Self::L, Self::R, Self::B]
    }

    #[inline]
    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            TrimEnd::L => Some(PossibleValue::new("l").alias("L")),
            TrimEnd::R => Some(PossibleValue::new("r").alias("R")),
            TrimEnd::B => Some(PossibleValue::new("b").alias("B")),
        }
    }
}

/// Ensures user has entered valid non-empty adapter or barcode literal for
/// trimming
fn validate_acgtn(value: &str) -> Result<Nucleotides, String> {
    if value.trim().is_empty() {
        // prevents panicking when `-A ""` is passed
        Err("Adapter (-A) or barcode (-B) cannot be empty!".to_string())
    } else if value.as_bytes().is_valid_dna(IsValidDNA::AcgtnNoGaps) {
        Ok(value.as_bytes().into())
    } else {
        Err("Adapter or barcode literal must only consist of canonical (ACGTN) bases".to_string())
    }
}

/// Custom validator for `kmer_length`
fn validate_kmer_length(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("`{value}` is not a valid integer."))?;
    if (2..=MAX_KMER_LENGTH).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!(
            "p-kmer-length must be between 2 and {MAX_KMER_LENGTH}, but {parsed} was provided."
        ))
    }
}

/// Validates barcode hamming distance to be between 0 and 3
fn validate_b_hdist(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("`{value}` is not a valid integer between 0 and 3."))?;
    if (0..=3).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!("b-hdist must be between 0 and 3, but `{parsed}` was provided."))
    }
}

/// For converting adapters and barcodes to uppercase (unless disabled) and
/// getting reverse complements
fn get_forward_reverse_sequence(mut adapter: Nucleotides, preserve_seq: bool) -> (Nucleotides, Nucleotides) {
    if !preserve_seq {
        adapter.as_mut_bytes().make_ascii_uppercase();
    }
    let reverse = adapter.to_reverse_complement();
    (adapter, reverse)
}

/// Reads a primer file and generates a k-mer set of all unique k-mers present
/// in the sequence and reverse complements.
///
/// If `fuzzy_kmer` is set, then k-mers with up to one mismatch are also
/// included. This mismatch can involve the introduction of an ambiguous base
/// 'N'.
///
/// ## Errors
///
/// `primer_path` must be successfully opened, and all lines must be parsed
/// without error.
///
/// ## Panics
///
/// `kmer_length` must be between 2 and [`MAX_KMER_LENGTH`], inclusive.
fn prepare_primer_kmers(
    primer_path: &PathBuf, kmer_length: usize, fuzzy_kmer: bool,
) -> std::io::Result<ThreeBitKmerSet<MAX_KMER_LENGTH, SeedableRandomState>> {
    let mut fasta_primer_reader = FastaReader::from_filename(primer_path)?;

    let mut unique_kmers =
        ThreeBitKmerSet::<MAX_KMER_LENGTH, _>::with_hasher(kmer_length, get_hasher()).expect("Expected valid kmer length");

    let insert_fn = if fuzzy_kmer {
        |kmer_set: &mut ThreeBitKmerSet<MAX_KMER_LENGTH, SeedableRandomState>, seq: &Nucleotides| {
            kmer_set.insert_from_sequence_with_variants::<1>(seq)
        }
    } else {
        |kmer_set: &mut ThreeBitKmerSet<MAX_KMER_LENGTH, SeedableRandomState>, seq: &Nucleotides| {
            kmer_set.insert_from_sequence(seq)
        }
    };

    fasta_primer_reader.try_for_each(|f| -> Result<(), std::io::Error> {
        let mut seq = Nucleotides::from_vec_unchecked(f?.sequence);
        insert_fn(&mut unique_kmers, &seq);
        seq.make_reverse_complement();
        insert_fn(&mut unique_kmers, &seq);
        Ok(())
    })?;

    Ok(unique_kmers)
}

/// Arguments specifying the types of clipping to be performed
#[derive(Debug)]
pub struct ParsedClippingArgs {
    pub preserve_bases:   bool,
    pub barcodes:         Option<(Nucleotides, Nucleotides)>,
    pub b_restrict_left:  Option<usize>,
    pub b_restrict_right: Option<usize>,
    pub b_hdist:          usize,
    pub adapters:         Option<(Nucleotides, Nucleotides)>,
    pub a_fuzzy:          bool,
    pub primer_kmers:     Option<ThreeBitKmerSet<MAX_KMER_LENGTH, SeedableRandomState>>,
    pub p_restrict_left:  Option<usize>,
    pub p_restrict_right: Option<usize>,
    pub polyg_left:       Option<usize>,
    pub polyg_right:      Option<usize>,
    pub hard_left:        usize,
    pub hard_right:       usize,
}

/// Parses all arguments related to clipping.
///
/// ## Errors
///
/// Any errors while processing the primers are propagated.
pub fn parse_clipping_args(args: ClippingArgs) -> std::io::Result<ParsedClippingArgs> {
    let ClippingArgs {
        preserve_bases,
        polyg_trim,
        g_polyg_end,
        g_polyg_left,
        g_polyg_right,
        barcode_trim,
        b_end,
        b_restrict,
        b_restrict_left,
        b_restrict_right,
        b_hdist,
        adapter_trim,
        a_fuzzy,
        primer_trim,
        p_fuzzy,
        p_kmer_length,
        p_end,
        p_restrict,
        p_restrict_left,
        p_restrict_right,
        hard_trim,
        h_left,
        h_right,
    } = args;

    let adapters = adapter_trim.map(|adapter| get_forward_reverse_sequence(adapter, preserve_bases));
    let barcodes = barcode_trim.map(|barcode| get_forward_reverse_sequence(barcode, preserve_bases));

    let primer_kmers = if let Some(primer_path) = &primer_trim {
        Some(prepare_primer_kmers(
            primer_path,
            // This is unreachable through clap due to being required
            p_kmer_length.expect("A kmer length must be provided for primer trimming"),
            p_fuzzy,
        )?)
    } else {
        None
    };

    // A value of None for left or right restricts will do full scan barcoding
    let default_b_restrict = b_restrict;
    let (b_restrict_left, b_restrict_right) = match b_end {
        TrimEnd::B => {
            // Apply b_restrict to both, unless explicitly overridden
            let left = b_restrict_left.or(default_b_restrict).map(NonZeroUsize::get);
            let right = b_restrict_right.or(default_b_restrict).map(NonZeroUsize::get);
            (left, right)
        }
        TrimEnd::L => {
            // Left gets b_left_restrict or b_restrict, right is forced to 0
            let left = b_restrict_left.or(default_b_restrict).map(NonZeroUsize::get);
            (left, Some(0))
        }
        TrimEnd::R => {
            // Right gets b_right_restrict or b_restrict, left is forced to 0
            let right = b_restrict_right.or(default_b_restrict).map(NonZeroUsize::get);
            (Some(0), right)
        }
    };

    let default_p_restrict = p_restrict;
    let (p_restrict_left, p_restrict_right) = match p_end {
        TrimEnd::B => {
            let left = p_restrict_left.unwrap_or(default_p_restrict);
            let right = p_restrict_right.unwrap_or(default_p_restrict);
            (Some(left.get()), Some(right.get()))
        }
        TrimEnd::L => {
            let left = p_restrict_left.unwrap_or(default_p_restrict);
            (Some(left.get()), None)
        }
        TrimEnd::R => {
            let right = p_restrict_right.unwrap_or(default_p_restrict);
            (None, Some(right.get()))
        }
    };

    let default_polyg = polyg_trim;
    let (polyg_left, polyg_right) = match g_polyg_end {
        TrimEnd::B => {
            let left = g_polyg_left.or(default_polyg).map(NonZeroUsize::get);
            let right = g_polyg_right.or(default_polyg).map(NonZeroUsize::get);
            (left, right)
        }
        TrimEnd::L => {
            let left = g_polyg_left.or(default_polyg).map(NonZeroUsize::get);
            (left, None)
        }
        TrimEnd::R => {
            let right = g_polyg_right.or(default_polyg).map(NonZeroUsize::get);
            (None, right)
        }
    };

    let default_hard_bases = hard_trim.unwrap_or(0);
    let hard_left = h_left.unwrap_or(default_hard_bases);
    let hard_right = h_right.unwrap_or(default_hard_bases);

    let parsed_args = ParsedClippingArgs {
        preserve_bases,
        barcodes,
        b_restrict_left,
        b_restrict_right,
        b_hdist,
        adapters,
        a_fuzzy,
        primer_kmers,
        p_restrict_left,
        p_restrict_right,
        polyg_left,
        polyg_right,
        hard_left,
        hard_right,
    };

    Ok(parsed_args)
}
