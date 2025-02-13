// Description:      Read FastQ files and trim with various options.
use clap::{builder::PossibleValue, Args, ValueEnum};
use either::Either;
use foldhash::fast::RandomState;
use std::{
    fs::{File, OpenOptions},
    io::{stdout, BufReader, BufWriter, Stdout, Write},
    num::NonZeroUsize,
    path::PathBuf,
};
use zoe::{data::types::nucleotides::reverse_complement, kmer::ThreeBitKmerSet, prelude::*};

use crate::qc::fastq::ReadTransforms;

const MAX_KMER_LENGTH: usize = 21;

#[derive(Args, Debug)]
pub struct TrimmerArgs {
    /// Path to .fastq file to be trimmed
    pub fastq_input_file1: PathBuf,

    /// Path to optional second .fastq file to be trimmed
    pub fastq_input_file2: Option<PathBuf>,

    #[arg(short = 'o', long)]
    /// Output path for trimmed files. Trimmed sequences print to stdout if not provided
    pub fastq_output_file: Option<PathBuf>,

    #[arg(short = 's', long)]
    /// Preserve original formatting: disables encoding sequences into uppercase canonical bases (A, C, T, G, N)
    pub preserve_fastq: bool,

    #[arg(short = 'm', long)]
    /// Perform masking with 'N' instead of clipping. Default behavior is clipping if not provided
    pub mask: bool,

    #[arg(short = 'n', long, default_value = "1")]
    /// Minimum sequence length required after trimming. Shorter sequences are filtered from output.
    pub min_length: NonZeroUsize,

    #[arg(short = 'B', long, value_parser = validate_non_empty, group = "adapter_vs_barcode")]
    /// Trim barcodes and their reverse complements from sequence using string matching. Requires literal barcode as argument
    pub barcode_trim: Option<String>,

    #[arg(long, value_enum, default_value = "b", requires = "barcode_trim")]
    /// Specifies the end of the sequence for barcode trimming : 'l' (left), 'r' (right), or 'b' (both)
    pub b_end: TrimEnd,

    #[arg(long, requires = "barcode_trim")]
    /// Restriction window size for barcode trimming on both ends of the sequence. If no restriction is provided, full scan is performed.
    pub b_restrict: Option<usize>,

    #[arg(long, requires = "barcode_trim")]
    /// Restriction window for trimming barcodes on the left end of the sequence. Overrides --b-restrict
    pub b_restrict_left: Option<usize>,

    #[arg(long, requires = "barcode_trim")]
    /// Restriction window for trimming barcodes on the right end of the sequence. Overrides --b_restrict
    pub b_restrict_right: Option<usize>,

    #[arg(long, value_parser = validate_barcode_hdist, default_value = "0", requires = "barcode_trim")]
    /// Accepted Hamming distance for fuzzy barcode matching and trimming, between 0 and 3
    pub b_hdist: usize,

    #[arg(short = 'A', long, value_parser = validate_non_empty, group = "adapter_vs_barcode")]
    /// Trim adapters and their reverse complements from sequence. Requires literal adapter as argument
    pub adapter_trim: Option<String>,

    #[arg(long, requires = "adapter_trim")]
    /// Allow up to one mismatch during adapter matching and trimming
    pub a_fuzzy_adapter: bool,

    #[arg(short = 'P', long)]
    /// Trim primers from sequence using kmer matching. Requires path to primer fasta file
    pub primer_trim: Option<PathBuf>,

    #[arg(long, requires = "primer_trim")]
    /// Enables fuzzy matching (one mismatch) for kmer searching of primers
    pub p_fuzzy_kmer: bool,

    #[arg(long, value_parser = validate_kmer_length, default_value = "17", requires = "primer_trim")]
    /// Length of k-mer used for matching primers.
    pub p_kmer_length: usize,

    #[arg(long, default_value = "b", requires = "primer_trim")]
    /// Specifies the end of the sequence for primer trimming : 'l' (left), 'r' (right), or 'b' (both)
    pub p_end: TrimEnd,

    #[arg(long, default_value = "30", requires = "primer_trim")]
    /// Restriction window size for primer trimming on both ends of the sequence
    pub p_restrict: usize,

    #[arg(long, requires = "primer_trim")]
    /// Restriction window for trimming primer on the left end of the sequence. Overrides --p_restrict
    pub p_restrict_left: Option<usize>,

    #[arg(long, requires = "primer_trim")]
    /// Restriction window for trimming barcodes on the right end of the sequence. Overrides --p_restrict
    pub p_restrict_right: Option<usize>,

    #[arg(short = 'H', long)]
    /// Hard trim from each end the specified number of bases
    pub hard_trim: Option<usize>,

    #[arg(long)]
    /// Hard trim range for only the left end of the sequence. Overrides hard-trim
    pub h_left: Option<usize>,

    #[arg(long)]
    /// Hard trim range for only the right end of the sequence. Overrides hard-trim
    pub h_right: Option<usize>,
}

/// Sub-program for trimming fastQ data.
pub fn trimmer_process(args: TrimmerArgs) -> Result<(), std::io::Error> {
    let mut args = parse_trim_args(args)?;
    let writer = &mut args.output_writer;

    let mut chained_reader = args.fastq_reader1.chain(args.fastq_reader2.into_iter().flatten());

    chained_reader.by_ref().try_for_each(|record| -> Result<(), std::io::Error> {
        let mut fq = record?;
        let mut fq_view = fq.as_view_mut();
        fq_view.to_canonical_bases(!args.preserve_seq);

        if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
            fq_view.transform_by_reverse_forward_search(args.a_fuzzy_adapter, !args.mask, reverse_adapter, forward_adapter);
        } else if let Some((barcode, reverse)) = &args.barcodes {
            fq_view.barcode_trim(
                barcode,
                reverse,
                args.b_hdist,
                args.mask,
                args.b_restrict_left,
                args.b_restrict_right,
            );
        }
        if let Some(ref kmers) = args.primer_kmers {
            if args.p_restrict_left > 0 {
                fq_view.left_primer_trim(args.p_restrict_left, kmers, args.mask);
            }
            if args.p_restrict_right > 0 {
                fq_view.right_primer_trim(args.p_restrict_right, kmers, args.mask);
            }
        }
        if args.hard_left > 0 || args.hard_right > 0 {
            fq_view.hard_trim(args.hard_left, args.hard_right, args.mask);
        }
        if args.mask {
            write!(writer, "{fq}")?;
        } else if fq_view.len() >= args.min_length {
            write!(writer, "{fq_view}")?;
        }
        Ok(())
    })?;
    writer.flush()?;

    Ok(())
}

// Enum for trimming end options
#[derive(Debug, Clone, Copy)]
pub enum TrimEnd {
    L, // Left
    R, // Right
    B, // Both
}

// Allows case insensitivity for trim ends
impl ValueEnum for TrimEnd {
    fn value_variants<'a>() -> &'a [Self] {
        &[Self::L, Self::R, Self::B]
    }

    fn to_possible_value(&self) -> Option<PossibleValue> {
        match self {
            TrimEnd::L => Some(PossibleValue::new("l").alias("L")),
            TrimEnd::R => Some(PossibleValue::new("r").alias("R")),
            TrimEnd::B => Some(PossibleValue::new("b").alias("B")),
        }
    }
}

// Custom validator for `kmer_length`.
fn validate_kmer_length(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("`{}` is not a valid integer.", value))?;
    if (2..=MAX_KMER_LENGTH).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!(
            "`kmer-length` must be between 2 and {MAX_KMER_LENGTH}, but {parsed} was provided."
        ))
    }
}

// validates barcode hamming distance to be between 0 and 3
fn validate_barcode_hdist(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("`{}` is not a valid integer between 0 and 3.", value))?;
    if (0..=3).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!(
            "`hamming-distance` must be between 0 and 3, but `{}` was provided.",
            parsed
        ))
    }
}

// ensures non-empty string passed to adapter and barcode (e.g. `-B ""...`)
fn validate_non_empty(value: &str) -> Result<String, String> {
    if value.trim().is_empty() {
        Err("Adapter (-A) or barcode (-B) cannot be empty!".to_string())
    } else {
        Ok(value.to_string())
    }
}

pub struct ParsedTrimmerArgs {
    pub fastq_reader1:    FastQReader<BufReader<File>>,
    pub fastq_reader2:    Option<FastQReader<BufReader<File>>>,
    pub output_writer:    BufWriter<Either<File, Stdout>>,
    pub preserve_seq:     bool,
    pub mask:             bool,
    pub min_length:       usize,
    pub barcodes:         Option<(Vec<u8>, Vec<u8>)>,
    pub b_restrict_left:  Option<usize>,
    pub b_restrict_right: Option<usize>,
    pub b_hdist:          usize,
    pub adapters:         Option<(Vec<u8>, Vec<u8>)>,
    pub a_fuzzy_adapter:  bool,
    pub primer_kmers:     Option<ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>>,
    pub p_restrict_left:  usize,
    pub p_restrict_right: usize,
    pub hard_left:        usize,
    pub hard_right:       usize,
}

pub fn parse_trim_args(args: TrimmerArgs) -> Result<ParsedTrimmerArgs, std::io::Error> {
    let fastq_reader1 = FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(&args.fastq_input_file1)?));

    let fastq_reader2 = if let Some(file2) = &args.fastq_input_file2 {
        Some(FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(file2)?)))
    } else {
        None
    };

    let output_writer = match &args.fastq_output_file {
        Some(path) => BufWriter::new(Either::Left(File::create(path)?)),
        None => BufWriter::new(Either::Right(stdout())),
    };

    let adapters = args
        .adapter_trim
        .map(|adapter| get_forward_reverse_sequence(&adapter, args.preserve_fastq));

    let barcodes = args
        .barcode_trim
        .map(|barcode| get_forward_reverse_sequence(&barcode, args.preserve_fastq));

    let primer_kmers = if let Some(primer_path) = &args.primer_trim {
        Some(prepare_primer_kmers(primer_path, args.p_kmer_length, args.p_fuzzy_kmer)?)
    } else {
        None
    };

    // A value of None for left or right restricts will do full scan barcoding
    let default_b_restrict = args.b_restrict;
    let (b_restrict_left, b_restrict_right) = match args.b_end {
        TrimEnd::B => {
            // Apply b_restrict to both, unless explicitly overridden
            let left = args.b_restrict_left.or(default_b_restrict);
            let right = args.b_restrict_right.or(default_b_restrict);
            (left, right)
        }
        TrimEnd::L => {
            // Left gets b_left_restrict or b_restrict, right is forced to 0
            let left = args.b_restrict_left.or(default_b_restrict);
            (left, Some(0))
        }
        TrimEnd::R => {
            // Right gets b_right_restrict or b_restrict, left is forced to 0
            let right = args.b_restrict_right.or(default_b_restrict);
            (Some(0), right)
        }
    };

    let default_p_restrict = args.p_restrict;
    let (p_restrict_left, p_restrict_right) = match args.p_end {
        TrimEnd::B => {
            let left = args.p_restrict_left.unwrap_or(default_p_restrict);
            let right = args.p_restrict_right.unwrap_or(default_p_restrict);
            (left, right)
        }
        TrimEnd::L => {
            let left = args.p_restrict_left.unwrap_or(default_p_restrict);
            (left, 0)
        }
        TrimEnd::R => {
            let right = args.p_restrict_right.unwrap_or(default_p_restrict);
            (0, right)
        }
    };

    let default_hard_bases = args.hard_trim.unwrap_or(0);
    let hard_left = args.h_left.unwrap_or(default_hard_bases);
    let hard_right = args.h_right.unwrap_or(default_hard_bases);

    let parsed_args = ParsedTrimmerArgs {
        fastq_reader1,
        fastq_reader2,
        output_writer,
        preserve_seq: args.preserve_fastq,
        mask: args.mask,
        min_length: args.min_length.get(),
        barcodes,
        b_restrict_left,
        b_restrict_right,
        b_hdist: args.b_hdist,
        adapters,
        a_fuzzy_adapter: args.a_fuzzy_adapter,
        primer_kmers,
        p_restrict_left,
        p_restrict_right,
        hard_left,
        hard_right,
    };
    Ok(parsed_args)
}

fn get_forward_reverse_sequence(adapter: &str, preserve_seq: bool) -> (Vec<u8>, Vec<u8>) {
    let mut forward = adapter.as_bytes().to_vec();
    // if the user doesn't specify to preserve the original fastq, we encode to uppercase canonical bases
    if !preserve_seq {
        forward.make_ascii_uppercase();
    }
    let reverse = reverse_complement(&forward);
    (forward, reverse)
}

/// # Panics
///
/// `kmer_length` must be between 2 and 21, inclusive.
fn prepare_primer_kmers(
    primer_path: &PathBuf, kmer_length: usize, fuzzy_kmer: bool,
) -> Result<ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, std::io::Error> {
    let mut fasta_primer_reader = FastaReader::from_filename(primer_path)?;

    let mut unique_kmers = ThreeBitKmerSet::<MAX_KMER_LENGTH, _>::with_hasher(kmer_length, RandomState::default())
        .expect("Expected valid kmer length");

    let insert_fn = if fuzzy_kmer {
        |kmer_set: &mut ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, seq: &Nucleotides| {
            kmer_set.insert_from_sequence_with_variants::<1>(seq)
        }
    } else {
        |kmer_set: &mut ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, seq: &Nucleotides| kmer_set.insert_from_sequence(seq)
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
