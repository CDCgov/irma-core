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
// An example run:
// cargo run -- trimmer 3003863475_N8KHVRSA_S384_R1_001.fastq -B CACAAAGACACCGACAACTTTCTT --b-end b -u -m -n 30 --b-restrict 36 -P swift_211206.fasta --p-end b --p-restrict 30 --p-kmer-length 17 -H 10 -o trimmed.fastq

/* Assumed args */

#[derive(Args, Debug)]
pub struct TrimmerArgs {
    pub fastq_input_file1: PathBuf,

    pub fastq_input_file2: Option<PathBuf>,

    #[arg(short = 'o', long)]
    /// Output path for trimmed files. Defaults to stdout if not provided
    pub fastq_output_file: Option<PathBuf>,

    #[arg(short = 'u', long)]
    /// Re-encode FASTQ sequence to only include canonical bases (A, C, T, G, N)
    pub canonical_bases: bool,

    #[arg(short = 'm', long)]
    /// Perform masking with 'N' instead of clipping. Default behavior is clipping if not provided
    pub mask: bool,

    #[arg(short = 'n', long, default_value = "1")]
    /// Minimum sequence length required after trimming. Shorter sequences are filtered from output. Must be greater than zero to ensure non-empty sequences
    pub min_length: NonZeroUsize,

    #[arg(short = 'H', long)]
    /// Hard trim from each end the specified number of bases
    pub hard_trim: Option<usize>,

    #[arg(long)]
    /// Hard trim range for only the left end of the sequence. Can be standalone, or override --hard_trim (-H)
    pub h_left: Option<usize>,

    #[arg(long)]
    /// Hard trim range for only the right end of the sequence. Can be standalone, or override --hard_trim (-H)
    pub h_right: Option<usize>,

    #[arg(short = 'B', long, value_parser = validate_non_empty, group = "adapter_vs_barcode")]
    /// Trim barcodes from sequence using fuzzy matching
    pub barcode_trim: Option<String>,

    #[arg(long, value_enum, default_value = "b")]
    /// Specifies the end of the sequence for barcode trimming : 'l' (left), 'r' (right), or 'b' (both)
    pub b_end: TrimEnd,

    // unsure if these are needed since we said we're doing full scan, but Illumina guarantees barcodes on a single side, so we might want this
    #[arg(long)]
    /// Restriction window size for barcode trimming on both ends of the sequence
    pub b_restrict: Option<usize>,

    #[arg(long)]
    /// Restriction window for trimming barcodes on the left end of the sequence. Can be standalone, or override --b_restrict
    pub b_restrict_left: Option<usize>,

    #[arg(long)]
    /// Restriction window for trimming barcodes on the right end of the sequence. Can be standalone, or override --b_restrict
    pub b_restrict_right: Option<usize>,

    #[arg(long, value_parser = validate_hamming_distance, default_value = "3")]
    /// Accepted Hamming distance for fuzzy barcode matching and trimming, between 0 and 3
    pub b_hdist: usize,

    #[arg(short = 'A', long, value_parser = validate_non_empty, group = "adapter_vs_barcode")]
    /// Trim a provided adapter from the sequence. Can be fuzzy (one mismatch) with `-f`
    pub adapter_trim: Option<String>,

    #[arg(long)]
    /// Allow up to one mismatch during adapter matching and trimming
    pub a_fuzzy_adapter: bool,

    #[arg(short = 'P', long)]
    /// Trim primers from provided primer fasta file
    pub primer_trim: Option<PathBuf>,

    #[arg(long, value_parser = validate_kmer_length, default_value = "17")]
    /// Length of k-mer used for matching primers.
    pub p_kmer_length: usize,

    #[arg(long, default_value = "b")]
    /// Specifies the end of the sequence for primer trimming : 'l' (left), 'r' (right), or 'b' (both)
    pub p_end: TrimEnd,

    #[arg(long)]
    /// Restriction window size for primer trimming on both ends of the sequence
    pub p_restrict: Option<usize>,

    #[arg(long)]
    /// Restriction window for trimming primer on the left end of the sequence. Can be standalone, or override --b_restrict
    pub p_restrict_left: Option<usize>,

    #[arg(long)]
    /// Restriction window for trimming barcodes on the right end of the sequence. Can be standalone, or override --b_restrict
    pub p_restrict_right: Option<usize>,
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

/// Custom validator for `kmer_length`.
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

fn validate_hamming_distance(value: &str) -> Result<usize, String> {
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
    pub canonical_bases:  bool,
    pub mask:             bool,
    pub min_length:       usize,
    pub barcode:          Option<Vec<u8>>,
    pub b_restrict_left:  usize,
    pub b_restrict_right: usize,
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
        .as_ref()
        .map(|adapter| get_forward_reverse_adapters(adapter));

    let barcode = args.barcode_trim.map(|b| b.into_bytes());

    let primer_kmers = if let Some(primer_path) = &args.primer_trim {
        Some(prepare_primer_kmers(primer_path, args.p_kmer_length)?)
    } else {
        None
    };

    // Default value for restriction based on size of ONT Barcode
    let default_b_restrict = args.b_restrict.unwrap_or(36);
    let (b_restrict_left, b_restrict_right) = match args.b_end {
        TrimEnd::B => {
            // Apply b_restrict to both, unless explicitly overridden
            let left = args.b_restrict_left.unwrap_or(default_b_restrict);
            let right = args.b_restrict_right.unwrap_or(default_b_restrict);
            (left, right)
        }
        TrimEnd::L => {
            // Left gets b_left_restrict or b_restrict, right is forced to 0
            let left = args.b_restrict_left.unwrap_or(default_b_restrict);
            (left, 0)
        }
        TrimEnd::R => {
            // Right gets b_right_restrict or b_restrict, left is forced to 0
            let right = args.b_restrict_right.unwrap_or(default_b_restrict);
            (0, right)
        }
    };

    let default_p_restrict = args.p_restrict.unwrap_or(30);
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
        canonical_bases: args.canonical_bases,
        mask: args.mask,
        min_length: args.min_length.get(),
        barcode,
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

//static MODULE: &str = module_path!();

/// Sub-program for trimming fastQ data.
pub fn trimmer_process(args: TrimmerArgs) -> Result<(), std::io::Error> {
    let mut args = parse_trim_args(args)?;
    let writer = &mut args.output_writer;

    let mut chained_reader = args.fastq_reader1.chain(args.fastq_reader2.into_iter().flatten());

    chained_reader.by_ref().try_for_each(|record| -> Result<(), std::io::Error> {
        let mut fq = record?;
        let mut fq_view = fq.as_view_mut();
        fq_view.to_canonical_bases(args.canonical_bases);

        if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
            fq_view.transform_by_reverse_forward_search(args.a_fuzzy_adapter, !args.mask, reverse_adapter, forward_adapter);
        } else if let Some(barcode) = &args.barcode {
            fq_view.barcode_trim(barcode, args.b_hdist, args.mask, args.b_restrict_left, args.b_restrict_right);
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

fn get_forward_reverse_adapters(adapter: &str) -> (Vec<u8>, Vec<u8>) {
    let forward = adapter.as_bytes().to_ascii_uppercase();
    let reverse = reverse_complement(&forward);
    (forward, reverse)
}

/// # Panics
///
/// `kmer_length` must be between 2 and 21, inclusive.
fn prepare_primer_kmers(
    primer_path: &PathBuf, kmer_length: usize,
) -> Result<ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, std::io::Error> {
    let fasta_primer_reader = FastaReader::from_filename(primer_path)?;

    let mut unique_kmers = ThreeBitKmerSet::<MAX_KMER_LENGTH, _>::with_hasher(kmer_length, RandomState::default())
        .expect("Expected valid kmer length");
    for f in fasta_primer_reader {
        let mut seq = Nucleotides::from_vec_unchecked(f?.sequence);
        if seq.len() > kmer_length {
            unique_kmers.insert_from_sequence_with_variants::<1>(&seq);
            seq.make_reverse_complement();
            unique_kmers.insert_from_sequence_with_variants::<1>(seq);
        }
    }
    Ok(unique_kmers)
}
