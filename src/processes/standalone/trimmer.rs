use crate::{
    qc::fastq::ReadTransforms,
    utils::{
        get_hasher, get_molecular_id_side,
        io::{ReadFileZip, WriteFileZipStdout, create_writer, open_fastq_file},
    },
};
use clap::{Args, ValueEnum, builder::PossibleValue};
use foldhash::fast::SeedableRandomState;
use std::{
    io::{Error as IOError, ErrorKind, Write},
    num::NonZeroUsize,
    path::PathBuf,
};
use zoe::{kmer::ThreeBitKmerSet, prelude::*};

pub const MAX_KMER_LENGTH: usize = 21;

#[derive(Args, Debug)]
pub struct TrimmerArgs {
    /// Path to .fastq or .fastq.gz file to be trimmed
    pub fastq_input_file: PathBuf,

    /// Path to optional second .fastq or .fastq.gz file to be trimmed
    pub fastq_input_file2: Option<PathBuf>,

    #[arg(short = '1', short_alias = 'o', long = "fastq-output")]
    /// Output filepath for trimmed reads. Trimmed reads print to STDOUT if not
    /// provided. May also use '-o'.
    pub fastq_output_file: Option<PathBuf>,

    #[arg(short = '2', long = "fastq-output2", requires = "fastq_input_file2")]
    /// Output path for secondary trimmed file if using paired reads. If this
    /// argument is omitted, output is interleaved.
    pub fastq_output_file2: Option<PathBuf>,

    #[arg(short = 'm', long)]
    /// Perform masking with 'N' instead of clipping. Default behavior is
    /// clipping if not provided
    pub mask: bool,

    #[arg(short = 'n', long, default_value = "1")]
    /// Minimum sequence length required after trimming. Shorter sequences are
    /// filtered from output.
    pub min_length: NonZeroUsize,

    #[arg(short = 'f', long)]
    /// Filter widowed reads
    pub filter_widows: bool,

    #[command(flatten)]
    pub clipping_args: ClippingArgs,
}

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

/// Sub-program for trimming FASTQ data.
pub fn trimmer_process(args: TrimmerArgs) -> Result<(), std::io::Error> {
    let ParsedTrimmerArgs {
        io_args:
            ParsedPairedIoArgs {
                mut reader1,
                reader2,
                mut writer1,
                writer2,
            },
        trimming_args,
    } = parse_trimmer_args(args)?;

    if let Some(mut reader2) = reader2 {
        match (writer2, trimming_args.filter_widows) {
            // Case 2: In 1, In 2, Out 1 (interleaved Illumina), no widow filtering
            (None, false) => {
                for r1 in reader1.by_ref() {
                    let mut read1 = r1?;
                    if let Some(trimmed1) = process_read(&mut read1, &trimming_args) {
                        write!(writer1, "{trimmed1}")?;
                    }

                    let Some(r2) = reader2.next() else {
                        break;
                    };
                    let mut read2 = r2?;

                    if let Some(trimmed2) = process_read(&mut read2, &trimming_args) {
                        write!(writer1, "{trimmed2}")?;
                    }
                }

                // If either reader has remaining sequences, this will trim them. (same as 1 in, 1 out)
                let mut remaining = reader1.chain(reader2);
                process_and_write(&mut remaining, &mut writer1, &trimming_args)?;

                writer1.flush()?;
            }

            // Case 3: In 1, In 2, Out 1, Filtering widows / orphan reads
            (None, true) => {
                for r1 in reader1 {
                    let mut read1 = r1?;
                    let Some(r2) = reader2.next() else {
                        return error_extra_read();
                    };
                    let mut read2 = r2?;
                    check_paired_headers(&read1, &read2)?;

                    if let Some((trimmed1, trimmed2)) = process_paired_reads(&mut read1, &mut read2, &trimming_args) {
                        write!(writer1, "{trimmed1}")?;
                        write!(writer1, "{trimmed2}")?;
                    }
                }

                writer1.flush()?;

                if reader2.next().is_some() {
                    return error_extra_read();
                }
            }

            // Case 4: In 1, In 2, Out 1, Out 2 (separated output Illumina), no filtering
            (Some(mut writer2), false) => {
                process_and_write(&mut reader1, &mut writer1, &trimming_args)?;
                process_and_write(&mut reader2, &mut writer2, &trimming_args)?;

                writer1.flush()?;
                writer2.flush()?;
            }

            // Case 5: In 1, In 2, Out 1, Out 2, filter widows
            (Some(mut writer2), true) => {
                for r1 in reader1 {
                    let mut read1 = r1?;
                    let Some(r2) = reader2.next() else {
                        return error_extra_read();
                    };
                    let mut read2 = r2?;
                    check_paired_headers(&read1, &read2)?;

                    if let Some((trimmed1, trimmed2)) = process_paired_reads(&mut read1, &mut read2, &trimming_args) {
                        write!(writer1, "{trimmed1}")?;
                        write!(writer2, "{trimmed2}")?;
                    }
                }

                writer1.flush()?;
                writer2.flush()?;

                if reader2.next().is_some() {
                    return error_extra_read();
                }
            }
        }
    } else {
        // Case 1: In 1, Out 1 (ONT, single-end, PacBio)
        process_and_write(&mut reader1, &mut writer1, &trimming_args)?;
        writer1.flush()?;
    }

    Ok(())
}

/// Processes an iterator of reads and writes it to output
fn process_and_write<I, W>(reads: &mut I, writer: &mut W, args: &ParsedTrimmingArgs) -> Result<(), std::io::Error>
where
    I: Iterator<Item = Result<FastQ, std::io::Error>>,
    W: Write, {
    reads.try_for_each(|record_res| {
        let mut record = record_res?;
        if let Some(trimmed) = process_read(&mut record, args) {
            write!(writer, "{trimmed}")?;
        }
        Ok(())
    })
}

/// Processes two reads together with filtering of widowed reads
fn process_paired_reads<'a>(
    record1: &'a mut FastQ, record2: &'a mut FastQ, args: &ParsedTrimmingArgs,
) -> Option<(FastQViewMut<'a>, FastQViewMut<'a>)> {
    let trimmed1 = process_read(record1, args)?;
    let trimmed2 = process_read(record2, args)?;
    Some((trimmed1, trimmed2))
}

/// Processes a single read and performs length filtering. Returns `None` if the
/// read is too short
fn process_read<'a>(record: &'a mut FastQ, args: &ParsedTrimmingArgs) -> Option<FastQViewMut<'a>> {
    if args.mask {
        let fq_view = record.as_view_mut();
        edit_read(fq_view, args.mask, &args.clipping_args);
        if record.len() >= args.min_length {
            return Some(record.as_view_mut());
        }
    } else {
        let fq_view = record.as_view_mut();
        let edited = edit_read(fq_view, args.mask, &args.clipping_args);
        if edited.len() >= args.min_length {
            return Some(edited);
        }
    }
    None
}

/// Trims or masks a read based on user provided arguments
pub fn edit_read<'a>(mut fq_view: FastQViewMut<'a>, mask: bool, args: &ParsedClippingArgs) -> FastQViewMut<'a> {
    fq_view.to_canonical_bases(!args.preserve_bases);

    fq_view.process_polyg(args.polyg_left, args.polyg_right, mask);

    if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
        fq_view.transform_by_reverse_forward_search(
            args.a_fuzzy,
            !mask,
            reverse_adapter.as_bytes(),
            forward_adapter.as_bytes(),
        );
    } else if let Some((barcode, reverse)) = &args.barcodes {
        fq_view.process_barcode(
            barcode.as_bytes(),
            reverse.as_bytes(),
            args.b_hdist,
            mask,
            args.b_restrict_left,
            args.b_restrict_right,
        );
    }

    if let Some(ref kmers) = args.primer_kmers {
        if let Some(p_restrict_left) = args.p_restrict_left {
            fq_view.process_left_primer(p_restrict_left, kmers, mask);
        }
        if let Some(p_restrict_right) = args.p_restrict_right {
            fq_view.process_right_primer(p_restrict_right, kmers, mask);
        }
    }

    if args.hard_left > 0 || args.hard_right > 0 {
        fq_view.hard_clip_or_mask(args.hard_left, args.hard_right, mask);
    }
    fq_view
}

/// Returns whether two reads have matching molecular IDs. Errors if the read
/// ID's don't match or can't be parsed.
pub fn check_paired_headers(read1: &FastQ, read2: &FastQ) -> Result<(), std::io::Error> {
    if let Some((id1, _)) = get_molecular_id_side(&read1.header, '0')
        && let Some((id2, _)) = get_molecular_id_side(&read2.header, '0')
    {
        if id1 == id2 {
            Ok(())
        } else {
            Err(IOError::new(
                ErrorKind::InvalidInput,
                format!(
                    "Paired read IDs out of sync:\n\t{h1}\n\t{h2}\n",
                    h1 = read1.header,
                    h2 = read2.header
                ),
            ))
        }
    } else {
        Err(IOError::new(ErrorKind::InvalidInput, "Could not parse the read IDs."))
    }
}

#[inline]
fn error_extra_read() -> Result<(), std::io::Error> {
    Err(std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        "Extra unpaired read(s) found at end of one of the input FASTQ files.",
    ))
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

/// Parsed arguments for the `trimmer` subprocess
struct ParsedTrimmerArgs {
    pub io_args:       ParsedPairedIoArgs,
    pub trimming_args: ParsedTrimmingArgs,
}

/// Parsed IO arguments for single or paired reads
struct ParsedPairedIoArgs {
    reader1: FastQReader<ReadFileZip>,
    reader2: Option<FastQReader<ReadFileZip>>,
    writer1: WriteFileZipStdout,
    writer2: Option<WriteFileZipStdout>,
}

/// Arguments related to clipping/masking reads, including length/widow
/// filtering
pub struct ParsedTrimmingArgs {
    pub mask:          bool,
    pub filter_widows: bool,
    pub min_length:    usize,
    pub clipping_args: ParsedClippingArgs,
}

/// Arguments specifying the types of clipping to be performed
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

fn parse_trimmer_args(args: TrimmerArgs) -> Result<ParsedTrimmerArgs, std::io::Error> {
    let TrimmerArgs {
        fastq_input_file,
        fastq_input_file2,
        fastq_output_file,
        fastq_output_file2,
        mask,
        filter_widows,
        min_length,
        clipping_args,
    } = args;

    let reader1 = open_fastq_file(fastq_input_file)?;
    let reader2 = match fastq_input_file2 {
        Some(file2) => Some(open_fastq_file(file2)?),
        None => None,
    };

    let writer1 = create_writer(fastq_output_file.clone())?;
    let writer2 = match fastq_output_file2 {
        Some(path) => Some(create_writer(Some(path))?),
        None => None,
    };

    let min_length = min_length.get();

    let clipping_args = parse_clipping_args(clipping_args)?;

    let parsed = ParsedTrimmerArgs {
        io_args:       ParsedPairedIoArgs {
            reader1,
            reader2,
            writer1,
            writer2,
        },
        trimming_args: ParsedTrimmingArgs {
            mask,
            filter_widows,
            min_length,
            clipping_args,
        },
    };

    Ok(parsed)
}

pub fn parse_clipping_args(args: ClippingArgs) -> Result<ParsedClippingArgs, std::io::Error> {
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

/// For converting adapters and barcodes to uppercase (unless disabled) and
/// getting reverse complements
pub fn get_forward_reverse_sequence(mut adapter: Nucleotides, preserve_seq: bool) -> (Nucleotides, Nucleotides) {
    if !preserve_seq {
        adapter.as_mut_bytes().make_ascii_uppercase();
    }
    let reverse = adapter.to_reverse_complement();
    (adapter, reverse)
}

/// # Panics
///
/// `kmer_length` must be between 2 and 21, inclusive.
pub fn prepare_primer_kmers(
    primer_path: &PathBuf, kmer_length: usize, fuzzy_kmer: bool,
) -> Result<ThreeBitKmerSet<MAX_KMER_LENGTH, SeedableRandomState>, std::io::Error> {
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
