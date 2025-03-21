// Description:      Read FastQ files and trim with various options.
use crate::{qc::fastq::ReadTransforms, utils::get_hasher, utils::get_molecular_id_side};
use clap::{Args, ValueEnum, builder::PossibleValue};
use either::Either;
use foldhash::fast::SeedableRandomState;
use std::{
    fs::{File, OpenOptions},
    io::{BufReader, BufWriter, Stdout, Write, stdout},
    num::NonZeroUsize,
    path::PathBuf,
};
use zoe::{data::nucleotides::CheckNucleotides, kmer::ThreeBitKmerSet, prelude::*};

const MAX_KMER_LENGTH: usize = 21;

#[derive(Args, Debug)]
pub struct TrimmerArgs {
    /// Path to FASTQ file to be trimmed
    pub fastq_input_file: PathBuf,

    /// Path to optional second FASTQ file to be trimmed
    pub fastq_input_file2: Option<PathBuf>,

    #[arg(short = 'o', long)]
    /// Output filepath for trimmed reads. Trimmed reads print to stdout if not
    /// provided
    pub fastq_output_file: Option<PathBuf>,

    #[arg(short = 'u', long, requires = "fastq_input_file2")]
    /// Output path for secondary trimmed file if using paired reads. Trimmed
    /// sequences from the second input fastq will be interleaved with the first
    /// if this argument is not given
    pub fastq_output_file2: Option<PathBuf>,

    #[arg(short = 's', long)]
    /// Preserve original formatting: disables encoding sequences into uppercase
    /// canonical bases (A, C, T, G, N)
    pub preserve_fastq: bool,

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
    let (mut io_args, trimmer_args) = parse_trim_args(args)?;
    let writer = &mut io_args.output_writer;

    match (io_args.fastq_reader2, io_args.output_writer2, trimmer_args.filter_widows) {
        // Case 1: In 1, Out 1 (ONT)
        (None, None, _) => {
            process_and_write(&mut io_args.fastq_reader1, writer, &trimmer_args)?;
            writer.flush()?;
        }

        // Case 2: In 1, In 2, Out 1 (interleaved Illumina), no widow filtering
        (Some(mut fastq_reader2), None, false) => {
            let mut zipped_reader = io_args.fastq_reader1.by_ref().zip(fastq_reader2.by_ref());
            zipped_reader.try_for_each(|(r1, r2)| -> Result<(), std::io::Error> {
                let mut read1 = r1?;
                let mut read2 = r2?;

                if let Some(trimmed1) = process_read(&mut read1, &trimmer_args) {
                    write!(writer, "{trimmed1}")?;
                }
                if let Some(trimmed2) = process_read(&mut read2, &trimmer_args) {
                    write!(writer, "{trimmed2}")?;
                }

                Ok(())
            })?;
            // if either reader has remaining sequences, this will trim them. (same as 1 in, 1 out)
            let mut remaining = io_args.fastq_reader1.chain(fastq_reader2);
            process_and_write(&mut remaining, writer, &trimmer_args)?;

            writer.flush()?;
        }

        // Case 3: In 1, In 2, Out 1, Filtering widows
        (Some(mut fastq_reader2), None, true) => {
            let mut zipped_reader = io_args.fastq_reader1.by_ref().zip(fastq_reader2.by_ref());
            zipped_reader.try_for_each(|(r1, r2)| -> Result<(), std::io::Error> {
                let mut read1 = r1?;
                let mut read2 = r2?;

                if !check_paired_headers(&read1, &read2) {
                    writer.flush()?;
                    panic!("A mismatch in read IDs was found. Terminating.");
                }

                if let Some((trimmed1, trimmed2)) = process_paired_reads(&mut read1, &mut read2, &trimmer_args) {
                    write!(writer, "{trimmed1}")?;
                    write!(writer, "{trimmed2}")?;
                }

                Ok(())
            })?;
            writer.flush()?;
            if io_args.fastq_reader1.next().is_some() || fastq_reader2.next().is_some() {
                panic!("Extra read found at end of one of the input FASTQ files. Terminating");
            }
        }

        // Case 4: In 1, In 2, Out 1, Out 2 (separated output Illumina), no filtering
        (Some(mut fastq_reader2), Some(output_writer2), false) => {
            let mut writer2 = output_writer2;
            process_and_write(&mut io_args.fastq_reader1, writer, &trimmer_args)?;
            process_and_write(&mut fastq_reader2, &mut writer2, &trimmer_args)?;

            writer.flush()?;
            writer2.flush()?;
        }

        // Case 5: In 1, In 2, Out 1, Out 2, filter widows
        (Some(mut fastq_reader2), Some(mut writer2), true) => {
            let mut zipped_reader = io_args.fastq_reader1.by_ref().zip(fastq_reader2.by_ref());
            zipped_reader.try_for_each(|(r1, r2)| -> Result<(), std::io::Error> {
                let mut read1 = r1?;
                let mut read2 = r2?;

                // terminates early if a mismatch is found
                if !check_paired_headers(&read1, &read2) {
                    writer.flush()?;
                    writer2.flush()?;
                    panic!("A mismatch in read IDs was found. Terminating.");
                }

                if let Some((trimmed1, trimmed2)) = process_paired_reads(&mut read1, &mut read2, &trimmer_args) {
                    write!(writer, "{trimmed1}")?;
                    write!(writer2, "{trimmed2}")?;
                }

                Ok(())
            })?;
            writer.flush()?;
            writer2.flush()?;
            if io_args.fastq_reader1.next().is_some() || fastq_reader2.next().is_some() {
                panic!("Extra read found at end of one of the input FASTQ files. Terminating");
            }
        }

        _ => {
            unreachable!();
        }
    }

    Ok(())
}

/// Processes an iterator of reads and writes it to output
pub fn process_and_write<I, W>(reads: &mut I, writer: &mut W, args: &ParsedTrimmerArgs) -> Result<(), std::io::Error>
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
pub fn process_paired_reads<'a>(
    record1: &'a mut FastQ, record2: &'a mut FastQ, args: &ParsedTrimmerArgs,
) -> Option<(FastQViewMut<'a>, FastQViewMut<'a>)> {
    let trimmed1 = process_read(record1, args)?;
    let trimmed2 = process_read(record2, args)?;
    Some((trimmed1, trimmed2))
}

/// Processes a single read and performs length filtering. Returns `None` if the
/// read is too short
pub fn process_read<'a>(record: &'a mut FastQ, args: &ParsedTrimmerArgs) -> Option<FastQViewMut<'a>> {
    if args.mask {
        let fq_view = record.as_view_mut();
        edit_read(fq_view, args);
        if record.len() >= args.min_length {
            return Some(record.as_view_mut());
        }
    } else {
        let fq_view = record.as_view_mut();
        let edited = edit_read(fq_view, args);
        if edited.len() >= args.min_length {
            return Some(edited);
        }
    }
    None
}

/// Trims or masks a read based on user provided arguments
pub fn edit_read<'a>(mut fq_view: FastQViewMut<'a>, args: &ParsedTrimmerArgs) -> FastQViewMut<'a> {
    fq_view.to_canonical_bases(!args.preserve_seq);

    fq_view.process_polyg(args.polyg_left, args.polyg_right, args.mask);

    if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
        fq_view.transform_by_reverse_forward_search(
            args.a_fuzzy,
            !args.mask,
            reverse_adapter.as_bytes(),
            forward_adapter.as_bytes(),
        );
    } else if let Some((barcode, reverse)) = &args.barcodes {
        fq_view.process_barcode(
            barcode.as_bytes(),
            reverse.as_bytes(),
            args.b_hdist,
            args.mask,
            args.b_restrict_left,
            args.b_restrict_right,
        );
    }

    if let Some(ref kmers) = args.primer_kmers {
        if let Some(p_restrict_left) = args.p_restrict_left {
            fq_view.process_left_primer(p_restrict_left, kmers, args.mask);
        }
        if let Some(p_restrict_right) = args.p_restrict_right {
            fq_view.process_right_primer(p_restrict_right, kmers, args.mask);
        }
    }

    if args.hard_left > 0 || args.hard_right > 0 {
        fq_view.hard_clip_or_mask(args.hard_left, args.hard_right, args.mask);
    }
    fq_view
}

/// Returns whether two reads have matching molecular IDs. If the IDs cannot be
/// parsed, then `false` is returned.
pub fn check_paired_headers(read1: &FastQ, read2: &FastQ) -> bool {
    let h1 = &read1.header;
    let h2 = &read2.header;
    if let (Some((id1, _)), Some((id2, _))) = (get_molecular_id_side(h1, '0'), get_molecular_id_side(h2, '1')) {
        return id1 == id2;
    }
    false
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
        .map_err(|_| format!("`{}` is not a valid integer.", value))?;
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
        .map_err(|_| format!("`{}` is not a valid integer between 0 and 3.", value))?;
    if (0..=3).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!("b-hdist must be between 0 and 3, but `{}` was provided.", parsed))
    }
}

/// Ensures user has entered valid non-empty adapter or barcode literal for
/// trimming
fn validate_acgtn(value: &str) -> Result<Nucleotides, String> {
    if value.trim().is_empty() {
        // prevents panicking when `-A ""` is passed
        Err("Adapter (-A) or barcode (-B) cannot be empty!".to_string())
    } else {
        let forward = Nucleotides::from(value.as_bytes());
        if forward.is_acgtn() {
            Ok(forward)
        } else {
            Err("Adapter or barcode literal must only consist of canonical (ACGTN) bases".to_string())
        }
    }
}

pub struct IOArgs {
    pub fastq_reader1:  FastQReader<BufReader<File>>,
    pub fastq_reader2:  Option<FastQReader<BufReader<File>>>,
    pub output_writer:  BufWriter<Either<File, Stdout>>,
    pub output_writer2: Option<BufWriter<File>>,
}

pub struct ParsedTrimmerArgs {
    pub filter_widows:    bool,
    pub preserve_seq:     bool,
    pub mask:             bool,
    pub min_length:       usize,
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

pub fn parse_trim_args(args: TrimmerArgs) -> Result<(IOArgs, ParsedTrimmerArgs), std::io::Error> {
    let fastq_reader1 = FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(&args.fastq_input_file)?));

    let fastq_reader2 = if let Some(file2) = &args.fastq_input_file2 {
        Some(FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(file2)?)))
    } else {
        None
    };

    let output_writer = match &args.fastq_output_file {
        Some(path) => BufWriter::new(Either::Left(File::create(path)?)),
        None => BufWriter::new(Either::Right(stdout())),
    };

    let output_writer2 = match &args.fastq_output_file2 {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };

    let adapters = args
        .adapter_trim
        .map(|adapter| get_forward_reverse_sequence(adapter, args.preserve_fastq));

    let barcodes = args
        .barcode_trim
        .map(|barcode| get_forward_reverse_sequence(barcode, args.preserve_fastq));

    let primer_kmers = if let Some(primer_path) = &args.primer_trim {
        Some(prepare_primer_kmers(
            primer_path,
            args.p_kmer_length
                // this is unreachable through clap due to being required
                .expect("A kmer length must be provided for primer trimming"),
            args.p_fuzzy,
        )?)
    } else {
        None
    };

    // A value of None for left or right restricts will do full scan barcoding
    let default_b_restrict = args.b_restrict;
    let (b_restrict_left, b_restrict_right) = match args.b_end {
        TrimEnd::B => {
            // Apply b_restrict to both, unless explicitly overridden
            let left = args.b_restrict_left.or(default_b_restrict).map(NonZeroUsize::get);
            let right = args.b_restrict_right.or(default_b_restrict).map(NonZeroUsize::get);
            (left, right)
        }
        TrimEnd::L => {
            // Left gets b_left_restrict or b_restrict, right is forced to 0
            let left = args.b_restrict_left.or(default_b_restrict).map(NonZeroUsize::get);
            (left, Some(0))
        }
        TrimEnd::R => {
            // Right gets b_right_restrict or b_restrict, left is forced to 0
            let right = args.b_restrict_right.or(default_b_restrict).map(NonZeroUsize::get);
            (Some(0), right)
        }
    };

    let default_p_restrict = args.p_restrict;
    let (p_restrict_left, p_restrict_right) = match args.p_end {
        TrimEnd::B => {
            let left = args.p_restrict_left.unwrap_or(default_p_restrict);
            let right = args.p_restrict_right.unwrap_or(default_p_restrict);
            (Some(left.get()), Some(right.get()))
        }
        TrimEnd::L => {
            let left = args.p_restrict_left.unwrap_or(default_p_restrict);
            (Some(left.get()), None)
        }
        TrimEnd::R => {
            let right = args.p_restrict_right.unwrap_or(default_p_restrict);
            (None, Some(right.get()))
        }
    };

    let default_polyg = args.polyg_trim;
    let (polyg_left, polyg_right) = match args.g_polyg_end {
        TrimEnd::B => {
            let left = args.g_polyg_left.or(default_polyg).map(NonZeroUsize::get);
            let right = args.g_polyg_right.or(default_polyg).map(NonZeroUsize::get);
            (left, right)
        }
        TrimEnd::L => {
            let left = args.g_polyg_left.or(default_polyg).map(NonZeroUsize::get);
            (left, None)
        }
        TrimEnd::R => {
            let right = args.g_polyg_right.or(default_polyg).map(NonZeroUsize::get);
            (None, right)
        }
    };

    let default_hard_bases = args.hard_trim.unwrap_or(0);
    let hard_left = args.h_left.unwrap_or(default_hard_bases);
    let hard_right = args.h_right.unwrap_or(default_hard_bases);

    let io_args = IOArgs {
        fastq_reader1,
        fastq_reader2,
        output_writer,
        output_writer2,
    };

    let parsed_args = ParsedTrimmerArgs {
        filter_widows: args.filter_widows,
        preserve_seq: args.preserve_fastq,
        mask: args.mask,
        min_length: args.min_length.get(),
        barcodes,
        b_restrict_left,
        b_restrict_right,
        b_hdist: args.b_hdist,
        adapters,
        a_fuzzy: args.a_fuzzy,
        primer_kmers,
        p_restrict_left,
        p_restrict_right,
        polyg_left,
        polyg_right,
        hard_left,
        hard_right,
    };
    Ok((io_args, parsed_args))
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

/// # Panics
///
/// `kmer_length` must be between 2 and 21, inclusive.
fn prepare_primer_kmers(
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
