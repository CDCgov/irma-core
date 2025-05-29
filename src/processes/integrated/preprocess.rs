// Description:      Read FastQ files, applies QC filtering (quality and
//                   length), adapter/primer/barcode trimming, and deflates into
//                   a custom XFL format + FASTA format.

use crate::{
    args::clipping::{ClippingArgs, ParsedClippingArgs, parse_clipping_args},
    io::{ReadFileZip, open_fastq_file},
    qc::{
        fastq::{ReadTransforms, fix_sra_format, keep_or_underscore_header},
        fastq_metadata::*,
    },
    utils::{
        get_hasher,
        paired_reads::{PairedReadFilterer, ReadSide},
        trimming::trim_read,
    },
};
use clap::{Args, ValueHint};
use foldhash::fast::SeedableRandomState;
use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{BufWriter, prelude::*},
    num::NonZeroUsize,
    path::PathBuf,
};
use zoe::prelude::*;

type DeflatedSequences = HashMap<Nucleotides, Vec<(String, QualityScores)>, SeedableRandomState>;

#[derive(Args, Debug)]
pub struct PreprocessArgs {
    /// Location to store the XFL file.
    table_file: PathBuf,

    /// Single-ended FASTQ or the R1 file.
    fastq_input_file1: PathBuf,

    /// The R2 paired-end FASTQ file.
    fastq_input_file2: Option<PathBuf>,

    #[arg(short = 'L', long, value_hint = ValueHint::FilePath)]
    /// Quality control log path and filename.
    log_file: Option<PathBuf>,

    #[arg(short = 'K', long)]
    /// Keep the fastq header as usual.
    keep_header: bool,

    #[arg(short = 'T', long, default_value_t = 0)]
    /// Specify the read quality threshold (geometric mean, median).
    min_read_quality: u8,

    #[arg(short = 'M', long)]
    /// Interprets the threshold (-T) as the median, not the geometric mean.
    use_median: bool,

    #[arg(short = 'n', long, default_value = "1")]
    /// Minimum length of sequence read data, filtered otherwise.
    min_length: NonZeroUsize,

    #[arg(short = 'E', long, requires = "min_length")]
    /// The minimum length threshold (-n) is enforced after all trimming.
    enforce_clipped_length: bool,

    #[arg(short = 'f', long, requires = "fastq_input_file2")]
    /// Filter widowed reads
    filter_widows: bool,

    #[command(flatten)]
    clipping_args: ClippingArgs,
}

const CLUSTER_PREFIX: &str = "C";
static MODULE: &str = module_path!();

/// # Panics
///
/// Sub-program for processing FASTQ data.
pub fn preprocess_process(args: PreprocessArgs) -> Result<(), std::io::Error> {
    let ParsedPreprocessArgs { mut io_args, options } = parse_preprocess_args(args)?;

    let paired_reads = io_args.reader2.is_some();

    let (metadata_by_sequence, metadata) = trim_and_deflate(&options, &mut io_args)?;

    let read_pattern_count_passing = if metadata.passed_qc_count == 0 {
        diagnose_none_passing(&metadata, paired_reads, &options);
        0
    } else {
        output_deflated_sequences(metadata_by_sequence, io_args.table_writer)?
    };

    if let Some(log_writer) = io_args.log_writer
        && let Some(log_file) = io_args.log_file
    {
        write_log(
            log_writer,
            &metadata,
            paired_reads,
            read_pattern_count_passing,
            &options,
            log_file,
        )?;
    }

    Ok(())
}

struct ParsedPreprocessIoArgs {
    table_writer: BufWriter<File>,
    reader1:      FastQReader<ReadFileZip>,
    reader2:      Option<FastQReader<ReadFileZip>>,
    log_writer:   Option<BufWriter<File>>,
    log_file:     Option<PathBuf>,
}

#[derive(Debug)]
struct ParsedPreprocessOptions {
    keep_header:            bool,
    min_read_quality:       u8,
    use_median:             bool,
    min_length:             usize,
    enforce_clipped_length: bool,
    filter_widows:          bool,
    clipping_args:          ParsedClippingArgs,
}

struct ParsedPreprocessArgs {
    io_args: ParsedPreprocessIoArgs,
    options: ParsedPreprocessOptions,
}

fn parse_preprocess_args(args: PreprocessArgs) -> Result<ParsedPreprocessArgs, std::io::Error> {
    let PreprocessArgs {
        table_file,
        fastq_input_file1,
        fastq_input_file2,
        log_file,
        keep_header,
        min_read_quality,
        use_median,
        min_length,
        enforce_clipped_length,
        filter_widows,
        clipping_args,
    } = args;

    let reader1 = open_fastq_file(fastq_input_file1)?;

    let reader2 = match fastq_input_file2 {
        Some(file2) => Some(open_fastq_file(file2)?),
        None => None,
    };

    let log_writer = match log_file {
        Some(ref file_path) => Some(BufWriter::new(
            OpenOptions::new().write(true).create(true).truncate(true).open(file_path)?,
        )),
        None => None,
    };

    let table_writer = BufWriter::new(OpenOptions::new().write(true).create(true).truncate(true).open(table_file)?);

    let min_length = min_length.get();

    let clipping_args = parse_clipping_args(clipping_args)?;

    let parsed = ParsedPreprocessArgs {
        io_args: ParsedPreprocessIoArgs {
            table_writer,
            reader1,
            reader2,
            log_writer,
            log_file,
        },
        options: ParsedPreprocessOptions {
            keep_header,
            min_read_quality,
            use_median,
            min_length,
            enforce_clipped_length,
            filter_widows,
            clipping_args,
        },
    };

    Ok(parsed)
}

/// Trims all sequences, applies quality filtering, and deflates the sequences.
/// Returns the deflated sequences and the log file metadata.
fn trim_and_deflate(
    options: &ParsedPreprocessOptions, io_args: &mut ParsedPreprocessIoArgs,
) -> std::io::Result<(DeflatedSequences, FastQMetadata)> {
    let reader1 = &mut io_args.reader1;
    let mut trimming_specs = Preprocessor::new(options);

    if let Some(reader2) = &mut io_args.reader2 {
        if options.filter_widows {
            trimming_specs.handle_paired_reads_with_filter_weak(
                reader1,
                reader2,
                header_mismatch_warning,
                extra_read_warning,
            )?;
        } else {
            trimming_specs.handle_paired_reads_no_filter(reader1, reader2)?;
        }
    } else {
        trimming_specs.handle_single_reads(reader1, ReadSide::Unpaired)?;
    };

    Ok(trimming_specs.finalize())
}

/// Writes the table file to `table_writer` and the XFL file to STDOUT. The
/// number of read patterns is returned.
fn output_deflated_sequences(
    metadata_by_sequence: DeflatedSequences, mut table_writer: BufWriter<File>,
) -> Result<usize, std::io::Error> {
    let mut stdout_writer = BufWriter::new(std::io::stdout());

    let mut read_pattern_number = 0;
    for (sequence, metadata) in metadata_by_sequence {
        let cluster_size = metadata.len();

        writeln!(
            stdout_writer,
            ">{CLUSTER_PREFIX}{read_pattern_number}%{cluster_size}\n{sequence}"
        )?;

        write!(table_writer, "{CLUSTER_PREFIX}{read_pattern_number}%{cluster_size}")?;
        for (header, quality_scores) in metadata {
            write!(table_writer, "\t{header}\t{quality_scores}")?;
        }
        writeln!(table_writer)?;
        read_pattern_number += 1;
    }

    table_writer.flush()?;
    stdout_writer.flush()?;

    Ok(read_pattern_number)
}

/// Writes the log file.
fn write_log(
    mut log_writer: BufWriter<File>, metadata: &FastQMetadata, paired_reads: bool, read_pattern_count_passing: usize,
    options: &ParsedPreprocessOptions, log_file: PathBuf,
) -> Result<(), std::io::Error> {
    let FastQMetadata {
        passed_qc_count,
        passed_len_count,
        observed_q_max,
        observed_raw_reads,
        observed_max_read_len,
        observed_max_clipped_read_len,
    } = metadata;

    writeln!(
        log_writer,
        "\n\
        NUMBER_INPUT_FILES\t{num_files}\n\
        OBSERVED_RAW_READS_OR_R1\t{r1_raw_reads}\n\
        OBSERVED_R2_READS\t{r2_raw_reads}\n\
        OBSERVED_MAX_READ_LEN\t{observed_max_read_len}\n\
        OBSERVED_MAX_CLIPPED_READ_LENGTH\t{observed_max_clipped_read_len}\n\
        OBSERVED_MAX_QUALITY\t{observed_q_max}\n\
        READ_COUNT_PASSING_ONLY_LENGTH_FILTER\t{passed_len_count}\n\
        READ_COUNT_PASSING_ALL_QUALITY_CONTROL_FILTERS\t{passed_qc_count}\n\
        READ_PATTERN_COUNT_PASSING\t{read_pattern_count_passing}\n\
        MIN_PHRED_QUALITY_THRESHOLD\t{min_read_quality}\n\
        MIN_READ_LENGTH_THRESHOLD\t{min_length}\n\
        QUALITY_MEASURE\t{center_type}\n\
        ",
        num_files = if paired_reads { 2 } else { 1 },
        r1_raw_reads = observed_raw_reads[0],
        r2_raw_reads = observed_raw_reads[1],
        observed_q_max = observed_q_max.map(|q| q.to_string()).unwrap_or_else(|| "NONE".to_string()),
        min_read_quality = options.min_read_quality,
        min_length = options.min_length,
        center_type = if options.use_median { "median" } else { "average" },
    )
    .unwrap_or_else(|e| {
        eprintln!("{MODULE} WARNING! Cannot write to {}. See: {e}", log_file.display());
    });

    Ok(())
}

/// Attempt to diagnose the problem when no reads pass all quality filters.
/// Warnings are printed to STDERR.
fn diagnose_none_passing(metadata: &FastQMetadata, paired_reads: bool, options: &ParsedPreprocessOptions) {
    match (metadata.observed_raw_reads[0], metadata.observed_raw_reads[1], paired_reads) {
        (0, _, false) => {
            eprintln!("{MODULE} WARNING! No reads were found in the input file.");
            return;
        }
        (0, 0, true) => {
            eprintln!("{MODULE} WARNING! No reads were found in either input file.");
            return;
        }
        (0, _, true) => {
            eprintln!("{MODULE} WARNING! No reads were found in the first input file.");
            return;
        }
        (_, 0, true) => {
            eprintln!("{MODULE} WARNING! No reads were found in the second input file.");
            return;
        }
        _ => {}
    }

    if let Some(obs_max) = metadata.observed_q_max
        && obs_max < f32::from(options.min_read_quality)
    {
        eprintln!(
            "{MODULE} WARNING! The observed max phred quality score ({obs_max}) is below the user specified threshold (QUAL_THRESHOLD = {}).",
            options.min_read_quality
        );
    }

    if metadata.observed_max_read_len < options.min_length {
        eprintln!(
            "{MODULE} WARNING! The observed max read length ({}) is below the user specified threshold (MIN_LEN = {}).",
            metadata.observed_max_read_len, options.min_length
        );
    }
}

#[inline]
#[must_use]
fn header_mismatch_warning(e: std::io::Error) -> String {
    format!(
        "{MODULE} WARNING! {e} `--filter-widows or -f` is being disabled for the remainder of the processing. Consider rerunning with corrected inputs."
    )
}

#[inline]
#[must_use]
fn extra_read_warning() -> String {
    format!(
        "{MODULE} WARNING! Extra unpaired read(s) found at end of first FASTQ file. `--filter-widows or -f` is being disabled for the remainder of the processing. Consider rerunning with corrected inputs."
    )
}

/// A [`PairedReadFilterer`] struct used by the preprocess process to deflate
/// and store trimmed/qc-filtered reads.
#[derive(Debug)]
struct Preprocessor<'a> {
    args:                 &'a ParsedPreprocessOptions,
    metadata_by_sequence: DeflatedSequences,
    metadata:             FastQMetadata,
}

impl<'a> Preprocessor<'a> {
    #[inline]
    fn new(args: &'a ParsedPreprocessOptions) -> Self {
        Self {
            args,
            metadata_by_sequence: HashMap::with_hasher(get_hasher()),
            metadata: FastQMetadata::default(),
        }
    }
}

struct TrimmedRead<'a> {
    sequence: NucleotidesViewMut<'a>,
    quality:  QualityScoresViewMut<'a>,
    header:   String,
}

impl<'a> TrimmedRead<'a> {
    #[inline]
    fn new(fastq: FastQViewMut<'a>) -> Self {
        Self {
            sequence: fastq.sequence,
            quality:  fastq.quality,
            header:   std::mem::take(fastq.header),
        }
    }
}

impl<'a> PairedReadFilterer for Preprocessor<'a> {
    const PRESERVE_ORDER: bool = false;
    type Processed<'b> = TrimmedRead<'b>;
    type Finalized = (DeflatedSequences, FastQMetadata);

    #[inline]
    fn process_read<'b>(&mut self, read: &'b mut FastQ, side: ReadSide) -> Option<Self::Processed<'b>> {
        self.metadata.observed_raw_reads += side.to_simd();
        self.metadata.observed_max_read_len = self.metadata.observed_max_read_len.max(read.sequence.len());
        if read.sequence.len() < self.args.min_length {
            return None;
        }

        let clipped = trim_read(read.as_view_mut(), false, &self.args.clipping_args);

        self.metadata.observed_max_clipped_read_len =
            self.metadata.observed_max_clipped_read_len.max(clipped.sequence.len());
        if (self.args.enforce_clipped_length && clipped.sequence.len() < self.args.min_length) || clipped.sequence.is_empty()
        {
            return None;
        }
        self.metadata.passed_len_count += 1;

        let read_q_center = clipped.get_q_center(self.args.use_median);
        self.metadata.observed_q_max = if read_q_center > self.metadata.observed_q_max {
            read_q_center
        } else {
            self.metadata.observed_q_max
        };
        if read_q_center < Some(f32::from(self.args.min_read_quality)) {
            return None;
        }

        self.metadata.passed_qc_count += 1;

        Some(TrimmedRead::new(clipped))
    }

    #[inline]
    fn output_read<'b>(&mut self, mut trimmed: Self::Processed<'b>, side: ReadSide) -> std::io::Result<()> {
        let header = {
            if let Some(side) = side.to_char() {
                trimmed.header = fix_sra_format(trimmed.header, side);
            }
            keep_or_underscore_header(trimmed.header, self.args.keep_header)
        };

        let sequence = trimmed.sequence.to_owned_data();
        let quality = trimmed.quality.to_owned_data();

        self.metadata_by_sequence.entry(sequence).or_default().push((header, quality));

        Ok(())
    }

    #[inline]
    fn finalize(&mut self) -> Self::Finalized {
        (
            std::mem::take(&mut self.metadata_by_sequence),
            std::mem::take(&mut self.metadata),
        )
    }
}
