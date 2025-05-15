// Description:      Read FastQ files, applies QC filtering (quality and length),
//                   adapter trimming, and deflates into a custom XFL format + FASTA format.

use crate::{
    qc::{fastq::*, fastq_metadata::*},
    trimmer::{ClippingArgs, ParsedClippingArgs, edit_read, parse_clipping_args},
    utils::get_hasher,
};
use clap::{Args, ValueHint};
use std::{
    collections::HashMap,
    fs::{File, OpenOptions},
    io::{BufReader, BufWriter, prelude::*},
    num::NonZeroUsize,
    path::PathBuf,
    simd::prelude::*,
};
use zoe::prelude::*;

#[derive(Args, Debug)]
pub struct PreprocessArgs {
    /// Location to store the XFL file.
    pub table_file: PathBuf,

    /// Single-ended FASTQ or the R1 file.
    pub fastq_input_file1: PathBuf,

    /// The R2 paired-end FASTQ file.
    pub fastq_input_file2: Option<PathBuf>,

    #[arg(short = 'L', long, value_hint = ValueHint::FilePath)]
    /// Quality control log path and filename.
    pub log_file: Option<PathBuf>,

    #[arg(short = 'K', long)]
    /// Keep the fastq header as usual.
    pub keep_header: bool,

    #[arg(short = 'T', long, default_value_t = 0)]
    /// Specify the read quality threshold (geometric mean, median).
    pub min_read_quality: u8,

    #[arg(short = 'M', long)]
    /// Interprets the threshold (-T) as the median, not the geometric mean.
    pub use_median: bool,

    #[arg(short = 'n', long, default_value = "1")]
    /// Minimum length of sequence read data, filtered otherwise.
    pub min_length: NonZeroUsize,

    #[arg(short = 'E', long)]
    /// The minimum length threshold (-n) is enforced after all trimming.
    pub enforce_clipped_length: bool,

    #[command(flatten)]
    pub clipping_args: ClippingArgs,
}

const CLUSTER_PREFIX: &str = "C";
static MODULE: &str = module_path!();

pub struct ParsedPreprocessIoArgs {
    pub table_writer: BufWriter<File>,
    pub reader1:      FastQReader<BufReader<File>>,
    pub reader2:      Option<FastQReader<BufReader<File>>>,
    pub log_writer:   Option<BufWriter<File>>,
    pub log_file:     Option<PathBuf>,
}

pub struct ParsedPreprocessArgs {
    pub io_args:                ParsedPreprocessIoArgs,
    pub keep_header:            bool,
    pub min_read_quality:       u8,
    pub use_median:             bool,
    pub min_length:             usize,
    pub enforce_clipped_length: bool,
    pub clipping_args:          ParsedClippingArgs,
}

pub fn parse_preprocess_args(args: PreprocessArgs) -> Result<ParsedPreprocessArgs, std::io::Error> {
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
        clipping_args,
    } = args;

    let reader1 = FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(fastq_input_file1)?));

    let reader2 = match fastq_input_file2 {
        Some(file2) => Some(FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(file2)?))),
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
        keep_header,
        min_read_quality,
        use_median,
        min_length,
        enforce_clipped_length,
        clipping_args,
    };

    Ok(parsed)
}

/// # Panics
///
/// Sub-program for processing fastQ data.
pub fn preprocess_process(args: PreprocessArgs) -> Result<(), std::io::Error> {
    let ParsedPreprocessArgs {
        io_args:
            ParsedPreprocessIoArgs {
                mut table_writer,
                reader1,
                reader2,
                mut log_writer,
                log_file,
            },
        keep_header,
        min_read_quality,
        use_median,
        min_length,
        enforce_clipped_length,
        clipping_args,
    } = parse_preprocess_args(args)?;

    let (r1, r2) = if reader2.is_some() {
        (Some('1'), Some('2'))
    } else {
        (None, None)
    };

    let mut stdout_writer = BufWriter::new(std::io::stdout());

    // try_for_each should be used on this iterator to avoid inefficiencies
    let mut chained_reader = reader1.map(|r| r.map(|q| (q, r1, Simd::from_array([1, 0])))).chain(
        reader2
            .into_iter()
            .flatten()
            .map(|r| r.map(|q| (q, r2, Simd::from_array([0, 1])))),
    );

    let mut metadata_by_sequence: HashMap<_, Vec<_>, _> = HashMap::with_hasher(get_hasher());

    let mut metadata = FastQMetadata::new();

    chained_reader.try_for_each(|record| {
        let (mut unclipped_fq, side, counts) = record?;

        metadata.observed_raw_reads += counts;

        metadata.observed_max_read_len = metadata.observed_max_read_len.max(unclipped_fq.sequence.len());
        if unclipped_fq.sequence.len() < min_length {
            return Ok::<_, std::io::Error>(());
        }

        let mut fq = edit_read(unclipped_fq.as_view_mut(), false, &clipping_args);

        metadata.observed_max_clipped_read_len = metadata.observed_max_clipped_read_len.max(fq.sequence.len());
        if enforce_clipped_length && fq.sequence.len() < min_length {
            return Ok(());
        }

        metadata.passed_len_count += 1;

        let read_q_center = fq.get_q_center(use_median);
        metadata.observed_q_max = if read_q_center > metadata.observed_q_max {
            read_q_center
        } else {
            metadata.observed_q_max
        };
        if read_q_center < Some(f32::from(min_read_quality)) {
            return Ok(());
        }

        metadata.passed_qc_count += 1;

        fq.fix_header(side).keep_or_underscore_header(keep_header);

        // Convert view to owned data, but without copying the header
        let sequence = fq.sequence.to_owned_data();
        let quality = fq.quality.to_owned_data();
        let header = unclipped_fq.header;

        metadata_by_sequence.entry(sequence).or_default().push((header, quality));

        Ok(())
    })?;

    let FastQMetadata {
        passed_qc_count,
        passed_len_count,
        observed_raw_reads,
        observed_q_max,
        observed_max_read_len,
        observed_max_clipped_read_len,
    } = metadata;

    let mut read_pattern_number = 0;
    if passed_qc_count == 0 {
        if let Some(obs_max) = observed_q_max
            && obs_max < f32::from(min_read_quality)
        {
            eprintln!(
                "WARNING: the observed max phred quality score ({obs_max}) is below the user specified threshold (QUAL_THRESHOLD = {min_read_quality})!",
            );
        }

        if observed_max_read_len < min_length {
            eprintln!(
                "WARNING: the observed max read length ({observed_max_read_len}) is below the user specified threshold (MIN_LEN = {min_length})!",
            );
        }
    } else {
        for (sequence, metadata) in metadata_by_sequence.into_iter() {
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
    }

    table_writer.flush()?;
    stdout_writer.flush()?;

    if let Some(ref mut w) = log_writer {
        writeln!(
            w,
            "\n\
            NUMBER_INPUT_FILES\t{num_files}\n\
            OBSERVED_RAW_READS_OR_R1\t{r1_raw_reads}\n\
            OBSERVED_R2_READS\t{r2_raw_reads}\n\
            OBSERVED_MAX_READ_LEN\t{observed_max_read_len}\n\
            OBSERVED_MAX_CLIPPED_READ_LENGTH\t{observed_max_clipped_read_len}\n\
            OBSERVED_MAX_QUALITY\t{observed_q_max}\n\
            READ_COUNT_PASSING_ONLY_LENGTH_FILTER\t{passed_len_count}\n\
            READ_COUNT_PASSING_ALL_QUALITY_CONTROL_FILTERS\t{passed_qc_count}\n\
            READ_PATTERN_COUNT_PASSING\t{read_pattern_number}\n\
            MIN_PHRED_QUALITY_THRESHOLD\t{min_read_quality}\n\
            MIN_READ_LENGTH_THRESHOLD\t{min_length}\n\
            QUALITY_MEASURE\t{center_type}\n\
            ",
            num_files = r2.is_some() as u8 + 1,
            r1_raw_reads = observed_raw_reads[0],
            r2_raw_reads = observed_raw_reads[1],
            observed_q_max = observed_q_max.map(|q| q.to_string()).unwrap_or_else(|| "NONE".to_string()),
            center_type = if use_median { "median" } else { "average" },
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "{MODULE} Warning! Cannot write to {}. See: {e}",
                log_file.as_ref().unwrap().display()
            );
        });
    }

    Ok(())
}
