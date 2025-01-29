// Description:      Read FastQ files, applies QC filtering (quality and length),
//                   adapter trimming, and deflates into a custom XFL format + FASTA format.

use crate::{
    qc::{fastq::*, fastq_metadata::*},
    utils::{get_seed, SeedableFoldHashMap},
};
use clap::{Args, ValueHint};
use indoc::writedoc;
use std::{
    fs::OpenOptions,
    io::{prelude::*, BufReader, BufWriter},
    path::PathBuf,
    simd::prelude::*,
};
use zoe::{data::types::nucleotides::reverse_complement, prelude::*};

#[derive(Args, Debug)]
pub struct QcTrimDeflateArgs {
    /// Location to store the XFL file.
    pub table_file: PathBuf,

    /// Single-ended FASTQ or the R1 file.
    pub fastq_input_file1: PathBuf,

    /// The R2 paired-end FASTQ file.
    pub fastq_input_file2: Option<PathBuf>,

    #[arg(short = 'H', long)]
    /// Keep the fastq header as usual.
    pub keep_header: bool,

    #[arg(short = 'T', long, default_value_t = 0)]
    /// Specify the read quality threshold (geometric mean, median).
    pub min_read_quality: u8,

    #[arg(short = 'M', long)]
    /// Interprets the threshold (-T) as the median, not the geometric mean.
    pub use_median: bool,

    #[arg(short = 'L', long, default_value_t = 0)]
    /// Minimum length of sequence read data, filtered otherwise.
    pub min_length: usize,

    #[arg(short = 'E', long)]
    /// The minimum length threshold (-L) is enforced when adapter clipped (-c).
    pub enforce_clipped_length: bool,

    #[arg(short = 'm', long)]
    /// Specify adapter sequence and mask when found in reads.
    pub mask_adapter: Option<String>,

    #[arg(short = 'c', long)]
    /// Specify adapter sequence and clip appropriate ends when found in reads.
    pub clip_adapter: Option<String>,

    #[arg(short = 'Z', long)]
    /// Allow up to one mismatch for adapter clipping (-c) or masking (-m).
    pub fuzzy_adapter: bool,

    #[arg(short = 'U', long)]
    /// Re-encode FASTQ sequence to expected input: A, C, T, G, N
    pub canonical_bases: bool,

    #[arg(short = 'G', long, value_hint = ValueHint::FilePath)]
    /// Quality control log path and filename.
    pub log_file: Option<PathBuf>,
}

const CLUSTER_PREFIX: &str = "C";
static MODULE: &str = module_path!();

/// # Panics
///
/// Sub-program for processing fastQ data.
pub fn qc_trim_deflate_process(args: QcTrimDeflateArgs) -> Result<(), std::io::Error> {
    let fastq_file_reader1 = FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(&args.fastq_input_file1)?));

    let (fastq_file_reader2, r1, r2) = if let Some(file2) = &args.fastq_input_file2 {
        (
            Some(FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(file2)?))),
            Some('1'),
            Some('2'),
        )
    } else {
        (None, None, None)
    };

    let mut log_file_writer = if let Some(ref file_path) = args.log_file {
        Some(BufWriter::new(
            OpenOptions::new().write(true).create(true).truncate(true).open(file_path)?,
        ))
    } else {
        None
    };

    let mut table_writer = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&args.table_file)?,
    );
    let mut stdout_writer = BufWriter::new(std::io::stdout());

    // TODO: What if user passes both mask_adapter and clip_adapter? Should throw error
    let (forward_adapter, reverse_adapter) = match (&args.mask_adapter, &args.clip_adapter) {
        (Some(ref a), _) | (_, Some(ref a)) => {
            let forward = a.as_bytes().to_ascii_uppercase();
            let reverse = reverse_complement(&forward);
            (forward, reverse)
        }
        _ => (Vec::new(), Vec::new()),
    };

    // Some uses of chain are considered inefficient.
    let mut chained_reader = fastq_file_reader1
        .map(|r| r.map(|q| (q, r1, Simd::from_array([1, 0]))))
        .chain(
            fastq_file_reader2
                .into_iter()
                .flatten()
                .map(|r| r.map(|q| (q, r2, Simd::from_array([0, 1])))),
        );

    let mut metadata_by_sequence: SeedableFoldHashMap<_, Vec<_>> = SeedableFoldHashMap::new(get_seed());
    let mut metadata = FastQMetadata::new();

    chained_reader.try_for_each(|record| {
        let (mut fq, side, counts) = record?;

        metadata.observed_raw_reads += counts;

        metadata.observed_max_read_len = metadata.observed_max_read_len.max(fq.sequence.len());
        if fq.sequence.len() < args.min_length {
            return Ok::<_, std::io::Error>(());
        }

        fq.to_canonical_bases(args.canonical_bases)
            .transform_by_reverse_forward_search(
                args.fuzzy_adapter,
                args.clip_adapter.is_some(),
                &reverse_adapter,
                &forward_adapter,
            );
        metadata.observed_max_clipped_read_len = metadata.observed_max_clipped_read_len.max(fq.sequence.len());
        if args.enforce_clipped_length && fq.sequence.len() < args.min_length {
            return Ok(());
        }

        metadata.passed_len_count += 1;

        let read_q_center = fq.get_q_center(args.use_median);
        metadata.observed_q_max = if read_q_center > metadata.observed_q_max {
            read_q_center
        } else {
            metadata.observed_q_max
        };
        if read_q_center < Some(f32::from(args.min_read_quality)) {
            return Ok(());
        }

        metadata.passed_qc_count += 1;

        fq.fix_header(side).keep_or_underscore_header(args.keep_header);

        metadata_by_sequence
            .entry(fq.sequence)
            .or_default()
            .push((fq.header, fq.quality));

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
            && obs_max < f32::from(args.min_read_quality)
        {
            eprintln!(
                "WARNING: the observed max phred quality score ({obs_max}) is below the user specified threshold (QUAL_THRESHOLD = {})!",
                args.min_read_quality
            );
        }

        if observed_max_read_len < args.min_length {
            eprintln!(
                "WARNING: the observed max read length ({observed_max_read_len}) is below the user specified threshold (MIN_LEN = {})!",
                args.min_length
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

    if let Some(ref mut w) = log_file_writer {
        writedoc!(
            w,
            "
            NUMBER_INPUT_FILES\t{num_files}
            OBSERVED_RAW_READS_OR_R1\t{r1_raw_reads}
            OBSERVED_R2_READS\t{r2_raw_reads}
            OBSERVED_MAX_READ_LEN\t{observed_max_read_len}
            OBSERVED_MAX_CLIPPED_READ_LENGTH\t{observed_max_clipped_read_len}
            OBSERVED_MAX_QUALITY\t{observed_q_max}
            READ_COUNT_PASSING_ONLY_LENGTH_FILTER\t{passed_len_count}
            READ_COUNT_PASSING_ALL_QUALITY_CONTROL_FILTERS\t{passed_qc_count}
            READ_PATTERN_COUNT_PASSING\t{read_pattern_number}
            MIN_PHRED_QUALITY_THRESHOLD\t{min_read_quality_threshold}
            MIN_READ_LENGTH_THRESHOLD\t{min_length_threshold}
            QUALITY_MEASURE\t{center_type}
            ",
            num_files = args.fastq_input_file2.is_some() as u8 + 1,
            r1_raw_reads = observed_raw_reads[0],
            r2_raw_reads = observed_raw_reads[1],
            observed_q_max = observed_q_max.map(|q| q.to_string()).unwrap_or_else(|| "NONE".to_string()),
            min_read_quality_threshold = args.min_read_quality,
            min_length_threshold = args.min_length,
            center_type = if args.use_median { "median" } else { "average" },
        )
        .unwrap_or_else(|e| {
            eprintln!(
                "{MODULE} Warning! Cannot write to {}. See: {e}",
                args.log_file.as_ref().unwrap().display()
            );
        });
    }

    Ok(())
}
