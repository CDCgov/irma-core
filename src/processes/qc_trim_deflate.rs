use clap::{Args, ValueHint};
use foldhash::{HashMap, HashMapExt};
use rayon::prelude::*;
use std::{
    borrow::Borrow,
    fs::OpenOptions,
    io::{prelude::*, BufReader, BufWriter},
    path::PathBuf,
};
use zoe::{
    data::{types::nucleotides::reverse_complement, types::phred::QualityScores},
    prelude::*,
};

use crate::qc::{fastq::*, fastq_metadata::*};

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

    #[arg(short = 'g', long)]
    /// Specify log ID tag (integer) for output collation.
    pub log_id: Option<usize>,
}

const CLUSTER_PREFIX: &str = "C";
static MODULE: &str = module_path!();

/// # Panics
///
/// Sub-program for processing fastQ data.
pub fn qc_trim_deflate_process(args: &QcTrimDeflateArgs) -> Result<(), std::io::Error> {
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

    let (mut log_file_writer, log_name_id) = if let Some(ref file_path) = args.log_file {
        let mut log_id = String::new();
        log_id.push_str(file_path.to_string_lossy().borrow());
        if let Some(id) = args.log_id {
            log_id.push(':');
            log_id.push_str(id.to_string().as_str());
        }

        (
            Some(BufWriter::new(OpenOptions::new().append(true).create(true).open(file_path)?)),
            log_id,
        )
    } else {
        (None, String::new())
    };

    let mut table_writer = BufWriter::new(
        OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&args.table_file)?,
    );
    let mut stdout_writer = BufWriter::new(std::io::stdout());

    let (forward_adapter, reverse_adapter) = match (&args.mask_adapter, &args.clip_adapter) {
        (Some(ref a), _) | (_, Some(ref a)) => {
            let forward = a.as_bytes().to_ascii_uppercase();
            let reverse = reverse_complement(&forward);
            (forward, reverse)
        }
        _ => (Vec::new(), Vec::new()),
    };

    // Some uses of chain are considered inefficient.
    let chained_reader = fastq_file_reader1
        .map(|r| r.map(|q| (q, r1)))
        .chain(fastq_file_reader2.into_iter().flatten().map(|r| r.map(|q| (q, r2))));

    type OkType = (FastQMetadata, HashMap<Nucleotides, Vec<(String, QualityScores)>>);
    let (metadata, metadata_by_sequence) = chained_reader
        .par_bridge()
        .try_fold(
            || (FastQMetadata::new(), HashMap::new()),
            |(set_metadata, mut xfl_partition), record| {
                let (mut fq, side) = record?;
                let mut item_metadata = FastQMetadata::new();

                if fq.sequence.len() < args.min_length {
                    return Ok::<OkType, std::io::Error>((set_metadata.merge(item_metadata), xfl_partition));
                }

                fq.to_canonical_bases(args.canonical_bases)
                    .transform_by_reverse_forward_search(
                        args.fuzzy_adapter,
                        args.clip_adapter.is_some(),
                        &reverse_adapter,
                        &forward_adapter,
                    );

                if args.enforce_clipped_length && fq.sequence.len() < args.min_length {
                    return Ok::<OkType, std::io::Error>((set_metadata.merge(item_metadata), xfl_partition));
                }

                item_metadata.dataset_max_clipped_read_len = Some(fq.sequence.len());
                item_metadata.passed_len_count = 1;

                let Some(read_q_center) = fq.get_q_center(args.use_median) else {
                    return Ok::<OkType, std::io::Error>((set_metadata.merge(item_metadata), xfl_partition));
                };

                item_metadata.dataset_q_max = Some(read_q_center);

                if read_q_center < f32::from(args.min_read_quality) {
                    return Ok::<OkType, std::io::Error>((set_metadata.merge(item_metadata), xfl_partition));
                }

                item_metadata.passed_qc_count = 1;

                fq.fix_header(side).keep_or_underscore_header(args.keep_header);

                xfl_partition.entry(fq.sequence).or_default().push((fq.header, fq.quality));

                Ok::<OkType, std::io::Error>((set_metadata.merge(item_metadata), xfl_partition))
            },
        )
        .try_reduce(
            || (FastQMetadata::new(), HashMap::new()),
            |(md1, mut hm1), (md2, hm2)| {
                for (seq, vec) in hm2 {
                    hm1.entry(seq)
                        .and_modify(|old_vec| {
                            old_vec.extend_from_slice(&vec);
                        })
                        .or_insert(vec);
                }
                Ok((md1.merge(md2), hm1))
            },
        )?;

    let FastQMetadata {
        passed_qc_count,
        passed_len_count,
        dataset_q_max,
        dataset_max_read_len,
        dataset_max_clipped_read_len,
    } = metadata;

    if let Some(ref mut w) = log_file_writer {
        writeln!(
            w,
            "{log_name_id}\t{passed_len_count}\t{passed_qc_count}\t{}\t{}\t{}",
            args.min_read_quality, args.min_length, args.use_median
        )
        .unwrap_or_else(|e| {
            eprintln!("{MODULE} Warning! Cannot write to {log_name_id}. See: {e}");
        });
    }

    if passed_qc_count == 0 {
        let center_type = if args.use_median { "median" } else { "average" };
        let dataset_q_max = if let Some(v) = dataset_q_max {
            v.to_string()
        } else {
            String::from("NONE")
        };
        let dataset_max_read_len = if let Some(v) = dataset_max_read_len {
            v.to_string()
        } else {
            String::from("NONE")
        };
        let dataset_max_clipped_read_len = if let Some(v) = dataset_max_clipped_read_len {
            v.to_string()
        } else {
            String::from("NONE")
        };

        eprintln!(
            "{MODULE} Warning! No reads passed QC.

            Configuration specified minimum {center_type} read quality score: {}
            Configuration specified minimum read length: {}
            Configuration minimum clipped read length enforced: {}

            Observed dataset read length max: {dataset_max_read_len}
            Observed dataset clipped read length max: {dataset_max_clipped_read_len}
            Observed dataset {center_type} read quality max: {dataset_q_max}\n",
            args.min_read_quality, args.min_length, args.enforce_clipped_length
        );
    }

    for (i, (sequence, metadata)) in metadata_by_sequence.into_iter().enumerate() {
        let cluster_size = metadata.len();

        stdout_writer.write_all(format!(">{CLUSTER_PREFIX}{i}%{cluster_size}\n{sequence}\n").as_bytes())?;

        table_writer.write_all(format!("{CLUSTER_PREFIX}{i}%{cluster_size}").as_bytes())?;
        for (header, quality_scores) in metadata {
            table_writer.write_all(format!("\t{header}\t{quality_scores}").as_bytes())?;
        }
        table_writer.write_all(b"\n")?;
    }

    table_writer.flush()?;
    stdout_writer.flush()?;

    Ok(())
}
