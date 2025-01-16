// Description:     Read FastQ files and deflates into a custom XFL format,
//                  converting to FASTA as well. Also can re-inflate back to FASTQ.

use crate::utils::{get_seed, SeedableFoldHashMap};
use clap::Parser;
use std::{
    fs::OpenOptions,
    io::{BufRead, BufReader, BufWriter, Write},
    path::{Path, PathBuf},
};
use zoe::{
    data::{
        fasta::FastaSeq,
        fastq::{FastQ, FastQReader},
        types::phred::QualityScores,
    },
    prelude::{FastaReader, Nucleotides},
};

const CLUSTER_PREFIX: &str = "C";

#[derive(Debug, Parser)]
#[command(version, about)]
pub struct XflateArgs {
    table_file: PathBuf,

    #[clap(required = true)]
    seq_files: Vec<PathBuf>,

    /// Inflate sequence files
    #[arg(short, long)]
    inflate: bool,
}

fn inflate(table_file: &Path, fasta_files: &Vec<PathBuf>) -> Result<(), std::io::Error> {
    let table_reader = BufReader::new(OpenOptions::new().read(true).open(table_file)?);
    let mut stdout_writer = BufWriter::new(std::io::stdout());

    let mut sequence_by_cluster = SeedableFoldHashMap::new(get_seed());

    for file in fasta_files {
        let reader = FastaReader::new(BufReader::new(OpenOptions::new().read(true).open(file)?));
        for record in reader {
            let FastaSeq { name, sequence } = record?;
            let mut sequence = Nucleotides::from_vec_unchecked(sequence);

            if let Some(name) = name.strip_prefix(CLUSTER_PREFIX)
                && let Some(cluster_id) = name.split('%').next()
                && let Ok(cluster_num) = cluster_id.parse::<usize>()
            {
                if name.ends_with("{c}") {
                    sequence.make_reverse_complement();
                }
                sequence_by_cluster.insert(cluster_num, sequence);
            } else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    "Invalid header in fasta file!",
                ));
            }
        }
    }

    for table_record in table_reader.lines() {
        let data = table_record?;
        let mut split = data.split('\t');
        if let Some(name) = split.next()
            && let Some(name) = name.strip_prefix(CLUSTER_PREFIX)
            && let Some(cluster_id) = name.split('%').next()
            && let Ok(cluster_num) = cluster_id.parse::<usize>()
        {
            if let Some(sequence) = sequence_by_cluster.get(&cluster_num) {
                while let (Some(header), Some(quality)) = (split.next(), split.next()) {
                    stdout_writer.write_all(format!("{header}\n{sequence}\n+\n{quality}\n").as_bytes())?;
                }
            }
        } else if !data.is_empty() {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid header in table file!",
            ));
        }
    }

    stdout_writer.flush()?;

    Ok(())
}

fn deflate(table_file: &Path, fastq_files: &Vec<PathBuf>) -> Result<(), std::io::Error> {
    let mut table_writer = BufWriter::new(OpenOptions::new().write(true).create(true).truncate(true).open(table_file)?);
    let mut stdout_writer = BufWriter::new(std::io::stdout());

    let mut metadata_by_sequence = SeedableFoldHashMap::<Nucleotides, Vec<(String, QualityScores)>>::new(get_seed());

    for file in fastq_files {
        let reader = FastQReader::new(BufReader::new(OpenOptions::new().read(true).open(file)?));
        for record in reader {
            let FastQ {
                header,
                sequence,
                quality,
            } = record?;

            metadata_by_sequence.entry(sequence).or_default().push((header, quality));
        }
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

pub fn xflate_process(args: &XflateArgs) -> Result<(), std::io::Error> {
    if args.inflate {
        inflate(&args.table_file, &args.seq_files)
    } else {
        deflate(&args.table_file, &args.seq_files)
    }
}
