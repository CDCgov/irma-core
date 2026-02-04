// Description:     Read FastQ files and deflates into a custom XFL format,
//                  converting to FASTA as well. Also can re-inflate back to FASTQ.

use crate::{
    io::{InputOptions, OutputOptions},
    utils::get_hasher,
};
use clap::Parser;
use std::{
    collections::HashMap,
    io::{BufRead, BufReader, Write},
    path::{Path, PathBuf},
};
use zoe::{
    data::{fasta::FastaSeq, fastq::FastQ, types::phred::QualityScores},
    prelude::Nucleotides,
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

/// ## Validity
///
/// This function returns an error intended to be displayed at the top-level. No
/// callers should add additional context other than calling a method in
/// [`OrFail`].
///
/// [`OrFail`]: zoe::data::err::OrFail
fn inflate(table_file: &Path, fasta_files: &Vec<PathBuf>) -> Result<(), std::io::Error> {
    let table_reader = BufReader::new(InputOptions::new_from_path(table_file).use_file().open()?);
    let mut stdout_writer = OutputOptions::new_stdout().open()?;

    let mut sequence_by_cluster = HashMap::with_hasher(get_hasher());

    for file in fasta_files {
        let reader = InputOptions::new_from_path(file).use_file().parse_fasta().open()?;

        for record in reader {
            let FastaSeq { name, sequence } = record?;
            let mut sequence = Nucleotides::from_vec_unchecked(sequence);

            let cluster_num = parse_cluster_num(&name, file)?;

            if name.ends_with("{c}") {
                sequence.make_reverse_complement();
            }
            sequence_by_cluster.insert(cluster_num, sequence);
        }
    }

    for table_record in table_reader.lines() {
        let data = table_record?;

        if data.is_empty() {
            continue;
        }

        let mut split = data.split('\t');

        let Some(name) = split.next() else {
            // This should be unreachable per the docs in split
            continue;
        };

        let cluster_num = parse_cluster_num(name, table_file)?;

        if let Some(sequence) = sequence_by_cluster.get(&cluster_num) {
            while let (Some(header), Some(quality)) = (split.next(), split.next()) {
                write!(stdout_writer, "@{header}\n{sequence}\n+\n{quality}\n")?;
            }
        }
    }

    stdout_writer.flush()?;

    Ok(())
}

fn deflate(table_file: &Path, fastq_files: &Vec<PathBuf>) -> Result<(), std::io::Error> {
    let mut table_writer = OutputOptions::new_from_path(table_file).use_file().open()?;
    let mut stdout_writer = OutputOptions::new_stdout().open()?;

    let mut metadata_by_sequence: HashMap<Nucleotides, Vec<(String, QualityScores)>, _> = HashMap::with_hasher(get_hasher());

    for file in fastq_files {
        let reader = InputOptions::new_from_path(file).use_file().parse_fastq().open()?;
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

        write!(stdout_writer, ">{CLUSTER_PREFIX}{i}%{cluster_size}\n{sequence}\n")?;

        write!(table_writer, "{CLUSTER_PREFIX}{i}%{cluster_size}")?;
        for (header, quality_scores) in metadata {
            write!(table_writer, "\t{header}\t{quality_scores}")?;
        }
        writeln!(table_writer)?;
    }

    table_writer.flush()?;
    stdout_writer.flush()?;

    Ok(())
}

pub fn xflate_process(args: XflateArgs) -> Result<(), std::io::Error> {
    if args.inflate {
        // Validity: No context is added to the result
        inflate(&args.table_file, &args.seq_files)
    } else {
        deflate(&args.table_file, &args.seq_files)
    }
}

/// Given a header containing the contents `name`, parse the cluster number from
/// it.
///
/// The header should be of the format `C<ID>%[REST]`, where `<ID>` is the
/// cluster number being parsed and `[REST]` is any additional optional
/// characters.
///
/// ## Errors
///
/// If `name` does not meet the required format, then an error is returned,
/// including `name` and `path` as context.
fn parse_cluster_num(name: &str, path: &Path) -> std::io::Result<usize> {
    if let Some(name) = name.strip_prefix(CLUSTER_PREFIX)
        && let Some(cluster_id) = name.split('%').next()
        && let Ok(cluster_num) = cluster_id.parse::<usize>()
    {
        Ok(cluster_num)
    } else {
        Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "Invalid header in file: {path}\nHeader: {name}\n\nInflation requires a header of the format: C<ID>%[REST], where:\n  - <ID> is a nonnegative integer\n  - [REST] is any additional optional characters",
                path = path.display()
            ),
        ))
    }
}
