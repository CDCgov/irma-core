use crate::{
    io::{FastXReader, RecordWriters, WriteFileZipStdout, WriteRecord, WriteRecords, get_paired_readers_and_writers},
    utils::paired_reads::{DeinterleavedPairedReadsExt, ZipPairedReadsExt},
};
use clap::Args;
use std::path::PathBuf;
use zoe::data::records::HeaderReadable;

#[derive(Args, Debug)]
pub struct XLeaveArgs {
    /// Path to FASTQ, FASTA, or .gz file to be interleaved/deinterleaved
    pub input_file1: PathBuf,

    /// Path to optional second FASTQ, FASTA, or .gz file to be deinterleaved
    pub input_file2: Option<PathBuf>,

    #[arg(short = '1', short_alias = 'o')]
    /// Output file path for interleaved/deinterleaved reads
    pub output_file1: Option<PathBuf>,

    #[arg(short = '2', requires = "output_file1")]
    /// Output path for a second sampled file if deinterleaving paired-end
    /// reads. If this argument is omitted, output is interleaved
    pub output_file2: Option<PathBuf>,
}

pub fn xleave_process(args: XLeaveArgs) -> Result<(), std::io::Error> {
    let (reader1, reader2, writer, threads) = get_paired_readers_and_writers(
        &args.input_file1,
        args.input_file2.as_ref(),
        args.output_file1,
        args.output_file2,
    )?;

    match (reader1, reader2) {
        (FastXReader::Fastq(reader), None) => match writer {
            RecordWriters::SingleEnd(_) => {
                return Err(std::io::Error::other(
                    "One input and one output were provided. No interleaving or de-interleaving can occur.",
                ));
            }
            RecordWriters::PairedEnd(_) => handle_single_input(reader, writer)?,
        },
        (FastXReader::Fastq(reader1), Some(FastXReader::Fastq(reader2))) => match writer {
            RecordWriters::SingleEnd(_) => reader1.zip_paired_reads(reader2).write_records(writer)?,
            RecordWriters::PairedEnd(_) => {
                return Err(std::io::Error::other(
                    "Two inputs and two outputs were provided. No interleaving or de-interleaving can occur.",
                ));
            }
        },
        (FastXReader::Fasta(reader), None) => match writer {
            RecordWriters::SingleEnd(_) => {
                return Err(std::io::Error::other(
                    "One input and one output were provided. No interleaving or de-interleaving can occur.",
                ));
            }
            RecordWriters::PairedEnd(_) => handle_single_input(reader, writer)?,
        },
        (FastXReader::Fasta(reader1), Some(FastXReader::Fasta(reader2))) => match writer {
            RecordWriters::SingleEnd(_) => reader1.zip_paired_reads(reader2).write_records(writer)?,
            RecordWriters::PairedEnd(_) => {
                return Err(std::io::Error::other(
                    "Two inputs and two outputs were provided. No interleaving or de-interleaving can occur.",
                ));
            }
        },
        (FastXReader::Fastq(_), Some(FastXReader::Fasta(_))) => {
            return Err(std::io::Error::other(
                "Paired read inputs must be both FASTQ or both FASTA. Found FASTQ for first input and FASTA for second input.",
            ));
        }
        (FastXReader::Fasta(_), Some(FastXReader::Fastq(_))) => {
            return Err(std::io::Error::other(
                "Paired read inputs must be both FASTQ or both FASTA. Found FASTA for first input and FASTQ for second input.",
            ));
        }
    }

    threads.finalize()?;
    Ok(())
}

fn handle_single_input<R1, A>(reader: R1, writer: RecordWriters<WriteFileZipStdout>) -> std::io::Result<()>
where
    R1: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable + WriteRecord<WriteFileZipStdout>,
    std::io::Result<A>: WriteRecord<WriteFileZipStdout>, {
    match writer {
        RecordWriters::SingleEnd(writer) => reader.write_records(writer),
        RecordWriters::PairedEnd(writer) => reader.deinterleave().write_records(writer),
    }
}
