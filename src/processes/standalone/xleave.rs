//! Interleaves or de-interleaves paired FastQ or FASTA files.

use crate::{
    io::{DispatchFastX, InputOptions, OutputOptions, RecordWriters, WriteRecords, check_distinct_files},
    utils::paired_reads::{DeinterleavedPairedReadsExt, ZipPairedReadsExt},
};
use clap::Args;
use std::path::PathBuf;

#[derive(Args, Debug)]
pub struct XleaveArgs {
    /// Path to FASTQ, FASTA, or .gz file to be interleaved/deinterleaved
    pub input_file1: PathBuf,

    /// Path to optional second FASTQ, FASTA, or .gz file to be deinterleaved
    pub input_file2: Option<PathBuf>,

    #[arg(short = '1', long, short_alias = 'o', aliases = ["output-file", "output-file1", "output1"])]
    /// Output file path for interleaved/deinterleaved reads
    pub output: Option<PathBuf>,

    #[arg(short = '2', long, requires = "output", alias = "output-file2")]
    /// Output path for a second sampled file if deinterleaving paired-end
    /// reads. If this argument is omitted, output is interleaved
    pub output2: Option<PathBuf>,
}

pub fn xleave_process(args: XleaveArgs) -> Result<(), std::io::Error> {
    check_distinct_files(
        &args.input_file1,
        args.input_file2.as_ref(),
        args.output.as_ref(),
        args.output2.as_ref(),
    )?;

    let readers = InputOptions::new_from_paths(&args.input_file1, args.input_file2.as_ref())
        .use_file_or_zip_threaded()
        .parse_fastx()
        .open()?;

    let writer = OutputOptions::new_from_opt_paths(args.output.as_ref(), args.output2.as_ref())
        .use_file_zip_or_stdout()
        .open()?;

    let reader1 = readers.reader1;
    let input_path1 = args.input_file1;

    if let Some((reader2, input_path2)) = readers.reader2.zip(args.input_file2) {
        let RecordWriters::SingleEnd(writer) = writer else {
            return Err(std::io::Error::other(
                "Two inputs and two outputs were provided. No interleaving or de-interleaving can occur.",
            ));
        };

        match (reader1.dispatch(), reader2.dispatch()) {
            (DispatchFastX::Fastq(reader1), DispatchFastX::Fastq(reader2)) => reader1
                .zip_paired_reads(reader2)
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1, &input_path2)))
                .write_records(writer)?,
            (DispatchFastX::Fasta(reader1), DispatchFastX::Fasta(reader2)) => reader1
                .zip_paired_reads(reader2)
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1, &input_path2)))
                .write_records(writer)?,
            (DispatchFastX::Fastq(_), DispatchFastX::Fasta(_)) => {
                return Err(std::io::Error::other(
                    "Paired read inputs must be both FASTQ or both FASTA. Found FASTQ for first input and FASTA for second input.",
                ));
            }
            (DispatchFastX::Fasta(_), DispatchFastX::Fastq(_)) => {
                return Err(std::io::Error::other(
                    "Paired read inputs must be both FASTQ or both FASTA. Found FASTA for first input and FASTQ for second input.",
                ));
            }
        }
    } else {
        let RecordWriters::PairedEnd(writer) = writer else {
            return Err(std::io::Error::other(
                "One input and one output were provided. No interleaving or de-interleaving can occur.",
            ));
        };

        match reader1.dispatch() {
            DispatchFastX::Fastq(reader) => reader
                .deinterleave()
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1)))
                .write_records(writer)?,
            DispatchFastX::Fasta(reader) => reader
                .deinterleave()
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1)))
                .write_records(writer)?,
        }
    }

    Ok(())
}
