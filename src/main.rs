#![allow(unreachable_patterns)]
#![feature(let_chains, round_char_boundary, portable_simd)]

use crate::processes::{fastq_converter::*, merge_sam_pairs::*, num_procs::*, qc_trim_deflate::*, trimmer::*, xflate::*};
use clap::{Parser, Subcommand};
use zoe::data::err::OrFail;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Performs all-in-one FastQ quality control, trimming, and deflation to
    /// XFL and FASTA formats.
    QcTrimDeflate(QcTrimDeflateArgs),
    /// Performs FastQ quality control, file conversions, and adapter trimming.
    FastqConverter(FastqConverterArgs),
    /// Merges Illumina paired-end reads with parsimonious error correction and
    /// detection.
    MergeSAM(MergeSAMArgs),
    /// Deflates FastQ files to deduplicated Fasta files, or reinflates
    /// deduplicated Fasta files to FastQ files
    Xflate(XflateArgs),
    /// Provides the physical or logical cores of a CPU portably.
    NumProcs(NumProcsArgs),
    /// Trims FASTQ files for genomic analysis with support for barcodes, adapters,
    /// primers, and hard trimming. Barcode and adapter trimming are mutually exclusive.
    /// If multiple trim operations are selected, trimming will proceed in
    /// barcode/adapter > primer > hard trim order.
    Trimmer(TrimmerArgs),
}

fn main() {
    let args = Cli::parse();
    let module = module_path!();

    match args.command {
        Commands::QcTrimDeflate(cmd_args) => {
            qc_trim_deflate_process(cmd_args).unwrap_or_die(&format!("{module}::QcTrimDeflate"))
        }
        Commands::FastqConverter(cmd_args) => fastqc_process(cmd_args).unwrap_or_die(&format!("{module}::FastqConverter")),
        Commands::MergeSAM(cmd_args) => merge_sam_pairs_process(cmd_args),
        Commands::Xflate(cmd_args) => xflate_process(cmd_args).unwrap_or_die(&format!("{module}::Xflate")),
        Commands::Trimmer(cmd_args) => trimmer_process(cmd_args).unwrap_or_die(&format!("{module}::Trimmer")),
        Commands::NumProcs(cmd_args) => num_procs_process(cmd_args).unwrap_or_die(&format!("{module}::NumProcs")),
        _ => {
            eprintln!("IRMA-CORE: unrecognized command {:?}", args.command);
            std::process::exit(1)
        }
    }
}

mod processes;

pub(crate) mod qc;
pub(crate) mod utils;

pub use crate::processes::*;
