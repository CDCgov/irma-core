#![allow(unreachable_patterns)]

use clap::{Parser, Subcommand};
use irma_core::{fastq_converter::*, merge_sam_pairs::*, xflate::*};
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
    /// Performs FastQ quality control, file conversions, and adapter trimming.
    FastqConverter(FastqConverterArgs),
    /// Merges Illumina paired-end reads with parsimonious error correction and
    /// detection.
    MergeSAM(MergeSAMArgs),
    /// Deflates FastQ files to deduplicated Fasta files, or reinflates
    /// deduplicated Fasta files to FastQ files
    Xflate(XflateArgs),
}

fn main() {
    let args = Cli::parse();
    let module = module_path!();

    match &args.command {
        Commands::FastqConverter(cmd_args) => fastqc_process(cmd_args).unwrap_or_die(&format!("{module}::FastqConverter")),
        Commands::MergeSAM(cmd_args) => merge_sam_pairs_process(cmd_args),
        Commands::Xflate(cmd_args) => xflate_process(cmd_args).unwrap_or_die(&format!("{module}::Xflate")),
        _ => {
            eprintln!("IRMA-CORE: unrecognized command {:?}", &args.command);
            std::process::exit(1)
        }
    }
}
