#![allow(unreachable_patterns)]
#![feature(portable_simd, try_trait_v2)]

use crate::processes::{aligner::*, merge_sam_pairs::*, num_procs::*, preprocess::*, trimmer::*, xflate::*, xleave::*};
use clap::{Parser, Subcommand};
use processes::sampler::{SamplerArgs, sampler_process};
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
    Preprocess(PreprocessArgs),
    /// Merges Illumina paired-end reads with parsimonious error correction and
    /// detection.
    MergeSAM(MergeSAMArgs),
    /// Deflates FastQ files to deduplicated Fasta files, or reinflates
    /// deduplicated Fasta files to FastQ files
    Xflate(XflateArgs),
    /// Provides the physical or logical cores of a CPU portably.
    NumProcs(NumProcsArgs),
    #[command(
        long_about = "Trims FASTQ files for genomic analysis with support for barcodes, adapters, primers, and hard trimming.\
         Barcode and adapter trimming are mutually exclusive. \
         If multiple trim operations are selected, trimming will proceed in barcode/adapter > primer > hard trim order."
    )]
    /// Read FastQ files and trim with various options.
    Trimmer(TrimmerArgs),
    /// Randomly downsamples FastQ or FASTA files. Deinterleaving supported.
    Sampler(SamplerArgs),
    /// Interleaves or De-interleaves paired FastQ or FASTA files.
    Xleave(XleaveArgs),
    /// Perform sequence alignment,
    Aligner(AlignerArgs),
}

fn main() {
    let args = Cli::parse();

    match args.command {
        Commands::Preprocess(cmd_args) => preprocess_process(cmd_args).unwrap_or_die("subcommand 'preprocess'"),
        Commands::MergeSAM(cmd_args) => merge_sam_pairs_process(cmd_args).unwrap_or_die("subcommand 'merge-sam'"),
        Commands::Xflate(cmd_args) => xflate_process(cmd_args).unwrap_or_die("subcommand 'xflate'"),
        Commands::Trimmer(cmd_args) => trimmer_process(cmd_args).unwrap_or_die("subcommand 'trimmer'"),
        Commands::Sampler(cmd_args) => sampler_process(cmd_args).unwrap_or_die("subcommand 'sampler'"),
        Commands::NumProcs(cmd_args) => num_procs_process(cmd_args).unwrap_or_die("subcommand 'num-procs'"),
        Commands::Xleave(cmd_args) => xleave_process(cmd_args).unwrap_or_die("subcommand 'xleave'"),
        Commands::Aligner(cmd_args) => aligner_process(cmd_args).unwrap_or_die("subcommand 'aligner'"),
        _ => {
            eprintln!("IRMA-CORE: unrecognized command {:?}", args.command);
            std::process::exit(1)
        }
    }
}

mod processes;

pub(crate) mod args;
pub(crate) mod io;
pub(crate) mod qc;
pub(crate) mod utils;

pub use crate::processes::*;
