#![doc = include_str!(concat!(env!("CARGO_MANIFEST_DIR"), "/../../README.md"))]
#![allow(unreachable_patterns)]
#![feature(portable_simd)]

use crate::processes::{
    aligner::*, merge_sam_pairs::*, num_procs::*, phase::*, preprocess::*, sampler::*, trimmer::*, xflate::*, xleave::*,
};
use clap::{CommandFactory, FromArgMatches, Parser, Subcommand};
use zoe::data::err::{ErrorWithContext, ResultWithErrorContext};

pub(crate) mod args;
mod processes;

pub(crate) use crate::processes::*;

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
    /// Reads in the variants table and associated distance matrix for a gene.
    /// Assigns a phase number to each variant using thresholded connected
    /// components and writes the results to the variants table.
    Phase(PhaseArgs),
    /// Provides the physical or logical cores of a CPU portably.
    NumProcs(NumProcsArgs),
    #[command(
        long_about = "Trims FASTQ files for genomic analysis with support for barcodes, adapters, primers, and hard trimming. \
         Barcode and adapter trimming are mutually exclusive. \
         If multiple trim operations are selected, trimming will proceed in barcode/adapter > primer > hard trim order."
    )]
    /// Read FastQ files and trim with various options.
    Trimmer(TrimmerArgs),
    /// Randomly downsamples FastQ or FASTA files. Deinterleaving supported.
    Sampler(SamplerArgs),
    /// Interleaves or De-interleaves paired FastQ or FASTA files.
    Xleave(XleaveArgs),
    /// Performs sequence alignment.
    Aligner(AlignerArgs),
}

fn main() -> Result<(), ErrorWithContext> {
    // Perform manual parsing so that we can get subcommand name automatically
    let mut matches = Cli::command().get_matches();
    let command_name = matches.subcommand_name().unwrap().to_owned();
    let args = Cli::from_arg_matches_mut(&mut matches).unwrap();

    let res = match args.command {
        Commands::Preprocess(cmd_args) => preprocess_process(cmd_args),
        Commands::MergeSAM(cmd_args) => merge_sam_pairs_process(cmd_args),
        Commands::Xflate(cmd_args) => xflate_process(cmd_args),
        Commands::Trimmer(cmd_args) => trimmer_process(cmd_args),
        Commands::Sampler(cmd_args) => sampler_process(cmd_args),
        Commands::NumProcs(cmd_args) => num_procs_process(cmd_args),
        Commands::Xleave(cmd_args) => xleave_process(cmd_args),
        Commands::Aligner(cmd_args) => aligner_process(cmd_args),
        Commands::Phase(cmd_args) => phase_process(cmd_args),
    };

    res.with_context(format!("Error in subcommand '{command_name}'"))
}
