#![allow(unreachable_patterns)]

use clap::{Parser, Subcommand};
use irma_core::fastq_converter::*;
use irma_core::merge_sam_pairs::*;
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
    FastqConverter(FastqConverterArgs),
    MergeSAM(MergeSAMArgs),
}

fn main() {
    let args = Cli::parse();
    let module = module_path!();

    match &args.command {
        Commands::FastqConverter(cmd_args) => fastqc_process(cmd_args).unwrap_or_die(&format!("{module}::FastqConverter")),
        Commands::MergeSAM(cmd_args) => merge_sam_pairs_process(cmd_args),
        _ => {
            eprintln!("IRMA-CORE: unrecognized command {:?}", &args.command);
            std::process::exit(1)
        }
    }
}
