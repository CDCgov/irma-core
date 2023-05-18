#![allow(unreachable_patterns)]

use clap::{Parser, Subcommand};
use irma_core::fastq_converter::*;
use zoe::data::err::OrFail;

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    FastqConverter(FastqConverterArgs),
}

fn main() {
    let args = Cli::parse();
    let module = module_path!();

    match &args.command {
        Commands::FastqConverter(cmd_args) => fastq_process(cmd_args).unwrap_or_die(&format!("{module}::FastqConverter")),
        _ => {
            eprintln!("Unrecognized command.");
            std::process::exit(1)
        }
    }
}
