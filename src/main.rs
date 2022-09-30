use clap::{Args, Parser, Subcommand};

#[derive(Parser)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    FastqConverter(FastqConverter),
}

#[derive(Args)]
struct FastqConverter {
    use_median: bool,
}

fn main() {
    let args = Cli::parse();
    println!("Hello, world!");
}
