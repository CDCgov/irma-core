//! Merges Illumina paired-end reads with parsimonious error correction and
//! detection.

use clap::Args;
use irma_records::{
    io::{Finish, InputOptions, OutputOptions, ValidatePaths},
    sam::{PairedMergeStats, pair_and_merge},
};
use std::{io::Write, path::PathBuf};
use zoe::data::sam::*;

#[derive(Args, Debug)]
pub struct MergeSAMArgs {
    /// Reference file used to generate the SAM.
    fasta_reference: PathBuf,

    /// SAM file to merge R1 and R2 pairs via alignment and parsimonious
    /// correction.
    sam_file: PathBuf,

    /// Output directory and prefix for merged SAM data.
    output_prefix: PathBuf,

    #[arg(short = 'S', long)]
    /// Serialize output observations for downstream analysis.
    store_stats: bool,

    #[arg(short = 'B', long)]
    /// SAM is in bowtie format.
    bowtie_format: bool,
}

struct ParsedMergeSamArgs {
    /// Reference file used to generate the SAM.
    fasta_reference: PathBuf,

    /// SAM file to merge R1 and R2 pairs via alignment and parsimonious
    /// correction.
    sam_file: PathBuf,

    /// The path for the output SAM file.
    merged_sam_file: PathBuf,

    /// If `Some`, the file to output observations for downstream analysis.
    paired_stats_file: Option<PathBuf>,

    /// SAM is in bowtie format.
    bowtie_format: bool,
}

fn parse_merge_sam_args(args: MergeSAMArgs) -> ParsedMergeSamArgs {
    if &args.output_prefix == "/dev/null" {
        ParsedMergeSamArgs {
            fasta_reference:   args.fasta_reference,
            sam_file:          args.sam_file,
            merged_sam_file:   args.output_prefix.clone(),
            paired_stats_file: args.store_stats.then_some(args.output_prefix),
            bowtie_format:     args.bowtie_format,
        }
    } else {
        ParsedMergeSamArgs {
            fasta_reference:   args.fasta_reference,
            sam_file:          args.sam_file,
            merged_sam_file:   args.output_prefix.with_extension("sam"),
            paired_stats_file: args.store_stats.then(|| args.output_prefix.with_extension("stats")),
            bowtie_format:     args.bowtie_format,
        }
    }
}

impl ValidatePaths for ParsedMergeSamArgs {
    fn inputs(&self) -> impl IntoIterator<Item = &PathBuf> {
        [&self.fasta_reference, &self.sam_file]
    }

    fn outputs(&self) -> impl IntoIterator<Item = &PathBuf> {
        let merged_sam_file = std::iter::once(&self.merged_sam_file);
        let paired_stats_file = self.paired_stats_file.iter();

        merged_sam_file.chain(paired_stats_file)
    }
}

pub fn merge_sam_pairs_process(args: MergeSAMArgs) -> Result<(), std::io::Error> {
    let args = parse_merge_sam_args(args);

    args.validate_paths()?;

    let mut ref_reader = InputOptions::new_from_path(&args.fasta_reference)
        .use_file()
        .parse_fasta()
        .open()?;

    let mut reference = ref_reader.next().ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("No valid fasta data in file: {file}", file = args.fasta_reference.display()),
        )
    })??;

    let Some(new_len) = reference.header.split_ascii_whitespace().next().map(str::len) else {
        return Err(std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!(
                "Empty header field in FASTA file: {file}",
                file = args.fasta_reference.display()
            ),
        ));
    };
    reference.header.truncate(new_len);

    const ONE_MB: usize = 2usize.pow(20);

    let mut sam_writer = OutputOptions::new_from_path(&args.merged_sam_file)
        .with_capacity(ONE_MB)
        .use_file()
        .open()?;

    let mut sam_rows = InputOptions::new_from_path(&args.sam_file)
        .use_file()
        .parse_sam()
        .open()?
        .peekable();

    // Write the header section before streaming the records to be merged
    while let Some(Ok(SamRow::Header(h))) = sam_rows.next_if(|row| matches!(row, Ok(SamRow::Header(_)))) {
        writeln!(sam_writer, "{h}")?;
    }

    let records = sam_rows.map(|row| match row? {
        SamRow::Data(record) => Ok(record),
        SamRow::Header(_) => Err(std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!(
                "SAM header lines must precede alignment records in: {file}",
                file = args.sam_file.display()
            ),
        )),
    });

    let mut stats = PairedMergeStats::default();

    for record in pair_and_merge(
        records,
        &reference.header,
        &reference.sequence,
        args.bowtie_format,
        &mut stats,
    ) {
        writeln!(sam_writer, "{}", record?)?;
    }

    if let Some(paired_stats_file) = args.paired_stats_file {
        let mut w = OutputOptions::new_from_path(&paired_stats_file).use_file().open()?;

        let PairedMergeStats {
            observations,
            true_variations,
            variant_errors,
            deletion_errors,
            insert_obs,
            insert_errors,
        } = stats;

        writeln!(
            &mut w,
            "{header}\tobs\t{observations}\n\
             {header}\ttmv\t{true_variations}\n\
             {header}\tfmv\t{variant_errors}\n\
             {header}\tdmv\t{deletion_errors}\n\
             {header}\tinsObs\t{insert_obs}\n\
             {header}\tinsErr\t{insert_errors}",
            header = reference.header
        )?;

        w.finish()?;
    }

    sam_writer.finish()
}
