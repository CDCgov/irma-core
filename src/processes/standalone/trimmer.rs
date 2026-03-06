//! Reads FastQ files and trims with various options.

use crate::{
    args::clipping::{ClippingArgs, ParsedClippingArgs, parse_clipping_args},
    io::{
        InputOptions, IterWithContext, OutputOptions, PairedWriters, ReadFileZipPipe, RecordWriters, WriteFileZipStdout,
        WriteRecord, check_distinct_files,
    },
    utils::{
        paired_reads::{DeinterleavedPairedReadsExt, ZipPairedReadsExt, ZipReadsError},
        trimming::{TrimmedCounts, trim_read},
    },
};
use clap::Args;
use core::fmt;
use std::{io::Write, num::NonZeroUsize, path::PathBuf};
use zoe::prelude::*;

#[derive(Args, Debug)]
pub struct TrimmerArgs {
    /// Path to .fastq or .fastq.gz file to be trimmed
    fastq_input: PathBuf,

    /// Path to optional second .fastq or .fastq.gz file to be trimmed
    fastq_input2: Option<PathBuf>,

    #[arg(short = '1', long, short_alias = 'o', aliases = ["output-file", "output-file1", "output1", "fastq-output", "fastq-output1"])]
    /// Output filepath for trimmed reads. Trimmed reads print to STDOUT if not
    /// provided. May also use '-o'.
    output: Option<PathBuf>,

    #[arg(short = '2', long, aliases = ["output-file2", "output2", "fastq-output2"])]
    /// Output path for secondary trimmed file if using paired reads. If this
    /// argument is omitted, output is interleaved.
    output2: Option<PathBuf>,

    #[arg(short = 'm', long)]
    /// Perform masking with 'N' instead of clipping. Default behavior is
    /// clipping if not provided
    mask: bool,

    #[arg(short = 'n', long, default_value = "1")]
    /// Minimum sequence length required after trimming. Shorter sequences are
    /// filtered from output.
    min_length: NonZeroUsize,

    #[arg(short = 'f', long)]
    /// Filter widowed reads
    filter_widows: bool,

    #[command(flatten)]
    clipping_args: ClippingArgs,

    #[arg(short = 'v', long)]
    /// Prints the number of records trimmed for each method to stderr
    verbose: bool,
}

/// Sub-program for trimming FASTQ data.
pub fn trimmer_process(args: TrimmerArgs) -> Result<(), std::io::Error> {
    let ParsedTrimmerArgs {
        io_args,
        strategy,
        trimming_args,
        primer_file,
    } = parse_trimmer_args(args)?;

    let mut counts = TrimmedCounts::default();

    match io_args {
        PairedIoArgs::OneInOneOutFilter {
            input_path1,
            reader1,
            mut writer,
        } => {
            reader1
                .deinterleave()
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1)))
                .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;
            writer.flush()?;
        }
        PairedIoArgs::TwoInOneOutFilter {
            input_path1,
            input_path2,
            reader1,
            reader2,
            mut writer,
        } => {
            reader1
                .zip_paired_reads(reader2)
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1, &input_path2)))
                .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;

            writer.flush()?;
        }
        PairedIoArgs::OneInTwoOutFilter {
            input_path1,
            reader1,
            mut writer,
        } => {
            reader1
                .deinterleave()
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1)))
                .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;
            writer.flush()?;
        }
        PairedIoArgs::TwoInTwoOutFilter {
            input_path1,
            input_path2,
            reader1,
            reader2,
            mut writer,
        } => {
            reader1
                .zip_paired_reads(reader2)
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1, &input_path2)))
                .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;
            writer.flush()?;
        }
        PairedIoArgs::OneInOneOutNoFilter { mut reader1, mut writer } => {
            reader1.try_for_each(|read| trim_and_write_seq(read?, &trimming_args, &mut writer, &mut counts))?;
            writer.flush()?;
        }
        PairedIoArgs::TwoInOneOutNoFilter {
            input_path1,
            input_path2,
            mut reader1,
            mut reader2,
            mut writer,
        } => {
            let result = reader1
                .by_ref()
                .zip_paired_reads_unchecked(reader2.by_ref())
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1, &input_path2)))
                .try_for_each(|pair| {
                    let [read1, read2] = pair?;
                    trim_and_write_seq(read1, &trimming_args, &mut writer, &mut counts)?;
                    trim_and_write_seq(read2, &trimming_args, &mut writer, &mut counts)?;
                    Ok(())
                });

            match result {
                Ok(()) => {}
                Err(ZipReadsError::IoError(e)) => return Err(e),
                Err(ZipReadsError::ExtraFirstRead(read1)) => {
                    std::iter::once(Ok(read1))
                        .chain(reader1)
                        .try_for_each(|read1| trim_and_write_seq(read1?, &trimming_args, &mut writer, &mut counts))?;
                }
                Err(ZipReadsError::ExtraSecondRead(read2)) => {
                    std::iter::once(Ok(read2))
                        .chain(reader2)
                        .try_for_each(|read2| trim_and_write_seq(read2?, &trimming_args, &mut writer, &mut counts))?;
                }
            }

            writer.flush()?;
        }
        PairedIoArgs::OneInTwoOutNoFilter {
            input_path1,
            reader1,
            mut writer,
        } => {
            reader1
                .deinterleave()
                .map(|res| res.map_err(|e| e.add_path_context(&input_path1)))
                .try_for_each(|pair| {
                    let [read1, read2] = pair?;
                    trim_and_write_seq(read1, &trimming_args, &mut writer.writer1, &mut counts)?;
                    trim_and_write_seq(read2, &trimming_args, &mut writer.writer2, &mut counts)
                })?;
            writer.flush()?;
        }
        PairedIoArgs::TwoInTwoOutNoFilter {
            mut reader1,
            mut reader2,
            mut writer,
        } => {
            let trimming_args = &trimming_args;
            let mut secondary_counts = TrimmedCounts::default();

            std::thread::scope(|s| {
                let handle = s.spawn(|| {
                    reader2.try_for_each(|read2| {
                        trim_and_write_seq(read2?, trimming_args, &mut writer.writer2, &mut secondary_counts)
                    })?;
                    writer.writer2.flush()
                });

                reader1.try_for_each(|read1| trim_and_write_seq(read1?, trimming_args, &mut writer.writer1, &mut counts))?;
                writer.writer1.flush()?;

                handle.join().unwrap()
            })?;

            counts = counts + secondary_counts;
        }
    }

    if trimming_args.verbose {
        counts.write_counts(&trimming_args.clipping_args, strategy, &trimming_args, primer_file);
    }
    Ok(())
}

/// Possible IO combinations for trimming inputs, outputs and widow filtering.
/// Input can be paired or single, output can be paired or single, and filtering
/// can be toggled, leaving 8 possible options.
///
/// This enum also holds the reader(s) and writer for each option.
#[allow(clippy::enum_variant_names)]
enum PairedIoArgs<I, W> {
    OneInOneOutFilter {
        input_path1: PathBuf,
        reader1:     I,
        writer:      W,
    },
    TwoInOneOutFilter {
        input_path1: PathBuf,
        input_path2: PathBuf,
        reader1:     I,
        reader2:     I,
        writer:      W,
    },
    OneInTwoOutFilter {
        input_path1: PathBuf,
        reader1:     I,
        writer:      PairedWriters<W>,
    },
    TwoInTwoOutFilter {
        input_path1: PathBuf,
        input_path2: PathBuf,
        reader1:     I,
        reader2:     I,
        writer:      PairedWriters<W>,
    },
    OneInOneOutNoFilter {
        reader1: I,
        writer:  W,
    },
    TwoInOneOutNoFilter {
        input_path1: PathBuf,
        input_path2: PathBuf,
        reader1:     I,
        reader2:     I,
        writer:      W,
    },
    OneInTwoOutNoFilter {
        input_path1: PathBuf,
        reader1:     I,
        writer:      PairedWriters<W>,
    },
    TwoInTwoOutNoFilter {
        reader1: I,
        reader2: I,
        writer:  PairedWriters<W>,
    },
}

/// Similar to [`PairedIoArgs`], but lacking state, used for display purposes.
#[derive(Debug)]
pub enum PairedIoStrategy {
    OneInOneOutFilter,
    TwoInOneOutFilter,
    OneInTwoOutFilter,
    TwoInTwoOutFilter,
    OneInOneOutNoFilter,
    TwoInOneOutNoFilter,
    OneInTwoOutNoFilter,
    TwoInTwoOutNoFilter,
}

impl fmt::Display for PairedIoStrategy {
    /// Display output for the different possible IO combinations, to be used in
    /// [`write_counts`](TrimmedCounts::write_counts).
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PairedIoStrategy::OneInOneOutFilter => write!(
                f,
                "one interleaved input and wrote to one interleaved output with filtering of widowed reads"
            ),
            PairedIoStrategy::TwoInOneOutFilter => write!(
                f,
                "two paired inputs and wrote to one interleaved output with filtering of widowed reads"
            ),

            PairedIoStrategy::OneInTwoOutFilter => write!(
                f,
                "one interleaved input and wrote to two paired outputs with filtering of widowed reads"
            ),
            PairedIoStrategy::TwoInTwoOutFilter => write!(
                f,
                "two paired inputs and wrote to two paired outputs with filtering of widowed reads"
            ),
            PairedIoStrategy::OneInOneOutNoFilter => write!(f, "one input and wrote to one output"),
            PairedIoStrategy::TwoInOneOutNoFilter => write!(
                f,
                "two paired inputs and wrote to one interleaved output without filtering of widowed read"
            ),
            PairedIoStrategy::OneInTwoOutNoFilter => write!(
                f,
                "one interleaved input and wrote to two paired outputs without filtering of widowed reads"
            ),
            PairedIoStrategy::TwoInTwoOutNoFilter => write!(
                f,
                "two paired inputs and wrote to two paired outputs without filtering of widowed reads"
            ),
        }
    }
}

/// Parsed arguments for the `trimmer` subprocess
struct ParsedTrimmerArgs {
    io_args:       PairedIoArgs<IterWithContext<FastQReader<ReadFileZipPipe>>, WriteFileZipStdout>,
    strategy:      PairedIoStrategy,
    trimming_args: ParsedTrimmerOptions,
    primer_file:   Option<PathBuf>,
}

/// Arguments related to clipping/masking reads, including length/widow
/// filtering
#[derive(Debug)]
struct ParsedTrimmerOptions {
    mask:          bool,
    min_length:    usize,
    verbose:       bool,
    clipping_args: ParsedClippingArgs,
}

/// Parses the trimmer arguments from the clap arguments
///
/// ## Errors
///
/// An error could occur when opening the first FASTQ file, the second FASTQ
/// file, creating the first or second writer, or processing the primer file.
/// Any errors generated will have customized error messages including
/// additional information.
fn parse_trimmer_args(args: TrimmerArgs) -> std::io::Result<ParsedTrimmerArgs> {
    let TrimmerArgs {
        fastq_input,
        fastq_input2,
        output,
        output2,
        mask,
        filter_widows,
        min_length,
        clipping_args,
        verbose,
    } = args;

    check_distinct_files(&fastq_input, fastq_input2.as_ref(), output.as_ref(), output2.as_ref())?;

    let readers = InputOptions::new_from_paths(&fastq_input, fastq_input2.as_ref())
        .use_file_or_zip_threaded()
        .parse_fastq()
        .open()?;

    let writer = OutputOptions::new_from_opt_paths(output.as_ref(), output2.as_ref())
        .use_file_zip_or_stdout()
        .open()?;

    let reader1 = readers.reader1;
    let input_path1 = fastq_input;

    let (io_args, strategy) = if let Some((reader2, input_path2)) = readers.reader2.zip(fastq_input2) {
        match (writer, filter_widows) {
            // Case 2: In 1, In 2, Out 1 (interleaved Illumina), no widow filtering
            (RecordWriters::SingleEnd(writer), false) => (
                PairedIoArgs::TwoInOneOutNoFilter {
                    input_path1,
                    input_path2,
                    reader1,
                    reader2,
                    writer,
                },
                PairedIoStrategy::TwoInOneOutNoFilter,
            ),

            // Case 3: In 1, In 2, Out 1, Filtering widows / orphan reads
            (RecordWriters::SingleEnd(writer), true) => (
                PairedIoArgs::TwoInOneOutFilter {
                    input_path1,
                    input_path2,
                    reader1,
                    reader2,
                    writer,
                },
                PairedIoStrategy::TwoInOneOutFilter,
            ),

            // Case 4: In 1, In 2, Out 1, Out 2 (separated output Illumina), no filtering
            (RecordWriters::PairedEnd(writer), false) => (
                PairedIoArgs::TwoInTwoOutNoFilter {
                    reader1,
                    reader2,
                    writer,
                },
                PairedIoStrategy::TwoInTwoOutNoFilter,
            ),

            // Case 5: In 1, In 2, Out 1, Out 2, filter widows
            (RecordWriters::PairedEnd(writer), true) => (
                PairedIoArgs::TwoInTwoOutFilter {
                    input_path1,
                    input_path2,
                    reader1,
                    reader2,
                    writer,
                },
                PairedIoStrategy::TwoInTwoOutFilter,
            ),
        }
    } else {
        match (writer, filter_widows) {
            (RecordWriters::SingleEnd(writer), false) => {
                // Case 1: In 1, Out 1 (ONT, single-end, PacBio)
                (
                    PairedIoArgs::OneInOneOutNoFilter { reader1, writer },
                    PairedIoStrategy::OneInOneOutNoFilter,
                )
            }

            (RecordWriters::SingleEnd(writer), true) => (
                PairedIoArgs::OneInOneOutFilter {
                    input_path1,
                    reader1,
                    writer,
                },
                PairedIoStrategy::OneInOneOutFilter,
            ),

            (RecordWriters::PairedEnd(writer), false) => (
                PairedIoArgs::OneInTwoOutNoFilter {
                    input_path1,
                    reader1,
                    writer,
                },
                PairedIoStrategy::OneInTwoOutNoFilter,
            ),

            (RecordWriters::PairedEnd(writer), true) => (
                PairedIoArgs::OneInTwoOutFilter {
                    input_path1,
                    reader1,
                    writer,
                },
                PairedIoStrategy::OneInTwoOutFilter,
            ),
        }
    };

    let min_length = min_length.get();

    let primer_file = clipping_args.primer_trim.clone();
    let clipping_args = parse_clipping_args(clipping_args)?;

    let parsed = ParsedTrimmerArgs {
        io_args,
        strategy,
        trimming_args: ParsedTrimmerOptions {
            mask,
            min_length,
            clipping_args,
            verbose,
        },
        primer_file,
    };

    Ok(parsed)
}

/// Trims a read (either with clipping or masking) and checks its length. `Some`
/// is returned if it passes the length filter.
fn trim_filter<'a>(
    read: &'a mut FastQ, args: &ParsedTrimmerOptions, trim_counts: &mut TrimmedCounts,
) -> Option<FastQViewMut<'a>> {
    if args.mask {
        let fq_view = read.as_view_mut();
        trim_read(fq_view, args.mask, &args.clipping_args, trim_counts, args.verbose);
        if read.len() >= args.min_length {
            Some(read.as_view_mut())
        } else {
            trim_counts.length_filtered += 1;
            None
        }
    } else {
        let fq_view = read.as_view_mut();
        let edited = trim_read(fq_view, args.mask, &args.clipping_args, trim_counts, args.verbose);
        if edited.len() >= args.min_length {
            Some(edited)
        } else {
            trim_counts.length_filtered += 1;
            None
        }
    }
}

/// Trims a read (either with clipping or masking) and writes it if it passes
/// the length filter.
fn trim_and_write_seq<W: Write>(
    mut read: FastQ, args: &ParsedTrimmerOptions, writer: &mut W, counts: &mut TrimmedCounts,
) -> std::io::Result<()> {
    counts.total_processed += 1;
    if let Some(trimmed) = trim_filter(&mut read, args, counts) {
        trimmed.write_record(writer)
    } else {
        Ok(())
    }
}

/// Trims a pair of reads (either with clipping or masking) and writes them if
/// both pass the length filter.
fn trim_and_write_pair<'a, W>(
    pair: [FastQ; 2], args: &ParsedTrimmerOptions, writer: &mut W, counts: &mut TrimmedCounts,
) -> std::io::Result<()>
where
    for<'b> [FastQViewMut<'b>; 2]: WriteRecord<W>, {
    counts.total_processed += 2;
    let [mut read1, mut read2] = pair;
    let Some(r1_trimmed) = trim_filter(&mut read1, args, counts) else {
        // Filtered first read, which we've counted as a length filter, so need to
        // count second read as being widow filtered
        counts.widow_filtered += 1;
        return Ok(());
    };
    let Some(r2_trimmed) = trim_filter(&mut read2, args, counts) else {
        counts.widow_filtered += 1;
        return Ok(());
    };
    [r1_trimmed, r2_trimmed].write_record(writer)
}

impl TrimmedCounts {
    /// Writes the trimmed and filtered counts for verbose output
    fn write_counts(
        self, args: &ParsedClippingArgs, strategy: PairedIoStrategy, options: &ParsedTrimmerOptions,
        primer_file: Option<PathBuf>,
    ) {
        let ParsedClippingArgs {
            preserve_bases: _,
            barcodes,
            b_restrict_left: _,
            b_restrict_right: _,
            b_hdist,
            adapters,
            a_fuzzy,
            primer_kmers,
            p_restrict_left: _,
            p_restrict_right: _,
            polyg_left,
            polyg_right,
            hard_left,
            hard_right,
        } = args;
        let ParsedTrimmerOptions {
            mask,
            min_length,
            clipping_args: _,
            verbose: _,
        } = options;

        let trim_mask = match mask {
            true => "masked",
            false => "trimmed",
        };

        eprintln!("IRMA-core trimmer processed reads from {strategy}");

        eprintln!("{:<20} {:>10} reads", "Input:", self.total_processed);

        if polyg_left.is_some() || polyg_right.is_some() {
            let polyg_left = polyg_left.unwrap_or(0);
            let polyg_right = polyg_right.unwrap_or(0);
            let percent = self.poly_g as f64 / self.total_processed as f64 * 100.0;
            let thresholds = if polyg_left == polyg_right {
                format!("with a PolyG threshold of {polyg_left} bases")
            } else {
                format!(
                    "with a left PolyG threshold of {polyg_left} bases and a right PolyG threshold of {polyg_right} bases"
                )
            };
            eprintln!(
                "{:<20} {:>10} reads ({percent:.2}%) {thresholds}",
                format!("PolyG {trim_mask}:"),
                self.poly_g
            );
        }
        if barcodes.is_some() {
            let percent = self.barcode as f64 / self.total_processed as f64 * 100.0;
            eprintln!(
                "{:<20} {:>10} reads ({percent:.2}%) with an allowable hamming distance of {b_hdist}",
                format!("Barcode {trim_mask}:"),
                self.barcode
            );
        }
        if adapters.is_some() {
            let percent = self.adapter as f64 / self.total_processed as f64 * 100.0;
            let fuzziness = match a_fuzzy {
                true => "fuzzy",
                false => "exact",
            };
            eprintln!(
                "{:<20} {:>10} reads ({percent:.2}%) with {fuzziness} matching",
                format!("Adapter {trim_mask}:"),
                self.adapter
            );
        }
        if primer_kmers.is_some() {
            let context = if let Some(path) = primer_file {
                format!("using primer set {}", path.display())
            } else {
                "".to_string()
            };
            let percent = self.primer as f64 / self.total_processed as f64 * 100.0;
            eprintln!(
                "{:<20} {:>10} reads ({percent:.2}%) {context}",
                format!("Primer {trim_mask}:"),
                self.primer
            );
        }
        if *hard_left > 0usize || *hard_right > 0usize {
            let percent = self.hard as f64 / self.total_processed as f64 * 100.0;
            let thresholds = if hard_left == hard_right {
                format!("with an amount of {hard_left} bases")
            } else {
                format!("with an amount of {hard_left} bases on the left and {hard_right} bases on the right")
            };
            eprintln!(
                "{:<20} {:>10} reads ({percent:.2}%) {thresholds}",
                format!("Hard {trim_mask}:"),
                self.hard
            );
        }

        let percent_trimmed = self.total_trimmed as f64 / self.total_processed as f64 * 100.0;
        eprintln!(
            "{:<20} {:>10} reads ({percent_trimmed:.2}%)",
            format!("Total {trim_mask}:"),
            self.total_trimmed
        );

        let percent_filtered = self.length_filtered as f64 / self.total_processed as f64 * 100.0;
        eprintln!(
            "{:<20} {:>10} reads ({percent_filtered:.2}%) for being shorter than the minimum post-trimming length of {min_length}",
            "Length filtered:", self.length_filtered,
        );

        if matches!(
            strategy,
            PairedIoStrategy::OneInOneOutFilter
                | PairedIoStrategy::OneInTwoOutFilter
                | PairedIoStrategy::TwoInOneOutFilter
                | PairedIoStrategy::TwoInTwoOutFilter
        ) {
            let percent_widowed = self.widow_filtered as f64 / self.total_processed as f64 * 100.;
            eprintln!(
                "{:<20} {:>10} reads ({percent_widowed:.2}%) for their paired read being shorter than the minimum post-trimming length of {min_length}",
                "Widow filtered:", self.widow_filtered
            )
        }
    }
}
