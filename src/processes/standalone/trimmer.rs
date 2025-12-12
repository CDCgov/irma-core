use crate::{
    args::clipping::{ClippingArgs, ParsedClippingArgs, parse_clipping_args},
    io::{
        InputOptions, IterWithContext, OutputOptions, ReadFileZipPipe, RecordReaders, RecordWriters, WriteFileZipStdout,
        WriteRecord, check_distinct_files,
    },
    utils::{
        paired_reads::{DeinterleavedPairedReadsExt, ZipPairedReadsExt, ZipReadsError},
        trimming::{TrimmedCounts, get_trimming_type, trim_read},
    },
};
use clap::Args;
use std::{io::Write, num::NonZeroUsize, path::PathBuf};
use zoe::prelude::*;

#[derive(Args, Debug)]
pub struct TrimmerArgs {
    /// Path to .fastq or .fastq.gz file to be trimmed
    fastq_input_file: PathBuf,

    /// Path to optional second .fastq or .fastq.gz file to be trimmed
    fastq_input_file2: Option<PathBuf>,

    #[arg(short = '1', short_alias = 'o', long = "fastq-output")]
    /// Output filepath for trimmed reads. Trimmed reads print to STDOUT if not
    /// provided. May also use '-o'.
    fastq_output_file: Option<PathBuf>,

    #[arg(short = '2', long = "fastq-output2")]
    /// Output path for secondary trimmed file if using paired reads. If this
    /// argument is omitted, output is interleaved.
    fastq_output_file2: Option<PathBuf>,

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
}

/// Sub-program for trimming FASTQ data.
pub fn trimmer_process(args: TrimmerArgs) -> Result<(), std::io::Error> {
    let ParsedTrimmerArgs {
        io_args: ParsedPairedIoArgs {
            mut reader1,
            reader2,
            writer,
        },
        trimming_args,
    } = parse_trimmer_args(args)?;

    let mut counts = TrimmedCounts::default();

    let paired_inputs = reader2.is_some();
    let single_out = matches!(writer, RecordWriters::SingleEnd(_));

    if let Some(mut reader2) = reader2 {
        match (writer, trimming_args.filter_widows) {
            // Case 2: In 1, In 2, Out 1 (interleaved Illumina), no widow filtering
            (RecordWriters::SingleEnd(mut writer), false) => {
                let result = reader1
                    .by_ref()
                    .zip_paired_reads_unchecked(reader2.by_ref())
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

            // Case 3: In 1, In 2, Out 1, Filtering widows / orphan reads
            (RecordWriters::SingleEnd(mut writer), true) => {
                reader1
                    .zip_paired_reads(reader2)
                    .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;

                writer.flush()?;
            }

            // Case 4: In 1, In 2, Out 1, Out 2 (separated output Illumina), no filtering
            (RecordWriters::PairedEnd(mut writer), false) => {
                let trimming_args = &trimming_args;
                let mut secondary_counts = TrimmedCounts::default();
                //let secondary_counts_ref = &mut secondary_counts;

                std::thread::scope(|s| {
                    let handle = s.spawn(|| {
                        reader2.try_for_each(|read2| {
                            trim_and_write_seq(read2?, trimming_args, &mut writer.writer2, &mut secondary_counts)
                        })?;
                        writer.writer2.flush()
                    });

                    reader1
                        .try_for_each(|read1| trim_and_write_seq(read1?, trimming_args, &mut writer.writer1, &mut counts))?;
                    writer.writer1.flush()?;

                    handle.join().unwrap()?;
                    std::io::Result::Ok(())
                })?;

                counts = counts + secondary_counts;
            }

            // Case 5: In 1, In 2, Out 1, Out 2, filter widows
            (RecordWriters::PairedEnd(mut writer), true) => {
                reader1
                    .zip_paired_reads(reader2)
                    .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;
                writer.flush()?;
            }
        }
    } else {
        match (writer, trimming_args.filter_widows) {
            (RecordWriters::SingleEnd(mut writer), false) => {
                // Case 1: In 1, Out 1 (ONT, single-end, PacBio)
                reader1.try_for_each(|read| trim_and_write_seq(read?, &trimming_args, &mut writer, &mut counts))?;
                writer.flush()?;
            }

            (RecordWriters::SingleEnd(mut writer), true) => {
                reader1
                    .deinterleave()
                    .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;
                writer.flush()?;
            }

            (RecordWriters::PairedEnd(mut writer), false) => {
                reader1.deinterleave().try_for_each(|pair| {
                    let [read1, read2] = pair?;
                    trim_and_write_seq(read1, &trimming_args, &mut writer.writer1, &mut counts)?;
                    trim_and_write_seq(read2, &trimming_args, &mut writer.writer2, &mut counts)?;
                    std::io::Result::Ok(())
                })?;
                writer.flush()?;
            }

            (RecordWriters::PairedEnd(mut writer), true) => {
                reader1
                    .deinterleave()
                    .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer, &mut counts))?;
                writer.flush()?;
            }
        }
    }

    if trimming_args.clipping_args.verbose {
        let trim_type = get_trimming_type(paired_inputs, single_out, trimming_args.filter_widows)?;
        counts.write_counts(&trimming_args.clipping_args, trim_type, trimming_args.mask);
    }
    Ok(())
}

/// Parsed arguments for the `trimmer` subprocess
pub struct ParsedTrimmerArgs {
    pub io_args:       ParsedPairedIoArgs,
    pub trimming_args: ParsedTrimmerOptions,
}

/// Parsed IO arguments for single or paired reads
pub struct ParsedPairedIoArgs {
    reader1: IterWithContext<FastQReader<ReadFileZipPipe>>,
    reader2: Option<IterWithContext<FastQReader<ReadFileZipPipe>>>,
    writer:  RecordWriters<WriteFileZipStdout>,
}

/// Arguments related to clipping/masking reads, including length/widow
/// filtering
#[derive(Debug)]
pub struct ParsedTrimmerOptions {
    mask:          bool,
    filter_widows: bool,
    min_length:    usize,
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
pub fn parse_trimmer_args(args: TrimmerArgs) -> std::io::Result<ParsedTrimmerArgs> {
    let TrimmerArgs {
        fastq_input_file,
        fastq_input_file2,
        fastq_output_file,
        fastq_output_file2,
        mask,
        filter_widows,
        min_length,
        clipping_args,
    } = args;

    check_distinct_files(
        &fastq_input_file,
        fastq_input_file2.as_ref(),
        fastq_output_file.as_ref(),
        fastq_output_file2.as_ref(),
    )?;

    let readers = InputOptions::new_from_paths(&fastq_input_file, fastq_input_file2.as_ref())
        .use_file_or_zip_threaded()
        .parse_fastq()
        .open()?;

    let writer = OutputOptions::new_from_opt_paths(fastq_output_file.as_ref(), fastq_output_file2.as_ref())
        .use_file_zip_or_stdout()
        .open()?;

    let RecordReaders { reader1, reader2 } = readers;
    let min_length = min_length.get();

    let clipping_args = parse_clipping_args(clipping_args)?;

    let parsed = ParsedTrimmerArgs {
        io_args:       ParsedPairedIoArgs {
            reader1,
            reader2,
            writer,
        },
        trimming_args: ParsedTrimmerOptions {
            mask,
            filter_widows,
            min_length,
            clipping_args,
        },
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
        trim_read(fq_view, args.mask, &args.clipping_args, trim_counts);
        if read.len() >= args.min_length {
            return Some(read.as_view_mut());
        } else {
            trim_counts.filtered += 1;
        }
    } else {
        let fq_view = read.as_view_mut();
        let edited = trim_read(fq_view, args.mask, &args.clipping_args, trim_counts);
        if edited.len() >= args.min_length {
            return Some(edited);
        } else {
            trim_counts.filtered += 1;
        }
    }
    None
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
        return Ok(());
    };
    let Some(r2_trimmed) = trim_filter(&mut read2, args, counts) else {
        return Ok(());
    };
    [r1_trimmed, r2_trimmed].write_record(writer)
}
