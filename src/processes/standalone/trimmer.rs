use crate::{
    args::clipping::{ClippingArgs, ParsedClippingArgs, parse_clipping_args},
    io::{
        WriteFileZipStdout, FlushWriter, ReadFileZipPipe, RecordReaders, RecordWriters, WriteRecord, check_distinct_files,
    },
    utils::{
        paired_reads::{ZipPairedReadsExt, ZipReadsError},
        trimming::trim_read,
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

    #[arg(short = '2', long = "fastq-output2", requires = "fastq_input_file2")]
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

    #[arg(short = 'f', long, requires = "fastq_input_file2")]
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

    if let Some(mut reader2) = reader2 {
        match (writer, trimming_args.filter_widows) {
            // Case 2: In 1, In 2, Out 1 (interleaved Illumina), no widow filtering
            (RecordWriters::SingleEnd(mut writer), false) => {
                let result = reader1
                    .by_ref()
                    .zip_paired_reads_unchecked(reader2.by_ref())
                    .try_for_each(|pair| {
                        let [read1, read2] = pair?;
                        trim_and_write_seq(read1, &trimming_args, &mut writer)?;
                        trim_and_write_seq(read2, &trimming_args, &mut writer)?;
                        Ok(())
                    });

                match result {
                    Ok(()) => {}
                    Err(ZipReadsError::IoError(e)) => return Err(e),
                    Err(ZipReadsError::ExtraFirstRead(read1)) => {
                        std::iter::once(Ok(read1))
                            .chain(reader1)
                            .try_for_each(|read1| trim_and_write_seq(read1?, &trimming_args, &mut writer))?;
                    }
                    Err(ZipReadsError::ExtraSecondRead(read2)) => {
                        std::iter::once(Ok(read2))
                            .chain(reader2)
                            .try_for_each(|read2| trim_and_write_seq(read2?, &trimming_args, &mut writer))?;
                    }
                }

                writer.flush()?;
            }

            // Case 3: In 1, In 2, Out 1, Filtering widows / orphan reads
            (RecordWriters::SingleEnd(mut writer), true) => {
                reader1
                    .zip_paired_reads(reader2)
                    .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer))?;
                writer.flush()?;
            }

            // Case 4: In 1, In 2, Out 1, Out 2 (separated output Illumina), no filtering
            (RecordWriters::PairedEnd(mut writer), false) => {
                let trimming_args = &trimming_args;

                std::thread::scope(|s| {
                    let handle = s.spawn(move || {
                        reader2.try_for_each(|read2| trim_and_write_seq(read2?, trimming_args, &mut writer.writer2))?;
                        writer.writer2.flush()
                    });

                    reader1.try_for_each(|read1| trim_and_write_seq(read1?, trimming_args, &mut writer.writer1))?;
                    writer.writer1.flush()?;

                    handle.join().unwrap()?;
                    std::io::Result::Ok(())
                })?;
            }

            // Case 5: In 1, In 2, Out 1, Out 2, filter widows
            (RecordWriters::PairedEnd(mut writer), true) => {
                reader1
                    .zip_paired_reads(reader2)
                    .try_for_each(|pair| trim_and_write_pair(pair?, &trimming_args, &mut writer))?;
                writer.flush_all()?;
            }
        }
    } else {
        match writer {
            RecordWriters::SingleEnd(mut writer) => {
                // Case 1: In 1, Out 1 (ONT, single-end, PacBio)
                reader1.try_for_each(|read| trim_and_write_seq(read?, &trimming_args, &mut writer))?;
                writer.flush()?;
            }
            RecordWriters::PairedEnd(_) => unreachable!("Validated by clap"),
        }
    }

    Ok(())
}

/// Parsed arguments for the `trimmer` subprocess
struct ParsedTrimmerArgs {
    io_args:       ParsedPairedIoArgs,
    trimming_args: ParsedTrimmerOptions,
}

/// Parsed IO arguments for single or paired reads
struct ParsedPairedIoArgs {
    reader1: FastQReader<ReadFileZipPipe>,
    reader2: Option<FastQReader<ReadFileZipPipe>>,
    writer:  RecordWriters<WriteFileZipStdout>,
}

/// Arguments related to clipping/masking reads, including length/widow
/// filtering
#[derive(Debug)]
struct ParsedTrimmerOptions {
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
fn parse_trimmer_args(args: TrimmerArgs) -> std::io::Result<ParsedTrimmerArgs> {
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
    let readers = RecordReaders::from_filenames(fastq_input_file, fastq_input_file2)?;
    let writer = RecordWriters::from_optional_filenames(fastq_output_file, fastq_output_file2)?;

    let RecordReaders { reader1, reader2 } = readers;

    let min_length = min_length.get();

    let primer_file = clipping_args.primer_trim.clone();
    let clipping_args = parse_clipping_args(clipping_args).map_err(|e| {
        std::io::Error::other(format!(
            "Failed to process the primers at path {path:#?} due to the error:\n{e}",
            path = primer_file.unwrap()
        ))
    })?;

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
fn trim_filter<'a>(read: &'a mut FastQ, args: &ParsedTrimmerOptions) -> Option<FastQViewMut<'a>> {
    if args.mask {
        let fq_view = read.as_view_mut();
        trim_read(fq_view, args.mask, &args.clipping_args);
        if read.len() >= args.min_length {
            return Some(read.as_view_mut());
        }
    } else {
        let fq_view = read.as_view_mut();
        let edited = trim_read(fq_view, args.mask, &args.clipping_args);
        if edited.len() >= args.min_length {
            return Some(edited);
        }
    }
    None
}

/// Trims a read (either with clipping or masking) and writes it if it passes
/// the length filter.
fn trim_and_write_seq<W: Write>(mut read: FastQ, args: &ParsedTrimmerOptions, writer: &mut W) -> std::io::Result<()> {
    if let Some(trimmed) = trim_filter(&mut read, args) {
        trimmed.write_record(writer)
    } else {
        Ok(())
    }
}

/// Trims a pair of reads (either with clipping or masking) and writes them if
/// both pass the length filter.
fn trim_and_write_pair<'a, W>(pair: [FastQ; 2], args: &ParsedTrimmerOptions, writer: &mut W) -> std::io::Result<()>
where
    for<'b> [FastQViewMut<'b>; 2]: WriteRecord<W>, {
    let [mut read1, mut read2] = pair;
    let Some(r1_trimmed) = trim_filter(&mut read1, args) else {
        return Ok(());
    };
    let Some(r2_trimmed) = trim_filter(&mut read2, args) else {
        return Ok(());
    };
    [r1_trimmed, r2_trimmed].write_record(writer)
}
