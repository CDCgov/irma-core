use crate::{
    args::clipping::{ClippingArgs, ParsedClippingArgs, parse_clipping_args},
    io::{ReadFileZip, WriteFileZipStdout, create_writer, open_fastq_file},
    utils::{
        paired_reads::{PairedReadFilterer, ReadSide},
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

    #[arg(short = 'f', long)]
    /// Filter widowed reads
    filter_widows: bool,

    #[command(flatten)]
    clipping_args: ClippingArgs,
}

/// Sub-program for trimming FASTQ data.
pub fn trimmer_process(args: TrimmerArgs) -> Result<(), std::io::Error> {
    let ParsedTrimmerArgs {
        io_args:
            ParsedPairedIoArgs {
                reader1,
                reader2,
                mut writer1,
                writer2,
            },
        trimming_args,
    } = parse_trimmer_args(args)?;

    if let Some(reader2) = reader2 {
        match (writer2, trimming_args.filter_widows) {
            // Case 2: In 1, In 2, Out 1 (interleaved Illumina), no widow filtering
            (None, false) => Trimmer::new(&trimming_args, [&mut writer1])
                .handle_paired_reads_no_filter(reader1, reader2)?
                .finalize(),

            // Case 3: In 1, In 2, Out 1, Filtering widows / orphan reads
            (None, true) => Trimmer::new(&trimming_args, [&mut writer1])
                .handle_paired_reads_with_filter_strict(reader1, reader2, |e| e, error_extra_read())?
                .finalize(),

            // Case 4: In 1, In 2, Out 1, Out 2 (separated output Illumina), no filtering
            (Some(mut writer2), false) => Trimmer::new(&trimming_args, [&mut writer1, &mut writer2])
                .handle_paired_reads_no_filter(reader1, reader2)?
                .finalize(),

            // Case 5: In 1, In 2, Out 1, Out 2, filter widows
            (Some(mut writer2), true) => Trimmer::new(&trimming_args, [&mut writer1, &mut writer2])
                .handle_paired_reads_with_filter_strict(reader1, reader2, |e| e, error_extra_read())?
                .finalize(),
        }
    } else {
        // Case 1: In 1, Out 1 (ONT, single-end, PacBio)
        let mut processor = Trimmer::new(&trimming_args, [&mut writer1]);
        processor.handle_single_reads(reader1, ReadSide::Unpaired)?.finalize()
    }
}

/// A [`PairedReadFilterer`] struct used by the trimmer process to perform its
/// core logic. Trimmed sequences are output to 1 (`N=1`) or 2 (`N=2`) files.
#[derive(Debug)]
struct Trimmer<'a, const N: usize> {
    args:    &'a ParsedTrimmerOptions,
    writers: [&'a mut WriteFileZipStdout; N],
}

impl<'a, const N: usize> Trimmer<'a, N> {
    #[inline]
    fn new(args: &'a ParsedTrimmerOptions, writers: [&'a mut WriteFileZipStdout; N]) -> Self {
        Self { args, writers }
    }
}

impl<'a, const N: usize> PairedReadFilterer for Trimmer<'a, N> {
    const PRESERVE_ORDER: bool = true;
    type Processed<'b> = FastQViewMut<'b>;
    type Finalized = std::io::Result<()>;

    #[inline]
    fn process_read<'b>(&mut self, read: &'b mut FastQ, _side: ReadSide) -> Option<FastQViewMut<'b>> {
        if self.args.mask {
            let fq_view = read.as_view_mut();
            trim_read(fq_view, self.args.mask, &self.args.clipping_args);
            if read.len() >= self.args.min_length {
                return Some(read.as_view_mut());
            }
        } else {
            let fq_view = read.as_view_mut();
            let edited = trim_read(fq_view, self.args.mask, &self.args.clipping_args);
            if edited.len() >= self.args.min_length {
                return Some(edited);
            }
        }
        None
    }

    #[inline]
    fn output_read<'b>(&mut self, trimmed: FastQViewMut<'b>, side: ReadSide) -> std::io::Result<()> {
        if N == 1 {
            write!(self.writers[0], "{trimmed}")?;
        } else {
            write!(self.writers[side.to_idx()], "{trimmed}")?;
        }

        Ok(())
    }

    #[inline]
    fn finalize(&mut self) -> Self::Finalized {
        for writer in &mut self.writers {
            writer.flush()?
        }
        Ok(())
    }
}

fn error_extra_read() -> std::io::Error {
    std::io::Error::new(
        std::io::ErrorKind::InvalidData,
        "Extra unpaired read(s) found at end of one of the input FASTQ files.",
    )
}

/// Parsed arguments for the `trimmer` subprocess
struct ParsedTrimmerArgs {
    io_args:       ParsedPairedIoArgs,
    trimming_args: ParsedTrimmerOptions,
}

/// Parsed IO arguments for single or paired reads
struct ParsedPairedIoArgs {
    reader1: FastQReader<ReadFileZip>,
    reader2: Option<FastQReader<ReadFileZip>>,
    writer1: WriteFileZipStdout,
    writer2: Option<WriteFileZipStdout>,
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

fn parse_trimmer_args(args: TrimmerArgs) -> Result<ParsedTrimmerArgs, std::io::Error> {
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

    let reader1 = open_fastq_file(fastq_input_file)?;
    let reader2 = match fastq_input_file2 {
        Some(file2) => Some(open_fastq_file(file2)?),
        None => None,
    };

    let writer1 = create_writer(fastq_output_file.clone())?;
    let writer2 = match fastq_output_file2 {
        Some(path) => Some(create_writer(Some(path))?),
        None => None,
    };

    let min_length = min_length.get();

    let clipping_args = parse_clipping_args(clipping_args)?;

    let parsed = ParsedTrimmerArgs {
        io_args:       ParsedPairedIoArgs {
            reader1,
            reader2,
            writer1,
            writer2,
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
