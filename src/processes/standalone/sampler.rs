//! Randomly downsamples FastQ or FASTA files. Deinterleaving supported.

use crate::{
    io::{
        DispatchFastX, FastXReader, InputOptions, IterWithContext, OutputOptions, ReadFileZipPipe, RecordReaders,
        RecordWriters, SequenceWriter, WriteFileZipStdout, WriteRecord, WriteRecordCompatibleItem, WriteRecords,
        check_distinct_files, is_gz,
    },
    utils::paired_reads::{DeinterleavedPairedReads, DeinterleavedPairedReadsExt, ZipPairedReadsExt},
};
use clap::Args;
use rand::{SeedableRng, make_rng};
use rand_xoshiro::Xoshiro256StarStar;
use std::{
    fmt::Debug,
    io::{BufRead, Read, Write},
    path::{Path, PathBuf},
};
use zoe::{
    data::records::HeaderReadable,
    iter_utils::{
        ProcessResultsExt,
        sampling::{DownsampleBernoulli, SkipSampler, downsample_reservoir},
    },
};

#[derive(Args, Debug)]
pub struct SamplerArgs {
    /// Path to FASTQ, FASTA, or .gz file to be sampled
    pub input_file: PathBuf,

    /// Path to optional second FASTQ, FASTA, or .gz file to be sampled
    pub input_file2: Option<PathBuf>,

    #[arg(short = '1', long, short_alias = 'o', aliases = ["output-file", "output-file1", "output1"])]
    /// Output file path for sampled reads. Sampled reads print to STDOUT if not
    /// provided. May also use '-o'
    pub output: Option<PathBuf>,

    #[arg(short = '2', long, requires = "output", alias = "output-file2")]
    /// Output path for a second sampled file if using paired-end reads. If this
    /// argument is omitted, output is interleaved
    pub output2: Option<PathBuf>,

    // this is for requiring either percent or subsample target, but not both
    #[command(flatten)]
    target: Target,

    #[arg(short = 's', long)]
    /// For reproducibility, provide an optional seed for the random number
    /// generator
    pub rng_seed: Option<u64>,

    #[arg(short = 'v', long)]
    /// Prints the original number of records and subsampled amount to stderr
    pub verbose: bool,
}

#[derive(Args, Debug)]
#[group(required = true, multiple = false)]
struct Target {
    #[arg(short = 't', long)]
    /// Target number of reads to be subsampled. Either a `subsample_target` or
    /// `percent_target` must be specified
    pub subsample_target: Option<usize>,

    #[arg(short = 'p', long, value_parser = validate_percent,)]
    /// Target percentage of reads to be sampled. Must be a positive integer in
    /// [0, 100]. Either a `subsample_target` or `percent_target` must be
    /// specified
    pub percent_target: Option<usize>,
}

/// Parses a percent (as a `usize`) from the command line
fn validate_percent(value: &str) -> Result<usize, String> {
    let parsed = value
        .parse::<usize>()
        .map_err(|_| format!("`{value}` is not a valid, positive, whole integer."))?;
    if (0..=100).contains(&parsed) {
        Ok(parsed)
    } else {
        Err(format!(
            "Sampling target percent must be between 0 and 100, but {parsed} was provided."
        ))
    }
}

/// main process getting called by irma-core main.rs
pub fn sampler_process(args: SamplerArgs) -> Result<(), std::io::Error> {
    let (io_args, rng, target, verbose) = parse_sampler_args(args)?;

    // Get the population sequence count from one of the files if possible
    let mut seq_count = get_paired_seq_count(&io_args)?;

    let is_single = io_args.reader2.is_none() && matches!(io_args.writer, RecordWriters::SingleEnd(_));

    // For de-interleaving, must divide sequence count by 2 to get number of
    // pairs
    if io_args.reader2.is_none() && matches!(io_args.writer, RecordWriters::PairedEnd(_)) {
        seq_count = seq_count.map(|seq_count| seq_count / 2)
    }

    // Update the target with the population sequence count
    let target = match (target, seq_count) {
        (SamplingTarget::Count(count), _) => SamplingTarget::Count(count),
        (SamplingTarget::Percent(percent), Some(seq_count)) => SamplingTarget::Count(seq_count * percent / 100),
        (SamplingTarget::Percent(percent), None) => SamplingTarget::Percent(percent),
    };

    let Reader {
        path: input_path1,
        iter: reader1,
    } = io_args.reader1;

    let (total_original, total_downsampled) = if let Some(reader2) = io_args.reader2 {
        let Reader {
            path: input_path2,
            iter: reader2,
        } = reader2;

        let input_paths = [input_path1, input_path2];

        match (reader1.dispatch(), reader2.dispatch()) {
            (DispatchFastX::Fastq(reader1), DispatchFastX::Fastq(reader2)) => {
                sample_paired_input(reader1, reader2, io_args.writer, target, seq_count, rng, input_paths)?
            }
            (DispatchFastX::Fasta(reader1), DispatchFastX::Fasta(reader2)) => {
                sample_paired_input(reader1, reader2, io_args.writer, target, seq_count, rng, input_paths)?
            }
            (DispatchFastX::Fastq(_), DispatchFastX::Fasta(_)) => {
                return Err(std::io::Error::other(
                    "Paired read inputs must be both FASTQ or both FASTA. Found FASTQ for first input and FASTA for second input.",
                ));
            }
            (DispatchFastX::Fasta(_), DispatchFastX::Fastq(_)) => {
                return Err(std::io::Error::other(
                    "Paired read inputs must be both FASTQ or both FASTA. Found FASTA for first input and FASTQ for second input.",
                ));
            }
        }
    } else {
        match reader1.dispatch() {
            DispatchFastX::Fastq(reader) => sample_single_input(reader, io_args.writer, target, seq_count, rng)?,
            DispatchFastX::Fasta(reader) => sample_single_input(reader, io_args.writer, target, seq_count, rng)?,
        }
    };

    if verbose {
        let single_paired = if is_single { "total records" } else { "pairs of records" };
        let percent = 100.0 * total_downsampled as f32 / total_original as f32;
        eprintln!("Downsampled {total_original} {single_paired} to {total_downsampled} ({percent:.02} %).");
    }
    Ok(())
}

/// Performs sampling for a single input file.
///
/// This may correspond to either single reads or interleaved paired reads,
/// depending on the number of writers.
///
/// This returns a tuple containing the original counts and downsampled counts.
/// For single end reads, the counts are the number of records. For paired end
/// reads, each pair counts once.
fn sample_single_input<R1, W, A>(
    reader: R1, writer: RecordWriters<W>, target: SamplingTarget, seq_count: Option<usize>, rng: Xoshiro256StarStar,
) -> std::io::Result<(usize, usize)>
where
    R1: Iterator<Item = std::io::Result<A>>,
    W: Write,
    A: HeaderReadable + WriteRecord<W>,
    std::io::Result<A>: WriteRecord<W>, {
    // Don't perform sampling if target is higher than population sequence count
    if let SamplingTarget::Count(target_count) = target
        && let Some(seq_count) = seq_count
        && target_count > seq_count
    {
        eprintln!(
            "Sampler Warning: Target sample size ({target_count}) was greater than population size ({seq_count}); no downsampling has occurred.",
        );
        match writer {
            RecordWriters::SingleEnd(writer) => {
                reader.write_records(writer)?;
            }
            RecordWriters::PairedEnd(writer) => {
                reader.deinterleave().write_records(writer)?;
            }
        }
        return Ok((seq_count, seq_count));
    }

    match writer {
        RecordWriters::SingleEnd(writer) => {
            let iterator = reader;
            sample_and_write_results(iterator, writer, target, seq_count, rng)
        }
        RecordWriters::PairedEnd(writer) => {
            let iterator: DeinterleavedPairedReads<R1, A> = reader.deinterleave();
            sample_and_write_results(iterator, writer, target, seq_count, rng)
        }
    }
}

/// Performs sampling for a pair of inputs.
///
/// This returns a tuple containing the original counts and downsampled counts.
/// Each pair of reads counts once.
fn sample_paired_input<R1, R2, W, A>(
    reader1: R1, reader2: R2, writer: RecordWriters<W>, target: SamplingTarget, seq_count: Option<usize>,
    rng: Xoshiro256StarStar, input_paths: [PathBuf; 2],
) -> std::io::Result<(usize, usize)>
where
    R1: Iterator<Item = std::io::Result<A>>,
    R2: Iterator<Item = std::io::Result<A>>,
    W: Write,
    A: HeaderReadable + WriteRecord<W> + Debug + Sync + Send + 'static, {
    // Zip the paired reads, and add context including the paths to any zipping
    // errors
    let iterator = reader1
        .zip_paired_reads(reader2)
        .map(|res| res.map_err(|e| e.add_path_context(&input_paths[0], &input_paths[1])));

    // Don't perform sampling if target is higher than population sequence count
    if let SamplingTarget::Count(target_count) = target
        && let Some(seq_count) = seq_count
        && target_count > seq_count
    {
        eprintln!(
            "Sampler Warning: Target sample size ({target_count}) was greater than population size ({seq_count}); no downsampling has occurred.",
        );
        iterator.write_records(writer)?;
        return Ok((seq_count, seq_count));
    }

    sample_and_write_results(iterator, writer, target, seq_count, rng)
}

/// Samples and writes all records from an iterator of results, propagating any
/// errors in the input.
///
/// See [`sample_and_write_records`] for more details. This function ensures
/// that any errors in the input (even those that are downsampled) are
/// propagated.
///
/// This returns a tuple containing the original counts and downsampled counts
/// from the iterator. For single end reads, the counts are the number of
/// records. For paired end reads, each pair counts once.
///
/// ## Assumptions
///
/// This function must be called on an iterator of results. The `Ok` variant of
/// the items should be a record such as [`FastQ`] or [`FastaSeq`], not a nested
/// result. Otherwise, this would permit logic errors such as silently ignoring
/// errors.
///
/// [`FastQ`]: zoe::data::records::fastq::FastQ
/// [`FastaSeq`]: zoe::data::records::fasta::FastaSeq
fn sample_and_write_results<I, W, A, E>(
    iterator: I, writer: W, target: SamplingTarget, seq_count: Option<usize>, rng: Xoshiro256StarStar,
) -> std::io::Result<(usize, usize)>
where
    I: Iterator<Item = Result<A, E>>,
    W: SequenceWriter,
    A: WriteRecordCompatibleItem<W>,
    std::io::Error: From<E>, {
    iterator.process_results(|mut iter| {
        let out = sample_and_write_records(&mut iter, writer, target, seq_count, rng);
        // Fully exhaust the input iterator, so that any errors are indeed
        // surfaced
        iter.last();
        out
    })?
}

/// Samples and writes all records from an iterator of records.
///
/// The method used is one of the following:
///
/// 1. Bernoulli sampling, if `target` is a [`Percent`] and the population
///    `seq_count` is [`None`].
/// 2. Method D sampling, if `target` is a [`Count`] and the population
///    `seq_count` is [`Some`].
/// 3. Resovoir sampling (method L), if `target` is a [`Count`] and the
///    population `seq_count` is [`None`].
///
/// This returns a tuple containing the original counts and downsampled counts
/// from the iterator. For single end reads, the counts are the number of
/// records. For paired end reads, each pair counts once.
///
/// ## Validity
///
/// This function should not be called on an iterator of results. Otherwise,
/// this would permit logic errors such as silently ignoring errors.
///
/// [`Percent`]: SamplingTarget::Percent
/// [`Count`]: SamplingTarget::Count
/// [`FastQ`]: zoe::data::records::fastq::FastQ
/// [`FastaSeq`]: zoe::data::records::fasta::FastaSeq
#[inline]
fn sample_and_write_records<I, W>(
    iterator: &mut I, writer: W, target: SamplingTarget, seq_count: Option<usize>, mut rng: Xoshiro256StarStar,
) -> std::io::Result<(usize, usize)>
where
    I: Iterator<Item: WriteRecordCompatibleItem<W>>,
    W: SequenceWriter, {
    let mut total_original = 0;
    let mut total_downsampled = 0;

    match target {
        SamplingTarget::Percent(percent) => {
            iterator
                .inspect(|_| total_original += 1)
                .downsample_bernoulli(percent as f32 / 100.0, &mut rng)
                .inspect(|_| total_downsampled += 1)
                .write_records(writer)?;
        }
        SamplingTarget::Count(target) => {
            if let Some(total_items) = seq_count {
                total_original = total_items;
                SkipSampler::new(iterator, target, total_items, &mut rng)?
                    .inspect(|_| total_downsampled += 1)
                    .write_records(writer)?;
            } else {
                let samples = downsample_reservoir(iterator.inspect(|_| total_original += 1), &mut rng, target);
                total_downsampled = samples.len();
                samples.into_iter().write_records(writer)?;
            }
        }
    }

    Ok((total_original, total_downsampled))
}

/// Gets the number of input sequences, using whichever paired input exists, is
/// a file, and is not zipped.
///
/// If neither meets these conditions, `None` is returned.
fn get_paired_seq_count(io_args: &IOArgs) -> std::io::Result<Option<usize>> {
    let IOArgs {
        reader1,
        reader2,
        writer: _,
    } = &io_args;

    if reader1.path.is_file() && !is_gz(&reader1.path) {
        Ok(Some(get_seq_count(&reader1.path, reader1.iter.inner_iter())?))
    } else if let Some(reader2) = reader2
        && reader2.path.is_file()
        && !is_gz(&reader2.path)
    {
        Ok(Some(get_seq_count(&reader2.path, reader2.iter.inner_iter())?))
    } else {
        Ok(None)
    }
}

/// The type sampler uses for input, along with the input path for error
/// context.
struct Reader {
    path: PathBuf,
    iter: IterWithContext<FastXReader<ReadFileZipPipe>>,
}

/// The IO arguments used by sampler, including up to two readers and writers.
struct IOArgs {
    reader1: Reader,
    reader2: Option<Reader>,
    writer:  RecordWriters<WriteFileZipStdout>,
}

/// The target number of sequences to sample
enum SamplingTarget {
    /// The target specified as a percent of the input number of sequences
    Percent(usize),
    /// The target as an exact count
    Count(usize),
}

fn parse_sampler_args(args: SamplerArgs) -> Result<(IOArgs, Xoshiro256StarStar, SamplingTarget, bool), std::io::Error> {
    let rng = if let Some(seed) = &args.rng_seed {
        Xoshiro256StarStar::seed_from_u64(*seed)
    } else {
        make_rng()
    };

    check_distinct_files(
        &args.input_file,
        args.input_file2.as_ref(),
        args.output.as_ref(),
        args.output2.as_ref(),
    )?;

    let readers = InputOptions::new_from_paths(&args.input_file, args.input_file2.as_ref())
        .use_file_or_zip_threaded()
        .parse_fastx()
        .open()?;

    let writer = OutputOptions::new_from_opt_paths(args.output.as_ref(), args.output2.as_ref())
        .use_file_zip_or_stdout()
        .open()?;

    let RecordReaders { reader1, reader2 } = readers;

    let reader1 = Reader {
        path: args.input_file,
        iter: reader1,
    };
    let reader2 = args.input_file2.zip(reader2).map(|(path, iter)| Reader { path, iter });

    let io_args = IOArgs {
        reader1,
        reader2,
        writer,
    };
    let target = if let Some(count) = args.target.subsample_target {
        SamplingTarget::Count(count)
    } else if let Some(percent) = args.target.percent_target {
        SamplingTarget::Percent(percent)
    } else {
        unreachable!("This can't be reached because clap requires a value for either count or percent")
    };
    Ok((io_args, rng, target, args.verbose))
}

/// Gets the count of the number of records in `input_file`.
///
/// This is achieved by counting the number of lines, and dividing it by the
/// proper amount (2 for FASTA, and 4 for FASTQ). The input file must exist, be
/// a file, and not be zipped.
fn get_seq_count<R: Read>(input_file: &Path, reader: &FastXReader<R>) -> std::io::Result<usize> {
    let input = InputOptions::new_from_path(input_file).use_file().open()?;
    let line_count = input.lines().process_results(|iter| iter.count())?;
    match reader {
        FastXReader::Fasta(_) => Ok(line_count / 2),
        FastXReader::Fastq(_) => Ok(line_count / 4),
    }
}
