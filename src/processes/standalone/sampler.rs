use crate::io::{
    FastXReader, IoThreads, ReadFileZip, RecordWriters, WriteFileZipStdout, WriteRecord, WriteRecords,
    get_paired_readers_and_writers, is_gz,
};
use crate::utils::paired_reads::{DeinterleavedPairedReads, DeinterleavedPairedReadsExt, ZipPairedReadsExt};
use clap::Args;
use rand::{Rng, SeedableRng, seq::IndexedMutRandom};
use rand_xoshiro::Xoshiro256StarStar;
use std::{
    fs::File,
    io::{BufRead, BufReader},
    path::PathBuf,
};
use zoe::data::records::HeaderReadable;
use zoe::unwrap_or_return_some_err;

#[derive(Args, Debug)]
pub struct SamplerArgs {
    /// Path to FASTQ, FASTA, or .gz file to be sampled
    pub input_file1: PathBuf,

    /// Path to optional second FASTQ, FASTA, or .gz file to be sampled
    pub input_file2: Option<PathBuf>,

    #[arg(short = '1', short_alias = 'o')]
    /// Output file path for sampled reads. Sampled reads print to STDOUT if not
    /// provided. May also use '-o'
    pub output_file1: Option<PathBuf>,

    #[arg(short = '2', requires = "output_file1")]
    /// Output path for a second sampled file if using paired-end reads. If this
    /// argument is omitted, output is interleaved
    pub output_file2: Option<PathBuf>,

    #[arg(short = 't', long, group = "target")]
    /// Target number of reads to be subsampled
    pub subsample_target: Option<usize>,

    #[arg(short = 'p', long, group = "target", value_parser = validate_percent,)]
    /// Target percentage of reads to be sampled. Must be a positive integer in
    /// [0, 100]
    pub percent_target: Option<usize>,

    #[arg(short = 's', long)]
    /// For reproducibility, provide an optional seed for the random number
    /// generator
    pub rng_seed: Option<u64>,
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
    let (io_args, rng, threads, target) = parse_sampler_args(args)?;

    // Get the population sequence count from one of the files if possible
    let mut seq_count = get_paired_seq_count(io_args.input_path1, io_args.input_path2, &io_args.reader1)?;

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

    match (io_args.reader1, io_args.reader2) {
        (FastXReader::Fastq(reader), None) => sample_single_input(reader, io_args.writer, target, seq_count, rng)?,
        (FastXReader::Fastq(reader1), Some(FastXReader::Fastq(reader2))) => {
            sample_paired_input(reader1, reader2, io_args.writer, target, seq_count, rng)?
        }
        (FastXReader::Fasta(reader), None) => sample_single_input(reader, io_args.writer, target, seq_count, rng)?,
        (FastXReader::Fasta(reader1), Some(FastXReader::Fasta(reader2))) => {
            sample_paired_input(reader1, reader2, io_args.writer, target, seq_count, rng)?
        }
        (FastXReader::Fastq(_), Some(FastXReader::Fasta(_))) => {
            return Err(std::io::Error::other(
                "Paired read inputs must be both FASTQ or both FASTA. Found FASTQ for first input and FASTA for second input.",
            ));
        }
        (FastXReader::Fasta(_), Some(FastXReader::Fastq(_))) => {
            return Err(std::io::Error::other(
                "Paired read inputs must be both FASTQ or both FASTA. Found FASTA for first input and FASTQ for second input.",
            ));
        }
    }

    threads.finalize()?;
    Ok(())
}

/// Performs sampling for a single input file.
fn sample_single_input<R1, A>(
    reader: R1, writer: RecordWriters<WriteFileZipStdout>, target: SamplingTarget, seq_count: Option<usize>,
    rng: Xoshiro256StarStar,
) -> std::io::Result<()>
where
    R1: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable + WriteRecord<WriteFileZipStdout>,
    std::io::Result<A>: WriteRecord<WriteFileZipStdout>, {
    // Don't perform sampling if target is higher than population sequence count
    if let SamplingTarget::Count(target_count) = target
        && let Some(seq_count) = seq_count
        && target_count > seq_count
    {
        eprintln!(
            "Sampler Warning: Target sample size ({target_count}) was greater than population size ({seq_count}); no downsampling has occurred.",
        );
        return match writer {
            RecordWriters::SingleEnd(writer) => reader.write_records(writer),
            RecordWriters::PairedEnd(writer) => reader.deinterleave().write_records(writer),
        };
    }

    match writer {
        RecordWriters::SingleEnd(writer) => {
            let iterator = reader;
            sample_and_write_records(iterator, writer, target, seq_count, rng)
        }
        RecordWriters::PairedEnd(writer) => {
            let iterator: DeinterleavedPairedReads<R1, A> = reader.deinterleave();
            sample_and_write_records(iterator, writer, target, seq_count, rng)
        }
    }
}

/// Performs sampling for a pair of inputs.
fn sample_paired_input<R1, R2, A>(
    reader1: R1, reader2: R2, writer: RecordWriters<WriteFileZipStdout>, target: SamplingTarget, seq_count: Option<usize>,
    rng: Xoshiro256StarStar,
) -> std::io::Result<()>
where
    R1: Iterator<Item = std::io::Result<A>>,
    R2: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable + WriteRecord<WriteFileZipStdout>, {
    // Don't perform sampling if target is higher than population sequence count
    if let SamplingTarget::Count(target_count) = target
        && let Some(seq_count) = seq_count
        && target_count > seq_count
    {
        eprintln!(
            "Sampler Warning: Target sample size ({target_count}) was greater than population size ({seq_count}); no downsampling has occurred.",
        );
        return reader1.zip_paired_reads(reader2).write_records(writer);
    }

    // Determine the proper iterator type for the inputs, then dispatch
    let iterator = reader1.zip_paired_reads(reader2);
    sample_and_write_records(iterator, writer, target, seq_count, rng)
}

/// Sample and write all records from an iterator. The method used is one of the
/// following:
/// 1. Bernoulli sampling, if `target` is a [`Percent`] and the population
///    `seq_count` is [`None`].
/// 2. Method D sampling, if `target` is a [`Count`] and the population
///    `seq_count` is [`Some`].
/// 3. Resovoir sampling (method L), if `target` is a [`Count`] and the
///    population `seq_count` is [`None`].
///
/// [`Percent`]: SamplingTarget::Percent
/// [`Count`]: SamplingTarget::Count
fn sample_and_write_records<I, W, A, E>(
    iterator: I, writer: W, target: SamplingTarget, seq_count: Option<usize>, mut rng: Xoshiro256StarStar,
) -> std::io::Result<()>
where
    for<'a> BernoulliSampler<'a, I>: WriteRecords<W>,
    for<'a> MethodDSampler<'a, I>: WriteRecords<W>,
    for<'a> std::vec::IntoIter<A>: WriteRecords<W>,
    I: Iterator<Item = Result<A, E>>,
    std::io::Error: From<E>, {
    match target {
        SamplingTarget::Percent(percent) => BernoulliSampler::new(iterator, percent, &mut rng).write_records(writer),
        SamplingTarget::Count(target) => {
            if let Some(total_items) = seq_count {
                MethodDSampler::new(iterator, target, total_items, &mut rng)?.write_records(writer)
            } else {
                method_l(iterator, &mut rng, target)?.into_iter().write_records(writer)
            }
        }
    }
}

fn get_paired_seq_count(
    fastq_path1: Option<PathBuf>, fastq_path2: Option<PathBuf>, reader1: &FastXReader<ReadFileZip>,
) -> std::io::Result<Option<usize>> {
    if let Some(path) = fastq_path1 {
        Ok(Some(get_seq_count(&path, reader1)?))
    } else if let Some(path) = fastq_path2 {
        Ok(Some(get_seq_count(&path, reader1)?))
    } else {
        Ok(None)
    }
}

struct IOArgs {
    /// This is only `Some` if the path corresponds to a file and is not zipped
    input_path1: Option<PathBuf>,
    /// This is only `Some` if paired ends are used, the path corresponds to a
    /// non-zipped file
    input_path2: Option<PathBuf>,
    reader1:     FastXReader<ReadFileZip>,
    reader2:     Option<FastXReader<ReadFileZip>>,
    writer:      RecordWriters<WriteFileZipStdout>,
}

/// The target number of sequences to sample
enum SamplingTarget {
    /// The target specified as a percent of the input number of sequences
    Percent(usize),
    /// The target as an exact count
    Count(usize),
}

fn parse_sampler_args(args: SamplerArgs) -> Result<(IOArgs, Xoshiro256StarStar, IoThreads, SamplingTarget), std::io::Error> {
    let rng = if let Some(seed) = &args.rng_seed {
        Xoshiro256StarStar::seed_from_u64(*seed)
    } else {
        Xoshiro256StarStar::from_os_rng()
    };

    let (fastq_reader1, fastq_reader2, writer, threads) = get_paired_readers_and_writers(
        &args.input_file1,
        args.input_file2.as_ref(),
        args.output_file1,
        args.output_file2,
    )?;

    let fastq_path1 = if args.input_file1.is_file() && !is_gz(&args.input_file1) {
        Some(args.input_file1)
    } else {
        None
    };
    let fastq_path2 = if let Some(path) = args.input_file2
        && path.is_file()
        && !is_gz(&path)
    {
        Some(path)
    } else {
        None
    };

    let io_args = IOArgs {
        input_path1: fastq_path1,
        input_path2: fastq_path2,
        reader1: fastq_reader1,
        reader2: fastq_reader2,
        writer,
    };
    let target = if let Some(count) = args.subsample_target {
        SamplingTarget::Count(count)
    } else if let Some(percent) = args.percent_target {
        SamplingTarget::Percent(percent)
    } else {
        unreachable!("This can't be reached because clap requires a value for either count or percent")
    };
    Ok((io_args, rng, threads, target))
}

fn get_seq_count(fastq: &PathBuf, reader: &FastXReader<ReadFileZip>) -> std::io::Result<usize> {
    let input = File::open(fastq)?;
    let buffered = BufReader::new(input);
    let line_count = buffered.lines().count();
    match reader {
        FastXReader::Fasta(_) => Ok(line_count / 2),
        FastXReader::Fastq(_) => Ok(line_count / 4),
    }
}

/// Method L (<https://dl.acm.org/doi/pdf/10.1145/198429.198435>)
fn method_l<I, T, E>(iter: I, rng: &mut Xoshiro256StarStar, target: usize) -> std::io::Result<Vec<T>>
where
    I: Iterator<Item = Result<T, E>>,
    std::io::Error: From<E>, {
    let mut reservoir = Vec::with_capacity(target);
    let n_sample = target;

    // Initial calculation of W and S
    let mut r = rng.random::<f32>();
    let mut w = (r.ln() / n_sample as f32).exp();
    r = rng.random::<f32>();
    let mut s = (r.ln() / (1.0 - w).ln()).floor() as usize;

    // Iterate through FastQ
    for (i, sample) in iter.enumerate() {
        // Initialize the reservoir
        if i < n_sample {
            reservoir.push(sample?);
        }
        // s=0 case, finished skipping, make swap
        else if s == 0 {
            // Insert into a random position in reservoir
            if let Some(slot) = reservoir.choose_mut(rng) {
                *slot = sample?;
            }

            // Recalculate S and W
            r = rng.random::<f32>();
            w *= (r.ln() / n_sample as f32).exp();
            r = rng.random::<f32>();
            s = (r.ln() / (1.0 - w).ln()).floor() as usize;
        }
        // s>0 case, keep skipping
        else {
            sample?;
            s -= 1;
        }
    }

    if n_sample > reservoir.len() {
        eprintln!(
            "Sampler Warning: Target sample size ({n_sample}) was greater than population size ({}); no downsampling has occurred.",
            reservoir.len()
        );
    }

    Ok(reservoir)
}

/// An iterator providing sampling using the Bernoulli method.
struct BernoulliSampler<'a, I: Iterator> {
    reader: I,
    prob:   f32,
    rng:    &'a mut Xoshiro256StarStar,
}

impl<'a, I: Iterator> BernoulliSampler<'a, I> {
    /// Creates an iterator which randomly downsamples from `reader` using the
    /// Bernoulli method, keeping approximately the specified percent of items.
    #[inline]
    fn new(reader: I, percent: usize, rng: &'a mut Xoshiro256StarStar) -> Self {
        Self {
            reader,
            prob: percent as f32 / 100.0,
            rng,
        }
    }
}

impl<'a, I, T, E> Iterator for BernoulliSampler<'a, I>
where
    I: Iterator<Item = Result<T, E>>,
{
    type Item = I::Item;

    fn next(&mut self) -> Option<Self::Item> {
        loop {
            let item = unwrap_or_return_some_err!(self.reader.next()?);
            if self.rng.random::<f32>() < self.prob {
                return Some(Ok(item));
            }
        }
    }
}

/// An iterator providing sampling using Algorithm D (Vitter,
/// <https://dl.acm.org/doi/pdf/10.1145/23002.23003>).
struct MethodDSampler<'a, I: Iterator> {
    reader:               I,
    rng:                  &'a mut Xoshiro256StarStar,
    remaining_samples:    usize,
    remaining_population: usize,
    /// $1/\alpha$, where $\alpha$ controls when to call method A.
    alphainv:             usize,
    /// Whether the iterator has fallen back to method A
    use_method_a:         bool,

    /// Next available random variate, computed as a uniform value [0,1) raised
    /// to `1.0 / remaining_samples`
    vprime:          f32,
    /// Cached value equal to `remaining_population - remaining_samples + 1`.
    /// Any valid skip must be below this value.
    max_start_index: usize,
    /// Cached value equal to `1.0 / remaining_samples`
    rem_samples_inv: f32,
    /// Cached value equal to $\alpha$ times `remaining_samples`.
    threshold:       usize,
}

impl<'a, I: Iterator> MethodDSampler<'a, I> {
    pub fn new(reader: I, target: usize, total_items: usize, rng: &'a mut Xoshiro256StarStar) -> std::io::Result<Self> {
        if target > total_items {
            return Err(std::io::Error::other(format!(
                "The sample target size {target} should not be larger than the population size {total_items}!"
            )));
        }

        const ALPHAINV: usize = 13;

        let remaining_samples = target;
        let remaining_population = total_items;
        let max_start_index = remaining_population - remaining_samples + 1;
        let rem_samples_inv = (remaining_samples as f32).recip();
        let vprime = rng.random::<f32>().powf(rem_samples_inv);

        Ok(MethodDSampler {
            reader,
            remaining_samples,
            remaining_population,
            max_start_index,
            rem_samples_inv,
            vprime,
            threshold: ALPHAINV * target,
            use_method_a: false,
            alphainv: ALPHAINV,
            rng,
        })
    }

    /// Computes a random number of reads to skip, after which the next one can
    /// be yielded by [`next`]. If `None` is returned, then downsampling is
    /// complete.
    ///
    /// This is accomplished via rejection sampling.
    ///
    /// [`next`]: MethodDSampler::next
    fn next_skip(&mut self) -> Option<usize> {
        if self.remaining_samples == 0 {
            return None;
        }

        // Once we start using method A, finish using method A
        if self.use_method_a {
            return self.next_skip_method_a();
        }

        // If remaining_samples*alpha >= remaining_population enter method A for
        // efficiency. If only one sample remains, enter method A to pick it
        // uniformly
        if self.threshold >= self.remaining_population || self.remaining_samples == 1 {
            self.use_method_a = true;
            return self.next_skip_method_a();
        }

        // The value of rem_samples_inv if remaining_samples is decremented by 1
        let rem_samples_min1_inv = ((self.remaining_samples - 1) as f32).recip();

        loop {
            let mut x = self.remaining_population as f32 * (1.0 - self.vprime);
            let mut skip = x as usize;

            // Randomly generate skips until a valid one is found (it must be
            // less than max_start_index)
            while skip >= self.max_start_index {
                // Get a new random variate to try again
                self.vprime = self.rng.random::<f32>().powf(self.rem_samples_inv);
                x = self.remaining_population as f32 * (1.0 - self.vprime);
                skip = x as usize;
            }

            let y1 = (self.rng.random::<f32>() * self.remaining_population as f32 / self.max_start_index as f32)
                .powf(rem_samples_min1_inv);
            self.vprime = y1
                * (1.0 - x / self.remaining_population as f32)
                * (self.max_start_index as f32 / (self.max_start_index - skip) as f32);

            if self.vprime <= 1.0 {
                return self.yield_skip(skip, rem_samples_min1_inv);
            }

            let mut y2 = 1.0;
            let (mut bottom, limit) = if self.remaining_samples - 1 > skip {
                (
                    self.remaining_population - self.remaining_samples,
                    self.remaining_population - skip,
                )
            } else {
                (self.remaining_population - skip - 1, self.max_start_index)
            };

            for top in (limit..self.remaining_population).rev() {
                y2 *= top as f32 / bottom as f32;
                bottom -= 1;
            }

            if self.remaining_population as f32 / (self.remaining_population as f32 - x)
                >= y1 * y2.powf(rem_samples_min1_inv)
            {
                // Need to use rem_samples_min1_inv, since we are computing for
                // the next iteration where remaining_samples is one smaller
                self.vprime = self.rng.random::<f32>().powf(rem_samples_min1_inv);
                return self.yield_skip(skip, rem_samples_min1_inv);
            }

            self.vprime = self.rng.random::<f32>().powf(self.rem_samples_inv);
        }
    }

    /// Given a skip to return, update the remaining population and sample, and
    /// update any cached values.
    ///
    /// `None` is never returned; an option is used for convenience in
    /// [`next_skip`].
    ///
    /// [`next_skip`]: MethodDSampler::next_skip
    fn yield_skip(&mut self, skip: usize, rem_samples_min1_inv: f32) -> Option<usize> {
        // Decrement due to skip, and due to the item yielded
        self.remaining_population -= skip + 1;
        self.remaining_samples -= 1;

        // Update cached values
        self.rem_samples_inv = rem_samples_min1_inv;
        self.max_start_index -= skip;
        self.threshold -= self.alphainv;
        Some(skip)
    }

    /// Similar to [`next_skip`], but uses Algorithm A which is faster when
    /// `remaining_samples` is large relative to `remaining_population`.
    ///
    /// Rather than using rejection sampling (Algorithm D) which has overhead,
    /// Algorithm A uses inverse transform sampling. This requires checking each
    /// value of `skips` from 0 until the final value returned, so when
    /// `remaining_samples` is large (and hence skips tend to be small), this
    /// may be more efficient.
    ///
    /// [`next_skip`]: MethodDSampler::next_skip
    fn next_skip_method_a(&mut self) -> Option<usize> {
        if self.remaining_samples > 1 {
            let mut top = self.remaining_population - self.remaining_samples;
            let mut skip = 0;
            let variate = self.rng.random::<f32>();

            // Find largest quotient of falling factorials (quot) greater than
            // random variate
            let mut quot = top as f32 / self.remaining_population as f32;
            while quot > variate {
                skip += 1;
                top -= 1;
                self.remaining_population -= 1;
                quot *= top as f32 / self.remaining_population as f32;
            }

            self.remaining_population -= 1;
            self.remaining_samples -= 1;
            Some(skip)
        } else {
            // The last remaining sample is picked uniformly from those
            // remaining
            let skip = self.rng.random_range(0..self.remaining_population);
            self.remaining_population -= skip + 1;
            self.remaining_samples = 0;
            Some(skip)
        }
    }
}

impl<T, E, I> Iterator for MethodDSampler<'_, I>
where
    I: Iterator<Item = Result<T, E>>,
{
    type Item = I::Item;

    /// Calculates the number of skips and then skips it for the iterator inside
    fn next(&mut self) -> Option<Self::Item> {
        // makes sure to validate all remaining reads even if the target number
        // of reads is reached
        let Some(skip) = self.next_skip() else {
            unwrap_or_return_some_err!(self.reader.try_for_each(|read| read.map(|_| ())));
            return None;
        };
        for _ in 0..skip {
            unwrap_or_return_some_err!(self.reader.next()?);
        }
        self.reader.next()
    }
}
