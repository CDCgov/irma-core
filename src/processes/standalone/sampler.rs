use crate::io::{IoThreads, ReadFileZip, WriteFileZipStdout, create_writer, open_fastq_files};
use clap::Args;
use rand::{Rng, SeedableRng, seq::IndexedMutRandom};
use rand_xoshiro::Xoshiro256StarStar;
use std::{
    fs::File,
    io::{BufRead, BufReader, Write},
    path::PathBuf,
};
use zoe::{prelude::*, unwrap_or_return_some_err};

#[derive(Args, Debug)]
pub struct SamplerArgs {
    /// Path to .fastq or .fastq.gz file to be trimmed
    pub fastq_input_file: PathBuf,

    /// Path to optional second .fastq or .fastq.gz file to be trimmed
    pub fastq_input_file2: Option<PathBuf>,

    #[arg(short = '1', short_alias = 'o', long = "fastq-output")]
    /// Output filepath for trimmed reads. Trimmed reads print to STDOUT if not
    /// provided. May also use '-o'.
    pub fastq_output_file: Option<PathBuf>,

    #[arg(short = '2', long = "fastq-output2", requires = "fastq_input_file2")]
    /// Output path for secondary trimmed file if using paired reads. If this
    /// argument is omitted, output is interleaved.
    pub fastq_output_file2: Option<PathBuf>,

    #[arg(short = 't', long, group = "target")]
    /// Target number of reads to be subsampled
    pub subsample_target: Option<usize>,

    #[arg(short = 'p', long, group = "target", value_parser = validate_percent,)]
    /// Target percentage of reads to be subsampled. Must be a positive integer [0, 100]
    pub percent_target: Option<usize>,

    #[arg(short = 's', long)]
    /// Rng seed
    pub rng_seed: Option<u64>,
}

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
    let (io_args, mut rng, threads, target) = parse_sampler_args(args)?;
    let mut writer = io_args.output_writer;

    match (io_args.fastq_reader2, io_args.output_writer2) {
        (None, None) => {
            if io_args.fastq_path1.is_file() {
                let seq_count = get_seq_count(&io_args.fastq_path1).unwrap();
                let target = {
                    match target {
                        SamplingTarget::Count(count) => count,
                        SamplingTarget::Percent(percent) => seq_count * percent / 100,
                    }
                };
                if target > seq_count {
                    eprintln!(
                        "Sampler Warning: Target sample size ({target}) was greater than population size ({seq_count}); no downsampling has occurred.",
                    );
                    for sample in io_args.fastq_reader1 {
                        write!(writer, "{}", sample?)?;
                    }
                } else {
                    let subsampled = method_d(io_args.fastq_reader1, &mut rng, target, seq_count)?;
                    for sample in subsampled {
                        write!(writer, "{}", sample?)?;
                    }
                }
            } else {
                match target {
                    SamplingTarget::Count(count) => {
                        let reservoir = method_l(io_args.fastq_reader1, &mut rng, count)?;
                        for sample in reservoir {
                            write!(writer, "{sample}")?;
                        }
                    }
                    SamplingTarget::Percent(percent) => {
                        let bernoulli_iterator = method_bernoulli(io_args.fastq_reader1, &mut rng, percent);
                        for sample in bernoulli_iterator {
                            write!(writer, "{}", sample?)?;
                        }
                    }
                }
            }
        }
        (Some(_reader2), None) => todo!("subsample and interleave"), // use the zipped thing,
        (Some(_reader2), Some(_writer2)) => todo!("subsample zipped"),
        (None, Some(_writer2)) => todo!("de-interleave"),
    };

    writer.flush()?;
    threads.finalize()?;
    Ok(())
}

struct IOArgs {
    fastq_path1:    PathBuf,
    fastq_reader1:  FastQReader<ReadFileZip>,
    fastq_reader2:  Option<FastQReader<ReadFileZip>>,
    output_writer:  WriteFileZipStdout,
    output_writer2: Option<WriteFileZipStdout>,
}

enum SamplingTarget {
    Percent(usize),
    Count(usize),
}

fn parse_sampler_args(args: SamplerArgs) -> Result<(IOArgs, Xoshiro256StarStar, IoThreads, SamplingTarget), std::io::Error> {
    let rng = if let Some(seed) = &args.rng_seed {
        Xoshiro256StarStar::seed_from_u64(*seed)
    } else {
        Xoshiro256StarStar::from_os_rng()
    };

    let (fastq_reader1, fastq_reader2, threads) = open_fastq_files(&args.fastq_input_file, args.fastq_input_file2.as_ref())?;

    let output_writer = create_writer(args.fastq_output_file.clone())?;
    let output_writer2 = match &args.fastq_output_file2 {
        Some(path) => Some(create_writer(Some(path))?),
        None => None,
    };
    let io_args = IOArgs {
        fastq_path1: args.fastq_input_file,
        fastq_reader1,
        fastq_reader2,
        output_writer,
        output_writer2,
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

fn get_seq_count(fastq: &PathBuf) -> std::io::Result<usize> {
    let input = File::open(fastq)?;
    let buffered = BufReader::new(input);
    let line_count = buffered.lines().count();
    Ok(line_count / 4)
}

/// Method L (https://dl.acm.org/doi/pdf/10.1145/198429.198435)
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

fn method_bernoulli<I, T, R>(iter: I, rng: &mut R, percent: usize) -> impl Iterator<Item = T>
where
    I: Iterator<Item = T>,
    R: Rng, {
    let percent = percent as f32 / 100.0;
    iter.filter(move |_| rng.random::<f32>() < percent)
}

/// Algorithm D (Vitter, https://dl.acm.org/doi/pdf/10.1145/23002.23003)
/// randomly produces how many reads to *skip* based on the total amount of
/// reads remaining
fn method_d<'a>(
    reader: FastQReader<ReadFileZip>, rng: &'a mut Xoshiro256StarStar, target: usize, seq_count: usize,
) -> std::io::Result<MethodDSampler<'a, FastQReader<ReadFileZip>>> {
    MethodDSampler::new(reader, target, seq_count, rng)
}

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
        let vprime = rand::random::<f32>().powf(rem_samples_inv);

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
        let skip = self.next_skip()?;
        for _ in 0..skip {
            unwrap_or_return_some_err!(self.reader.next()?);
        }
        self.reader.next()
    }
}
