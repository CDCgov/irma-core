// Description:      Read FastQ files and trim with various options.
use clap::Args;
use either::Either;
use foldhash::fast::RandomState;
use std::{
    fs::OpenOptions,
    io::{stdin, BufReader, BufWriter, Write},
    path::PathBuf,
};
use zoe::{data::types::nucleotides::reverse_complement, kmer::ThreeBitOneMismatchKmerSet, prelude::KmerSet, prelude::*};

// use crate::qc::fastq::ReadTransforms;

/* Assumed args */

#[derive(Args, Debug)]
pub struct TrimmerArgs {
    pub fastq_input_file: Option<PathBuf>,

    #[arg(short = 'H', long)]
    /// Hard trim from each end the specified length.
    pub hard_trim: Option<usize>,

    #[arg(short = 'm', long)]
    /// Specify adapter sequence and mask when found in reads.
    pub mask_adapter: Option<String>,

    #[arg(short = 'c', long)]
    /// Specify adapter sequence and clip appropriate ends when found in reads.
    pub clip_adapter: Option<String>,

    #[arg(short = 'Z', long)]
    /// Allow up to one mismatch for adapter clipping (-c) or masking (-m).
    pub fuzzy_adapter: bool,

    #[arg(short = 'U', long)]
    /// Re-encode FASTQ sequence to expected input: A, C, T, G, N
    pub canonical_bases: bool,

    // primer trimming args to go here
    #[arg(short = 'k', long)]
    /// primers
    pub fasta_primer_file: PathBuf,
}

//static MODULE: &str = module_path!();

/// # Panics
///
/// Sub-program for trimming fastQ data.
pub fn trimmer_process(args: &TrimmerArgs) -> Result<(), std::io::Error> {
    let kmer_length = 17;

    let fastq_file_reader = if let Some(ref file_path) = args.fastq_input_file {
        FastQReader::new(BufReader::new(Either::Left(OpenOptions::new().read(true).open(file_path)?)))
    } else {
        FastQReader::new(BufReader::new(Either::Right(stdin())))
    };
    let mut stdout_writer = BufWriter::new(std::io::stdout());

    let fasta_primer_reader = FastaReader::new(BufReader::new(OpenOptions::new().read(true).open(&args.fasta_primer_file)?));

    let restrict_left = 30;
    let mut unique_kmers = ThreeBitOneMismatchKmerSet::<21, _>::with_hasher(kmer_length, RandomState::default()).unwrap();

    fasta_primer_reader.into_iter().for_each(|f| {
        let seq = f.unwrap().sequence;
        if seq.len() > kmer_length {
            let rev_comp = reverse_complement(&seq);
            unique_kmers.insert_from_sequence(seq);
            unique_kmers.insert_from_sequence(rev_comp);
        }
    });

    /*let (forward_adapter, reverse_adapter) = match (&args.mask_adapter, &args.clip_adapter) {
        (Some(ref a), _) | (_, Some(ref a)) => {
            let forward = a.as_bytes().to_ascii_uppercase();
            let reverse = reverse_complement(&forward);
            (forward, reverse)
        }
        _ => (Vec::new(), Vec::new()),
    };*/
    let mut i = 0;
    let mut chopped = 0;
    for record in fastq_file_reader {
        let mut fq = record?;

        /*fq.to_canonical_bases(args.canonical_bases)
            .transform_by_reverse_forward_search(
                args.fuzzy_adapter,
                args.clip_adapter.is_some(),
                &reverse_adapter,
                &forward_adapter,
            )
            .hard_trim(args.hard_trim);
        */
        if let Some(range) = unique_kmers.find_kmers_rev(&fq.sequence[0..restrict_left]) {
            fq.sequence.cut_to_start(range.end);
            fq.quality.cut_to_start(range.end);
            chopped += 1;
        }
        i += 1;
        write!(stdout_writer, "{fq}")?;
    }

    eprintln!("Processed {i} reads.");
    let percent = chopped as f64 / i as f64;
    eprintln!("Chopped {} reads ({}%)", chopped, percent);

    Ok(())
}
