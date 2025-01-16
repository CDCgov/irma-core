// Description:      Read FastQ files and trim with various options.
use clap::Args;
use either::Either;
use foldhash::fast::RandomState;
use std::{
    fs::OpenOptions,
    io::{stdin, BufReader, BufWriter, Write},
    path::PathBuf,
};
use zoe::{data::types::nucleotides::reverse_complement, kmer::ThreeBitKmerSet, prelude::*};

use crate::qc::fastq::ReadTransforms;

// An example run:
// cargo run -- trimmer 3003863475_N8KHVRSA_S384_R1_001.fastq -k swift_211206.fasta -Z -c AGATGTGTATAAGAGACAG -U -H 4 > /dev/null

// Benchmarking:
// hyperfine -r 30 -w 1 -N \
// './irma-core-view trimmer 3003863475_N8KHVRSA_S384_R1_001.fastq -k swift_211206.fasta -Z -c AGATGTGTATAAGAGACAG -U -H 4' \
// './irma-core-orig trimmer 3003863475_N8KHVRSA_S384_R1_001.fastq -k swift_211206.fasta -Z -c AGATGTGTATAAGAGACAG -U -H 4'

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
    const KMER_LENGTH: usize = 17;
    const PRIMER_RESTRICT_LEFT: usize = 30;
    const PRIMER_RESTRICT_RIGHT: usize = 30;

    let fastq_file_reader = if let Some(ref file_path) = args.fastq_input_file {
        FastQReader::new(BufReader::new(Either::Left(OpenOptions::new().read(true).open(file_path)?)))
    } else {
        FastQReader::new(BufReader::new(Either::Right(stdin())))
    };
    let mut stdout_writer = BufWriter::new(std::io::stdout());

    let fasta_primer_reader = FastaReader::new(BufReader::new(OpenOptions::new().read(true).open(&args.fasta_primer_file)?));

    let mut unique_kmers = ThreeBitKmerSet::<21, _>::with_hasher(KMER_LENGTH, RandomState::default()).unwrap();

    fasta_primer_reader.into_iter().for_each(|f| {
        let seq = f.unwrap().sequence;
        if seq.len() > KMER_LENGTH {
            let rev_comp = reverse_complement(&seq);
            unique_kmers.insert_from_sequence_one_mismatch(seq);
            unique_kmers.insert_from_sequence_one_mismatch(rev_comp);
        }
    });

    let (forward_adapter, reverse_adapter) = match (&args.mask_adapter, &args.clip_adapter) {
        (Some(ref a), _) | (_, Some(ref a)) => {
            let forward = a.as_bytes().to_ascii_uppercase();
            let reverse = reverse_complement(&forward);
            (forward, reverse)
        }
        _ => (Vec::new(), Vec::new()),
    };

    let mut i = 0;
    let mut chopped_left = 0;
    let mut chopped_right = 0;
    for record in fastq_file_reader {
        let mut fq = record?;
        let mut fq = fq.as_view_mut();

        fq.to_canonical_bases(args.canonical_bases)
            .transform_by_reverse_forward_search(
                args.fuzzy_adapter,
                args.clip_adapter.is_some(),
                &reverse_adapter,
                &forward_adapter,
            );

        let left_len = fq.sequence.len().min(PRIMER_RESTRICT_LEFT);
        let left_end = &fq.sequence[..left_len];
        if let Some(range) = unique_kmers.find_kmers_rev(left_end) {
            fq.restrict(range.end..);
            chopped_left += 1;
        }

        let right_len = fq.sequence.len().min(PRIMER_RESTRICT_RIGHT);
        let right_starting_idx = fq.sequence.len() - right_len;
        let right_end = &fq.sequence[right_starting_idx..];
        if let Some(range) = unique_kmers.find_kmers(right_end) {
            let new_start = right_starting_idx + range.start;
            fq.restrict(..new_start);
            chopped_right += 1;
        }

        fq.hard_trim(args.hard_trim);

        i += 1;
        write!(stdout_writer, "{fq}")?;
    }

    eprintln!("Processed {i} reads.");
    let percent_left = chopped_left as f64 / i as f64;
    let percent_right = chopped_right as f64 / i as f64;
    eprintln!("Chopped {} reads on left ({}%)", chopped_left, percent_left);
    eprintln!("Chopped {} reads on right ({}%)", chopped_right, percent_right);

    Ok(())
}
