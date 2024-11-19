// Description:      Read FastQ files and trim with various options.
#![allow(clippy::needless_range_loop)]
use clap::Args;
use either::Either;
use foldhash::HashSet;
use std::{
    fs::OpenOptions,
    io::{stdin, BufReader, BufWriter, Write},
    path::PathBuf,
};
use zoe::{data::to_dna_profile_index, data::types::nucleotides::reverse_complement, prelude::*};

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

const BITS_PER_BASE: usize = 3;
const MAX_KMER_SIZE: usize = 21;

// pre-computing a table of set masks to quickly change a base within a pre-encoded kmer
const SET_MASKS: [[u64; 5]; MAX_KMER_SIZE] = {
    let mut set_masks = [[0u64; 5]; MAX_KMER_SIZE];
    let mut index = 0;
    while index < MAX_KMER_SIZE {
        let mut base = 0;
        while base < 5 {
            set_masks[index][base] = (base << (3 * index)) as u64; // take the integer value of the base and move it over 3*i times to create a set mask for the given index, base pair
            base += 1;
        }
        index += 1;
    }
    set_masks
};

// pre-computing a table of clear masks to quickly clear a base position within a pre-encoded kmer
const CLEAR_MASKS: [u64; MAX_KMER_SIZE] = {
    let mut clear_masks = [0u64; 21];
    let mut index = 0;
    while index < MAX_KMER_SIZE {
        clear_masks[index] = !(0b111 << (3 * index)); // 111, we need to shift this 3*i times to get to the correct index, then negate it to put a 000, or empty base at that position
        index += 1;
    }
    clear_masks
};

/// Encodes a kmer into a u64 with 3-bit representations
fn encode_kmer(seq: &[u8]) -> u64 {
    let mask = (1u64 << (BITS_PER_BASE * seq.len())) - 1;

    let mut kmer: u64 = 0;

    for &base in seq {
        let base_val = to_dna_profile_index(base) as u64;
        kmer = ((kmer << BITS_PER_BASE) | base_val) & mask;
    }

    kmer
}

/// Finds all possible 1 hamming distance variants of a given kmer and inserts them into a hash set
/// Panics if primer.len() < kmer_length
fn mutate_primer_add_kmers(primer: &[u8], kmer_length: usize, hash: &mut HashSet<u64>) {
    let mut current_kmer = encode_kmer(&primer[primer.len() - kmer_length..]);
    for i in (0..=primer.len() - kmer_length).rev() {
        for index in 0..kmer_length {
            for base in 0..b"ACGTN".len() {
                let variant = (current_kmer & CLEAR_MASKS[index]) | SET_MASKS[index][base]; // clear the position of the changing base with &, then put the new base in with |
                hash.insert(variant);
            }
        }
        if i > 0 {
            let next_base = primer[i - 1];
            current_kmer = (current_kmer >> BITS_PER_BASE) | SET_MASKS[kmer_length - 1][to_dna_profile_index(next_base)];
            // slide the window to include the next kmer
        }
    }
}

/// moves right to left within a range of a slice, checking if each encoded kmer is contained within the hash set
/// if not, mutates the kmer using bitwise operations to slide the window
fn find_kmer(sequence: &[u8], kmer_length: usize, unique_kmers: &HashSet<u64>) -> Option<usize> {
    if sequence.len() < kmer_length {
        return None;
    }

    let start_index = sequence.len() - kmer_length;
    let mut current_kmer = encode_kmer(&sequence[start_index..]);

    if unique_kmers.contains(&current_kmer) {
        return Some(start_index);
    }

    //let mask = (1 << (3 * (kmer_length - 1))) - 1;
    //not needed in a right-to-left implementation

    for i in (0..start_index).rev() {
        let next_base = sequence[i];
        current_kmer = (current_kmer >> BITS_PER_BASE) | SET_MASKS[kmer_length - 1][to_dna_profile_index(next_base)];

        if unique_kmers.contains(&current_kmer) {
            return Some(i);
        }
    }

    None
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
    let mut unique_kmers = HashSet::default();

    fasta_primer_reader.into_iter().for_each(|f| {
        let seq = f.unwrap().sequence;
        if seq.len() > kmer_length {
            let rev_comp = reverse_complement(&seq);
            mutate_primer_add_kmers(&seq, kmer_length, &mut unique_kmers);
            mutate_primer_add_kmers(&rev_comp, kmer_length, &mut unique_kmers);
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

        if let Some(max_index) = find_kmer(&fq.sequence[0..restrict_left], kmer_length, &unique_kmers) {
            fq.sequence.cut_to_start(max_index + kmer_length);
            fq.quality.cut_to_start(max_index + kmer_length);
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
