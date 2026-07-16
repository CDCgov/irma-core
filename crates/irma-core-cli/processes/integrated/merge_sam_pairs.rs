//! Merges Illumina paired-end reads with parsimonious error correction and
//! detection.

use clap::Args;
use irma_records::{
    hashing::get_hasher,
    io::{InputOptions, OutputOptions, ValidatePaths},
    paired::get_molecular_id_side,
    sam::{PairedMergeStats, SamMergeablePairs},
};
use std::{collections::HashMap, io::Write, path::PathBuf};
use zoe::data::{sam::*, views::Len};

#[derive(Args, Debug)]
pub struct MergeSAMArgs {
    /// Reference file used to generate the SAM.
    fasta_reference: PathBuf,

    /// SAM file to merge R1 and R2 pairs via alignment and parsimonious
    /// correction.
    sam_file: PathBuf,

    /// Output directory and prefix for merged SAM data.
    output_prefix: PathBuf,

    #[arg(short = 'S', long)]
    /// Serialize output observations for downstream analysis.
    store_stats: bool,

    #[arg(short = 'B', long)]
    /// SAM is in bowtie format.
    bowtie_format: bool,
}

struct ParsedMergeSamArgs {
    /// Reference file used to generate the SAM.
    fasta_reference: PathBuf,

    /// SAM file to merge R1 and R2 pairs via alignment and parsimonious
    /// correction.
    sam_file: PathBuf,

    /// The path for the output SAM file.
    merged_sam_file: PathBuf,

    /// If `Some`, the file to output observations for downstream analysis.
    paired_stats_file: Option<PathBuf>,

    /// SAM is in bowtie format.
    bowtie_format: bool,
}

fn parse_merge_sam_args(args: MergeSAMArgs) -> ParsedMergeSamArgs {
    ParsedMergeSamArgs {
        fasta_reference:   args.fasta_reference,
        sam_file:          args.sam_file,
        merged_sam_file:   args.output_prefix.with_extension("sam"),
        paired_stats_file: args.store_stats.then(|| args.output_prefix.with_extension("stats")),
        bowtie_format:     args.bowtie_format,
    }
}

impl ValidatePaths for ParsedMergeSamArgs {
    fn inputs(&self) -> impl IntoIterator<Item = &PathBuf> {
        [&self.fasta_reference, &self.sam_file]
    }

    fn outputs(&self) -> impl IntoIterator<Item = &PathBuf> {
        let merged_sam_file = std::iter::once(&self.merged_sam_file);
        let paired_stats_file = self.paired_stats_file.iter();

        merged_sam_file.chain(paired_stats_file)
    }
}

pub fn merge_sam_pairs_process(args: MergeSAMArgs) -> Result<(), std::io::Error> {
    let args = parse_merge_sam_args(args);

    args.validate_paths()?;

    let mut ref_reader = InputOptions::new_from_path(&args.fasta_reference)
        .use_file()
        .parse_fasta()
        .open()?;

    let reference = ref_reader.next().ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("No valid fasta data in file: {file}", file = args.fasta_reference.display()),
        )
    })??;

    const ONE_MB: usize = 2usize.pow(20);

    let mut sam_writer = OutputOptions::new_from_path(&args.merged_sam_file)
        .with_capacity(ONE_MB)
        .use_file()
        .open()?;

    let mut sam_data: Vec<SamData> = Vec::new();
    let mut pairs: HashMap<String, IndexPair, _> = HashMap::with_hasher(get_hasher());
    let mut index = 0;

    let sam_records = InputOptions::new_from_path(&args.sam_file).use_file().parse_sam().open()?;

    for sam_row in sam_records {
        let row = match sam_row? {
            SamRow::Data(d) => {
                if d.rname != reference.name {
                    continue;
                }
                d
            }
            SamRow::Header(h) => {
                writeln!(sam_writer, "{h}")?;
                continue;
            }
        };

        match get_molecular_id_side(&row.qname, '0') {
            Some((mol_name_id, read_side @ ('1' | '2'))) => {
                let new_pair = read_side.into_indexpair(index);
                pairs
                    .entry(mol_name_id.to_string())
                    .and_modify(|old_pair| old_pair.merge(&new_pair))
                    .or_insert(new_pair);
            }
            // If the header is parseable but does not encode an explicit read
            // side, fall back to order-based pairing instead of panicking.
            Some((mol_name_id, _)) => {
                pairs
                    .entry(mol_name_id.to_string())
                    .and_modify(|pair| pair.merge(&IndexPair::new_r2(index)))
                    .or_insert(IndexPair::new_r1(index));
            }
            None => {
                pairs
                    .entry(row.qname.clone())
                    .and_modify(|pair| pair.merge(&IndexPair::new_r2(index)))
                    .or_insert(IndexPair::new_r1(index));
            }
        }

        sam_data.push(row);
        index += 1;
    }

    // Store statistics: Observations, deletion minor variants, true SNV, false
    // SNV, insertion observations, insertion discrepancy.
    let mut paired_merging_stats = PairedMergeStats::default();

    for pair in pairs.values() {
        match (pair.r1, pair.r2) {
            (Some(pair_index1), Some(pair_index2)) => {
                let (sam1, sam2) = (&sam_data[pair_index1], &sam_data[pair_index2]);

                // IRMA does not define read-pair merging yet for the empty quality score case.
                // TODO: in v0.0.32 Zoe will only require checking for empty
                if !sam1.qual.is_empty()
                    && !sam2.qual.is_empty()
                    && sam1.qual.as_bytes() != b"*"
                    && sam2.qual.as_bytes() != b"*"
                {
                    let (s, stats) = sam1.merge_pair_using_reference(sam2, &reference.sequence, args.bowtie_format);
                    paired_merging_stats += stats;

                    writeln!(sam_writer, "{s}")?;
                } else {
                    writeln!(sam_writer, "{sam1}")?;
                    writeln!(sam_writer, "{sam2}")?;
                }
            }
            (Some(index), None) | (None, Some(index)) => {
                writeln!(sam_writer, "{}", sam_data[index])?;
            }
            _ => continue,
        }
    }

    if let Some(paired_stats_file) = args.paired_stats_file {
        let mut w = OutputOptions::new_from_path(&paired_stats_file).use_file().open()?;

        let PairedMergeStats {
            observations,
            true_variations,
            variant_errors,
            deletion_errors,
            insert_obs,
            insert_errors,
        } = paired_merging_stats;

        writeln!(
            &mut w,
            "{name}\tobs\t{observations}\n\
             {name}\ttmv\t{true_variations}\n\
             {name}\tfmv\t{variant_errors}\n\
             {name}\tdmv\t{deletion_errors}\n\
             {name}\tinsObs\t{insert_obs}\n\
             {name}\tinsErr\t{insert_errors}",
            name = reference.name
        )?;

        w.flush()?;
    }

    sam_writer.flush()
}

#[derive(Debug)]
struct IndexPair {
    r1: Option<usize>,
    r2: Option<usize>,
}

impl IndexPair {
    fn merge(&mut self, other: &IndexPair) {
        if let (None, Some(p)) = (self.r1, other.r1) {
            self.r1 = Some(p);
        }

        if let (None, Some(p)) = (self.r2, other.r2) {
            self.r2 = Some(p);
        }
    }

    fn new_r1(index: usize) -> Self {
        Self {
            r1: Some(index),
            r2: None,
        }
    }

    fn new_r2(index: usize) -> Self {
        Self {
            r1: None,
            r2: Some(index),
        }
    }
}

trait IntoIndexPair {
    fn into_indexpair(self, index: usize) -> IndexPair;
}

impl IntoIndexPair for char {
    fn into_indexpair(self, index: usize) -> IndexPair {
        match self {
            '1' => IndexPair {
                r1: Some(index),
                r2: None,
            },
            '2' => IndexPair {
                r1: None,
                r2: Some(index),
            },
            _ => panic!("Unexpected value for SAM pair side for read-pair merging: {self}"),
        }
    }
}
