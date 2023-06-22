#![allow(non_snake_case, unused_mut, dead_code, clippy::wrong_self_convention)]

const PROGRAM: &str = "mergeSAMpairs";

use clap::Args;
use indoc::writedoc;
use regex::Regex;
use std::collections::HashMap;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use zoe::data::{err::OrFail, fasta::FastaReader, sam::*};

#[derive(Args, Debug)]
pub struct MergeSAMArgs {
    /// Reference file used to generate the SAM.
    fasta_reference: PathBuf,

    /// SAM file to merge R1 and R2 pairs via alignment and parsimonious correction.
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

static MODULE: &str = module_path!();

pub fn merge_sam_pairs_process(args: &MergeSAMArgs) {
    let reference = FastaReader::from_filename(&args.fasta_reference)
        .unwrap_or_die(&format!("{MODULE}: cannot open reference file"))
        .filter_map(|f| f.ok())
        .next()
        .unwrap_or_else(|| panic!("{MODULE}: no valid fasta data"));

    const ONE_MB: usize = 2usize.pow(20);
    let merged_sam_file = args.output_prefix.with_extension("sam");
    let mut sam_writer = BufWriter::with_capacity(
        ONE_MB,
        std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&merged_sam_file)
            .unwrap_or_die(&format!("Cannot open output SAM file: {}", merged_sam_file.display())),
    );

    // Make lazy static for a lib function
    let MOL_NAME_ID = Regex::new(r"(.+?)[_ ]([12]):.+").expect("REGEX molecular ID didn't compile.");

    let mut sam_data: Vec<SamData> = Vec::new();
    let mut pairs: HashMap<String, IndexPair> = HashMap::new();
    let mut index = 0;

    // TO-DO: could be more functional
    // TO-DO: push down a reference
    let sam_records = SAMReader::from_filename(&args.sam_file).unwrap_or_die(&format!("{MODULE}: cannot read SAM file."));
    for sam_row in sam_records {
        let row = match sam_row {
            Ok(SamRow::Data(d)) => {
                // TO-DO: inconsistent, fix types
                if d.rname.as_bytes() != reference.name {
                    continue;
                }
                *d
            }
            Ok(SamRow::Header(h)) => {
                writeln!(sam_writer, "{h}").unwrap_or_die(&format!(
                    "Failed to write header to merged sam file: {}",
                    merged_sam_file.display()
                ));
                continue;
            }
            Err(e) => Err(e).unwrap_or_die(&format!("{MODULE}: parsing SAM row failed")),
        };

        // bowtie is order based
        if let Some(cap) = MOL_NAME_ID.captures_iter(&row.qname).next() {
            let new_pair = cap[2].into_indexpair(index);

            pairs
                .entry(cap[1].to_string())
                .and_modify(|mut old_pair| old_pair.merge(&new_pair))
                .or_insert(new_pair);
        } else {
            pairs
                .entry(row.qname.clone())
                .and_modify(|mut r1| r1.merge(&IndexPair::new_r2(index)))
                .or_insert(IndexPair::new_r1(index));
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
                let (s, stats) = sam_data[pair_index1].merge_pair_using_reference(
                    &sam_data[pair_index2],
                    &reference.sequence,
                    args.bowtie_format,
                );
                paired_merging_stats += stats;

                writeln!(sam_writer, "{s}")
                    .unwrap_or_die(&format!("Failed to write to paired sam file: {}", merged_sam_file.display()));
            }
            (Some(index), None) | (None, Some(index)) => {
                writeln!(sam_writer, "{}", sam_data[index])
                    .unwrap_or_die(&format!("Failed to write to merged sam file: {}", merged_sam_file.display()));
            }
            _ => continue,
        }
    }

    if args.store_stats {
        let paired_stats_file = args.output_prefix.with_extension("stats");
        let mut w = BufWriter::new(
            std::fs::OpenOptions::new()
                .write(true)
                .create(true)
                .truncate(true)
                .open(&paired_stats_file)
                .unwrap_or_die(&format!("Cannot open stats file: {}", paired_stats_file.display())),
        );

        let name = String::from_utf8_lossy(&reference.name);
        let PairedMergeStats {
            observations,
            true_variations,
            variant_errors,
            deletion_errors,
            insert_obs,
            insert_errors,
        } = paired_merging_stats;

        writedoc!(
            &mut w,
            "{name}\tobs\t{observations}
             {name}\ttmv\t{true_variations}
             {name}\tfmv\t{variant_errors}
             {name}\tdmv\t{deletion_errors}
             {name}\tinsObs\t{insert_obs}
             {name}\tinsErr\t{insert_errors}\n",
        )
        .unwrap_or_die(&format!("Failed to write paired stats file: {}", paired_stats_file.display()));
    }
}

#[derive(Debug)]
pub struct IndexPair {
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

    fn is_widowed(&self) -> bool {
        self.r1.is_some() ^ self.r2.is_some()
    }

    fn get_either(&self) -> Option<usize> {
        self.r1.or(self.r2)
    }
}

trait IntoIndexPair {
    fn into_indexpair(self, index: usize) -> IndexPair;
}

impl IntoIndexPair for &str {
    fn into_indexpair(self, index: usize) -> IndexPair {
        match self {
            "1" => IndexPair {
                r1: Some(index),
                r2: None,
            },
            "2" => IndexPair {
                r1: None,
                r2: Some(index),
            },
            _ => panic!("{MODULE}: unexpected value for SAM pair"),
        }
    }
}
