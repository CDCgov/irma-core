// Filename:         merge_sam_pairs
//
// Description:      Merges Illumina paired-end reads with parsimonious error
//                   correction and detection.
//
// Date dedicated:   2023-08-22
// Author:           Samuel S. Shepard, Centers for Disease Control and Prevention
//
// Citation:         Shepard SS, Meno S, Bahl J, Wilson MM, Barnes J, Neuhaus E.
//                   Viral deep sequencing needs an adaptive approach: IRMA, the
//                   iterative refinement meta-assembler. BMC Genomics.
//                   2016;17(1). doi:10.1186/s12864-016-3030-6
//
// =============================================================================
//
//                            PUBLIC DOMAIN NOTICE
//
//  This source code file or script constitutes a work of the United States
//  Government and is not subject to domestic copyright protection under 17 USC §
//  105. This file is in the public domain within the United States, and
//  copyright and related rights in the work worldwide are waived through the CC0
//  1.0 Universal public domain dedication:
//  https://creativecommons.org/publicdomain/zero/1.0/
//
//  The material embodied in this software is provided to you "as-is" and without
//  warranty of any kind, express, implied or otherwise, including without
//  limitation, any warranty of fitness for a particular purpose. In no event
//  shall the Centers for Disease Control and Prevention (CDC) or the United
//  States (U.S.) government be liable to you or anyone else for any direct,
//  special, incidental, indirect or consequential damages of any kind, or any
//  damages whatsoever, including without limitation, loss of profit, loss of
//  use, savings or revenue, or the claims of third parties, whether or not CDC
//  or the U.S. government has been advised of the possibility of such loss,
//  however caused and on any theory of liability, arising out of or in
//  connection with the possession, use or performance of this software.
//
//  Please provide appropriate attribution in any work or product based on this
//  material.
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
        .unwrap_or_die(&format!(
            "cannot open reference file '{}'\n  In: {MODULE}",
            &args.fasta_reference.display()
        ))
        .filter_map(|f| f.ok())
        .next()
        .unwrap_or_else(|| {
            panic!(
                "Error: no valid fasta data in file '{}'\n  In: {MODULE}",
                &args.fasta_reference.display()
            )
        });

    const ONE_MB: usize = 2usize.pow(20);
    let merged_sam_file = args.output_prefix.with_extension("sam");
    let mut sam_writer = BufWriter::with_capacity(
        ONE_MB,
        std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .open(&merged_sam_file)
            .unwrap_or_die(&format!(
                "cannot open output SAM file '{}'\n  In: {MODULE}",
                merged_sam_file.display()
            )),
    );

    // Make lazy static for a lib function
    let mol_name_id = Regex::new(r"(.+?)[_ ]([12]):.+")
        .unwrap_or_else(|e| panic!("REGEX molecular ID didn't compile: {e}\n  In: {MODULE}"));

    let mut sam_data: Vec<SamData> = Vec::new();
    let mut pairs: HashMap<String, IndexPair> = HashMap::new();
    let mut index = 0;

    // TO-DO: could be more functional
    // TO-DO: push down a reference
    let sam_records = SAMReader::from_filename(&args.sam_file).unwrap_or_die(&format!(
        "cannot read SAM file '{}'\n  In: {MODULE}",
        &args.sam_file.display()
    ));
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
                    "failed to write header to merged sam file: {}\n  In: {MODULE}",
                    merged_sam_file.display()
                ));
                continue;
            }
            Err(e) => Err(e).unwrap_or_die(&format!("parsing SAM row failed\n  In: {MODULE}")),
        };

        // bowtie is order based
        if let Some(cap) = mol_name_id.captures_iter(&row.qname).next() {
            let new_pair = cap[2].into_indexpair(index);

            pairs
                .entry(cap[1].to_string())
                .and_modify(|old_pair| old_pair.merge(&new_pair))
                .or_insert(new_pair);
        } else {
            pairs
                .entry(row.qname.clone())
                .and_modify(|r1| r1.merge(&IndexPair::new_r2(index)))
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

                writeln!(sam_writer, "{s}").unwrap_or_die(&format!(
                    "failed to write to paired sam file: {}\n  In: {MODULE}",
                    merged_sam_file.display()
                ));
            }
            (Some(index), None) | (None, Some(index)) => {
                writeln!(sam_writer, "{}", sam_data[index]).unwrap_or_die(&format!(
                    "failed to write to merged sam file: {}\n  In: {MODULE}",
                    merged_sam_file.display()
                ));
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
                .unwrap_or_die(&format!(
                    "cannot open stats file: {}\n  In: {MODULE}",
                    paired_stats_file.display()
                )),
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
        .unwrap_or_die(&format!(
            "failed to write paired stats file: {}\n  In: {MODULE}",
            paired_stats_file.display()
        ));
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
            _ => panic!("unexpected value for SAM pair\n  In: {MODULE}"),
        }
    }
}
