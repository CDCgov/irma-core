// Description:      Merges Illumina paired-end reads with parsimonious error
//                   correction and detection.
//

use clap::Args;
use std::collections::HashMap;
use std::io::Write;
use std::path::PathBuf;
use zoe::data::err::ResultWithErrorContext;
use zoe::data::sam::*;

use crate::io::{InputOptions, OutputOptions};
use crate::utils::paired_reads::get_molecular_id_side;

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

pub fn merge_sam_pairs_process(args: MergeSAMArgs) -> Result<(), std::io::Error> {
    let mut ref_reader = InputOptions::new_from_path(&args.fasta_reference)
        .use_file()
        .parse_fasta()
        .open()?;

    let reference = ref_reader.next().ok_or_else(|| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidInput,
            format!("No valid fasta data in file: {file}", file = &args.fasta_reference.display()),
        )
    })??;

    const ONE_MB: usize = 2usize.pow(20);
    let merged_sam_file = args.output_prefix.with_extension("sam");

    let mut sam_writer = OutputOptions::new_from_path(&merged_sam_file)
        .with_capacity(ONE_MB)
        .use_file()
        .open()?;

    let mut sam_data: Vec<SamData> = Vec::new();
    let mut pairs: HashMap<String, IndexPair> = HashMap::new();
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
                writeln!(sam_writer, "{h}")
                    .with_file_context("Failed to write header to merged SAM file", &merged_sam_file)?;
                continue;
            }
        };

        // bowtie is order based
        if let Some((mol_name_id, read_side)) = get_molecular_id_side(&row.qname, '0') {
            let new_pair = read_side.into_indexpair(index);

            pairs
                .entry(mol_name_id.to_string())
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

                writeln!(sam_writer, "{s}").with_file_context("Failed to write to merged SAM file", &merged_sam_file)?;
            }
            (Some(index), None) | (None, Some(index)) => {
                writeln!(sam_writer, "{}", sam_data[index])
                    .with_file_context("Failed to write to merged SAM file", &merged_sam_file)?;
            }
            _ => continue,
        }
    }

    if args.store_stats {
        let paired_stats_file = args.output_prefix.with_extension("stats");

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
        )
        .with_file_context("Failed to write paired stats file", &paired_stats_file)?;
    }

    Ok(())
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

#[cfg(test)]
mod test {
    use crate::merge_sam_pairs::get_molecular_id_side;

    static QNAMES: [&str; 26] = [
        "SRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "SRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "SRR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "SRR26182418.1.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "A00350:691:HCKYLDSX3:2:2119:23863:2456/2",
        "A00350:691:HCKYLDSX3:2:2119:23863:2456/1",
        "M02989:9:000000000-L4PJL:1:2112:9890:15606 1:N:0:AACGCACGAG+GCCTCGGATA",
        "M02989:9:000000000-L4PJL:1:2112:9890:15606 2:N:0:AACGCACGAG+GCCTCGGATA",
        "NS500500:69:HKJFLAFX5:1:11204:14878:14643 1:N:0:TTCTCGTGCA+CTCTGTGTAT",
        "NS500500:69:HKJFLAFX5:1:11204:14878:14643 2:N:0:TTCTCGTGCA+CTCTGTGTAT",
        "A01000:249:HJFFWDRX2:1:2107:24605:18082 1:N:0:TAGGCATG+ATAGCCTT",
        "A01000:249:HJFFWDRX2:1:2107:24605:18082 2:N:0:TAGGCATG+ATAGCCTT",
        "M02989:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
        "M02989:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
        "M02989_1:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
        "M02989_1:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
        "SRR26182418.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
        "SRR26182418.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
        "SRR26182418.1.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
        "SRR26182418.1.2_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
        "SRR26182418.1 1:N:18:NULL",
        "SRR26182418.1.1 1:N:18:NULL",
        "ERR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "DRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "ERR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "DRR26182418.2.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
    ];

    #[test]
    fn test_get_mol_id() {
        let extracted = [
            ("SRR26182418.1", '0'),
            ("SRR26182418.1", '0'),
            ("SRR26182418.1", '1'),
            ("SRR26182418.1", '2'),
            ("A00350:691:HCKYLDSX3:2:2119:23863:2456", '2'),
            ("A00350:691:HCKYLDSX3:2:2119:23863:2456", '1'),
            ("M02989:9:000000000-L4PJL:1:2112:9890:15606", '1'),
            ("M02989:9:000000000-L4PJL:1:2112:9890:15606", '2'),
            ("NS500500:69:HKJFLAFX5:1:11204:14878:14643", '1'),
            ("NS500500:69:HKJFLAFX5:1:11204:14878:14643", '2'),
            ("A01000:249:HJFFWDRX2:1:2107:24605:18082", '1'),
            ("A01000:249:HJFFWDRX2:1:2107:24605:18082", '2'),
            ("M02989:9:000000000-L4PJL:1:2114:17393:19614", '1'),
            ("M02989:9:000000000-L4PJL:1:2114:17393:19614", '2'),
            ("M02989_1:9:000000000-L4PJL:1:2114:17393:19614", '1'),
            ("M02989_1:9:000000000-L4PJL:1:2114:17393:19614", '2'),
            ("SRR26182418.1", '0'),
            ("SRR26182418.1", '0'),
            ("SRR26182418.1", '1'),
            ("SRR26182418.1", '2'),
            ("SRR26182418.1", '0'),
            ("SRR26182418.1", '1'),
            ("ERR26182418.1", '0'),
            ("DRR26182418.1", '0'),
            ("ERR26182418.1", '1'),
            ("DRR26182418.2", '2'),
        ];

        for (i, o) in QNAMES.iter().enumerate() {
            assert_eq!(get_molecular_id_side(o, '0'), Some(extracted[i]), "'{o}'");
        }
    }
}
