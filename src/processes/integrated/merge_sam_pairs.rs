// Description:      Merges Illumina paired-end reads with parsimonious error
//                   correction and detection.
//

use clap::Args;
use indoc::writedoc;
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

pub fn merge_sam_pairs_process(args: MergeSAMArgs) {
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
                if d.rname != reference.name {
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
            name = reference.name
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
            _ => panic!("unexpected value for SAM pair\n  In: {MODULE}"),
        }
    }
}

fn get_molecular_id_side(s: &str, default_side: char) -> Option<(&str, char)> {
    let (the_id, the_side) = if s.contains(' ') {
        let mut pieces = s.split(' ');
        let id = pieces.next().unwrap_or_default();

        if !(id.starts_with("SRR") || id.starts_with("DRR") || id.starts_with("ERR")) || !id.contains('.') {
            // Illumina format
            (
                id,
                pieces
                    .next()
                    .unwrap_or_default()
                    .split(':')
                    .next()
                    .unwrap_or_default()
                    .chars()
                    .next(),
            )
        } else if let Some(index) = id.match_indices('.').nth(1).map(|(i, _)| i) {
            // SRA format, read side included
            let (new_id, side) = id.split_at(index);
            (new_id, side.chars().nth(1))
        } else {
            // SRA format, no read side
            (id, Some(default_side))
        }
    } else if let Some(index) = s.find('/') {
        // Legacy Illumina
        let (new_id, side) = s.split_at(index);
        (new_id, side.chars().nth(1))
    } else if (s.starts_with("SRR") || s.starts_with("DRR") || s.starts_with("ERR")) && s.contains('.') {
        let mut pieces = s.split('_');
        let id = pieces.next().unwrap_or_default();

        if let Some(index) = id.match_indices('.').nth(1).map(|(i, _)| i) {
            // SRA with read side
            let (new_id, side) = id.split_at(index);
            (new_id, side.chars().nth(1))
        } else {
            // SRA, no read side
            (id, Some(default_side))
        }
    } else {
        // IRMA Illumina legacy output
        let mut indices = s.match_indices(':');
        let (left, right) = (indices.nth(5), indices.next());
        if let (Some((start, _)), Some((stop, _))) = (left, right)
            && let Some(us) = s[start..stop].find('_')
        {
            let underscore_index = start + us;
            (&s[..underscore_index], s[..stop].chars().next_back())
        } else {
            return None;
        }
    };

    if let (id, Some(side @ '0'..='3')) = (the_id, the_side) {
        Some((id, side))
    } else {
        Some((the_id, default_side))
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
