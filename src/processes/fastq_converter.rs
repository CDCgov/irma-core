// Filename:         fastq_converter
//
// Description:      Read FastQ files, applies QC filtering (quality and length),
//                   adapter trimming, and format conversion as requested.
//
// Date dedicated:   2022-09-30
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
#![allow(
    clippy::struct_excessive_bools,
    clippy::module_name_repetitions,
    clippy::too_many_lines,
    clippy::missing_errors_doc
)]

use clap::{Args, ValueHint};
use either::Either;
use regex::bytes::Regex;
use std::borrow::Borrow;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::{stdin, BufReader, BufWriter, Error as IOError, ErrorKind};
use std::path::PathBuf;
use zoe::data::fastq::{FastQ, FastQReader};
use zoe::data::types::nucleotides::reverse_complement;

#[derive(Args, Debug)]
pub struct FastqConverterArgs {
    fastq_input_file: Option<PathBuf>,

    #[arg(short = 'Q', long)]
    /// Outputs fastQ instead of fastA format.
    fastq_output:   bool,
    #[arg(short = 'C', long)]
    /// Take the reverse complement and add to data.
    complement_add: bool,
    #[arg(short = 'H', long)]
    /// Keep the fastq header as usual.
    keep_header:    bool,

    #[arg(short = 'T', long, default_value_t = 0)]
    /// Specify the read quality threshold (geometric mean, median).
    min_read_quality: u8,
    #[arg(short = 'M', long)]
    /// Interprets the threshold (-T) as the median, not the geometric mean.
    use_median:       bool,

    #[arg(short = 'L', long, default_value_t = 0)]
    /// Minimum length of sequence read data, filtered otherwise.
    min_length:             usize,
    #[arg(short = 'E', long)]
    /// The minimum length threshold (-L) is enforced when adapter clipped (-c).
    enforce_clipped_length: bool,

    #[arg(short = 'm', long)]
    /// Specify adapter sequence and mask when found in reads.
    mask_adapter:  Option<String>,
    #[arg(short = 'c', long)]
    /// Specify adapter sequence and clip appropriate ends when found in reads.
    clip_adapter:  Option<String>,
    #[arg(short = 'Z', long)]
    /// Allow up to one mismatch for adapter clipping (-c) or masking (-m).
    fuzzy_adapter: bool,

    #[arg(short = 'U', long)]
    /// Re-encode FASTQ sequence to expected input: A, C, T, G, N
    canonical_bases: bool,

    #[arg(short = 'S', long, value_hint = ValueHint::FilePath)]
    /// Save quality file for back-mapping.
    save_quality_file: Option<PathBuf>,

    #[arg(short = 'G', long, value_hint = ValueHint::FilePath)]
    /// Quality control log path and filename.
    log_file: Option<PathBuf>,
    #[arg(short = 'g', long)]
    /// Specify log ID tag (integer) for output collation.
    log_id:   Option<usize>,

    #[arg(short = 'O', long)]
    /// Replace header with monotonically increasing ordinal headers.
    ordinal_headers: bool,
    #[arg(short = 'F', long)]
    /// File id tag for ordinal header mode (-O).
    ordinal_file_id: Option<String>,

    #[arg(short = 'A', long, value_hint = ValueHint::FilePath)]
    /// Save quality vs. length statistics file for analysis.
    save_stats:     Option<PathBuf>,
    #[arg(short = 'K', long)]
    /// Do not output data FASTA/FASTQ data (assumes -A).
    skip_remaining: bool,

    /// If FASTQ header is SRA format and missing a read identifier, alter the header to add it.
    #[arg(short = 'R', long)]
    read_side: Option<char>,
}

static MODULE: &str = module_path!();

// Create a fuzzy pattern from a simple nucleotide sequence (one mismatch)
const N_PATTERN: &[u8; 7] = b"[ATCGN]";
fn create_fuzzy_pattern(pattern: &[u8]) -> Vec<u8> {
    let mut fuzzy_pattern = Vec::with_capacity(pattern.len() + N_PATTERN.len() + 1);
    let mut new_combined_pattern = Vec::new();

    for i in 0..pattern.len() {
        fuzzy_pattern.extend_from_slice(&pattern[0..i]);
        fuzzy_pattern.extend_from_slice(N_PATTERN);
        fuzzy_pattern.extend_from_slice(&pattern[(i + 1)..]);

        new_combined_pattern.extend_from_slice(&fuzzy_pattern);
        new_combined_pattern.extend_from_slice(b"|");
        fuzzy_pattern.clear();
    }
    new_combined_pattern.pop();

    new_combined_pattern
}

// Helper function to compile regex
fn compile_byte_regex(b: &[u8]) -> Result<Regex, IOError> {
    match Regex::new(String::from_utf8_lossy(b).borrow()) {
        Ok(r) => Ok(r),
        Err(e) => Err(IOError::new(ErrorKind::InvalidData, e)),
    }
}

const BAM_QNAME_LIMIT: usize = 254;

/// # Panics
///
/// Sub-program for processing fastQ data.
pub fn fastqc_process(args: &FastqConverterArgs) -> Result<(), std::io::Error> {
    let fastq_file_reader = if let Some(ref file_path) = args.fastq_input_file {
        FastQReader::new(BufReader::new(Either::Left(OpenOptions::new().read(true).open(file_path)?)))
    } else {
        FastQReader::new(BufReader::new(Either::Right(stdin())))
    };

    let (mut log_file_writer, log_name_id) = if let Some(ref file_path) = args.log_file {
        let mut log_id = String::new();
        log_id.push_str(file_path.to_string_lossy().borrow());
        if let Some(id) = args.log_id {
            log_id.push(':');
            log_id.push_str(id.to_string().as_str());
        }

        (
            Some(BufWriter::new(OpenOptions::new().append(true).create(true).open(file_path)?)),
            log_id,
        )
    } else {
        (None, String::new())
    };

    let mut save_stats_writer = if let Some(ref file_path) = args.save_stats {
        Some(BufWriter::new(
            OpenOptions::new().write(true).create(true).truncate(true).open(file_path)?,
        ))
    } else {
        None
    };

    let mut quality_file_writer = if let Some(ref file_path) = args.save_quality_file {
        Some(BufWriter::new(
            OpenOptions::new().write(true).create(true).truncate(true).open(file_path)?,
        ))
    } else {
        None
    };

    let file_id = if let Some(ref id) = args.ordinal_file_id {
        if args.complement_add {
            "|".to_owned() + id + "|"
        } else {
            id.clone() + "|"
        }
    } else {
        String::new()
    };

    let (forward_adapter, reverse_adapter) = if let Some(ref a) = args.mask_adapter {
        let forward = a.as_bytes().to_ascii_uppercase();
        let reverse = reverse_complement(&forward);

        (forward, reverse)
    } else if let Some(ref a) = args.clip_adapter {
        let forward = a.as_bytes().to_ascii_uppercase();
        let reverse = reverse_complement(&forward);

        (forward, reverse)
    } else {
        (Vec::new(), Vec::new())
    };

    let (forward_fuzzy_regex, reverse_fuzzy_regex) =
        if args.fuzzy_adapter && (args.mask_adapter.is_some() || args.clip_adapter.is_some()) {
            let fuzzy_adapter_forward = create_fuzzy_pattern(&forward_adapter);
            let fuzzy_adapter_reverse = create_fuzzy_pattern(&reverse_adapter);
            (
                Some(compile_byte_regex(&fuzzy_adapter_forward)?),
                Some(compile_byte_regex(&fuzzy_adapter_reverse)?),
            )
        } else {
            (None, None)
        };

    let (forward_regex, reverse_regex) = if args.mask_adapter.is_some() || args.clip_adapter.is_some() {
        (
            Some(compile_byte_regex(&forward_adapter)?),
            Some(compile_byte_regex(&reverse_adapter)?),
        )
    } else {
        (None, None)
    };

    let mut reads_passing_qc: u32 = 0;
    let mut ordinal_id: u32 = 0;

    let (pp, nn, dnp) = if args.complement_add { ("P", "N", "_") } else { ("", "", "") };

    // Maximum dataset median or average read quality score.
    let mut dataset_q_max = None;

    let mut give_warning_for_long_fastq = true;

    // Maximum read lengths at different points in the pipeline.
    let (mut dataset_max_read_len, mut dataset_max_clipped_read_len) = (None, None);

    for result in fastq_file_reader {
        let FastQ {
            mut header,
            mut sequence,
            mut quality,
        } = result?;

        if let Some(rside) = args.read_side {
            header = fix_sra_format(header, rside);
        }

        if sequence.len() > dataset_max_read_len.unwrap_or_default() {
            dataset_max_read_len = Some(sequence.len());
        }

        if sequence.len() < args.min_length {
            continue;
        }

        if args.canonical_bases {
            sequence.recode_iupac_to_actgn();
        }

        if let Some(ref r) = reverse_regex
            && let Some(c) = r.captures(sequence.as_bytes())
            && let Some(m) = c.get(0)
        {
            if args.clip_adapter.is_some() {
                // Chop 3' end of sequence data
                let new_length = m.start();
                sequence.as_mut_vec().truncate(new_length);
                quality.as_mut_vec().truncate(new_length);
            } else {
                // Mask data
                for i in m.start()..m.end() {
                    sequence[i] = b'N';
                }
            }
        } else if let Some(ref r) = forward_regex
            && let Some(c) = r.captures(sequence.as_bytes())
            && let Some(m) = c.get(0)
        {
            if args.clip_adapter.is_some() {
                // Remove the 5' and clone back in
                let new_start = m.end();
                sequence = sequence.as_mut_vec().drain(new_start..).collect();
                quality = quality.as_mut_vec().drain(new_start..).collect();
            } else {
                // Mask data
                for i in m.start()..m.end() {
                    sequence[i] = b'N';
                }
            }
        } else if let Some(ref r) = reverse_fuzzy_regex
            && let Some(c) = r.captures(sequence.as_bytes())
            && let Some(m) = c.get(0)
        {
            if args.clip_adapter.is_some() {
                // Chop 3' end of sequence data
                let new_length = m.start();
                sequence.as_mut_vec().truncate(new_length);
                quality.as_mut_vec().truncate(new_length);
            } else {
                // Mask data
                for i in m.start()..m.end() {
                    sequence[i] = b'N';
                }
            }
        } else if let Some(ref r) = forward_fuzzy_regex
            && let Some(c) = r.captures(sequence.as_bytes())
            && let Some(m) = c.get(0)
        {
            if args.clip_adapter.is_some() {
                // Remove the 5' and clone back in
                let new_start = m.end();
                sequence = sequence.as_mut_vec().drain(new_start..).collect();
                quality = quality.as_mut_vec().drain(new_start..).collect();
            } else {
                // Mask data
                for i in m.start()..m.end() {
                    sequence[i] = b'N';
                }
            }
        }

        let n = sequence.len();
        if args.enforce_clipped_length && n < args.min_length {
            continue;
        }

        if sequence.len() > dataset_max_clipped_read_len.unwrap_or_default() {
            dataset_max_clipped_read_len = Some(sequence.len());
        }

        ordinal_id += 1;

        let read_q_center = if args.use_median {
            quality.median()
        } else {
            quality.geometric_mean()
        };

        let read_q_center = if let Some(center) = read_q_center {
            center.as_f32()
        } else {
            continue;
        };

        if read_q_center > dataset_q_max.unwrap_or_default() {
            dataset_q_max = Some(read_q_center);
        }

        if let Some(ref mut w) = save_stats_writer {
            if let Err(e) = writeln!(w, "{ordinal_id}\t{read_q_center}\t{n}") {
                eprintln!(
                    "{MODULE} Warning! Failed to write to '{}' for id '{ordinal_id}'. See: {e}",
                    args.save_stats.as_ref().unwrap().display()
                );
            }

            if args.skip_remaining {
                continue;
            }
        }

        if read_q_center >= f32::from(args.min_read_quality) {
            reads_passing_qc += 1;
            if !args.keep_header {
                if header.len() > BAM_QNAME_LIMIT {
                    // Cannot panic given use of `floor_char_boundary`
                    header.truncate(header.floor_char_boundary(BAM_QNAME_LIMIT));
                    if give_warning_for_long_fastq {
                        give_warning_for_long_fastq = false;
                        eprintln!(
                            "{MODULE} WARNING: FASTQ headers truncated, downstream BAM format expects no more than 254 bytes!"
                        )
                    }
                }

                header = header.replace(' ', "_");
            }

            if args.ordinal_headers {
                if args.fastq_output {
                    print!("@{pp}{file_id}{ordinal_id}\n{sequence}\n+\n{quality}\n");
                } else {
                    print!(">{pp}{file_id}{ordinal_id}\n{sequence}\n");
                }
            // Regular headers
            } else if args.fastq_output {
                print!("{header}{dnp}{pp}\n{sequence}\n+\n{quality}\n");
            } else {
                let header2: String = header.chars().skip(1).collect();
                print!(">{header2}{dnp}{pp}|{read_q_center:.1}|{n}\n{sequence}\n");
                if let Some(ref mut w) = quality_file_writer {
                    writeln!(w, "{header2}{dnp}{pp}\t{quality}").unwrap_or_else(|e| {
                        eprintln!(
                            "{MODULE} Warning! Cannot write to '{}'. See: {e}",
                            args.save_quality_file.as_ref().unwrap().display()
                        );
                    });
                }
            }

            if args.complement_add {
                let reverse_sequence = sequence.reverse_complement();
                let reverse_quality = quality.reverse();

                if args.ordinal_headers {
                    if args.fastq_output {
                        print!("@{nn}{file_id}{ordinal_id}\n{reverse_sequence}\n+\n{reverse_quality}\n");
                    } else {
                        print!(">{nn}{file_id}{ordinal_id}\n{reverse_sequence}\n");
                    }
                // Regular headers
                } else if args.fastq_output {
                    print!("{header}{dnp}{nn}\n{reverse_sequence}\n+\n{reverse_quality}\n");
                } else {
                    let header2: String = header.chars().skip(1).collect();
                    print!(">{header2}{dnp}{nn}|{read_q_center:.1}|{n}\n{reverse_sequence}\n");
                    if let Some(ref mut w) = quality_file_writer {
                        writeln!(w, "{header2}{dnp}{nn}\t{quality}").unwrap_or_else(|e| {
                            eprintln!(
                                "{MODULE} Warning! Cannot write to '{}'. See: {e}",
                                args.save_quality_file.as_ref().unwrap().display()
                            );
                        });
                    }
                }
            }
        }
    }

    if let Some(ref mut w) = log_file_writer {
        let reads_passing_length = ordinal_id;
        writeln!(
            w,
            "{log_name_id}\t{reads_passing_length}\t{reads_passing_qc}\t{}\t{}\t{}",
            args.min_read_quality, args.min_length, args.use_median
        )
        .unwrap_or_else(|e| {
            eprintln!("{MODULE} Warning! Cannot write to {log_name_id}. See: {e}");
        });
    }

    if reads_passing_qc == 0 && !args.skip_remaining {
        let center_type = if args.use_median { "median" } else { "average" };
        let dataset_q_max = if let Some(v) = dataset_q_max {
            v.to_string()
        } else {
            String::from("NONE")
        };
        let dataset_max_read_len = if let Some(v) = dataset_max_read_len {
            v.to_string()
        } else {
            String::from("NONE")
        };
        let dataset_max_clipped_read_len = if let Some(v) = dataset_max_clipped_read_len {
            v.to_string()
        } else {
            String::from("NONE")
        };

        eprintln!(
            "{MODULE} Warning! No reads passed QC.

            Configuration specified minimum {center_type} read quality score: {}
            Configuration specified minimum read length: {}
            Configuration minimum clipped read length enforced: {}

            Observed dataset read length max: {dataset_max_read_len}
            Observed dataset clipped read length max: {dataset_max_clipped_read_len}
            Observed dataset {center_type} read quality max: {dataset_q_max}\n",
            args.min_read_quality, args.min_length, args.enforce_clipped_length
        );
    }

    Ok(())
}

fn fix_sra_format(header: String, read_side: char) -> String {
    let delim = if header.contains(' ') { ' ' } else { '_' };
    let mut pieces = header.split(delim);
    let maybe_id = pieces.next();

    if let Some(id) = maybe_id
        && (id.starts_with("@SRR") || id.starts_with("@DRR") || id.starts_with("@ERR"))
        && id.chars().filter(|&c| c == '.').count() == 1
    {
        let mut updated = String::with_capacity(header.len() + 2);
        updated.push_str(id);
        updated.push('.');
        updated.push(read_side);
        for p in pieces {
            updated.push(delim);
            updated.push_str(p);
        }

        updated
    } else {
        header
    }
}

#[cfg(test)]
mod test {
    use crate::fastq_converter::fix_sra_format;

    static QNAMES: [&str; 26] = [
        "@SRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "@SRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "@SRR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "@SRR26182418.1.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "@A00350:691:HCKYLDSX3:2:2119:23863:2456/2",
        "@A00350:691:HCKYLDSX3:2:2119:23863:2456/1",
        "@M02989:9:000000000-L4PJL:1:2112:9890:15606 1:N:0:AACGCACGAG+GCCTCGGATA",
        "@M02989:9:000000000-L4PJL:1:2112:9890:15606 2:N:0:AACGCACGAG+GCCTCGGATA",
        "@NS500500:69:HKJFLAFX5:1:11204:14878:14643 1:N:0:TTCTCGTGCA+CTCTGTGTAT",
        "@NS500500:69:HKJFLAFX5:1:11204:14878:14643 2:N:0:TTCTCGTGCA+CTCTGTGTAT",
        "@A01000:249:HJFFWDRX2:1:2107:24605:18082 1:N:0:TAGGCATG+ATAGCCTT",
        "@A01000:249:HJFFWDRX2:1:2107:24605:18082 2:N:0:TAGGCATG+ATAGCCTT",
        "@M02989:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
        "@M02989:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
        "@M02989_1:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
        "@M02989_1:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
        "@SRR26182418.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
        "@SRR26182418.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
        "@SRR26182418.1.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
        "@SRR26182418.1.2_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
        "@SRR26182418.1 1:N:18:NULL",
        "@SRR26182418.1.1 1:N:18:NULL",
        "@ERR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "@DRR26182418.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        "@ERR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
        "@DRR26182418.2.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
    ];

    #[test]
    fn test_fix_sra_names() {
        let fixed = [
            "@SRR26182418.1.0 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
            "@SRR26182418.1.0 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
            "@SRR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
            "@SRR26182418.1.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
            "@A00350:691:HCKYLDSX3:2:2119:23863:2456/2",
            "@A00350:691:HCKYLDSX3:2:2119:23863:2456/1",
            "@M02989:9:000000000-L4PJL:1:2112:9890:15606 1:N:0:AACGCACGAG+GCCTCGGATA",
            "@M02989:9:000000000-L4PJL:1:2112:9890:15606 2:N:0:AACGCACGAG+GCCTCGGATA",
            "@NS500500:69:HKJFLAFX5:1:11204:14878:14643 1:N:0:TTCTCGTGCA+CTCTGTGTAT",
            "@NS500500:69:HKJFLAFX5:1:11204:14878:14643 2:N:0:TTCTCGTGCA+CTCTGTGTAT",
            "@A01000:249:HJFFWDRX2:1:2107:24605:18082 1:N:0:TAGGCATG+ATAGCCTT",
            "@A01000:249:HJFFWDRX2:1:2107:24605:18082 2:N:0:TAGGCATG+ATAGCCTT",
            "@M02989:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
            "@M02989:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
            "@M02989_1:9:000000000-L4PJL:1:2114:17393:19614_1:N:0:CTCTGCAGCG+GATGGATGTA",
            "@M02989_1:9:000000000-L4PJL:1:2114:17393:19614_2:N:0:CTCTGCAGCG+GATGGATGTA",
            "@SRR26182418.1.0_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
            "@SRR26182418.1.0_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
            "@SRR26182418.1.1_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=147",
            "@SRR26182418.1.2_M07901:28:000000000-KP3NB:1:1101:10138:2117_length=301",
            "@SRR26182418.1.0 1:N:18:NULL",
            "@SRR26182418.1.1 1:N:18:NULL",
            "@ERR26182418.1.0 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
            "@DRR26182418.1.0 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
            "@ERR26182418.1.1 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=147",
            "@DRR26182418.2.2 M07901:28:000000000-KP3NB:1:1101:10138:2117 length=301",
        ];

        for (i, &o) in QNAMES.iter().enumerate() {
            assert_eq!(fix_sra_format(o.to_string(), '0'), fixed[i]);
        }
    }
}
