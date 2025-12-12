use std::ops::Add;

use crate::{args::clipping::ParsedClippingArgs, qc::fastq::ReadTransforms};
use zoe::prelude::{FastQViewMut, Len};

/// Trims or masks a read based on user provided arguments. This edits the
/// underlying FASTQ data for masking and recoding.
pub fn trim_read<'a>(
    mut fq_view: FastQViewMut<'a>, mask: bool, args: &ParsedClippingArgs, counts: &mut TrimmedCounts,
) -> FastQViewMut<'a> {
    fq_view.to_canonical_bases(!args.preserve_bases);

    counts.last_read_len = fq_view.sequence.len();
    let mut original_len = fq_view.sequence.len();

    fq_view.process_polyg(args.polyg_left, args.polyg_right, mask);
    update_trimmed_counts_field(&mut counts.poly_g, &fq_view, &mut counts.last_read_len, args);

    if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
        fq_view.transform_by_reverse_forward_search(
            args.a_fuzzy,
            !mask,
            reverse_adapter.as_bytes(),
            forward_adapter.as_bytes(),
        );
        update_trimmed_counts_field(&mut counts.adapter, &fq_view, &mut counts.last_read_len, args);
    } else if let Some((barcode, reverse)) = &args.barcodes {
        fq_view.process_barcode(
            barcode.as_bytes(),
            reverse.as_bytes(),
            args.b_hdist,
            mask,
            args.b_restrict_left,
            args.b_restrict_right,
        );
        update_trimmed_counts_field(&mut counts.barcode, &fq_view, &mut counts.last_read_len, args);
    }

    if let Some(ref kmers) = args.primer_kmers {
        if let Some(p_restrict_left) = args.p_restrict_left {
            fq_view.process_left_primer(p_restrict_left, kmers, mask);
        }
        if let Some(p_restrict_right) = args.p_restrict_right {
            fq_view.process_right_primer(p_restrict_right, kmers, mask);
        }
        update_trimmed_counts_field(&mut counts.primer, &fq_view, &mut counts.last_read_len, args);
    }

    if args.hard_left > 0 || args.hard_right > 0 {
        fq_view.hard_clip_or_mask(args.hard_left, args.hard_right, mask);
        update_trimmed_counts_field(&mut counts.hard, &fq_view, &mut counts.last_read_len, args);
    }
    update_trimmed_counts_field(&mut counts.total_trimmed, &fq_view, &mut original_len, args);
    fq_view
}

pub struct TrimmedCounts {
    pub last_read_len:   usize,
    pub hard:            usize,
    pub poly_g:          usize,
    pub adapter:         usize,
    pub barcode:         usize,
    pub primer:          usize,
    pub filtered:        usize,
    pub total_trimmed:   usize,
    pub total_processed: usize,
}

impl Default for TrimmedCounts {
    fn default() -> Self {
        Self {
            last_read_len:   0,
            hard:            0,
            poly_g:          0,
            adapter:         0,
            barcode:         0,
            primer:          0,
            filtered:        0,
            total_trimmed:   0,
            total_processed: 0,
        }
    }
}

pub fn get_trimming_type(paired_input: bool, single_out: bool, filter_widows: bool) -> Result<&'static str, std::io::Error> {
    if paired_input {
        match (single_out, filter_widows) {
            // Case 2: In 1, In 2, Out 1 (interleaved Illumina), no widow filtering
            (true, false) => {
                return Ok("Paired, Interleaved output, no widow filtering");
            }

            // Case 3: In 1, In 2, Out 1, Filtering widows / orphan reads
            (true, true) => {
                return Ok("Paired, interleaved output, widow filtering");
            }

            // Case 4: In 1, In 2, Out 1, Out 2 (separated output Illumina), no filtering
            (false, false) => {
                return Ok("Paired input and output, no filtering");
            }

            // Case 5: In 1, In 2, Out 1, Out 2, filter widows
            (false, true) => {
                return Ok("paired input and output, filtering");
            }
        }
    } else {
        match (single_out, filter_widows) {
            (true, false) => {
                return Ok("single input and output");
            }

            (true, true) => {
                return Ok("interleaved input and output, filtering");
            }

            (false, false) => {
                return Ok("interleaved input, two output, no filtering");
            }

            (false, true) => {
                return Ok("interleaved input, two output, no filtering");
            }
        }
    }
}

impl TrimmedCounts {
    pub fn write_counts(self, args: &ParsedClippingArgs, trim_type: &str, masking: bool) {
        let ParsedClippingArgs {
            preserve_bases: _,
            barcodes,
            b_restrict_left: _,
            b_restrict_right: _,
            b_hdist: _,
            adapters,
            a_fuzzy: _,
            primer_kmers,
            p_restrict_left: _,
            p_restrict_right: _,
            polyg_left,
            polyg_right,
            hard_left,
            hard_right,
            verbose: _,
        } = args;

        let trim_mask = match masking {
            true => "masked",
            false => "trimmed",
        };

        eprintln!("Trimmer processed {} total reads with {trim_type}", self.total_processed);

        if polyg_left.is_some() || polyg_right.is_some() {
            eprintln!("{} reads PolyG {trim_mask}", self.poly_g);
        }
        if barcodes.is_some() {
            eprintln!("{} reads barcode {trim_mask}", self.barcode);
        }
        if adapters.is_some() {
            eprintln!("{} reads adapter {trim_mask}", self.adapter);
        }
        if primer_kmers.is_some() {
            eprintln!("{} reads primer {trim_mask}", self.primer);
        }
        if *hard_left > 0usize || *hard_right > 0usize {
            eprintln!("{} reads hard {trim_mask}", self.hard);
        }
        eprintln!("{} total reads {trim_mask}", self.total_trimmed);
    }
}

impl Add for TrimmedCounts {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        TrimmedCounts {
            last_read_len:   self.last_read_len,
            hard:            self.hard + other.hard,
            poly_g:          self.poly_g + other.poly_g,
            adapter:         self.adapter + other.adapter,
            barcode:         self.barcode + other.barcode,
            primer:          self.primer + other.primer,
            filtered:        self.filtered + other.filtered,
            total_trimmed:   self.total_trimmed + other.total_trimmed,
            total_processed: self.total_processed + other.total_processed,
        }
    }
}

pub fn update_trimmed_counts_field(
    field: &mut usize, read: &FastQViewMut<'_>, last_read_len: &mut usize, args: &ParsedClippingArgs,
) {
    if args.verbose {
        if read.len() < *last_read_len {
            *field += 1;
        }
        *last_read_len = read.len();
    }
}
