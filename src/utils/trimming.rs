use crate::{args::clipping::ParsedClippingArgs, qc::fastq::ReadTransforms};
use std::ops::Add;
use zoe::prelude::{FastQViewMut, Len};

/// Trims or masks a read based on user provided arguments. This edits the
/// underlying FASTQ data for masking and recoding.
pub fn trim_read<'a>(
    mut fq_view: FastQViewMut<'a>, mask: bool, args: &ParsedClippingArgs, counts: &mut TrimmedCounts, verbose: bool,
) -> FastQViewMut<'a> {
    fq_view.to_canonical_bases(!args.preserve_bases);

    counts.last_read_len = fq_view.sequence.len();
    let mut original_len = fq_view.sequence.len();

    fq_view.process_polyg(args.polyg_left, args.polyg_right, mask);
    update_trimmed_counts_field(&mut counts.poly_g, &fq_view, &mut counts.last_read_len, verbose);

    if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
        fq_view.transform_by_reverse_forward_search(
            args.a_fuzzy,
            !mask,
            reverse_adapter.as_bytes(),
            forward_adapter.as_bytes(),
        );
        update_trimmed_counts_field(&mut counts.adapter, &fq_view, &mut counts.last_read_len, verbose);
    } else if let Some((barcode, reverse)) = &args.barcodes {
        fq_view.process_barcode(
            barcode.as_bytes(),
            reverse.as_bytes(),
            args.b_hdist,
            mask,
            args.b_restrict_left,
            args.b_restrict_right,
        );
        update_trimmed_counts_field(&mut counts.barcode, &fq_view, &mut counts.last_read_len, verbose);
    }

    if let Some(ref kmers) = args.primer_kmers {
        if let Some(p_restrict_left) = args.p_restrict_left {
            fq_view.process_left_primer(p_restrict_left, kmers, mask);
        }
        if let Some(p_restrict_right) = args.p_restrict_right {
            fq_view.process_right_primer(p_restrict_right, kmers, mask);
        }
        update_trimmed_counts_field(&mut counts.primer, &fq_view, &mut counts.last_read_len, verbose);
    }

    if args.hard_left > 0 || args.hard_right > 0 {
        fq_view.hard_clip_or_mask(args.hard_left, args.hard_right, mask);
        update_trimmed_counts_field(&mut counts.hard, &fq_view, &mut counts.last_read_len, verbose);
    }
    update_trimmed_counts_field(&mut counts.total_trimmed, &fq_view, &mut original_len, verbose);
    fq_view
}

#[derive(Default, Debug)]
pub struct TrimmedCounts {
    pub last_read_len:   usize,
    pub hard:            usize,
    pub poly_g:          usize,
    pub adapter:         usize,
    pub barcode:         usize,
    pub primer:          usize,
    pub length_filtered: usize,
    pub widow_filtered:  usize,
    pub total_trimmed:   usize,
    pub total_processed: usize,
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
            length_filtered: self.length_filtered + other.length_filtered,
            widow_filtered:  self.widow_filtered + other.widow_filtered,
            total_trimmed:   self.total_trimmed + other.total_trimmed,
            total_processed: self.total_processed + other.total_processed,
        }
    }
}

pub fn update_trimmed_counts_field(field: &mut usize, read: &FastQViewMut<'_>, last_read_len: &mut usize, verbose: bool) {
    if verbose {
        if read.len() < *last_read_len {
            *field += 1;
        }
        *last_read_len = read.len();
    }
}
