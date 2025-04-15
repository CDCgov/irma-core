use crate::{args::clipping::ParsedClippingArgs, qc::fastq::ReadTransforms};
use zoe::prelude::FastQViewMut;

/// Trims or masks a read based on user provided arguments. This edits the
/// underlying FASTQ data for masking and recoding.
pub fn trim_read<'a>(mut fq_view: FastQViewMut<'a>, mask: bool, args: &ParsedClippingArgs) -> FastQViewMut<'a> {
    fq_view.to_canonical_bases(!args.preserve_bases);

    fq_view.process_polyg(args.polyg_left, args.polyg_right, mask);

    if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
        fq_view.transform_by_reverse_forward_search(
            args.a_fuzzy,
            !mask,
            reverse_adapter.as_bytes(),
            forward_adapter.as_bytes(),
        );
    } else if let Some((barcode, reverse)) = &args.barcodes {
        fq_view.process_barcode(
            barcode.as_bytes(),
            reverse.as_bytes(),
            args.b_hdist,
            mask,
            args.b_restrict_left,
            args.b_restrict_right,
        );
    }

    if let Some(ref kmers) = args.primer_kmers {
        if let Some(p_restrict_left) = args.p_restrict_left {
            fq_view.process_left_primer(p_restrict_left, kmers, mask);
        }
        if let Some(p_restrict_right) = args.p_restrict_right {
            fq_view.process_right_primer(p_restrict_right, kmers, mask);
        }
    }

    if args.hard_left > 0 || args.hard_right > 0 {
        fq_view.hard_clip_or_mask(args.hard_left, args.hard_right, mask);
    }
    fq_view
}
