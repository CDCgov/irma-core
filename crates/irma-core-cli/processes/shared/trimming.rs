use crate::args::clipping::ParsedClippingArgs;
use irma_records::fastq::ReadTransforms;
use std::ops::{Add, Index, IndexMut};
use zoe::prelude::{FastQViewMut, Len};

/// Trims or masks a read based on user provided arguments. This edits the
/// underlying FASTQ data for masking and recoding.
///
/// This updates [`TrimmedCounts::pretrim_lengths`], as well as the count fields
/// if `verbose` is set.
pub fn trim_read<'a>(
    mut fq_view: FastQViewMut<'a>, mask: bool, args: &ParsedClippingArgs, counts: &mut TrimmedCounts, verbose: bool,
) -> FastQViewMut<'a> {
    fq_view.to_canonical_bases(!args.preserve_bases);

    let original_len = fq_view.sequence.len();
    counts.pretrim_lengths[original_len] += 1;

    // Poly-G trimming
    let polyg_trimmed = fq_view.process_polyg(args.polyg_left, args.polyg_right, mask);

    // Adapter or barcode trimming
    let (adapter_trimmed, barcode_trimmed) = if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
        let adapter_trimmed =
            fq_view.process_adapter(reverse_adapter.as_bytes(), forward_adapter.as_bytes(), args.a_fuzzy, mask);
        (adapter_trimmed, false)
    } else if let Some((barcode, reverse)) = &args.barcodes {
        let barcode_trimmed = fq_view.process_barcode(
            barcode.as_bytes(),
            reverse.as_bytes(),
            args.b_hdist,
            mask,
            args.b_restrict_left,
            args.b_restrict_right,
        );
        (false, barcode_trimmed)
    } else {
        (false, false)
    };

    // Primer trimming
    let primer_trimmed = if let Some(ref kmers) = args.primer_kmers {
        let mut primer_trimmed = false;
        if let Some(p_restrict_left) = args.p_restrict_left {
            primer_trimmed |= fq_view.process_left_primer(p_restrict_left, kmers, mask);
        }
        if let Some(p_restrict_right) = args.p_restrict_right {
            primer_trimmed |= fq_view.process_right_primer(p_restrict_right, kmers, mask);
        }
        primer_trimmed
    } else {
        false
    };

    // Hard trimming
    let hard_trimmed = if args.hard_left > 0 || args.hard_right > 0 {
        fq_view.hard_clip_or_mask(args.hard_left, args.hard_right, mask)
    } else {
        false
    };

    // Update counts if in verbose mode
    if verbose {
        counts.poly_g += usize::from(polyg_trimmed);
        counts.adapter += usize::from(adapter_trimmed);
        counts.barcode += usize::from(barcode_trimmed);
        counts.primer += usize::from(primer_trimmed);
        counts.hard += usize::from(hard_trimmed);

        let any_trimmed = polyg_trimmed | adapter_trimmed | barcode_trimmed | primer_trimmed | hard_trimmed;
        counts.total_trimmed += usize::from(any_trimmed);
    }

    fq_view
}

/// Counts observed read lengths, using the read length as the index.
#[derive(Debug, Default)]
pub struct LengthCounts(Vec<u32>);

impl LengthCounts {
    pub fn iter(&self) -> impl Iterator<Item = &u32> {
        self.0.iter()
    }
}

impl Index<usize> for LengthCounts {
    type Output = u32;

    fn index(&self, length: usize) -> &Self::Output {
        self.0.get(length).unwrap_or(&0)
    }
}

impl IndexMut<usize> for LengthCounts {
    fn index_mut(&mut self, length: usize) -> &mut Self::Output {
        if length >= self.0.len() {
            self.0.resize(length + 1, 0);
        }

        &mut self.0[length]
    }
}

impl Add for LengthCounts {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        let (mut longer, shorter) = if self.0.len() >= other.0.len() {
            (self, other)
        } else {
            (other, self)
        };

        for (total, count) in longer.0.iter_mut().zip(shorter.0) {
            *total += count;
        }

        longer
    }
}

#[derive(Debug, Default)]
pub struct TrimmedCounts {
    pub hard:             usize,
    pub poly_g:           usize,
    pub adapter:          usize,
    pub barcode:          usize,
    pub primer:           usize,
    pub length_filtered:  usize,
    pub widow_filtered:   usize,
    pub total_trimmed:    usize,
    pub total_processed:  usize,
    pub pretrim_lengths:  LengthCounts,
    pub posttrim_lengths: LengthCounts,
}

impl Add for TrimmedCounts {
    type Output = Self;

    fn add(self, other: Self) -> Self::Output {
        TrimmedCounts {
            hard:             self.hard + other.hard,
            poly_g:           self.poly_g + other.poly_g,
            adapter:          self.adapter + other.adapter,
            barcode:          self.barcode + other.barcode,
            primer:           self.primer + other.primer,
            length_filtered:  self.length_filtered + other.length_filtered,
            widow_filtered:   self.widow_filtered + other.widow_filtered,
            total_trimmed:    self.total_trimmed + other.total_trimmed,
            total_processed:  self.total_processed + other.total_processed,
            pretrim_lengths:  self.pretrim_lengths + other.pretrim_lengths,
            posttrim_lengths: self.posttrim_lengths + other.posttrim_lengths,
        }
    }
}

#[cfg(test)]
mod test {
    use crate::shared::trimming::LengthCounts;

    // tests for the issue that simply zipping and adding the two vectors containing
    // counts of seen lengths could possibly lose data at the end if the two vecs
    // are of different lengths
    #[test]
    fn adds_length_counts() {
        assert_eq!(
            (LengthCounts(vec![1, 2]) + LengthCounts(vec![3, 4, 5, 6])).0,
            vec![4, 6, 5, 6]
        );
        assert_eq!(
            (LengthCounts(vec![3, 4, 5, 6]) + LengthCounts(vec![1, 2])).0,
            vec![4, 6, 5, 6]
        );
        assert_eq!(
            (LengthCounts(vec![3, 4, 5, 6]) + LengthCounts(vec![3, 4, 5, 6])).0,
            vec![6, 8, 10, 12],
        );
    }
}
