use crate::args::clipping::ParsedClippingArgs;
use irma_records::fastq::ReadTransforms;
use std::ops::{Add, Index, IndexMut};
use zoe::prelude::{FastQViewMut, Len};

/// Trims or masks a read based on user provided arguments. This edits the
/// underlying FASTQ data for masking and recoding.
pub fn trim_read<'a>(
    mut fq_view: FastQViewMut<'a>, mask: bool, args: &ParsedClippingArgs, counts: &mut TrimmedCounts, verbose: bool,
) -> FastQViewMut<'a> {
    fq_view.to_canonical_bases(!args.preserve_bases);

    counts.last_read_len = fq_view.sequence.len();

    let mut original_len = fq_view.sequence.len();
    counts.pretrim_lengths[original_len] += 1;

    fq_view.process_polyg(args.polyg_left, args.polyg_right, mask);
    update_trimmed_counts_field(&mut counts.poly_g, &fq_view, &mut counts.last_read_len, verbose);

    if let Some((ref forward_adapter, ref reverse_adapter)) = args.adapters {
        fq_view.process_adapter(reverse_adapter.as_bytes(), forward_adapter.as_bytes(), args.a_fuzzy, mask);
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

/// Counts observed read lengths, using the read length as the index.
#[derive(Default, Debug)]
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

#[derive(Default, Debug)]
pub struct TrimmedCounts {
    pub last_read_len:    usize,
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
            last_read_len:    self.last_read_len,
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

pub fn update_trimmed_counts_field(field: &mut usize, read: &FastQViewMut<'_>, last_read_len: &mut usize, verbose: bool) {
    if verbose {
        if read.len() < *last_read_len {
            *field += 1;
        }
        *last_read_len = read.len();
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
