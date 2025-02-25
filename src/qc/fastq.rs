use foldhash::fast::RandomState;
use std::sync::LazyLock;
use zoe::{
    data::fastq::FastQ,
    kmer::ThreeBitKmerSet,
    prelude::*,
    search::{RangeSearch, ToRangeSearch},
};

static GIVE_WARNING_FOR_LONG_FASTQ: LazyLock<()> =
    LazyLock::new(|| eprintln!("WARNING: FASTQ headers truncated, downstream BAM format expects no more than 254 bytes!"));

const BAM_QNAME_LIMIT: usize = 254;
const MAX_KMER_LENGTH: usize = 21;

pub(crate) fn fix_sra_format(header: String, read_side: char) -> String {
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

pub(crate) trait ReadTransforms {
    fn hard_clip(&mut self, left_bases: usize, right_bases: usize) -> &mut Self;
    fn hard_mask(&mut self, left_bases: usize, right_bases: usize) -> &mut Self;

    #[inline]
    fn hard_clip_or_mask(&mut self, left_bases: usize, right_bases: usize, masking: bool) -> &mut Self {
        if masking {
            self.hard_mask(left_bases, right_bases)
        } else {
            self.hard_clip(left_bases, right_bases)
        }
    }

    fn process_left_primer(
        &mut self, restrict_left: usize, kmer_set: &ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, masking: bool,
    ) -> &mut Self;
    fn process_right_primer(
        &mut self, restrict_right: usize, kmer_set: &ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, masking: bool,
    ) -> &mut Self;
    fn process_barcode(
        &mut self, barcode: &[u8], reverse: &[u8], hdist: usize, masking: bool, b_restrict_left: Option<usize>,
        b_restrict_right: Option<usize>,
    ) -> &mut Self;
    fn process_polyg(&mut self, polyg_literals: &(Option<Vec<u8>>, Option<Vec<u8>>), masking: bool) -> &mut Self;
    fn process_left_polyg(&mut self, left_literal: &[u8], masking: bool) -> &mut Self;
    fn process_right_polyg(&mut self, right_literal: &[u8], masking: bool) -> &mut Self;
    fn fix_header(&mut self, read_side: Option<char>) -> &mut Self;
    fn clip_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn clip_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn mask_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn mask_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn to_canonical_bases(&mut self, recode: bool) -> &mut Self;

    #[inline]
    fn transform_by_reverse_forward_search(
        &mut self, is_fuzzy: bool, is_clipping: bool, reverse: &[u8], forward: &[u8],
    ) -> &mut Self {
        match (is_fuzzy, is_clipping) {
            (true, true) => self.clip_exact_or_fuzzy(reverse, forward),
            (true, false) => self.mask_exact_or_fuzzy(reverse, forward),
            (false, true) => self.clip_exact(reverse, forward),
            (false, false) => self.mask_exact(reverse, forward),
        }
    }

    fn get_q_center(&self, use_median: bool) -> Option<f32>;
    fn keep_or_underscore_header(&mut self, keep_header: bool) -> &mut Self;
}

impl ReadTransforms for FastQ {
    #[inline]
    fn hard_clip(&mut self, left_bases: usize, right_bases: usize) -> &mut Self {
        if self.sequence.len() > left_bases + right_bases {
            self.sequence.cut_to_start(left_bases);
            self.sequence.shorten_to(self.sequence.len() - right_bases);
            self.quality.cut_to_start(left_bases);
            self.quality.shorten_to(self.quality.len() - right_bases);
        } else {
            self.sequence.shorten_to(0);
            self.quality.shorten_to(0);
        }
        self
    }

    #[inline]
    fn hard_mask(&mut self, left_bases: usize, right_bases: usize) -> &mut Self {
        if self.sequence.len() > left_bases + right_bases {
            self.sequence.mask_if_exists(..left_bases);
            self.sequence.mask_if_exists(self.sequence.len() - right_bases..);
        } else {
            self.sequence.mask_if_exists(..);
        }
        self
    }

    #[inline]
    fn process_left_primer(
        &mut self, restrict_left: usize, kmer_set: &ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, masking: bool,
    ) -> &mut Self {
        if !masking {
            if let Some(range) = self.sequence.search_in_first(restrict_left).find_kmers_rev(kmer_set) {
                self.sequence.cut_to_start(range.end);
                self.quality.cut_to_start(range.end);
            }
        } else {
            let mut ranges = self.sequence.search_in_first(restrict_left).find_all_kmers_rev(kmer_set);
            if let Some(mut masking_range) = ranges.next() {
                for other_range in ranges {
                    if other_range.start <= masking_range.end {
                        masking_range.end = other_range.end;
                    } else {
                        break;
                    }
                }
                self.mask_if_exists(masking_range);
            }
        }
        self
    }

    #[inline]
    fn process_right_primer(
        &mut self, restrict_right: usize, kmer_set: &ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, masking: bool,
    ) -> &mut Self {
        if !masking {
            if let Some(range) = self.sequence.search_in_last(restrict_right).find_kmers(kmer_set) {
                self.sequence.shorten_to(range.start);
                self.quality.shorten_to(range.start);
            }
        } else {
            let mut ranges = self.sequence.search_in_last(restrict_right).find_all_kmers(kmer_set);
            if let Some(mut masking_range) = ranges.next() {
                for other_range in ranges {
                    if other_range.start <= masking_range.end {
                        masking_range.end = other_range.end;
                    } else {
                        break;
                    }
                }
                self.mask_if_exists(masking_range);
            }
        }
        self
    }

    #[inline]
    fn process_barcode(
        &mut self, barcode: &[u8], reverse: &[u8], hdist: usize, masking: bool, b_restrict_left: Option<usize>,
        b_restrict_right: Option<usize>,
    ) -> &mut Self {
        let restricted_substring_fn = match hdist {
            0 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_substring(needle),
            1 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_fuzzy_substring::<1>(needle),
            2 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_fuzzy_substring::<2>(needle),
            3 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_fuzzy_substring::<3>(needle),
            _ => unreachable!(),
        };

        let substring_fn = match hdist {
            0 => |needle: &[u8], seq: &Nucleotides| seq.find_substring(needle),
            1 => |needle: &[u8], seq: &Nucleotides| seq.find_fuzzy_substring::<1>(needle),
            2 => |needle: &[u8], seq: &Nucleotides| seq.find_fuzzy_substring::<2>(needle),
            3 => |needle: &[u8], seq: &Nucleotides| seq.find_fuzzy_substring::<3>(needle),
            _ => unreachable!(),
        };

        let left_barcode_pos = match b_restrict_left {
            Some(b_restrict_left) => restricted_substring_fn(barcode, &self.sequence.search_in_first(b_restrict_left)),
            None => substring_fn(barcode, &self.sequence),
        };

        if let Some(left_range) = left_barcode_pos {
            if masking {
                self.sequence.mask_if_exists(left_range);
            } else {
                self.sequence.cut_to_start(left_range.end);
                self.quality.cut_to_start(left_range.end);
            }
        }

        let right_barcode_pos = match b_restrict_right {
            Some(b_restrict_right) => restricted_substring_fn(reverse, &self.sequence.search_in_last(b_restrict_right)),
            None => substring_fn(reverse, &self.sequence),
        };

        if let Some(right_range) = right_barcode_pos {
            if masking {
                self.sequence.mask_if_exists(right_range);
            } else {
                self.sequence.shorten_to(right_range.start);
                self.quality.shorten_to(right_range.start);
            }
        }
        self
    }

    fn process_polyg(&mut self, polyg_literals: &(Option<Vec<u8>>, Option<Vec<u8>>), masking: bool) -> &mut Self {
        if let Some(left_literal) = &polyg_literals.0 {
            self.process_left_polyg(left_literal, masking);
        }
        if let Some(right_literal) = &polyg_literals.1 {
            self.process_right_polyg(right_literal, masking);
        }

        self
    }

    #[inline]
    fn process_left_polyg(&mut self, left_literal: &[u8], masking: bool) -> &mut Self {
        if self.sequence.as_bytes().starts_with(left_literal) {
            let polyg_end = self
                .sequence
                .search_in(left_literal.len()..)
                .position(|base| base != b'G')
                .unwrap_or(self.sequence.len());
            if masking {
                self.sequence.mask_if_exists(..polyg_end);
            } else {
                self.sequence.cut_to_start(polyg_end);
                self.quality.cut_to_start(polyg_end);
            }
        }
        self
    }

    #[inline]
    fn process_right_polyg(&mut self, right_literal: &[u8], masking: bool) -> &mut Self {
        if self.sequence.as_bytes().ends_with(right_literal) {
            let polyg_start = self
                .sequence
                .search_in(..self.sequence.len() - right_literal.len())
                .rposition(|base| base != b'G')
                .map_or(0, |pos| pos + 1);
            if masking {
                self.sequence.mask_if_exists(polyg_start..);
            } else {
                self.sequence.shorten_to(polyg_start);
                self.quality.shorten_to(polyg_start);
            }
        }

        self
    }

    #[inline]
    fn fix_header(&mut self, read_side: Option<char>) -> &mut Self {
        if let Some(read_side) = read_side {
            self.header = fix_sra_format(std::mem::take(&mut self.header), read_side);
        }
        self
    }

    #[inline]
    fn clip_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self.sequence.find_substring(reverse) {
            // Chop 3' end of sequence data
            self.sequence.shorten_to(r.start);
            self.quality.shorten_to(r.start);
        } else if let Some(r) = self.sequence.find_substring(forward) {
            // Remove the 5' and clone back in
            self.sequence.cut_to_start(r.end);
            self.quality.cut_to_start(r.end);
        }
        self
    }

    #[inline]
    fn clip_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self.sequence.find_substring(reverse) {
            // Chop 3' end of sequence data
            self.sequence.shorten_to(r.start);
            self.quality.shorten_to(r.start);
        } else if let Some(r) = self.sequence.find_substring(forward) {
            // Remove the 5' and clone back in
            self.sequence.cut_to_start(r.end);
            self.quality.cut_to_start(r.end);
        } else if let Some(r) = self.sequence.find_fuzzy_substring::<1>(reverse) {
            // Chop 3' end of sequence data
            self.sequence.shorten_to(r.start);
            self.quality.shorten_to(r.start);
        } else if let Some(r) = self.sequence.find_fuzzy_substring::<1>(forward) {
            // Remove the 5' and clone back in
            self.sequence.cut_to_start(r.end);
            self.quality.cut_to_start(r.end);
        }
        self
    }

    #[inline]
    fn mask_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self
            .sequence
            .find_substring(reverse)
            .or_else(|| self.sequence.find_substring(forward))
        {
            self.sequence.mask_if_exists(r);
        }
        self
    }

    #[inline]
    fn mask_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self
            .sequence
            .find_substring(reverse)
            .or_else(|| self.sequence.find_substring(forward))
            .or_else(|| self.sequence.find_fuzzy_substring::<1>(reverse))
            .or_else(|| self.sequence.find_fuzzy_substring::<1>(forward))
        {
            self.sequence.mask_if_exists(r);
        }
        self
    }

    #[inline]
    fn to_canonical_bases(&mut self, recode: bool) -> &mut Self {
        if recode {
            self.sequence.recode_any_to_acgtn_uc();
        }
        self
    }

    #[inline]
    fn get_q_center(&self, use_median: bool) -> Option<f32> {
        if use_median {
            self.quality.median()
        } else {
            self.quality.geometric_mean()
        }
        .map(|q| q.as_f32())
    }

    #[inline]
    fn keep_or_underscore_header(&mut self, keep_header: bool) -> &mut Self {
        if !keep_header {
            if self.header.len() > BAM_QNAME_LIMIT {
                // Cannot panic given use of `floor_char_boundary`
                self.header.truncate(self.header.floor_char_boundary(BAM_QNAME_LIMIT));
                // Print warning once
                *GIVE_WARNING_FOR_LONG_FASTQ;
            }

            self.header = self.header.replace(' ', "_");
        }
        self
    }
}

impl ReadTransforms for FastQViewMut<'_> {
    #[inline]
    fn hard_clip(&mut self, left_bases: usize, right_bases: usize) -> &mut Self {
        if self.sequence.len() > left_bases + right_bases {
            self.restrict(left_bases..self.len() - right_bases);
        } else {
            self.restrict(..0);
        }
        self
    }

    #[inline]
    fn hard_mask(&mut self, left_bases: usize, right_bases: usize) -> &mut Self {
        if self.sequence.len() > left_bases + right_bases {
            self.mask_if_exists(..left_bases);
            self.mask_if_exists(self.sequence.len() - right_bases..);
            self.restrict(left_bases..self.len() - right_bases);
        } else {
            self.mask_if_exists(..);
            self.restrict(..0);
        }
        self
    }

    #[inline]
    fn process_left_primer(
        &mut self, restrict_left: usize, kmer_set: &ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, masking: bool,
    ) -> &mut Self {
        if !masking {
            if let Some(range) = self.sequence.search_in_first(restrict_left).find_kmers_rev(kmer_set) {
                self.restrict(range.end..);
            }
        } else {
            let mut ranges = self.sequence.search_in_first(restrict_left).find_all_kmers_rev(kmer_set);
            if let Some(mut masking_range) = ranges.next() {
                for other_range in ranges {
                    if other_range.start <= masking_range.end {
                        masking_range.end = other_range.end;
                    } else {
                        break;
                    }
                }
                self.mask_if_exists(masking_range.clone());
                self.restrict(masking_range.end..);
            }
        }
        self
    }

    #[inline]
    fn process_right_primer(
        &mut self, restrict_right: usize, kmer_set: &ThreeBitKmerSet<MAX_KMER_LENGTH, RandomState>, masking: bool,
    ) -> &mut Self {
        if !masking {
            if let Some(range) = self.sequence.search_in_last(restrict_right).find_kmers(kmer_set) {
                self.restrict(..range.start);
            }
        } else {
            let mut ranges = self.sequence.search_in_last(restrict_right).find_all_kmers(kmer_set);
            if let Some(mut masking_range) = ranges.next() {
                for other_range in ranges {
                    if other_range.start <= masking_range.end {
                        masking_range.end = other_range.end;
                    } else {
                        break;
                    }
                }
                self.mask_if_exists(masking_range.clone());
                self.restrict(..masking_range.start);
            }
        }
        self
    }

    #[inline]
    fn process_barcode(
        &mut self, barcode: &[u8], reverse: &[u8], hdist: usize, masking: bool, b_restrict_left: Option<usize>,
        b_restrict_right: Option<usize>,
    ) -> &mut Self {
        let substring_fn = match hdist {
            0 => |needle: &[u8], seq: &NucleotidesViewMut<'_>| seq.find_substring(needle),
            1 => |needle: &[u8], seq: &NucleotidesViewMut<'_>| seq.find_fuzzy_substring::<1>(needle),
            2 => |needle: &[u8], seq: &NucleotidesViewMut<'_>| seq.find_fuzzy_substring::<2>(needle),
            3 => |needle: &[u8], seq: &NucleotidesViewMut<'_>| seq.find_fuzzy_substring::<3>(needle),
            _ => unreachable!(),
        };

        let restricted_substring_fn = match hdist {
            0 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_substring(needle),
            1 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_fuzzy_substring::<1>(needle),
            2 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_fuzzy_substring::<2>(needle),
            3 => |needle: &[u8], seq: &RangeSearch<'_>| seq.find_fuzzy_substring::<3>(needle),
            _ => unreachable!(),
        };

        let left_barcode_pos = match b_restrict_left {
            Some(b_restrict_left) => restricted_substring_fn(barcode, &self.sequence.search_in_first(b_restrict_left)),
            None => substring_fn(barcode, &self.sequence),
        };

        if let Some(left_range) = left_barcode_pos {
            if masking {
                self.sequence.mask_if_exists(left_range.clone());
            }
            self.restrict(left_range.end..);
        }

        let right_barcode_pos = match b_restrict_right {
            Some(b_restrict_right) => restricted_substring_fn(reverse, &self.sequence.search_in_last(b_restrict_right)),
            None => substring_fn(reverse, &self.sequence),
        };

        if let Some(right_range) = right_barcode_pos {
            if masking {
                self.sequence.mask_if_exists(right_range.clone());
            }
            self.restrict(..right_range.start);
        }
        self
    }

    fn process_polyg(&mut self, polyg_literals: &(Option<Vec<u8>>, Option<Vec<u8>>), masking: bool) -> &mut Self {
        if let Some(left_literal) = &polyg_literals.0 {
            self.process_left_polyg(left_literal, masking);
        }
        if let Some(right_literal) = &polyg_literals.1 {
            self.process_right_polyg(right_literal, masking);
        }
        self
    }

    fn process_left_polyg(&mut self, left_literal: &[u8], masking: bool) -> &mut Self {
        if self.sequence.as_bytes().starts_with(left_literal) {
            let polyg_end = self
                .sequence
                .search_in(left_literal.len()..)
                .position(|base| base != b'G')
                .unwrap_or(self.sequence.len());
            if masking {
                self.sequence.mask_if_exists(..polyg_end)
            }
            self.sequence.restrict(polyg_end..);
        }
        self
    }

    fn process_right_polyg(&mut self, right_literal: &[u8], masking: bool) -> &mut Self {
        if self.sequence.as_bytes().ends_with(right_literal) {
            let polyg_start = self
                .sequence
                .search_in(..self.sequence.len() - right_literal.len())
                .rposition(|base| base != b'G')
                .map_or(0, |pos| pos + 1);
            if masking {
                self.sequence.mask_if_exists(polyg_start..);
            }
            self.sequence.restrict(..polyg_start);
        }
        self
    }

    #[inline]
    fn fix_header(&mut self, read_side: Option<char>) -> &mut Self {
        if let Some(read_side) = read_side {
            *self.header = fix_sra_format(std::mem::take(self.header), read_side);
        }
        self
    }

    #[inline]
    fn clip_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self.sequence.find_substring(reverse) {
            // Chop 3' end of sequence data
            self.restrict(..r.start);
        } else if let Some(r) = self.sequence.find_substring(forward) {
            // Remove the 5' end
            self.restrict(r.end..);
        }
        self
    }

    #[inline]
    fn clip_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self.sequence.find_substring(reverse) {
            // Chop 3' end of sequence data
            self.restrict(..r.start);
        } else if let Some(r) = self.sequence.find_substring(forward) {
            // Remove the 5' end
            self.restrict(r.end..);
        } else if let Some(r) = self.sequence.find_fuzzy_substring::<1>(reverse) {
            // Chop 3' end of sequence data
            self.restrict(..r.start);
        } else if let Some(r) = self.sequence.find_fuzzy_substring::<1>(forward) {
            // Remove the 5' end
            self.restrict(r.end..);
        }
        self
    }

    #[inline]
    fn mask_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self
            .sequence
            .find_substring(reverse)
            .or_else(|| self.sequence.find_substring(forward))
        {
            self.sequence.mask_if_exists(r);
        }
        self
    }

    #[inline]
    fn mask_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self {
        if let Some(r) = self
            .sequence
            .find_substring(reverse)
            .or_else(|| self.sequence.find_substring(forward))
            .or_else(|| self.sequence.find_fuzzy_substring::<1>(reverse))
            .or_else(|| self.sequence.find_fuzzy_substring::<1>(forward))
        {
            self.sequence.mask_if_exists(r);
        }
        self
    }

    #[inline]
    fn to_canonical_bases(&mut self, recode: bool) -> &mut Self {
        if recode {
            self.sequence.recode_any_to_acgtn_uc();
        }
        self
    }

    #[inline]
    fn get_q_center(&self, use_median: bool) -> Option<f32> {
        if use_median {
            self.quality.median()
        } else {
            self.quality.geometric_mean()
        }
        .map(|q| q.as_f32())
    }

    #[inline]
    fn keep_or_underscore_header(&mut self, keep_header: bool) -> &mut Self {
        if !keep_header {
            if self.header.len() > BAM_QNAME_LIMIT {
                // Cannot panic given use of `floor_char_boundary`
                self.header.truncate(self.header.floor_char_boundary(BAM_QNAME_LIMIT));
                // Print warning once
                *GIVE_WARNING_FOR_LONG_FASTQ;
            }

            *self.header = self.header.replace(' ', "_");
        }
        self
    }
}
