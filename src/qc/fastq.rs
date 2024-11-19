use std::sync::LazyLock;
use zoe::{data::fastq::FastQ, prelude::*};

static GIVE_WARNING_FOR_LONG_FASTQ: LazyLock<()> =
    LazyLock::new(|| eprintln!("WARNING: FASTQ headers truncated, downstream BAM format expects no more than 254 bytes!"));

const BAM_QNAME_LIMIT: usize = 254;

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
    #[allow(dead_code)]
    fn hard_trim(&mut self, bases: Option<usize>) -> &mut Self;
    fn fix_header(&mut self, read_side: Option<char>) -> &mut Self;
    fn clip_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn clip_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn mask_exact(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn mask_exact_or_fuzzy(&mut self, reverse: &[u8], forward: &[u8]) -> &mut Self;
    fn to_canonical_bases(&mut self, recode: bool) -> &mut Self;
    fn transform_by_reverse_forward_search(
        &mut self, is_fuzzy: bool, is_clipping: bool, reverse: &[u8], forward: &[u8],
    ) -> &mut Self;
    fn get_q_center(&self, use_median: bool) -> Option<f32>;
    fn keep_or_underscore_header(&mut self, keep_header: bool) -> &mut Self;
}

impl ReadTransforms for FastQ {
    #[inline]
    fn hard_trim(&mut self, bases: Option<usize>) -> &mut Self {
        if let Some(n) = bases
            && self.sequence.len() > 2 * n
        {
            self.sequence.cut_to_start(n);
            self.sequence.shorten_to(self.sequence.len() - n);
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
            self.sequence.mask_if_exists(r, b'N');
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
            self.sequence.mask_if_exists(r, b'N');
        }
        self
    }

    #[inline]
    fn to_canonical_bases(&mut self, recode: bool) -> &mut Self {
        if recode {
            self.sequence.recode_iupac_to_actgn();
        }
        self
    }

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
