use crate::paired::PairedHeaderError;
use std::simd::Simd;
use zoe::data::records::HeaderReadable;

/// An enum representing the read side for paired or unpaired reads.
#[derive(Copy, Clone)]
pub enum ReadSide {
    R1,
    R2,
    Unpaired,
}

impl ReadSide {
    /// Convert to a char. [`ReadSide::R1`] is `Some('1')`, [`ReadSide::R2`] is
    /// `Some('2')`, and `Unpaired` is `None`.
    #[inline]
    pub fn to_char(self) -> Option<char> {
        match self {
            ReadSide::R1 => Some('1'),
            ReadSide::R2 => Some('2'),
            ReadSide::Unpaired => None,
        }
    }

    /// Convert to a SIMD vector for generating count metadata. [`ReadSide::R1`]
    /// and `Unpaired` are both `[1, 0]` and [`ReadSide::R2`] is `[0, 1]`.
    #[inline]
    pub fn to_simd(self) -> Simd<usize, 2> {
        match self {
            ReadSide::R1 | ReadSide::Unpaired => Simd::from_array([1, 0]),
            ReadSide::R2 => Simd::from_array([0, 1]),
        }
    }
}

/// Takes a FASTQ header and returns the molecular ID and side (for paired
/// reads)
pub fn get_molecular_id_side(s: &str, default_side: char) -> Option<(&str, char)> {
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

/// Returns whether two reads have matching molecular IDs. Errors if the read
/// ID's don't match or can't be parsed.
pub fn check_paired_headers<A: HeaderReadable, B: HeaderReadable>(read1: &A, read2: &B) -> Result<(), PairedHeaderError> {
    let Some((id1, _)) = get_molecular_id_side(read1.header(), '0') else {
        return Err(PairedHeaderError::ParsingErrorFirst);
    };

    let Some((id2, _)) = get_molecular_id_side(read2.header(), '1') else {
        return Err(PairedHeaderError::ParsingErrorSecond);
    };

    if id1 == id2 {
        Ok(())
    } else {
        Err(PairedHeaderError::Mismatch)
    }
}
