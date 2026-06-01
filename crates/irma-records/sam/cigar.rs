use core::fmt::NumBuffer;
use zoe::data::{cigar::Ciglet, types::cigar::Cigar};

/// A CIGAR string as an "expanded" opcode-only byte string, e.g. `3M` → `MMM`.
#[derive(Clone, PartialEq)]
pub(crate) struct ExpandedCigar(Vec<u8>);

impl ExpandedCigar {
    /// Condenses the [`ExpandedCigar`] back to its standard form,  e.g. `MMM` →
    /// `3M`.
    #[inline]
    #[must_use]
    pub(crate) fn condense_to_cigar(self) -> Cigar {
        let mut condensed: Vec<u8> = Vec::new();
        let mut buff = NumBuffer::new();

        let mut cigars = self.0.iter().copied().filter(|op| is_valid_op(*op));

        let Some(mut previous) = cigars.next() else {
            return Cigar::from_vec_unchecked(condensed);
        };

        let mut count: usize = 1;

        for op in cigars {
            if previous == op {
                count += 1;
            } else {
                condensed.extend_from_slice(count.format_into(&mut buff).as_bytes());
                condensed.push(previous);
                previous = op;
                count = 1;
            }
        }

        condensed.extend_from_slice(count.format_into(&mut buff).as_bytes());
        condensed.push(previous);

        Cigar::from_vec_unchecked(condensed)
    }
}

impl std::fmt::Display for ExpandedCigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.0))
    }
}

impl std::fmt::Debug for ExpandedCigar {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&String::from_utf8_lossy(&self.0))
    }
}

impl From<Cigar> for ExpandedCigar {
    #[inline]
    fn from(c: Cigar) -> ExpandedCigar {
        let mut expanded = Vec::new();

        for Ciglet { inc, op } in &c {
            expanded.extend(std::iter::repeat_n(op, inc));
        }

        ExpandedCigar(expanded)
    }
}

impl From<Vec<u8>> for ExpandedCigar {
    #[inline]
    fn from(vec: Vec<u8>) -> Self {
        ExpandedCigar(vec)
    }
}

impl From<&str> for ExpandedCigar {
    #[inline]
    fn from(s: &str) -> Self {
        ExpandedCigar(s.as_bytes().to_owned())
    }
}

impl From<&[u8]> for ExpandedCigar {
    #[inline]
    fn from(v: &[u8]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}

impl<const N: usize> From<&[u8; N]> for ExpandedCigar {
    #[inline]
    fn from(v: &[u8; N]) -> Self {
        ExpandedCigar(v.to_vec())
    }
}

/// Returns whether a byte represents a valid CIGAR operation.
///
/// The following operations are valid: `MIDNSHPX=`.
//
// Copied/redundant with Zoe.
#[inline]
pub(crate) const fn is_valid_op(op: u8) -> bool {
    matches!(op, b'M' | b'I' | b'D' | b'N' | b'S' | b'H' | b'P' | b'X' | b'=')
}
