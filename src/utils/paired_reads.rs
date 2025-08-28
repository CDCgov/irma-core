use std::{
    fmt::Debug,
    io::{Error as IOError, ErrorKind},
    simd::Simd,
};
use zoe::{data::fasta::FastaSeq, prelude::FastQ, unwrap_or_return_some_err};

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
pub fn check_paired_headers<A: RecordWithHeader, B: RecordWithHeader>(read1: &A, read2: &B) -> Result<(), std::io::Error> {
    // Strip the @ from the header before calling get_molecular_id_side as needed (transitional)
    if let Some((id1, _)) = get_molecular_id_side(read1.header().strip_prefix('@').unwrap_or(read1.header()), '0')
        && let Some((id2, _)) = get_molecular_id_side(read2.header().strip_prefix('@').unwrap_or(read2.header()), '0')
    {
        if id1 == id2 {
            Ok(())
        } else {
            Err(IOError::new(
                ErrorKind::InvalidInput,
                format!(
                    "Paired read IDs out of sync:\n\t{h1}\n\t{h2}\n",
                    h1 = read1.header(),
                    h2 = read2.header()
                ),
            ))
        }
    } else {
        Err(IOError::new(ErrorKind::InvalidInput, "Could not parse the read IDs."))
    }
}

/// An enum representing the read side for paired or unpaired reads.
#[derive(Copy, Clone)]
pub enum ReadSide {
    R1,
    R2,
    Unpaired,
}

impl ReadSide {
    #[inline]
    pub fn to_idx(self) -> usize {
        match self {
            ReadSide::R1 => 0,
            ReadSide::R2 => 1,
            ReadSide::Unpaired => 0,
        }
    }

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

/// A trait for unifying the different ways to process paired reads and
/// potentially filter them. Currently used in `trimmer` and `preprocess`.
///
/// A single read is first processed and potentially filtered using
/// [`process_read`]. If filtering does not occur, then [`output_read`] is used
/// to perform any final processing and then output/store the result. Together,
/// these two actions can be performed with [`handle_read`].
///
/// [`process_read`]: PairedReadFilterer::process_read
/// [`output_read`]: PairedReadFilterer::output_read
/// [`handle_read`]: PairedReadFilterer::handle_read
pub trait PairedReadFilterer: Debug {
    /// Controls whether the order of processing must be preserved. When `true`,
    /// reads will be processed in alternating order between R1 and R2. When
    /// `false`, the order may differ.
    const PRESERVE_ORDER: bool;

    /// The type of the trimmed output.
    type Processed<'a>;

    /// The type returned by [`finalize`]. For example, this may be
    /// `std::io::Result` if writers are being flushed, or it could be any
    /// metadata tallied.
    ///
    /// [`finalize`]: Self::finalize
    type Finalized;

    /// Processes the read (e.g., trimming), performs any filtering (in which
    /// case `None` is returned) and tallies any relevant metadata (by mutating
    /// `self`).
    fn process_read<'a>(&mut self, read: &'a mut FastQ, side: ReadSide) -> Option<Self::Processed<'a>>;

    /// For any read that passed filtering with [`process_read`], perform any
    /// final processing and output/store it (either mutating `self` or using IO
    /// operations)
    ///
    /// [`process_read`]: PairedReadFilterer::process_read
    fn output_read<'a>(&mut self, trimmed: Self::Processed<'a>, side: ReadSide) -> std::io::Result<()>;

    /// Finalize the [`PairedReadFilterer`] (e.g., flushing buffers or returning
    /// saved state).
    fn finalize(&mut self) -> Self::Finalized;

    /// Processes a read using [`process_read`] followed by [`output_read`].
    ///
    /// [`process_read`]: Self::process_read
    /// [`output_read`]: Self::output_read
    fn handle_read(&mut self, mut read: FastQ, side: ReadSide) -> std::io::Result<&mut Self> {
        let Some(trimmed) = self.process_read(&mut read, side) else {
            return Ok(self);
        };
        self.output_read(trimmed, side)?;
        Ok(self)
    }

    /// Handles both reads in a read pair using [`process_read`], and then
    /// outputs/stores them with [`output_read`] only if neither is filtered.
    ///
    /// [`process_read`]: Self::process_read
    /// [`output_read`]: Self::output_read
    fn handle_read_pair(&mut self, mut read1: FastQ, mut read2: FastQ) -> std::io::Result<&mut Self> {
        let Some(trimmed1) = self.process_read(&mut read1, ReadSide::R1) else {
            return Ok(self);
        };
        let Some(trimmed2) = self.process_read(&mut read2, ReadSide::R2) else {
            return Ok(self);
        };

        self.output_read(trimmed1, ReadSide::R1)?;
        self.output_read(trimmed2, ReadSide::R2)?;
        Ok(self)
    }

    /// Handles all the reads for the given read `side` using [`process_read`].
    ///
    /// [`process_read`]: Self::process_read
    fn handle_single_reads(
        &mut self, mut reader: impl Iterator<Item = Result<FastQ, std::io::Error>>, side: ReadSide,
    ) -> std::io::Result<&mut Self> {
        reader.try_for_each(|read| self.handle_read(read?, side).map(|_| ()))?;
        Ok(self)
    }

    /// Handles all paired reads without filtering widows.
    fn handle_paired_reads_no_filter(
        &mut self, mut reader1: impl Iterator<Item = Result<FastQ, std::io::Error>>,
        mut reader2: impl Iterator<Item = Result<FastQ, std::io::Error>>,
    ) -> std::io::Result<&mut Self> {
        if Self::PRESERVE_ORDER {
            for r1 in reader1.by_ref() {
                self.handle_read(r1?, ReadSide::R1)?;
                let Some(r2) = reader2.next() else {
                    break;
                };
                self.handle_read(r2?, ReadSide::R2)?;
            }
        }

        self.handle_single_reads(reader1, ReadSide::R1)?;
        self.handle_single_reads(reader2, ReadSide::R2)
    }

    /// Handles all reads with widow filtering. Strict checking is performed: if
    /// any mismatch is found between the headers, an error is issued.
    fn handle_paired_reads_with_filter_strict<I, M>(
        &mut self, reader1: I, mut reader2: I, header_mismatch: M, extra_read: std::io::Error,
    ) -> std::io::Result<&mut Self>
    where
        I: Iterator<Item = Result<FastQ, std::io::Error>>,
        M: FnOnce(std::io::Error) -> std::io::Error, {
        for r1 in reader1 {
            let Some(r2) = reader2.next() else {
                return Err(extra_read);
            };
            let r1 = r1?;
            let r2 = r2?;
            if let Err(e) = check_paired_headers(&r1, &r2) {
                return Err(header_mismatch(e));
            }
            self.handle_read_pair(r1, r2)?;
        }

        if reader2.next().is_some() {
            return Err(extra_read);
        }

        Ok(self)
    }

    /// Handles all reads with widow filtering. Weak checking is performed: if
    /// any mismatch is found between the headers, a warning is issued and
    /// processing continues without widow filtering.
    fn handle_paired_reads_with_filter_weak<I, M, X>(
        &mut self, mut reader1: I, mut reader2: I, header_mismatch_warning: M, extra_read_warning: X,
    ) -> std::io::Result<&mut Self>
    where
        I: Iterator<Item = Result<FastQ, std::io::Error>>,
        M: FnOnce(std::io::Error) -> String,
        X: FnOnce() -> String, {
        for r1 in reader1.by_ref() {
            let Some(r2) = reader2.next() else {
                eprintln!("{}", extra_read_warning());
                return self.handle_single_reads(std::iter::once(r1).chain(reader1), ReadSide::R1);
            };
            let r1 = r1?;
            let r2 = r2?;
            if let Err(e) = check_paired_headers(&r1, &r2) {
                eprintln!("{}", header_mismatch_warning(e));
                return self.handle_paired_reads_no_filter(
                    std::iter::once(Ok(r1)).chain(reader1),
                    std::iter::once(Ok(r2)).chain(reader2),
                );
            }
            self.handle_read_pair(r1, r2)?;
        }

        if let Some(r2) = reader2.next() {
            eprintln!("{}", extra_read_warning());
            self.handle_single_reads(std::iter::once(r2).chain(reader2), ReadSide::R2)?;
        }

        Ok(self)
    }
}

/// A trait unifying the record types that have a header field.
pub trait RecordWithHeader {
    /// Gets the header from the record.
    #[must_use]
    fn header(&self) -> &str;
}

impl RecordWithHeader for FastQ {
    #[inline]
    fn header(&self) -> &str {
        &self.header
    }
}

impl RecordWithHeader for FastaSeq {
    #[inline]
    fn header(&self) -> &str {
        &self.name
    }
}

/// An iterator over paired reads. This behaves similarly to a zipped iterator,
/// but ensures the input iterators are the same length and provides appropriate
/// error messages.
pub struct ZipPairedReads<I, J, A, B>
where
    I: Iterator<Item = std::io::Result<A>>,
    J: Iterator<Item = std::io::Result<B>>,
    A: RecordWithHeader,
    B: RecordWithHeader, {
    reads1: I,
    reads2: J,
}

impl<I, J, A, B> Iterator for ZipPairedReads<I, J, A, B>
where
    I: Iterator<Item = std::io::Result<A>>,
    J: Iterator<Item = std::io::Result<B>>,
    A: RecordWithHeader,
    B: RecordWithHeader,
{
    type Item = std::io::Result<(A, B)>;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.reads1.next(), self.reads2.next()) {
            (Some(a), Some(b)) => {
                let (a, b) = (unwrap_or_return_some_err!(a), unwrap_or_return_some_err!(b));
                if check_paired_headers(&a, &b).is_ok() {
                    Some(Ok((a, b)))
                } else {
                    Some(Err(IOError::new(
                        ErrorKind::InvalidInput,
                        format!(
                            "Paired read IDs out of sync:\n\t{h1}\n\t{h2}\n",
                            h1 = a.header(),
                            h2 = b.header()
                        ),
                    )))
                }
            }
            (Some(read1), None) => {
                let read1 = unwrap_or_return_some_err!(read1);
                Some(Err(std::io::Error::other(format!(
                    "Zip reads failed! An extra read in the first file was found with header {header1}",
                    header1 = read1.header()
                ))))
            }
            (None, Some(read2)) => {
                let read2 = unwrap_or_return_some_err!(read2);
                Some(Err(std::io::Error::other(format!(
                    "Zip reads failed! An extra read in the second file was found with header {header2}",
                    header2 = read2.header()
                ))))
            }
            (None, None) => None,
        }
    }
}

/// An extension trait for forming a [`ZipPairedReads`] iterator.
pub trait ZipPairedReadsExt<A>: Sized + Iterator<Item = std::io::Result<A>>
where
    A: RecordWithHeader, {
    /// Zips the iterator of reads with their paired reads. Any IO error is
    /// propagated, and a mismatch in iterator lengths will yield an error
    /// containing the header of the extra item.
    #[inline]
    #[must_use]
    fn zip_paired_reads<J, B>(self, other: J) -> ZipPairedReads<Self, J::IntoIter, A, B>
    where
        J: IntoIterator<Item = std::io::Result<B>>,
        B: RecordWithHeader, {
        ZipPairedReads {
            reads1: self,
            reads2: other.into_iter(),
        }
    }
}

impl<I, A> ZipPairedReadsExt<A> for I
where
    I: Iterator<Item = std::io::Result<A>>,
    A: RecordWithHeader,
{
}

pub struct DeinterleavedPairedReads<I, A>(I)
where
    I: Iterator<Item = std::io::Result<A>>,
    A: RecordWithHeader;

impl<I, A> Iterator for DeinterleavedPairedReads<I, A>
where
    I: Iterator<Item = std::io::Result<A>>,
    A: RecordWithHeader,
{
    type Item = std::io::Result<[A; 2]>;

    fn next(&mut self) -> Option<Self::Item> {
        let read1 = unwrap_or_return_some_err!(self.0.next()?);
        if let Some(read2) = self.0.next() {
            let read2 = unwrap_or_return_some_err!(read2);
            if check_paired_headers(&read1, &read2).is_ok() {
                Some(Ok([read1, read2]))
            } else {
                Some(Err(IOError::new(
                    ErrorKind::InvalidInput,
                    format!(
                        "Paired read IDs out of sync:\n\t{h1}\n\t{h2}\n",
                        h1 = read1.header(),
                        h2 = read2.header()
                    ),
                )))
            }
        } else {
            Some(Err(std::io::Error::other(format!(
                "An odd number of reads was found while de-interleaving: {header1}",
                header1 = read1.header()
            ))))
        }
    }
}

pub trait DeinterleavedPairedReadsExt<A>: Sized + Iterator<Item = std::io::Result<A>>
where
    A: RecordWithHeader, {
    #[inline]
    #[must_use]
    fn deinterleave(self) -> DeinterleavedPairedReads<Self, A> {
        DeinterleavedPairedReads(self)
    }
}

impl<I, A> DeinterleavedPairedReadsExt<A> for I
where
    I: Iterator<Item = std::io::Result<A>>,
    A: RecordWithHeader,
{
}
