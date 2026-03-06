use std::{
    error::Error,
    fmt::{Debug, Display},
    io::{Error as IOError, ErrorKind},
    marker::PhantomData,
    path::Path,
    simd::Simd,
};
use zoe::{
    data::{
        err::{GetCode, WithErrorContext},
        records::HeaderReadable,
    },
    unwrap_or_return_some_err,
};

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
pub fn check_paired_headers<A: HeaderReadable, B: HeaderReadable>(read1: &A, read2: &B) -> Result<(), std::io::Error> {
    if let Some((id1, _)) = get_molecular_id_side(read1.header(), '0')
        && let Some((id2, _)) = get_molecular_id_side(read2.header(), '0')
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

/// An iterator over paired reads. This behaves similarly to a zipped iterator,
/// but ensures the input iterators are the same length, and provides
/// appropriate error messages.
///
/// If `CHECKED` is true, paired header checking is performed.
pub struct ZipReads<I, J, A, C>
where
    I: Iterator<Item = std::io::Result<A>>,
    J: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable, {
    reads1:  I,
    reads2:  J,
    phantom: PhantomData<C>,
}

/// The error type when [`ZipReads`] is used without paired header checking.
#[derive(Debug)]
pub enum ZipReadsError<A> {
    /// An IO error from one of the readers
    IoError(std::io::Error),
    /// An extra read in the first iterator
    ExtraFirstRead(A),
    /// An extra read in the second iterator
    ExtraSecondRead(A),
}

/// The error type when [`ZipReads`] is used with paired header checking.
#[derive(Debug)]
pub enum ZipPairedReadsError<A> {
    /// An IO error from one of the readers
    IoError(std::io::Error),
    /// A mismatch in the header IDs of the paired reads
    MismatchedHeaders([A; 2]),
    /// An extra read in the first iterator
    ExtraFirstRead(A),
    /// An extra read in the second iterator
    ExtraSecondRead(A),
}

impl<A> ZipReadsError<A>
where
    A: HeaderReadable,
{
    /// Maps the error to include messages and context with the provided paths.
    ///
    /// It is assumed that IO errors already have path context included.
    pub fn add_path_context(self, path1: &Path, path2: &Path) -> std::io::Error {
        match self {
            ZipReadsError::IoError(e) => e,
            ZipReadsError::ExtraFirstRead(r1) => {
                std::io::Error::other(format!("Unexpected read found with header: {header1}", header1 = r1.header()))
                    .with_file_context("An extra read was found in file", path1)
                    .into()
            }
            ZipReadsError::ExtraSecondRead(r2) => {
                std::io::Error::other(format!("Unexpected read found with header: {header2}", header2 = r2.header()))
                    .with_file_context("An extra read was found in file", path2)
                    .into()
            }
        }
    }
}

impl<A> ZipPairedReadsError<A>
where
    A: HeaderReadable,
{
    /// Maps the error to include messages and context with the provided paths.
    ///
    /// It is assumed that IO errors already have path context included.
    pub fn add_path_context(self, path1: &Path, path2: &Path) -> std::io::Error {
        match self {
            ZipPairedReadsError::IoError(e) => e,
            ZipPairedReadsError::MismatchedHeaders([r1, r2]) => std::io::Error::other(format!(
                "Did not find corresponding paired reads in {path1} and {path2}.\n\nHeader 1: {header1}\nHeader 2: {header2}",
                path1 = path1.display(),
                path2 = path2.display(),
                header1 = r1.header(),
                header2 = r2.header(),
            )),
            ZipPairedReadsError::ExtraFirstRead(r1) => {
                std::io::Error::other(format!("Unexpected read found with header: {header1}", header1 = r1.header()))
                    .with_file_context("An extra read was found in file", path1)
                    .into()
            }
            ZipPairedReadsError::ExtraSecondRead(r2) => {
                std::io::Error::other(format!("Unexpected read found with header: {header2}", header2 = r2.header()))
                    .with_file_context("An extra read was found in file", path2)
                    .into()
            }
        }
    }
}

impl<A> From<std::io::Error> for ZipReadsError<A> {
    #[inline]
    fn from(value: std::io::Error) -> Self {
        Self::IoError(value)
    }
}

impl<A> From<std::io::Error> for ZipPairedReadsError<A> {
    #[inline]
    fn from(value: std::io::Error) -> Self {
        Self::IoError(value)
    }
}

impl<A: HeaderReadable> From<ZipReadsError<A>> for std::io::Error {
    #[inline]
    fn from(value: ZipReadsError<A>) -> Self {
        match value {
            ZipReadsError::IoError(e) => e,
            other => std::io::Error::other(other.to_string()),
        }
    }
}

impl<A: HeaderReadable> From<ZipPairedReadsError<A>> for std::io::Error {
    #[inline]
    fn from(value: ZipPairedReadsError<A>) -> Self {
        match value {
            ZipPairedReadsError::IoError(e) => e,
            other => std::io::Error::other(other.to_string()),
        }
    }
}

/// A trait marking the checking level in [`ZipReads`], as well as unifying
/// [`ZipReadsError`] and [`ZipPairedReadsError`] by providing necessary
/// constructors and a method for zipping a pair of reads (which may or may not
/// perform checking).
pub trait ZipReadsCheckingLevel {
    type Err<A>;

    /// Constructs an error variant corresponding to an IO error from one of the
    /// readers.
    #[must_use]
    fn new_io_error<A>(err: std::io::Error) -> Self::Err<A>;

    /// Constructs an error variant corresponding to an extra read in the first
    /// iterator.
    #[must_use]
    fn new_extra_first_read<A>(r1: A) -> Self::Err<A>;

    /// Constructs an error variant corresponding to an extra read in the second
    /// iterator.
    #[must_use]
    fn new_extra_second_read<A>(r2: A) -> Self::Err<A>;

    /// Groups two reads together, potentially performing paired header
    /// checking.
    fn zip_pair<A: HeaderReadable>(r1: A, r2: A) -> Result<[A; 2], Self::Err<A>>;
}

/// A marker for skipping header checking in [`ZipReads`].
pub struct UncheckedHeaders;

/// A marker for performing header checking in [`ZipReads`].
pub struct CheckedHeaders;

impl ZipReadsCheckingLevel for UncheckedHeaders {
    type Err<A> = ZipReadsError<A>;

    #[inline]
    fn new_io_error<A>(err: std::io::Error) -> ZipReadsError<A> {
        ZipReadsError::IoError(err)
    }

    #[inline]
    fn new_extra_first_read<A>(r1: A) -> ZipReadsError<A> {
        ZipReadsError::ExtraFirstRead(r1)
    }

    #[inline]
    fn new_extra_second_read<A>(r2: A) -> ZipReadsError<A> {
        ZipReadsError::ExtraSecondRead(r2)
    }

    #[inline]
    fn zip_pair<A>(r1: A, r2: A) -> Result<[A; 2], ZipReadsError<A>> {
        Ok([r1, r2])
    }
}

impl ZipReadsCheckingLevel for CheckedHeaders {
    type Err<A> = ZipPairedReadsError<A>;

    #[inline]
    fn new_io_error<A>(err: std::io::Error) -> ZipPairedReadsError<A> {
        ZipPairedReadsError::IoError(err)
    }

    #[inline]
    fn new_extra_first_read<A>(r1: A) -> ZipPairedReadsError<A> {
        ZipPairedReadsError::ExtraFirstRead(r1)
    }

    #[inline]
    fn new_extra_second_read<A>(r2: A) -> ZipPairedReadsError<A> {
        ZipPairedReadsError::ExtraSecondRead(r2)
    }

    #[inline]
    fn zip_pair<A: HeaderReadable>(r1: A, r2: A) -> Result<[A; 2], ZipPairedReadsError<A>> {
        if check_paired_headers(&r1, &r2).is_ok() {
            Ok([r1, r2])
        } else {
            Err(ZipPairedReadsError::MismatchedHeaders([r1, r2]))
        }
    }
}

impl<A: HeaderReadable> Display for ZipReadsError<A> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ZipReadsError::IoError(e) => write!(f, "{e}"),
            ZipReadsError::ExtraFirstRead(r1) => write!(
                f,
                "An extra read in the first file was found with header {header1}",
                header1 = r1.header()
            ),
            ZipReadsError::ExtraSecondRead(r2) => write!(
                f,
                "An extra read in the second file was found with header {header2}",
                header2 = r2.header()
            ),
        }
    }
}

impl<A: HeaderReadable> Display for ZipPairedReadsError<A> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ZipPairedReadsError::IoError(e) => write!(f, "{e}"),
            ZipPairedReadsError::MismatchedHeaders([r1, r2]) => write!(
                f,
                "Paired read IDs out of sync:\n\t{h1}\n\t{h2}\n",
                h1 = r1.header(),
                h2 = r2.header()
            ),
            ZipPairedReadsError::ExtraFirstRead(r1) => write!(
                f,
                "An extra read in the first file was found with header {header1}",
                header1 = r1.header()
            ),
            ZipPairedReadsError::ExtraSecondRead(r2) => write!(
                f,
                "An extra read in the second file was found with header {header2}",
                header2 = r2.header()
            ),
        }
    }
}

impl<A: HeaderReadable + Debug> Error for ZipReadsError<A> {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ZipReadsError::IoError(e) => e.source(),
            _ => None,
        }
    }
}

impl<A: HeaderReadable> GetCode for ZipReadsError<A> {
    fn get_code(&self) -> i32 {
        match self {
            ZipReadsError::IoError(e) => e.get_code(),
            _ => 1,
        }
    }
}

impl<A: HeaderReadable + Debug> Error for ZipPairedReadsError<A> {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ZipPairedReadsError::IoError(e) => e.source(),
            _ => None,
        }
    }
}

impl<A: HeaderReadable> GetCode for ZipPairedReadsError<A> {
    fn get_code(&self) -> i32 {
        match self {
            ZipPairedReadsError::IoError(e) => e.get_code(),
            _ => 1,
        }
    }
}

impl<I, J, A, C> Iterator for ZipReads<I, J, A, C>
where
    I: Iterator<Item = std::io::Result<A>>,
    J: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable,
    C: ZipReadsCheckingLevel,
{
    type Item = Result<[A; 2], C::Err<A>>;

    fn next(&mut self) -> Option<Self::Item> {
        match (self.reads1.next(), self.reads2.next()) {
            (Some(read1), Some(read2)) => {
                let read1 = unwrap_or_return_some_err!(read1.map_err(C::new_io_error));
                let read2 = unwrap_or_return_some_err!(read2.map_err(C::new_io_error));
                Some(C::zip_pair(read1, read2))
            }
            (Some(read1), None) => {
                let read1 = unwrap_or_return_some_err!(read1.map_err(C::new_io_error));
                Some(Err(C::new_extra_first_read(read1)))
            }
            (None, Some(read2)) => {
                let read2 = unwrap_or_return_some_err!(read2.map_err(C::new_io_error));
                Some(Err(C::new_extra_second_read(read2)))
            }
            (None, None) => None,
        }
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        let (i_lower, i_upper) = self.reads1.size_hint();
        let (j_lower, j_upper) = self.reads2.size_hint();
        let lower = i_lower.min(j_lower);
        let upper = match (i_upper, j_upper) {
            (Some(i_upper), Some(j_upper)) if i_upper == j_upper => Some(i_upper.min(j_upper)),
            // Potentially unequal iterator lengths, so add 1 for error message
            (Some(i_upper), Some(j_upper)) => Some(i_upper.min(j_upper) + 1),
            // Potentially unequal iterator lengths, so add 1 for error message
            (Some(upper), None) | (None, Some(upper)) => Some(upper + 1),
            (None, None) => None,
        };
        (lower, upper)
    }

    fn try_fold<B, F, R>(&mut self, init: B, mut f: F) -> R
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> R,
        R: std::ops::Try<Output = B>, {
        let accum = self.reads1.try_fold(init, |accum, r1| {
            let r1 = match r1 {
                Ok(r1) => r1,
                Err(err) => return f(accum, Err(C::new_io_error(err))),
            };

            let r2 = match self.reads2.next() {
                None => return f(accum, Err(C::new_extra_first_read(r1))),
                Some(Err(err)) => return f(accum, Err(C::new_io_error(err))),
                Some(Ok(r2)) => r2,
            };

            f(accum, C::zip_pair(r1, r2))
        })?;

        match self.reads2.next() {
            None => R::from_output(accum),
            Some(Ok(r2)) => f(accum, Err(C::new_extra_second_read(r2))),
            Some(Err(e)) => f(accum, Err(C::new_io_error(e))),
        }
    }
}

/// An extension trait for forming a [`ZipReads`] iterator.
pub trait ZipPairedReadsExt<A>: Sized + Iterator<Item = std::io::Result<A>>
where
    A: HeaderReadable, {
    /// Zips the iterator of reads with their paired reads, while checking
    /// compatible headers.
    ///
    /// Any IO error is propagated, and a mismatch in iterator lengths will
    /// yield an error containing the header of the extra item.
    #[inline]
    #[must_use]
    fn zip_paired_reads<J>(self, other: J) -> ZipReads<Self, J::IntoIter, A, CheckedHeaders>
    where
        J: IntoIterator<Item = std::io::Result<A>>, {
        ZipReads {
            reads1:  self,
            reads2:  other.into_iter(),
            phantom: PhantomData,
        }
    }

    /// Zips the iterator of reads with their paired reads, without checking
    /// compatible headers.
    ///
    /// Any IO error is propagated, and a mismatch in iterator lengths will
    /// yield an error containing the header of the extra item.
    #[inline]
    #[must_use]
    fn zip_paired_reads_unchecked<J>(self, other: J) -> ZipReads<Self, J::IntoIter, A, UncheckedHeaders>
    where
        J: IntoIterator<Item = std::io::Result<A>>, {
        ZipReads {
            reads1:  self,
            reads2:  other.into_iter(),
            phantom: PhantomData,
        }
    }
}

impl<I, A> ZipPairedReadsExt<A> for I
where
    I: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable,
{
}

/// The error type for [`DeinterleavedPairedReads`].
#[derive(Debug)]
pub enum DeinterleaveError<A> {
    /// An IO error from the reader
    IoError(std::io::Error),
    /// A mismatch in the header IDs of the paired subsequent reads
    MismatchedHeaders([A; 2]),
    /// An odd number of reads in the iterator
    OddNumberOfReads(A),
}

impl<A> From<std::io::Error> for DeinterleaveError<A> {
    #[inline]
    fn from(value: std::io::Error) -> Self {
        Self::IoError(value)
    }
}

impl<A: HeaderReadable> From<DeinterleaveError<A>> for std::io::Error {
    #[inline]
    fn from(value: DeinterleaveError<A>) -> Self {
        match value {
            DeinterleaveError::IoError(e) => e,
            other => std::io::Error::other(other.to_string()),
        }
    }
}

impl<A: HeaderReadable> Display for DeinterleaveError<A> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            DeinterleaveError::IoError(e) => write!(f, "{e}"),
            DeinterleaveError::MismatchedHeaders([r1, r2]) => write!(
                f,
                "Paired read IDs out of sync:\n\t{h1}\n\t{h2}\n",
                h1 = r1.header(),
                h2 = r2.header()
            ),
            DeinterleaveError::OddNumberOfReads(r1) => write!(
                f,
                "An odd number of reads was found while de-interleaving. See header: {header1}",
                header1 = r1.header()
            ),
        }
    }
}

impl<A: HeaderReadable + Debug> Error for DeinterleaveError<A> {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            DeinterleaveError::IoError(e) => e.source(),
            _ => None,
        }
    }
}

impl<A: HeaderReadable> GetCode for DeinterleaveError<A> {
    fn get_code(&self) -> i32 {
        match self {
            DeinterleaveError::IoError(e) => e.get_code(),
            _ => 1,
        }
    }
}

impl<A> DeinterleaveError<A>
where
    A: HeaderReadable,
{
    /// Maps the error to include messages and context with the provided path.
    ///
    /// It is assumed that IO errors already have path context included.
    pub fn add_path_context(self, path: &Path) -> std::io::Error {
        match self {
            DeinterleaveError::IoError(e) => e,
            DeinterleaveError::MismatchedHeaders([r1, r2]) => std::io::Error::other(format!(
                "Paired read IDs out of sync:\n| Header 1: {header1}\n| Header 2: {header2}",
                header1 = r1.header(),
                header2 = r2.header(),
            ))
            .with_file_context("Failed to deinterleave the reads in file", path)
            .into(),
            e @ DeinterleaveError::OddNumberOfReads(_) => std::io::Error::from(e)
                .with_file_context("Failed to deinterleave the reads in file", path)
                .into(),
        }
    }
}

pub struct DeinterleavedPairedReads<I, A>(I)
where
    I: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable;

impl<I, A> Iterator for DeinterleavedPairedReads<I, A>
where
    I: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable,
{
    type Item = Result<[A; 2], DeinterleaveError<A>>;

    fn next(&mut self) -> Option<Self::Item> {
        let read1 = unwrap_or_return_some_err!(self.0.next()?.map_err(DeinterleaveError::IoError));

        if let Some(read2) = self.0.next() {
            let read2 = unwrap_or_return_some_err!(read2.map_err(DeinterleaveError::IoError));

            if check_paired_headers(&read1, &read2).is_ok() {
                Some(Ok([read1, read2]))
            } else {
                Some(Err(DeinterleaveError::MismatchedHeaders([read1, read2])))
            }
        } else {
            Some(Err(DeinterleaveError::OddNumberOfReads(read1)))
        }
    }
}

pub trait DeinterleavedPairedReadsExt<A>: Sized + Iterator<Item = std::io::Result<A>>
where
    A: HeaderReadable, {
    #[inline]
    #[must_use]
    fn deinterleave(self) -> DeinterleavedPairedReads<Self, A> {
        DeinterleavedPairedReads(self)
    }
}

impl<I, A> DeinterleavedPairedReadsExt<A> for I
where
    I: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable,
{
}
