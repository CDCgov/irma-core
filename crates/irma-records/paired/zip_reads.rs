use crate::paired::PairedHeaderError;
use std::{
    error::Error,
    fmt::{Debug, Display},
    marker::PhantomData,
    path::Path,
};
use zoe::{
    data::{
        err::{ErrorWithContext, GetCode, WithErrorContext, WithSubitem},
        records::HeaderReadable,
    },
    unwrap_or_return_some_err,
};

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
    /// Errors with the headers, such as failure to parse or being mismatched
    BadHeaders { records: [A; 2], source: PairedHeaderError },
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
                    .with_path_context("An extra read was found in file", path1)
                    .into()
            }
            ZipReadsError::ExtraSecondRead(r2) => {
                std::io::Error::other(format!("Unexpected read found with header: {header2}", header2 = r2.header()))
                    .with_path_context("An extra read was found in file", path2)
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
            ZipPairedReadsError::BadHeaders {
                records: [r1, r2],
                source,
            } => source
                .with_context(format!(
                    "Did not find corresponding paired reads in {path1} and {path2}",
                    path1 = path1.display(),
                    path2 = path2.display(),
                ))
                .with_subitem(format!("Header 1: {header1}", header1 = r1.header()))
                .with_subitem(format!("Header 2: {header2}", header2 = r2.header()))
                .into(),
            ZipPairedReadsError::ExtraFirstRead(r1) => {
                std::io::Error::other(format!("Unexpected read found with header: {header1}", header1 = r1.header()))
                    .with_path_context("An extra read was found in file", path1)
                    .into()
            }
            ZipPairedReadsError::ExtraSecondRead(r2) => {
                std::io::Error::other(format!("Unexpected read found with header: {header2}", header2 = r2.header()))
                    .with_path_context("An extra read was found in file", path2)
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

impl<A> From<ZipPairedReadsError<A>> for std::io::Error
where
    A: HeaderReadable + Debug + Sync + Send + 'static,
{
    #[inline]
    fn from(value: ZipPairedReadsError<A>) -> Self {
        match value {
            ZipPairedReadsError::IoError(e) => e,
            other => std::io::Error::other(other),
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
        if let Err(source) = super::check_paired_headers(&r1, &r2) {
            Err(ZipPairedReadsError::BadHeaders {
                records: [r1, r2],
                source,
            })
        } else {
            Ok([r1, r2])
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
            ZipPairedReadsError::BadHeaders { records: [r1, r2], .. } => {
                let err = ErrorWithContext::new("Paired read IDs out of sync:")
                    .with_subitem(format!("Header 1: {h1}", h1 = r1.header()))
                    .with_subitem(format!("Header 2: {h2}", h2 = r2.header()));

                write!(f, "{err}")
            }
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
            ZipPairedReadsError::BadHeaders { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl<A: HeaderReadable> GetCode for ZipPairedReadsError<A> {
    fn get_code(&self) -> i32 {
        match self {
            ZipPairedReadsError::IoError(e) => e.get_code(),
            ZipPairedReadsError::BadHeaders { source, .. } => source.get_code(),
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
