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
/// If `C` is [`CheckedHeaders`], then paired header checking is performed.
/// Otherwise if it is [`UncheckedHeaders`], then no header checking is
/// performed.
///
/// ## Behavior After Errors
///
/// On an error, the mismatching or unmatched read(s) will be included in the
/// payload of the error. In most realistic cases, the application will likely
/// not call [`next`] on [`ZipReads`] again, either because the error is
/// immediately propagated, or because special handling on the two iterators
/// individually is performed.
///
/// However, in the case that items continue to be accessed from [`ZipReads`]
/// after an error, the behavior is:
///
/// - A [`next`] call always attempts to consume an item from each input
///   iterator. This means that if either iterator yields an IO error, the
///   result from the other will be discarded. IO errors are considered
///   irrecoverable by IRMA-core.
/// - For all other errors related to mismatched headers or unpaired remaining
///   reads, the reads causing the errors are contained within the errors as
///   payloads, and can be accessed and manually handled by the application. The
///   reads are consumed from the iterators, but remaining elements of the
///   iterators may remain, and hence [`ZipReads`] may return additional `Err`
///   or `Ok` variants depending on the error variant, whether the input
///   iterators were fused, etc.
pub struct ZipReads<I, J, A, C>
where
    I: Iterator<Item = std::io::Result<A>>,
    J: Iterator<Item = std::io::Result<A>>,
    A: HeaderReadable, {
    /// The iterator of left reads.
    reads1:  I,
    /// The iterator of right reads.
    reads2:  J,
    /// A phantom field to track the checking level ([`ZipReadsCheckingLevel`]).
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
        // We always attempt to consume an item from each iterator, even if one
        // of them is an IO error. Also, IO errors have precedence before
        // pairing errors, and read1 has precendence before read2
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
        let (i_lower, _) = self.reads1.size_hint();
        let (j_lower, _) = self.reads2.size_hint();

        (i_lower.max(j_lower), None)
    }

    fn try_fold<B, F, R>(&mut self, init: B, mut f: F) -> R
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> R,
        R: std::ops::Try<Output = B>, {
        let mut accum = self.reads1.try_fold(init, |accum, r1| {
            // We always attempt to consume an item from each iterator, even if
            // r1 is an error
            let r2 = self.reads2.next();

            // r1 has precedence before r2, and IO errors have precendence
            // before pairing errors
            let r1 = match r1 {
                Ok(r1) => r1,
                Err(err) => return f(accum, Err(C::new_io_error(err))),
            };

            let r2 = match r2 {
                None => return f(accum, Err(C::new_extra_first_read(r1))),
                Some(Err(err)) => return f(accum, Err(C::new_io_error(err))),
                Some(Ok(r2)) => r2,
            };

            f(accum, C::zip_pair(r1, r2))
        })?;

        // reads1 had an extra next call containing `None`, which caused
        // `try_fold` to end. Call next on reads2 to balance it again.

        accum = match self.reads2.next() {
            // If reads2 is also None, then next would have yielded None, which
            // terminates try_fold
            None => return R::from_output(accum),
            Some(Ok(r2)) => f(accum, Err(C::new_extra_second_read(r2)))?,
            Some(Err(e)) => f(accum, Err(C::new_io_error(e)))?,
        };

        // Assuming the user is aborting after the first error, this will not be
        // reached, so we defer to default impl
        for item in self {
            accum = f(accum, item)?;
        }
        R::from_output(accum)
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
