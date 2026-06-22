use crate::paired::PairedHeaderError;
use std::{
    error::Error,
    fmt::{Debug, Display},
    path::Path,
};
use zoe::{
    data::{
        err::{ErrorWithContext, GetCode, WithErrorContext, WithSubitem},
        records::HeaderReadable,
    },
    unwrap_or_return_some_err,
};

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

            if let Err(source) = super::check_paired_headers(&read1, &read2) {
                Some(Err(DeinterleaveError::BadHeaders {
                    records: [read1, read2],
                    source,
                }))
            } else {
                Some(Ok([read1, read2]))
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

/// The error type for [`DeinterleavedPairedReads`].
#[derive(Debug)]
pub enum DeinterleaveError<A> {
    /// An IO error from the reader
    IoError(std::io::Error),
    /// Errors with the headers of subsequent reads, such as failure to parse or
    /// being mismatched
    BadHeaders { records: [A; 2], source: PairedHeaderError },
    /// An odd number of reads in the iterator
    OddNumberOfReads(A),
}

impl<A> From<std::io::Error> for DeinterleaveError<A> {
    #[inline]
    fn from(value: std::io::Error) -> Self {
        Self::IoError(value)
    }
}

impl<A> From<DeinterleaveError<A>> for std::io::Error
where
    A: HeaderReadable + Debug + Sync + Send + 'static,
{
    #[inline]
    fn from(value: DeinterleaveError<A>) -> Self {
        match value {
            DeinterleaveError::IoError(e) => e,
            other => std::io::Error::other(other),
        }
    }
}

impl<A: HeaderReadable> Display for DeinterleaveError<A> {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            DeinterleaveError::IoError(e) => write!(f, "{e}"),
            DeinterleaveError::BadHeaders { records: [r1, r2], .. } => {
                let err = ErrorWithContext::new("Paired read IDs out of sync:")
                    .with_subitem(format!("Header 1: {h1}", h1 = r1.header()))
                    .with_subitem(format!("Header 2: {h2}", h2 = r2.header()));

                write!(f, "{err}")
            }
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
            DeinterleaveError::BadHeaders { source, .. } => Some(source),
            _ => None,
        }
    }
}

impl<A: HeaderReadable> GetCode for DeinterleaveError<A> {
    fn get_code(&self) -> i32 {
        match self {
            DeinterleaveError::IoError(e) => e.get_code(),
            DeinterleaveError::BadHeaders { source, .. } => source.get_code(),
            _ => 1,
        }
    }
}

impl<A> DeinterleaveError<A>
where
    A: HeaderReadable + Debug + Sync + Send + 'static,
{
    /// Maps the error to include messages and context with the provided path.
    ///
    /// It is assumed that IO errors already have path context included.
    pub fn add_path_context(self, path: &Path) -> std::io::Error {
        match self {
            DeinterleaveError::IoError(e) => e,
            DeinterleaveError::BadHeaders {
                records: [r1, r2],
                source,
            } => source
                .with_context("Paired read IDs out of sync:")
                .with_subitem(format!("Header 1: {header1}", header1 = r1.header()))
                .with_subitem(format!("Header 2: {header2}", header2 = r2.header()))
                .with_path_context("Failed to deinterleave the reads in file", path)
                .into(),
            e @ DeinterleaveError::OddNumberOfReads(_) => std::io::Error::from(e)
                .with_path_context("Failed to deinterleave the reads in file", path)
                .into(),
        }
    }
}
