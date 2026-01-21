use std::{fmt::Display, io::Read, path::Path};
use zoe::{
    data::err::{ResultWithErrorContext, WithErrorContext},
    prelude::{FastQReader, FastaReader},
};

mod fastx;
mod readers;
mod write_records;
mod writers;

pub use fastx::*;
pub(crate) use readers::*;
pub use write_records::*;
pub use writers::*;

/// A wrapper around a fallible iterator, which adds context to any errors in
/// the items.
pub struct IterWithContext<I> {
    /// The inner, fallible iterator.
    iter:        I,
    /// The context to add for any errors.
    description: String,
}

impl<I> IterWithContext<I> {
    /// Retrieves a reference to the inner iterator.
    pub fn inner_iter(&self) -> &I {
        &self.iter
    }
}

impl<I, V, E> Iterator for IterWithContext<I>
where
    I: Iterator<Item = Result<V, E>>,
    E: WithErrorContext,
{
    type Item = std::io::Result<V>;

    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|val| Ok(val.with_context(&self.description)?))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.iter.size_hint()
    }

    fn count(self) -> usize
    where
        Self: Sized, {
        self.iter.count()
    }

    fn last(self) -> Option<Self::Item>
    where
        Self: Sized, {
        self.iter.last().map(|val| Ok(val.with_context(&self.description)?))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.iter.nth(n).map(|val| Ok(val.with_context(&self.description)?))
    }

    fn try_fold<B, F, R>(&mut self, init: B, mut f: F) -> R
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> R,
        R: std::ops::Try<Output = B>, {
        self.iter.try_fold(init, |accum, val| {
            f(accum, val.with_context(&self.description).map_err(Into::into))
        })
    }

    fn fold<B, F>(self, init: B, mut f: F) -> B
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> B, {
        self.iter.fold(init, |accum, val| {
            f(accum, val.with_context(&self.description).map_err(Into::into))
        })
    }
}

/// An extension trait for fallible iterators allowing additional context to be
/// added to each [`Err`] variant that gets yielded (via a
/// [`ErrorWithContext`]).
///
/// [`ErrorWithContext`]: zoe::data::err::ErrorWithContext
pub trait IterWithErrorContext: Sized {
    /// Wraps any errors that get yielded in an [`ErrorWithContext`] with the
    /// given description.
    ///
    /// The `description` field may be anything implementing `Into<String>`.
    /// Passing an owned `String` avoids an extra allocation.
    ///
    /// [`ErrorWithContext`]: zoe::data::err::ErrorWithContext
    fn iter_with_context(self, description: impl Into<String>) -> IterWithContext<Self>;

    /// Convenience function for adding file context to an yielded errors.
    ///
    /// The context will be formatted as `msg: file`. The `msg` field may be
    /// anything implementing [`Display`].
    fn iter_with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> IterWithContext<Self>;
}

impl<I, V, E> IterWithErrorContext for I
where
    I: Iterator<Item = Result<V, E>>,
    E: WithErrorContext,
{
    fn iter_with_context(self, description: impl Into<String>) -> IterWithContext<Self> {
        IterWithContext {
            iter:        self,
            description: description.into(),
        }
    }

    fn iter_with_file_context(self, msg: impl Display, file: impl AsRef<Path>) -> IterWithContext<Self> {
        Self::iter_with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

pub enum DispatchFastX<R>
where
    R: Read, {
    Fastq(IterWithContext<FastQReader<R>>),
    Fasta(IterWithContext<FastaReader<R>>),
}

impl<R: Read> IterWithContext<FastXReader<R>> {
    pub fn dispatch(self) -> DispatchFastX<R> {
        let IterWithContext { iter, description } = self;

        match iter {
            FastXReader::Fastq(iter) => DispatchFastX::Fastq(iter.iter_with_context(description)),
            FastXReader::Fasta(iter) => DispatchFastX::Fasta(iter.iter_with_context(description)),
        }
    }
}

/// Checks that all the provided paths are distinct from each other.
///
/// ## Errors
///
/// If any of the paths are equal, then an appropriate error message is
/// provided.
pub(crate) fn check_distinct_files(
    input1: impl AsRef<Path>, input2: Option<impl AsRef<Path>>, output1: Option<impl AsRef<Path>>,
    output2: Option<impl AsRef<Path>>,
) -> std::io::Result<()> {
    let input1 = Some(input1.as_ref());
    let input2 = input2.as_ref().map(AsRef::as_ref);
    let output1 = output1.as_ref().map(AsRef::as_ref);
    let output2 = output2.as_ref().map(AsRef::as_ref);

    if input1 == input2 {
        Err(std::io::Error::other("The two input files are the same"))
    } else if input1 == output1 {
        Err(std::io::Error::other(
            "The first input file is the same as the first output file",
        ))
    } else if input1 == output2 {
        Err(std::io::Error::other(
            "The first input file is the same as the second output file",
        ))
    } else if input2 == output1 && input2.is_some() {
        Err(std::io::Error::other(
            "The second input file is the same as the first output file",
        ))
    } else if input2 == output2 && input2.is_some() {
        Err(std::io::Error::other(
            "The second input file is the same as the second output file",
        ))
    } else if output1 == output2 && output1.is_some() {
        Err(std::io::Error::other("The two output files are the same"))
    } else {
        Ok(())
    }
}

/// Checks whether a file is a [gzip
/// file](https://www.rfc-editor.org/rfc/rfc1952#page-5).
///
/// This is currently done naively by seeing if it ends with a `gz` extension.
#[inline]
pub(crate) fn is_gz<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref().extension().is_some_and(|ext| ext == "gz")
}

pub enum FromFilenameType {
    Reader,
    Writer,
}

/// A trait unifying file-like structs or enums that can be created from a
/// filename. This is designed for types implementing [`Read`] or [`Write`].
pub(crate) trait FromFilename
where
    Self: Sized, {
    const TYPE: FromFilenameType;

    /// Creates the reader/writer from the path, without adding any context to
    /// error messages.
    fn from_filename_no_context<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>;

    /// Creates the reader/writer from the path.
    ///
    /// ## Errors
    ///
    /// Context including the path is added to the error.
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        let msg = match Self::TYPE {
            FromFilenameType::Reader => "Failed to open file",
            FromFilenameType::Writer => "Failed to create file",
        };

        Ok(Self::from_filename_no_context(&path).with_file_context(msg, path)?)
    }
}

pub(crate) trait FromOptionalFilename: FromFilename {
    const DEFAULT_NAME: &str;

    fn on_none() -> Self;

    /// Creates the reader/writer from an optional path, using [`default`] if it
    /// is not provided, without adding any context to error messages.
    ///
    /// [`default`]: Default::default
    fn from_optional_filename_no_context<P>(path: Option<P>) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        match path {
            Some(path) => Self::from_filename_no_context(path),
            None => Ok(Self::on_none()),
        }
    }

    /// Creates the reader/writer from an optional path, using [`default`] if it
    /// is not provided.
    ///
    /// ## Errors
    ///
    /// Context including the path is added to the error.
    ///
    /// [`default`]: Default::default
    fn from_optional_filename<P>(path: Option<P>) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        match path {
            Some(path) => Self::from_filename(path),
            None => Ok(Self::on_none()),
        }
    }
}

pub(crate) trait IterFromFilename: IterWithErrorContext
where
    Self: Sized, {
    fn get_record_name(&self) -> &'static str;

    fn from_filename_no_context<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>;

    fn from_filename<P>(path: P) -> std::io::Result<IterWithContext<Self>>
    where
        P: AsRef<Path>, {
        let record_reader = Self::from_filename_no_context(&path).with_file_context("Failed to open file", &path)?;

        let record_name = record_reader.get_record_name();

        Ok(record_reader.iter_with_file_context(format!("Invalid {record_name} record while reading file"), path))
    }
}

pub(crate) trait IterFromOptionalFilename: IterFromFilename
where
    Self: Sized, {
    fn on_none_no_context() -> std::io::Result<Self>;

    fn on_none() -> std::io::Result<IterWithContext<Self>>;

    fn from_optional_filename_no_context<P>(path: Option<P>) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        match path {
            Some(path) => Self::from_filename_no_context(path),
            None => Self::on_none_no_context(),
        }
    }

    fn from_optional_filename<P>(path: Option<P>) -> std::io::Result<IterWithContext<Self>>
    where
        P: AsRef<Path>, {
        match path {
            Some(path) => Self::from_filename(path),
            None => Self::on_none(),
        }
    }
}
