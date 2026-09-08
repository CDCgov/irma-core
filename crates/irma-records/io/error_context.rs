use crate::io::{FastXReader, FastXType, Finish};
use std::{
    fmt::Display,
    io::{BufRead, Read, Write},
    path::Path,
};
use zoe::data::{
    err::{ResultWithErrorContext, WithErrorContext},
    fasta::FastaReader,
    fastq::FastQReader,
};

/// A wrapper around a struct with fallible methods/traits, which adds context to
/// those errors. For example, this could wrap a fallible iterator, writer, or
/// reader.
#[derive(Clone, Eq, PartialEq, Hash, Debug)]
pub struct WithContext<T> {
    /// The inner type with fallible operations.
    inner:       T,
    /// The context to add for any errors.
    description: String,
}

impl<T> WithContext<T> {
    /// Manually constructs a new [`WithContext`].
    ///
    /// For support with [`Write`], [`Read`], or fallible iterators, consider
    /// using the `with_context` extension methods (from the
    /// [`WriterWithContext`], [`ReaderWithContext`], and [`IterWithContext`]
    /// extension traits).
    pub fn new(inner: T, description: impl Into<String>) -> Self {
        Self {
            inner,
            description: description.into(),
        }
    }

    /// Runs a custom fallible function on the inner struct, with context added
    /// if it fails.
    ///
    /// Many traits such as [`Iterator`], [`Write`], and [`Read`] are exposed on
    /// [`WithContext`], but for custom methods, this can be used.
    pub fn run_with_context<F, R>(&mut self, f: F) -> std::io::Result<R>
    where
        F: FnOnce(&mut T) -> std::io::Result<R>, {
        Ok(f(&mut self.inner).with_context(&self.description)?)
    }

    /// Maps the inner struct while retaining the same error context.
    pub fn map<U, F>(self, f: F) -> WithContext<U>
    where
        F: FnOnce(T) -> U, {
        WithContext {
            inner:       f(self.inner),
            description: self.description,
        }
    }
}

impl<I, V, E> Iterator for WithContext<I>
where
    I: Iterator<Item = Result<V, E>>,
    E: WithErrorContext,
{
    type Item = std::io::Result<V>;

    fn next(&mut self) -> Option<Self::Item> {
        self.inner.next().map(|val| Ok(val.with_context(&self.description)?))
    }

    fn size_hint(&self) -> (usize, Option<usize>) {
        self.inner.size_hint()
    }

    fn count(self) -> usize
    where
        Self: Sized, {
        self.inner.count()
    }

    fn last(self) -> Option<Self::Item>
    where
        Self: Sized, {
        self.inner.last().map(|val| Ok(val.with_context(&self.description)?))
    }

    fn nth(&mut self, n: usize) -> Option<Self::Item> {
        self.inner.nth(n).map(|val| Ok(val.with_context(&self.description)?))
    }

    fn try_fold<B, F, R>(&mut self, init: B, mut f: F) -> R
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> R,
        R: std::ops::Try<Output = B>, {
        self.inner.try_fold(init, |accum, val| {
            f(accum, val.with_context(&self.description).map_err(Into::into))
        })
    }

    fn fold<B, F>(self, init: B, mut f: F) -> B
    where
        Self: Sized,
        F: FnMut(B, Self::Item) -> B, {
        self.inner.fold(init, |accum, val| {
            f(accum, val.with_context(&self.description).map_err(Into::into))
        })
    }
}

impl<W> Write for WithContext<W>
where
    W: Write,
{
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        self.run_with_context(|writer| writer.write(buf))
    }

    fn flush(&mut self) -> std::io::Result<()> {
        self.run_with_context(|writer| writer.flush())
    }

    fn write_vectored(&mut self, bufs: &[std::io::IoSlice<'_>]) -> std::io::Result<usize> {
        self.run_with_context(|writer| writer.write_vectored(bufs))
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.run_with_context(|writer| writer.write_all(buf))
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        self.run_with_context(|writer| writer.write_fmt(args))
    }
}

impl<W: Finish> Finish for WithContext<W> {
    fn finish(self) -> std::io::Result<()> {
        self.inner.finish().with_context(self.description)?;
        Ok(())
    }
}

impl<R> Read for WithContext<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        Ok(self.inner.read(buf).with_context(&self.description)?)
    }

    fn read_vectored(&mut self, bufs: &mut [std::io::IoSliceMut<'_>]) -> std::io::Result<usize> {
        Ok(self.inner.read_vectored(bufs).with_context(&self.description)?)
    }

    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        Ok(self.inner.read_to_end(buf).with_context(&self.description)?)
    }

    fn read_to_string(&mut self, buf: &mut String) -> std::io::Result<usize> {
        Ok(self.inner.read_to_string(buf).with_context(&self.description)?)
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> std::io::Result<()> {
        Ok(self.inner.read_exact(buf).with_context(&self.description)?)
    }
}

impl<R> BufRead for WithContext<R>
where
    R: BufRead,
{
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        Ok(self.inner.fill_buf().with_context(&self.description)?)
    }

    fn consume(&mut self, amount: usize) {
        self.inner.consume(amount);
    }

    fn read_until(&mut self, byte: u8, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        Ok(self.inner.read_until(byte, buf).with_context(&self.description)?)
    }

    fn skip_until(&mut self, byte: u8) -> std::io::Result<usize> {
        Ok(self.inner.skip_until(byte).with_context(&self.description)?)
    }

    fn read_line(&mut self, buf: &mut String) -> std::io::Result<usize> {
        Ok(self.inner.read_line(buf).with_context(&self.description)?)
    }
}

/// An extension trait for fallible iterators allowing additional context to be
/// added to each [`Err`] variant that gets yielded (via [`WithContext`]).
pub trait IterWithContext: Sized {
    /// Wraps any errors that get yielded in an [`ErrorWithContext`] with the
    /// given description.
    ///
    /// The `description` field may be anything implementing `Into<String>`.
    /// Passing an owned `String` avoids an extra allocation.
    ///
    /// [`ErrorWithContext`]: zoe::data::err::ErrorWithContext
    fn with_context(self, description: impl Into<String>) -> WithContext<Self>;

    /// Convenience function for adding path context to any yielded errors.
    ///
    /// The context will be formatted as `msg: path`. The `msg` field may be
    /// anything implementing [`Display`].
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WithContext<Self>;
}

impl<I, V, E> IterWithContext for I
where
    I: Iterator<Item = Result<V, E>>,
    E: WithErrorContext,
{
    fn with_context(self, description: impl Into<String>) -> WithContext<Self> {
        WithContext {
            inner:       self,
            description: description.into(),
        }
    }

    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WithContext<Self> {
        Self::with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

/// An extension trait for [`Write`] allowing additional context to be added to
/// each failed write (via [`WithContext`]).
pub trait WriterWithContext: Sized {
    /// Wraps any errors that get produced during writing in an
    /// [`ErrorWithContext`] with the given description.
    ///
    /// The `description` field may be anything implementing `Into<String>`.
    /// Passing an owned `String` avoids an extra allocation.
    ///
    /// [`ErrorWithContext`]: zoe::data::err::ErrorWithContext
    fn with_context(self, description: impl Into<String>) -> WithContext<Self>;

    /// Convenience function for adding path context to any produced errors.
    ///
    /// The context will be formatted as `msg: path`. The `msg` field may be
    /// anything implementing [`Display`].
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WithContext<Self>;
}

impl<W> WriterWithContext for W
where
    W: Write,
{
    fn with_context(self, description: impl Into<String>) -> WithContext<Self> {
        WithContext {
            inner:       self,
            description: description.into(),
        }
    }

    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WithContext<Self> {
        Self::with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

/// An extension trait for [`Read`] allowing additional context to be added to
/// each failed read (via [`WithContext`]).
pub trait ReaderWithContext: Sized {
    /// Wraps any errors that get produced during reading in an
    /// [`ErrorWithContext`] with the given description.
    ///
    /// The `description` field may be anything implementing `Into<String>`.
    /// Passing an owned `String` avoids an extra allocation.
    ///
    /// [`ErrorWithContext`]: zoe::data::err::ErrorWithContext
    fn with_context(self, description: impl Into<String>) -> WithContext<Self>;

    /// Convenience function for adding path context to any produced errors.
    ///
    /// The context will be formatted as `msg: path`. The `msg` field may be
    /// anything implementing [`Display`].
    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WithContext<Self>;
}

impl<R> ReaderWithContext for R
where
    R: Read,
{
    fn with_context(self, description: impl Into<String>) -> WithContext<Self> {
        WithContext {
            inner:       self,
            description: description.into(),
        }
    }

    fn with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WithContext<Self> {
        Self::with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

/// A dispatch-ready version [`FastXReader`] where each variant is wrapped with
/// context.
///
/// This is a transposed version of `IterWithContext<FastXReader<...>>`, where
/// each variant contains the [`IterWithContext`] so that the value can be
/// matched on.
pub enum DispatchFastX<R>
where
    R: Read, {
    Fastq(WithContext<FastQReader<R>>),
    Fasta(WithContext<FastaReader<R>>),
}

impl<R: Read> WithContext<FastXReader<R>> {
    /// Returns the type of record that the [`FastXReader`] is parsing.
    pub fn record_type(&self) -> FastXType {
        self.inner.record_type()
    }

    /// Moves the context inside each variant so that the reader type can be
    /// matched on.
    pub fn dispatch(self) -> DispatchFastX<R> {
        let WithContext { inner, description } = self;

        match inner {
            FastXReader::Fastq(inner) => DispatchFastX::Fastq(inner.with_context(description)),
            FastXReader::Fasta(inner) => DispatchFastX::Fasta(inner.with_context(description)),
        }
    }
}
