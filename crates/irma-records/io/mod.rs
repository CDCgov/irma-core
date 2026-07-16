use std::{
    fmt::Display,
    io::{BufRead, Read, Write},
    path::{Path, PathBuf},
};
use zoe::{
    data::err::{ResultWithErrorContext, WithErrorContext},
    prelude::{FastQReader, FastaReader},
};

mod fastx;
mod open_options;
mod readers;
mod write_records;
mod writers;

pub use fastx::*;
pub use open_options::*;
pub use readers::*;
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

    /// Convenience function for adding path context to any yielded errors.
    ///
    /// The context will be formatted as `msg: path`. The `msg` field may be
    /// anything implementing [`Display`].
    fn iter_with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> IterWithContext<Self>;
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

    fn iter_with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> IterWithContext<Self> {
        Self::iter_with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
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
    Fastq(IterWithContext<FastQReader<R>>),
    Fasta(IterWithContext<FastaReader<R>>),
}

impl<R: Read> IterWithContext<FastXReader<R>> {
    /// Moves the context inside each variant so that the reader type can be
    /// matched on.
    pub fn dispatch(self) -> DispatchFastX<R> {
        let IterWithContext { iter, description } = self;

        match iter {
            FastXReader::Fastq(iter) => DispatchFastX::Fastq(iter.iter_with_context(description)),
            FastXReader::Fasta(iter) => DispatchFastX::Fasta(iter.iter_with_context(description)),
        }
    }
}

/// Checks whether a path represents a device file, such as `/dev/null`.
fn is_linux_device(path: &Path) -> bool {
    path.starts_with("/dev/")
}

/// A trait for validating that input and output paths do not have conflicts.
pub trait ValidatePaths {
    /// Returns the paths that will be read from by the process.
    fn inputs(&self) -> impl IntoIterator<Item = &PathBuf>;

    /// Returns the paths that will be written to by the process.
    fn outputs(&self) -> impl IntoIterator<Item = &PathBuf>;

    /// Validates that no path is both an input and an output, and that all
    /// output paths are distinct.
    ///
    /// Device files (paths beginning with `/dev`/) are ignored.
    ///
    /// ## Errors
    ///
    /// All input paths must exist, and the parent directories of the output
    /// paths must exist. The paths must be successfully canonicalized. All
    /// output paths must be distinct and cannot also be input paths.
    fn validate_paths(&self) -> std::io::Result<()> {
        let inputs = self
            .inputs()
            .into_iter()
            .filter(|path| !is_linux_device(path))
            .map(|path| std::fs::canonicalize(path).with_path_context("Failed to canonicalize path", path));

        let outputs = self
            .outputs()
            .into_iter()
            .filter(|path| !is_linux_device(path))
            .map(|path| {
                // If the output path already exists (including as a symlink),
                // canonicalize it directly so that aliases such as symlinks
                // pointing at an input file are resolved to their real path.
                if path.exists() {
                    return Ok(std::fs::canonicalize(path).with_path_context("Failed to canonicalize path", path)?);
                }

                let filename = path.file_name().ok_or_else(|| {
                    std::io::Error::other(format!("Failed to find filename of path: {path}", path = path.display()))
                })?;
                let parent = match path.parent() {
                    Some(parent) if !parent.as_os_str().is_empty() => parent,
                    _ => Path::new("."),
                };
                let canonical_parent =
                    std::fs::canonicalize(parent).with_path_context("Failed to canonicalize parent path", parent)?;
                Ok(canonical_parent.join(filename))
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        for input1 in inputs {
            let input1 = input1?;

            for output in &outputs {
                if &input1 == output {
                    return Err(std::io::Error::other(format!(
                        "Found a file as both an input and an output: {input1}",
                        input1 = input1.display()
                    )));
                }
            }
        }

        for (i, output1) in outputs.iter().enumerate() {
            let rest = &outputs[i + 1..];
            for output2 in rest {
                if output1 == output2 {
                    return Err(std::io::Error::other(format!(
                        "Two output files were the same: {output1}",
                        output1 = output1.display()
                    )));
                }
            }
        }

        Ok(())
    }
}

/// Checks whether a file is a [gzip
/// file](https://www.rfc-editor.org/rfc/rfc1952#page-5).
///
/// This is currently done naively by seeing if it ends with a `gz` extension.
#[inline]
pub fn is_gz<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref().extension().is_some_and(|ext| ext == "gz")
}

/// A wrapper around a writer of type `W` such that error context is added to
/// any failed writes.
#[derive(Debug)]
pub struct WriterWithContext<W> {
    /// The inner writer.
    writer:      W,
    /// The context to add to any failed writes.
    description: String,
}

impl<W> Write for WriterWithContext<W>
where
    W: Write,
{
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        Ok(self.writer.write(buf).with_context(&self.description)?)
    }

    fn flush(&mut self) -> std::io::Result<()> {
        Ok(self.writer.flush().with_context(&self.description)?)
    }

    fn write_vectored(&mut self, bufs: &[std::io::IoSlice<'_>]) -> std::io::Result<usize> {
        Ok(self.writer.write_vectored(bufs).with_context(&self.description)?)
    }

    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        Ok(self.writer.write_all(buf).with_context(&self.description)?)
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        Ok(self.writer.write_fmt(args).with_context(&self.description)?)
    }
}

/// An extension trait for [`Write`] allowing additional context to be added to
/// each failed write (via [`WriterWithContext`]).
pub trait WriterWithErrorContext: Sized {
    /// Wraps any errors that get produced during writing in an
    /// [`ErrorWithContext`] with the given description.
    ///
    /// The `description` field may be anything implementing `Into<String>`.
    /// Passing an owned `String` avoids an extra allocation.
    ///
    /// [`ErrorWithContext`]: zoe::data::err::ErrorWithContext
    fn writer_with_context(self, description: impl Into<String>) -> WriterWithContext<Self>;

    /// Convenience function for adding path context to any produced errors.
    ///
    /// The context will be formatted as `msg: path`. The `msg` field may be
    /// anything implementing [`Display`].
    fn writer_with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WriterWithContext<Self>;
}

impl<W> WriterWithErrorContext for W
where
    W: Write,
{
    fn writer_with_context(self, description: impl Into<String>) -> WriterWithContext<Self> {
        WriterWithContext {
            writer:      self,
            description: description.into(),
        }
    }

    fn writer_with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> WriterWithContext<Self> {
        Self::writer_with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}

/// A wrapper around a reader of type `R` such that error context is added to
/// any failed reads.
pub struct ReaderWithContext<R> {
    /// The inner reader.
    reader:      R,
    /// The context to add to any failed reads.
    description: String,
}

impl<R> Read for ReaderWithContext<R>
where
    R: Read,
{
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        Ok(self.reader.read(buf).with_context(&self.description)?)
    }

    fn read_vectored(&mut self, bufs: &mut [std::io::IoSliceMut<'_>]) -> std::io::Result<usize> {
        Ok(self.reader.read_vectored(bufs).with_context(&self.description)?)
    }

    fn read_to_end(&mut self, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        Ok(self.reader.read_to_end(buf).with_context(&self.description)?)
    }

    fn read_to_string(&mut self, buf: &mut String) -> std::io::Result<usize> {
        Ok(self.reader.read_to_string(buf).with_context(&self.description)?)
    }

    fn read_exact(&mut self, buf: &mut [u8]) -> std::io::Result<()> {
        Ok(self.reader.read_exact(buf).with_context(&self.description)?)
    }
}

impl<R> BufRead for ReaderWithContext<R>
where
    R: BufRead,
{
    fn fill_buf(&mut self) -> std::io::Result<&[u8]> {
        Ok(self.reader.fill_buf().with_context(&self.description)?)
    }

    fn consume(&mut self, amount: usize) {
        self.reader.consume(amount);
    }

    fn read_until(&mut self, byte: u8, buf: &mut Vec<u8>) -> std::io::Result<usize> {
        Ok(self.reader.read_until(byte, buf).with_context(&self.description)?)
    }

    fn skip_until(&mut self, byte: u8) -> std::io::Result<usize> {
        Ok(self.reader.skip_until(byte).with_context(&self.description)?)
    }

    fn read_line(&mut self, buf: &mut String) -> std::io::Result<usize> {
        Ok(self.reader.read_line(buf).with_context(&self.description)?)
    }
}

/// An extension trait for [`Read`] allowing additional context to be added to
/// each failed read (via [`ReaderWithContext`]).
pub trait ReaderWithErrorContext: Sized {
    /// Wraps any errors that get produced during reading in an
    /// [`ErrorWithContext`] with the given description.
    ///
    /// The `description` field may be anything implementing `Into<String>`.
    /// Passing an owned `String` avoids an extra allocation.
    ///
    /// [`ErrorWithContext`]: zoe::data::err::ErrorWithContext
    fn reader_with_context(self, description: impl Into<String>) -> ReaderWithContext<Self>;

    /// Convenience function for adding path context to any produced errors.
    ///
    /// The context will be formatted as `msg: path`. The `msg` field may be
    /// anything implementing [`Display`].
    fn reader_with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> ReaderWithContext<Self>;
}

impl<R> ReaderWithErrorContext for R
where
    R: Read,
{
    fn reader_with_context(self, description: impl Into<String>) -> ReaderWithContext<Self> {
        ReaderWithContext {
            reader:      self,
            description: description.into(),
        }
    }

    fn reader_with_path_context(self, msg: impl Display, file: impl AsRef<Path>) -> ReaderWithContext<Self> {
        Self::reader_with_context(self, format!("{msg}: '{path}'", path = file.as_ref().display()))
    }
}
