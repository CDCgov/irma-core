use std::{fmt::Display, io::Read, path::Path};
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
    fn identical_path_input_output(path: &Path) -> std::io::Error {
        std::io::Error::other(format!(
            "An identical path was found for an input file and output file: {}",
            path.display()
        ))
    }

    let input1 = input1.as_ref();
    let input2 = input2.as_ref().map(AsRef::as_ref);
    let output1 = output1.as_ref().map(AsRef::as_ref);
    let output2 = output2.as_ref().map(AsRef::as_ref);

    if let Some(input2) = input2
        && input1 == input2
    {
        Err(std::io::Error::other(format!(
            "An identical path was found for the two input files: {}",
            input1.display()
        )))
    } else if let Some(output1) = output1
        && input1 == output1
    {
        Err(identical_path_input_output(input1))
    } else if let Some(output2) = output2
        && input1 == output2
    {
        Err(identical_path_input_output(input1))
    } else if let Some(input2) = input2
        && let Some(output1) = output1
        && input2 == output1
    {
        Err(identical_path_input_output(input1))
    } else if let Some(input2) = input2
        && let Some(output2) = output2
        && input2 == output2
    {
        Err(identical_path_input_output(input1))
    } else if let Some(output1) = output1
        && let Some(output2) = output2
        && output1 == output2
    {
        Err(std::io::Error::other(format!(
            "An identical path was found for the two output files: {}",
            output1.display(),
        )))
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
