//! Builder patterns for opening input files and creating output files in
//! IRMA-core.
//!
//! These builder patterns support the enum types in [`readers`] and
//! [`writers`], support potentially paired inputs/outputs with
//! [`RecordReaders`] and [`RecordWriters`], and automatically add error
//! context.
//!
//! See [`InputOptions`] and [`OutputOptions`] for more details and step-by-step
//! instructions.
//!
//! [`readers`]: crate::io::readers
//! [`writers`]: crate::io::writers

mod context;
mod input_options;
mod output_options;

pub use context::*;
pub use input_options::*;
pub use output_options::*;

use crate::io::{PairedWriters, RecordReaders, RecordWriters};
use std::path::Path;

/// A struct holding two optional paired paths.
///
/// This is used by [`InputOptions::new_from_opt_paths`] and
/// [`OutputOptions::new_from_opt_paths`].
pub struct OptionalPaths<'a> {
    /// The optional first path. If [`None`], this corresponds to
    /// [`InputType::Stdin`] or [`OutputType::Stdout`].
    path1: Option<&'a Path>,
    /// The optional second path. If [`None`], the means `path1` is unpaired.
    path2: Option<&'a Path>,
}

impl<'a> OptionalPaths<'a> {
    /// Attempts to map the two paths to [`RecordReaders`] via a fallible
    /// closure.
    ///
    /// This is similar to [`PairedStruct::try_map`], but it changes the type of
    /// the paired data to [`RecordReaders`].
    #[allow(dead_code)]
    fn try_map_readers<U, F>(self, op: F) -> Result<RecordReaders<U>, PairedErrors>
    where
        F: Fn(Option<&'a Path>) -> std::io::Result<U>, {
        let reader1 = op(self.path1).map_err(PairedErrors::Err1)?;
        let reader2 = self
            .path2
            .map(|path2| op(Some(path2)).map_err(PairedErrors::Err2))
            .transpose()?;

        Ok(RecordReaders { reader1, reader2 })
    }

    /// Attempts to map the two paths to [`RecordWriters`] via a fallible
    /// closure.
    ///
    /// This is similar to [`PairedStruct::try_map`], but it changes the type of
    /// the paired data to [`RecordWriters`].
    fn try_map_writers<U, F>(self, op: F) -> Result<RecordWriters<U>, PairedErrors>
    where
        F: Fn(Option<&'a Path>) -> std::io::Result<U>, {
        let writer1 = op(self.path1).map_err(PairedErrors::Err1)?;
        let writer2 = self
            .path2
            .map(|path2| op(Some(path2)).map_err(PairedErrors::Err2))
            .transpose()?;

        Ok(RecordWriters::new(writer1, writer2))
    }
}

/// A trait providing basic mapping operations for structs containing
/// potentially paired data.
trait PairedStruct<T> {
    /// The same type as `Self`, but with a different generic.
    type Paired<U>;

    /// Maps the paired data (potentially to a different type) via a closure.
    #[allow(dead_code)]
    fn map<U, F>(self, op: F) -> Self::Paired<U>
    where
        F: Fn(T) -> U;

    /// Attempts to map the paired data (potentially to a different type) via a
    /// fallible closure.
    ///
    /// An error type of [`PairedErrors`] is used to indicate whether the first
    /// or second data field failed to map. This function short-circuits if the
    /// first field fails to map.
    fn try_map<U, F>(self, op: F) -> Result<Self::Paired<U>, PairedErrors>
    where
        F: Fn(T) -> std::io::Result<U>;
}

impl<R> PairedStruct<R> for RecordReaders<R> {
    type Paired<U> = RecordReaders<U>;

    fn map<U, F>(self, op: F) -> Self::Paired<U>
    where
        F: Fn(R) -> U, {
        RecordReaders {
            reader1: op(self.reader1),
            reader2: self.reader2.map(op),
        }
    }

    fn try_map<U, F>(self, op: F) -> Result<Self::Paired<U>, PairedErrors>
    where
        F: Fn(R) -> std::io::Result<U>, {
        let reader1 = op(self.reader1).map_err(PairedErrors::Err1)?;
        let reader2 = self
            .reader2
            .map(|reader2| op(reader2).map_err(PairedErrors::Err2))
            .transpose()?;

        Ok(RecordReaders { reader1, reader2 })
    }
}

impl<W> PairedStruct<W> for PairedWriters<W> {
    type Paired<U> = PairedWriters<U>;

    fn map<U, F>(self, op: F) -> Self::Paired<U>
    where
        F: Fn(W) -> U, {
        PairedWriters {
            writer1: op(self.writer1),
            writer2: op(self.writer2),
        }
    }

    fn try_map<U, F>(self, op: F) -> Result<Self::Paired<U>, PairedErrors>
    where
        F: Fn(W) -> std::io::Result<U>, {
        Ok(PairedWriters {
            writer1: op(self.writer1).map_err(PairedErrors::Err1)?,
            writer2: op(self.writer2).map_err(PairedErrors::Err2)?,
        })
    }
}

impl<W> PairedStruct<W> for RecordWriters<W> {
    type Paired<U> = RecordWriters<U>;

    fn map<U, F>(self, op: F) -> Self::Paired<U>
    where
        F: Fn(W) -> U, {
        match self {
            RecordWriters::SingleEnd(w) => RecordWriters::SingleEnd(op(w)),
            RecordWriters::PairedEnd(w) => RecordWriters::PairedEnd(w.map(op)),
        }
    }

    fn try_map<U, F>(self, op: F) -> Result<Self::Paired<U>, PairedErrors>
    where
        F: Fn(W) -> std::io::Result<U>, {
        match self {
            RecordWriters::SingleEnd(w) => op(w).map_err(PairedErrors::Err1).map(RecordWriters::SingleEnd),
            RecordWriters::PairedEnd(w) => w.try_map(op).map(RecordWriters::PairedEnd),
        }
    }
}
