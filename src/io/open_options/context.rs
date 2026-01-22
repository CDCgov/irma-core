use crate::io::{IterWithContext, IterWithErrorContext};
use std::{error::Error, fmt::Display, path::Path};
use zoe::data::err::{ErrorWithContext, WithErrorContext};

/// An enum to represent the possible reader types for [`InputOptions`], for the
/// purpose of storing context in case of an error.
///
/// [`InputOptions`]: crate::io::open_options::InputOptions
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum ReaderType {
    /// A FASTQ record, as used in [`FastQReader`] or [`FastXReader`].
    ///
    /// [`FastQReader`]: zoe::prelude::FastQReader
    /// [`FastXReader`]: crate::io::FastXReader
    FastQ,
    /// A FASTA record, as used in [`FastaReader`] or [`FastXReader`].
    ///
    /// [`FastaReader`]: zoe::prelude::FastaReader
    /// [`FastXReader`]: crate::io::FastXReader
    Fasta,
    /// A record that is either FASTQ or FASTA, used in [`FastXReader`] when the
    /// proper record type cannot be determined.
    ///
    /// [`FastXReader`]: crate::io::FastXReader
    FastX,
    /// A SAM record, as used in [`SAMReader`].
    ///
    /// [`SAMReader`]: zoe::data::records::sam::SAMReader
    Sam,
}

/// An enum to represent the possible input source for [`InputOptions`], for the
/// purpose of storing context in case of an error.
///
/// [`InputOptions`]: crate::io::InputOptions
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum InputType<'a> {
    /// The content is being read from a file with the specified path.
    ///
    /// This may also include pipes or zipped files.
    File(&'a Path),
    /// The content is being read from stdin.
    Stdin,
}

impl<'a> InputType<'a> {
    /// Creates a new [`InputType`] from an optional path. If present,
    /// [`InputType::File`] is used, otherwise [`InputType::Stdin`] is assumed.
    fn new<P>(path: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        match path {
            Some(path) => InputType::File(path.as_ref()),
            None => InputType::Stdin,
        }
    }
}

/// An enum to represent the possible output source for [`OutputOptions`], for
/// the purpose of storing context in case of an error.
///
/// [`OutputOptions`]: crate::io::OutputOptions
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum OutputType<'a> {
    /// The content is being written to a file with the specified path.
    File(&'a Path),
    /// The content is being written to stdout.
    Stdout,
}

impl<'a> OutputType<'a> {
    /// Creates a new [`OutputType`] from an optional path. If present,
    /// [`OutputType::File`] is used, otherwise [`OutputType::Stdout`] is
    /// assumed.
    fn new<P>(path: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        match path {
            Some(path) => OutputType::File(path.as_ref()),
            None => OutputType::Stdout,
        }
    }
}

/// The complete context for [`InputOptions`] necessary for displaying any
/// errors.
///
/// [`InputOptions`]: crate::io::InputOptions
pub struct InputContext<'a> {
    /// The [`ReaderType`] used by the first input, if that input is being
    /// parsed as a reader.
    pub reader1: Option<ReaderType>,
    /// The [`ReaderType`] used by the second input, if that input is being
    /// parsed as a reader. If there is no second input, this is [`None`], and
    /// this will not influence the error context displayed.
    pub reader2: Option<ReaderType>,
    /// The [`InputType`] used by the first input.
    pub input1:  InputType<'a>,
    /// The [`InputType`] used by the second input. If there is no second input,
    /// this defaults to [`InputType::Stdin`], but it will not influence the
    /// error context displayed.
    pub input2:  InputType<'a>,
}

impl<'a> InputContext<'a> {
    /// Creates a new [`InputContext`] from two optional paths.
    pub fn new<P>(path1: Option<&'a P>, path2: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        Self {
            reader1: None,
            reader2: None,
            input1:  InputType::new(path1),
            input2:  InputType::new(path2),
        }
    }

    /// Attaches a [`ReaderType`] to the context for the first input.
    pub fn with_reader1(mut self, reader: ReaderType) -> Self {
        self.reader1 = Some(reader);
        self
    }

    /// Attaches a [`ReaderType`] to the context for the second input.
    pub fn with_reader2(mut self, reader: ReaderType) -> Self {
        self.reader1 = Some(reader);
        self
    }
}

impl InputContext<'_> {
    /// Adds context to a [`PairedErrors`] from an [`InputContext`].
    ///
    /// The context will include the path if available, and the record type if a
    /// `parse` method was called and the error originated during parsing.
    pub fn add_context(&self, e: PairedErrors) -> ErrorWithContext {
        match e {
            PairedErrors::Err1(e) => Self::add_context_helper(e, self.reader1, self.input1),
            PairedErrors::Err2(e) => Self::add_context_helper(e, self.reader2, self.input2),
        }
    }

    /// Given one of the `reader` and `input` fields for an [`InputContext`],
    /// add the context to a provided [`std::io::Error`].
    ///
    /// This is a helper method for [`InputContext::add_context`].
    fn add_context_helper(e: std::io::Error, reader: Option<ReaderType>, input: InputType) -> ErrorWithContext {
        match (reader, input) {
            (None, InputType::File(path)) => e.with_file_context("Failed to open file", path),
            (None, InputType::Stdin) => e.with_context("Failed to read from stdin"),
            (Some(ReaderType::FastQ), InputType::File(path)) => {
                e.with_file_context("Failed to read FASTQ records from file", path)
            }
            (Some(ReaderType::Fasta), InputType::File(path)) => {
                e.with_file_context("Failed to read FASTA records from file", path)
            }
            (Some(ReaderType::FastX), InputType::File(path)) => {
                e.with_file_context("Failed to read records from file", path)
            }
            (Some(ReaderType::Sam), InputType::File(path)) => {
                e.with_file_context("Failed to read SAM records from file", path)
            }
            (Some(ReaderType::FastQ), InputType::Stdin) => e.with_context("Failed to read FASTQ records from stdin"),
            (Some(ReaderType::Fasta), InputType::Stdin) => e.with_context("Failed to read FASTA records from stdin"),
            (Some(ReaderType::FastX), InputType::Stdin) => e.with_context("Failed to read records from stdin"),
            (Some(ReaderType::Sam), InputType::Stdin) => e.with_context("Failed to read SAM records from stdin"),
        }
    }

    /// Adds context to the items in a fallible iterator via
    /// [`IterWithContext`].
    ///
    /// `reader` and `input` should be corresponding fields in an
    /// [`InputContext`] struct. The context will include the path if available
    /// and the record type.
    pub fn add_iter_context<I>(iter: I, reader: Option<ReaderType>, input: InputType) -> IterWithContext<I>
    where
        I: IterWithErrorContext, {
        match (reader, input) {
            (None, InputType::File(path)) => iter.iter_with_file_context("Invalid record in file", path),
            (None, InputType::Stdin) => iter.iter_with_context("Invalid record in stdin"),
            (Some(ReaderType::FastQ), InputType::File(path)) => {
                iter.iter_with_file_context("Invalid FASTQ record in file", path)
            }
            (Some(ReaderType::Fasta), InputType::File(path)) => {
                iter.iter_with_file_context("Invalid FASTA record in file", path)
            }
            (Some(ReaderType::FastX), InputType::File(path)) => iter.iter_with_file_context("Invalid record in file", path),
            (Some(ReaderType::Sam), InputType::File(path)) => {
                iter.iter_with_file_context("Invalid SAM record in file", path)
            }
            (Some(ReaderType::FastQ), InputType::Stdin) => iter.iter_with_context("Invalid FASTQ record from stdin"),
            (Some(ReaderType::Fasta), InputType::Stdin) => iter.iter_with_context("Invalid FASTA record from stdin"),
            (Some(ReaderType::FastX), InputType::Stdin) => iter.iter_with_context("Invalid record from stdin"),
            (Some(ReaderType::Sam), InputType::Stdin) => iter.iter_with_context("Invalid SAM record from stdin"),
        }
    }
}

/// The complete context for [`OutputOptions`] necessary for displaying any
/// errors.
///
/// [`OutputOptions`]: crate::io::open_options::OutputOptions
pub struct OutputContext<'a> {
    /// The [`OutputType`] used by the first output.
    output1: OutputType<'a>,
    /// The [`OutputType`] used by the second output. If there is no second
    /// output, this defaults to [`OutputType::Stdout`], but it will not
    /// influence the error context displayed.
    output2: OutputType<'a>,
}

impl<'a> OutputContext<'a> {
    /// Creates a new [`OutputContext`] from two optional paths.
    pub fn new<P>(path1: Option<&'a P>, path2: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        Self {
            output1: OutputType::new(path1),
            output2: OutputType::new(path2),
        }
    }
}

impl OutputContext<'_> {
    /// Adds context to a [`PairedErrors`] from an [`OutputContext`].
    pub fn add_context(&self, e: PairedErrors) -> ErrorWithContext {
        match e {
            PairedErrors::Err1(e) => Self::add_context_helper(e, self.output1),
            PairedErrors::Err2(e) => Self::add_context_helper(e, self.output2),
        }
    }

    /// Given one of the `output` fields for an [`OutputContext`], add the
    /// context to a provided [`std::io::Error`].
    ///
    /// This is a helper method for [`OutputContext::add_context`].
    fn add_context_helper(e: std::io::Error, output: OutputType) -> ErrorWithContext {
        match output {
            OutputType::File(path) => e.with_file_context("Failed to open file for writing", path),
            OutputType::Stdout => e.with_context("Failed to write to stdout"),
        }
    }
}

/// An error occurring while working with potentially-paired inputs or outputs,
/// with two variants to hold whether the error occurred for the first
/// input/output or the second.
///
/// If paired inputs/outputs are not being used, then [`PairedErrors::Err1`] is
/// used.
///
/// This enum does not alter the [`Display`] or [`Error::source`]
/// implementations, and is instead used to facilitate adding proper context at
/// the end of [`InputOptions`] or [`OutputOptions`] when an `open` method is
/// called.
///
/// [`InputOptions`]: crate::io::open_options::InputOptions
/// [`OutputOptions`]: crate::io::open_options::OutputOptions
#[derive(Debug)]
pub enum PairedErrors {
    /// An error that occurred while working with the first input/output (or an
    /// unpaired input/output).
    Err1(std::io::Error),
    /// An error that occurred while working with the second input/output.
    Err2(std::io::Error),
}

impl Display for PairedErrors {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PairedErrors::Err1(e) => write!(f, "{e}"),
            PairedErrors::Err2(e) => write!(f, "{e}"),
        }
    }
}

impl Error for PairedErrors {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            PairedErrors::Err1(e) => Some(e),
            PairedErrors::Err2(e) => Some(e),
        }
    }
}
