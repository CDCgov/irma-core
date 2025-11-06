use std::{
    error::Error,
    path::{Path, PathBuf},
};

mod fastx;
mod readers;
mod write_records;
mod writers;

pub use fastx::*;
pub(crate) use readers::*;
pub use write_records::*;
pub use writers::*;

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

/// A trait unifying readers/writers that can be created from filenames.
///
/// This allows for composibility, such as an implementation of [`FromFilename`]
/// on `FastQReader<ReadFileZip>`.
pub(crate) trait FromFilename
where
    Self: Sized, {
    /// Creates the reader/writer from the path.
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>;

    /// Creates the reader/writer from an optional path, using [`default`] if it
    /// is not provided.
    ///
    /// [`default`]: Default::default
    fn from_optional_filename<P>(path: Option<P>) -> std::io::Result<Self>
    where
        Self: Default,
        P: AsRef<Path>, {
        match path {
            Some(path) => Self::from_filename(path),
            None => Ok(Self::default()),
        }
    }
}

/// A wrapper around [`std::io::Error`] used to indicate whether an error
/// occurred with the first file or second file specified.
#[non_exhaustive]
#[derive(Debug)]
pub enum ReaderFromFileError {
    File1 { path: PathBuf, source: std::io::Error },
    File2 { path: PathBuf, source: std::io::Error },
}

impl std::fmt::Display for ReaderFromFileError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            ReaderFromFileError::File1 { path, source } => {
                write!(
                    f,
                    "Failed to read the data in first file {path:#?} due to the error:\n{source}"
                )
            }
            ReaderFromFileError::File2 { path, source } => {
                write!(
                    f,
                    "Failed to read the data in second file {path:#?} due to the error:\n{source}"
                )
            }
        }
    }
}

impl Error for ReaderFromFileError {
    #[inline]
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ReaderFromFileError::File1 { source, .. } => Some(source),
            ReaderFromFileError::File2 { source, .. } => Some(source),
        }
    }
}

impl From<ReaderFromFileError> for std::io::Error {
    #[inline]
    fn from(error: ReaderFromFileError) -> Self {
        std::io::Error::other(error)
    }
}

impl ReaderFromFileError {
    #[inline]
    fn file1(path: impl AsRef<Path>, source: std::io::Error) -> Self {
        Self::File1 {
            path: path.as_ref().to_path_buf(),
            source,
        }
    }

    #[inline]
    fn file2(path: impl AsRef<Path>, source: std::io::Error) -> Self {
        Self::File2 {
            path: path.as_ref().to_path_buf(),
            source,
        }
    }
}

pub trait MapFailedWriteExt<T> {
    fn map_failed_write<P: AsRef<Path>>(self, path: P) -> std::io::Result<T>;
}

impl<T> MapFailedWriteExt<T> for std::io::Result<T> {
    fn map_failed_write<P: AsRef<Path>>(self, path: P) -> std::io::Result<T> {
        self.map_err(|e| {
            std::io::Error::other(format!(
                "Failed to open {path} for writing due to the error:\n{e}",
                path = path.as_ref().display()
            ))
        })
    }
}
