use std::path::Path;

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
