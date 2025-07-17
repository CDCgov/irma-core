use crate::utils::whichever::define_whichever;
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use std::{
    error::Error,
    fs::File,
    io::{BufWriter, PipeReader, Stdin, Stdout, stdout},
    path::{Path, PathBuf},
    thread::{self, JoinHandle},
};
use zoe::{data::err::GetCode, prelude::FastQReader};

define_whichever! {
    #[allow(clippy::large_enum_variant)]
    #[doc="An enum for the different acceptable input types"]
    pub(crate) enum ReadFileZip {
        #[doc="A reader for a regular uncompressed file"]
        File(File),
        #[doc="A reader for a gzip compressed file, using a thread and an anonymous pipe for decoding"]
        Zipped(PipeReader),
    }

    impl Read for ReadFileZip {}
}

/// A type alias for the `FastQReader` used by IRMA-core's trimmer and
/// preprocess
type FastQReaderIc = FastQReader<ReadFileZip>;

define_whichever! {
    #[doc="An enum for the different acceptable output types"]
    #[derive(Debug)]
    pub(crate) enum  WriteFileZipStdout {
        #[doc="A writer for a regular uncompressed file"]
        File(BufWriter<File>),
        #[doc="A writer for a gzip compressed file"]
        Zipped(GzEncoder<BufWriter<File>>),
        #[doc="A writer for uncompressed data to stdout"]
        Stdout(BufWriter<Stdout>),
    }

    impl Write for WriteFileZipStdout {}
}

define_whichever! {
    pub(crate) enum ReadFileStdin {
        File(File),
        Stdin(Stdin),
    }

    impl Read for ReadFileStdin {}
}

/// Opens a single FASTQ file.
///
/// If the filename ends in `gz`, a thread is spawned with [`spawn_decoder`] to
/// decode the input. The decoded lines are sent via a pipe to a
/// [`FastQReader`]. The second return value is the handle to the thread.
///
/// If the filename does not end in `gz`, the [`FastQReader`] is backed directly
/// by the file, and the [`IoThread`] return is `None`.
///
/// ## Errors
///
/// `path` must exist and contain FASTQ data, and if the file is zipped, then
/// creation of the pipe must succeed.
#[inline]
pub(crate) fn open_fastq_file<P: AsRef<Path>>(path: P) -> std::io::Result<(FastQReaderIc, Option<IoThread>)> {
    let file = File::open(&path)?;

    let is_gz = path.as_ref().extension().is_some_and(|ext| ext == "gz");

    if is_gz {
        let (pipe, thread) = spawn_decoder(path)?;
        Ok((FastQReader::from_readable(ReadFileZip::Zipped(pipe))?, Some(thread)))
    } else {
        Ok((FastQReader::from_readable(ReadFileZip::File(file))?, None))
    }
}

/// Open one or two FASTQ files. The type of data backing the [`FastQReader`]
/// depends on the output of [`open_fastq_file`]. The thread handles are grouped
/// together in [`IoThreads`], if either file was zipped.
///
/// ## Errors
///
/// `path1` and `path2` must both exist and contain FASTQ data (if not `None`).
/// Also, if either is zipped, the creation of the pipe must succeed.
#[inline]
pub(crate) fn open_fastq_files<P: AsRef<Path>>(
    path1: P, path2: Option<P>,
) -> Result<(FastQReaderIc, Option<FastQReaderIc>, IoThreads), OpenFastqError> {
    let Some(path2) = path2 else {
        let (reader, thread) = open_fastq_file(path1).map_err(OpenFastqError::File1)?;
        let threads = IoThreads(thread, None);

        return Ok((reader, None, threads));
    };

    let (reader1, thread1) = open_fastq_file(path1).map_err(OpenFastqError::File1)?;
    let (reader2, thread2) = open_fastq_file(path2).map_err(OpenFastqError::File2)?;
    let threads = IoThreads(thread1, thread2);
    Ok((reader1, Some(reader2), threads))
}

/// Creates a [`WriteFileZipStdout`], using `path` to determine whether a
/// regular file, zipped file, or stdout should be used.
///
/// ## Errors
///
/// Creation of `path` must be successful, if a path is specified.
#[inline]
pub(crate) fn create_writer<P: AsRef<Path>>(path: Option<P>) -> std::io::Result<WriteFileZipStdout> {
    let writer = match path {
        Some(ref p) => {
            let is_gz = p.as_ref().extension().is_some_and(|ext| ext == "gz");
            let file = File::create(p)?;
            let buf_writer = BufWriter::new(file);

            if is_gz {
                WriteFileZipStdout::Zipped(GzEncoder::new(buf_writer, Compression::default()))
            } else {
                WriteFileZipStdout::File(buf_writer)
            }
        }
        None => WriteFileZipStdout::Stdout(BufWriter::new(stdout())),
    };

    Ok(writer)
}

/// Spawns a thread that decodes the input file using [`MultiGzDecoder`].
/// Returns a [`PipeReader`] for receiving the decoded data and an [`IoThread`]
/// handle for handling the thread and propagating errors.
///
/// ## Errors
///
/// `file_path` must exist, and the creation of the pipe must succeed.
#[inline]
fn spawn_decoder(file_path: impl AsRef<Path>) -> std::io::Result<(std::io::PipeReader, IoThread)> {
    let (reader, mut writer) = std::io::pipe()?;

    let mut decoder = MultiGzDecoder::new(File::open(file_path)?);

    let thread = thread::spawn(move || -> std::io::Result<_> {
        std::io::copy(&mut decoder, &mut writer)?;
        Ok(())
    });

    Ok((reader, thread))
}

/// The handle for a thread used for IO.
type IoThread = JoinHandle<std::io::Result<()>>;

/// A struct holding two optional [`IoThread`] handles.
pub(crate) struct IoThreads(Option<IoThread>, Option<IoThread>);

impl IoThreads {
    /// Calls `join` on the underlying threads and propagate any errors.
    #[inline]
    pub(crate) fn finalize(self) -> std::io::Result<()> {
        if let Some(thread1) = self.0 {
            thread1.join().unwrap()?;
        };
        if let Some(thread2) = self.1 {
            thread2.join().unwrap()?;
        };
        Ok(())
    }
}

/// A wrapper around [`std::io::Error`] used to indicate whether an error
/// occurred with the first file or second file specified.
#[non_exhaustive]
#[derive(Debug)]
pub enum OpenFastqError {
    File1(std::io::Error),
    File2(std::io::Error),
}

impl std::fmt::Display for OpenFastqError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            OpenFastqError::File1(error) => write!(f, "{error}"),
            OpenFastqError::File2(error) => write!(f, "{error}"),
        }
    }
}

impl Error for OpenFastqError {}
impl GetCode for OpenFastqError {}

pub trait MapFailedOpenExt<T> {
    fn map_failed_open(self, path1: &Path, path2: Option<&PathBuf>) -> std::io::Result<T>;
}

impl<T> MapFailedOpenExt<T> for Result<T, OpenFastqError> {
    fn map_failed_open(self, path1: &Path, path2: Option<&PathBuf>) -> std::io::Result<T> {
        self.map_err(|e| match e {
            OpenFastqError::File1(error) => std::io::Error::other(format!(
                "Failed to read the data in file {path1:#?} due to the error:\n{error}"
            )),
            OpenFastqError::File2(error) => std::io::Error::other(format!(
                "Failed to read the data in file {path:#?} due to the error:\n{error}",
                path = path2.unwrap()
            )),
        })
    }
}

pub trait MapFailedWriteExt<T> {
    fn map_failed_write(self, path: Option<&PathBuf>) -> std::io::Result<T>;
}

impl<T> MapFailedWriteExt<T> for std::io::Result<T> {
    fn map_failed_write(self, path: Option<&PathBuf>) -> std::io::Result<T> {
        self.map_err(|e| {
            std::io::Error::other(format!(
                "Failed to open {path:#?} for writing due to the error:\n{e}",
                path = path.unwrap()
            ))
        })
    }
}
