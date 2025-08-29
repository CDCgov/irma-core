use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use std::{
    error::Error,
    fs::File,
    io::{BufWriter, ErrorKind, PipeReader, Read, Stdin, Stdout, Write, stdout},
    path::{Path, PathBuf},
    thread::{self, JoinHandle},
};
use zoe::{define_whichever, prelude::FastQReader};

mod fastx;
mod write_records;

pub use fastx::*;
pub use write_records::*;

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
    // TODO: Implement reading from stdin for select processes
    #[allow(dead_code)]
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
/// Checks whether the file is a gz zipped file.
///
/// IRMA-core's current strategy for this is by checking whether the extension
/// is `.gz`.
pub(crate) fn is_gz<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref().extension().is_some_and(|ext| ext == "gz")
}

/// Open a single FASTQ file.
///
/// If it ends in `gz`, return a [`FastQReader`] backed by
/// [`ReadFileZip::Zipped`], and return the thread handle for error propagation.
/// Otherwise, return a [`FastQReader`] backed by [`ReadFileZip::File`].
#[inline]
pub(crate) fn open_fastq_file<P: AsRef<Path>>(path: P) -> std::io::Result<(FastQReaderIc, Option<IoThread>)> {
    let file = File::open(&path)?;

    if is_gz(&path) {
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
/// Open a single FASTQ/FASTA file.
///
/// If it ends in `gz`, return a [`FastQReader`] backed by
/// [`ReadFileZip::Zipped`], and return the thread handle for error propagation.
/// Otherwise, return a [`FastQReader`] backed by [`ReadFileZip::File`].
#[inline]
pub(crate) fn open_fastx_file<P: AsRef<Path>>(path: P) -> std::io::Result<(FastXReader<ReadFileZip>, Option<IoThread>)> {
    let file = File::open(&path)?;

    if is_gz(&path) {
        let (pipe, thread) = spawn_decoder(path)?;
        Ok((FastXReader::from_readable(ReadFileZip::Zipped(pipe))?, Some(thread)))
    } else {
        Ok((FastXReader::from_readable(ReadFileZip::File(file))?, None))
    }
}

/// Open one or two FASTQ files using the strategy given by [`open_fastq_file`].
/// The thread handles are grouped together in [`IoThreads`].
#[inline]
pub(crate) fn open_fastq_files<P: AsRef<Path>>(
    path1: P, path2: Option<P>,
) -> Result<(FastQReaderIc, Option<FastQReaderIc>, IoThreads), OpenFastqError> {
    let Some(path2) = path2 else {
        let (reader, thread) = open_fastq_file(&path1).map_err(|e| OpenFastqError::file1(&path1, e))?;
        let threads = IoThreads(thread, None);

        return Ok((reader, None, threads));
    };

    let (reader1, thread1) = open_fastq_file(&path1).map_err(|e| OpenFastqError::file1(&path1, e))?;
    let (reader2, thread2) = open_fastq_file(&path2).map_err(|e| OpenFastqError::file2(&path2, e))?;
    let threads = IoThreads(thread1, thread2);
    Ok((reader1, Some(reader2), threads))
}

/// Creates a [`WriteFileZipStdout`], using `path` to determine whether a
/// regular file, zipped file, or stdout should be used.
///
/// ## Errors
///
/// Creation of `path` must be successful, if a path is specified.
/// Open one or two FASTQ/FASTA files using the strategy given by
/// [`open_fastx_file`]. The thread handles are grouped together in
/// [`IoThreads`].
#[inline]
pub(crate) fn open_fastx_files<P: AsRef<Path>>(
    path1: P, path2: Option<P>,
) -> std::io::Result<(FastXReader<ReadFileZip>, Option<FastXReader<ReadFileZip>>, IoThreads)> {
    let Some(path2) = path2 else {
        let (reader, thread) = open_fastx_file(path1)?;
        let threads = IoThreads(thread, None);

        return Ok((reader, None, threads));
    };

    let (reader1, thread1) = open_fastx_file(path1)?;
    let (reader2, thread2) = open_fastx_file(path2)?;
    let threads = IoThreads(thread1, thread2);
    Ok((reader1, Some(reader2), threads))
}

/// Creates a [`WriteFileZipStdout`], using `path` to determine whether a regular
/// file, zipped file, or stdout should be used.
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
        match std::io::copy(&mut decoder, &mut writer) {
            // A broken pipe signals that the reader has been dropped, so the
            // writer can end early (e.g., in sampler). `MultiGzDecoder` will
            // never yield such an error since it is reading from a file.
            Err(e) if e.kind() == ErrorKind::BrokenPipe => Ok(()),
            Ok(_) => Ok(()),
            Err(e) => Err(e),
        }
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
    File1 { path: PathBuf, source: std::io::Error },
    File2 { path: PathBuf, source: std::io::Error },
}

impl std::fmt::Display for OpenFastqError {
    #[inline]
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            OpenFastqError::File1 { path, source } => {
                write!(
                    f,
                    "Failed to read the data in first file {path:#?} due to the error:\n{source}"
                )
            }
            OpenFastqError::File2 { path, source } => {
                write!(
                    f,
                    "Failed to read the data in second file {path:#?} due to the error:\n{source}"
                )
            }
        }
    }
}
/// A struct containing two writers for paired reads: one for left reads and one
/// for the right.
pub struct PairedWriters<W: Write> {
    writer1: W,
    writer2: W,
}

impl<W: Write> PairedWriters<W> {
    /// Creates a new [`PairedWriters`] from two writers.
    #[inline]
    fn new(writer1: W, writer2: W) -> Self {
        Self { writer1, writer2 }
    }
}

/// A trait providing the ability to flush all writers in a struct.
pub trait FlushWriter {
    /// Flushes all the writers.
    fn flush_all(&mut self) -> std::io::Result<()>;
}

impl<W: Write> FlushWriter for W {
    #[inline]
    fn flush_all(&mut self) -> std::io::Result<()> {
        self.flush()
    }
}

impl<W: Write> FlushWriter for PairedWriters<W> {
    #[inline]
    fn flush_all(&mut self) -> std::io::Result<()> {
        self.writer1.flush()?;
        self.writer2.flush()
    }
}

/// An enum for holding either a single writer (single reads) or
/// [`PairedWriters`] (for paired reads).
pub enum RecordWriters<W: Write> {
    SingleEnd(W),
    PairedEnd(PairedWriters<W>),
}

impl<W: Write> RecordWriters<W> {
    /// Creates a new [`RecordWriters`] object to represent either a single
    /// writer (single reads) or two writers (paired reads). This is used for
    /// parsing clap arguments.
    #[inline]
    pub fn new(writer1: W, writer2: Option<W>) -> Self {
        match writer2 {
            Some(writer2) => Self::PairedEnd(PairedWriters::new(writer1, writer2)),
            None => Self::SingleEnd(writer1),
        }
    }
}

impl RecordWriters<WriteFileZipStdout> {
    /// Creates a [`RecordWriters`] from two files, interpretted as
    /// [`WriteFileZipStdout`] writers.
    pub fn from_filename<P: AsRef<Path>>(writer1: Option<P>, writer2: Option<P>) -> std::io::Result<Self> {
        let writer1 = create_writer(writer1)?;
        let writer2 = match writer2 {
            Some(path) => Some(create_writer(Some(path))?),
            None => None,
        };
        Ok(Self::new(writer1, writer2))
    }
}

impl Error for OpenFastqError {
    #[inline]
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            OpenFastqError::File1 { source, .. } => Some(source),
            OpenFastqError::File2 { source, .. } => Some(source),
        }
    }
}

impl From<OpenFastqError> for std::io::Error {
    fn from(error: OpenFastqError) -> Self {
        std::io::Error::other(error)
    }
}

impl OpenFastqError {
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
