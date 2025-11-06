use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use std::{
    error::Error,
    fs::File,
    io::{BufWriter, PipeReader, Read, Stdin, Stdout, Write, stdout},
    path::{Path, PathBuf},
    thread::{self, JoinHandle},
};
use zoe::{define_whichever, prelude::FastQReader};

mod fastx;
mod write_records;

pub use fastx::*;
pub use write_records::*;

/// A reader for a [gzip file](https://www.rfc-editor.org/rfc/rfc1952#page-5)
/// that may have multiple members, spawning a separate thread for the unzipping
/// process.
///
/// A separate thread for unzipping to minimize IO overhead, and uses an
/// anonymous pipe to communicate the unzipped data.
///
/// ## Limitations
///
/// This is designed for scenarios where the file is read in its entirety. Any
/// pipe failures are yielded only when EOF is reached (otherwise they may be
/// ignored). Furthermore, dropping this reader may not necessarily terminate
/// the thread.
pub struct GzipReaderPiped {
    reader: PipeReader,
    thread: Option<JoinHandle<std::io::Result<()>>>,
}

impl GzipReaderPiped {
    /// Creates a new [`GzipReaderPiped`] from a type implementing [`Read`].
    ///
    /// `readable` should contain
    /// [gzip](https://www.rfc-editor.org/rfc/rfc1952#page-5) encoded data.
    pub fn from_readable<R>(readable: R) -> std::io::Result<Self>
    where
        R: Read + Send + 'static, {
        let (reader, mut writer) = std::io::pipe()?;

        let mut decoder = MultiGzDecoder::new(readable);

        let thread = thread::spawn(move || -> std::io::Result<_> {
            // This thread may throw a broken pipe error if the reader is
            // dropped early, but in that case, the thread will never be joined
            // either
            std::io::copy(&mut decoder, &mut writer)?;
            Ok(())
        });

        Ok(Self {
            reader,
            thread: Some(thread),
        })
    }

    /// Creates a new [`GzipReaderPiped`] over a file.
    ///
    /// The file should contain
    /// [gzip](https://www.rfc-editor.org/rfc/rfc1952#page-5) encoded data.
    pub fn from_filename<P>(filename: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        let file = File::open(filename)?;
        Self::from_readable(file)
    }
}

impl Read for GzipReaderPiped {
    #[inline]
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let bytes_read = self.reader.read(buf)?;

        // Check if at EOF
        if bytes_read == 0
            && !buf.is_empty()
            && let Some(thread) = std::mem::take(&mut self.thread)
        {
            thread.join().unwrap()?;
        }
        Ok(bytes_read)
    }
}

define_whichever! {
    #[allow(clippy::large_enum_variant)]
    #[doc="An enum for the different acceptable input types"]
    pub(crate) enum ReadFileZip {
        #[doc="A reader for a regular uncompressed file"]
        File(File),
        #[doc="A reader for a gzip compressed file, using a thread and an anonymous pipe for decoding"]
        Zipped(GzipReaderPiped),
    }

    impl Read for ReadFileZip {}
}

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

/// Uses [`readers_from_files`] and [`RecordWriters::from_filename`] to generate
/// readers from the input paths and writers for the output paths.
///
/// ## Errors
///
/// In addition to any IO errors returned by the above functions, if any of the
/// input paths equal any of the output paths, then an error is returned.
#[inline]
pub(crate) fn get_paired_readers_and_writers<R: RecordReader<ReadFileZip>>(
    input1: impl AsRef<Path>, input2: Option<impl AsRef<Path>>, output1: Option<impl AsRef<Path>>,
    output2: Option<impl AsRef<Path>>,
) -> std::io::Result<(R, Option<R>, RecordWriters<WriteFileZipStdout>)> {
    fn inner<R: RecordReader<ReadFileZip>>(
        input1: &Path, input2: Option<&Path>, output1: Option<&Path>, output2: Option<&Path>,
    ) -> std::io::Result<(R, Option<R>, RecordWriters<WriteFileZipStdout>)> {
        if let Some(output1) = output1 {
            if input1 == output1 {
                return Err(std::io::Error::other(
                    "The first input file is the same as the first output file",
                ));
            } else if let Some(input2) = input2
                && input2 == output1
            {
                return Err(std::io::Error::other(
                    "The second input file is the same as the first output file",
                ));
            }
        }

        if let Some(output2) = output2 {
            if input1 == output2 {
                return Err(std::io::Error::other(
                    "The first input file is the same as the second output file",
                ));
            } else if let Some(input2) = input2
                && input2 == output2
            {
                return Err(std::io::Error::other(
                    "The second input file is the same as the second output file",
                ));
            }
        }

        let (reader1, reader2) = readers_from_files::<R, _>(input1, input2).map_err(Into::<std::io::Error>::into)?;
        let writer = RecordWriters::from_filename(output1, output2)?;
        Ok((reader1, reader2, writer))
    }

    inner(
        input1.as_ref(),
        input2.as_ref().map(AsRef::as_ref),
        output1.as_ref().map(AsRef::as_ref),
        output2.as_ref().map(AsRef::as_ref),
    )
}

/// Checks whether a file is a [gzip
/// file](https://www.rfc-editor.org/rfc/rfc1952#page-5).
///
/// This is currently done naively by seeing if it ends with a `gz` extension.
#[inline]
pub(crate) fn is_gz<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref().extension().is_some_and(|ext| ext == "gz")
}

/// Creates a reader from a single file.
///
/// If it ends in `gz`, return a reader backed by [`ReadFileZip::Zipped`] (which
/// uses a thread for unzipping). Otherwise, returns a reader backed by
/// [`ReadFileZip::File`].
#[inline]
pub(crate) fn reader_from_file<R: RecordReader<ReadFileZip>, P: AsRef<Path>>(path: P) -> std::io::Result<R> {
    let file = File::open(&path)?;

    if is_gz(&path) {
        R::from_readable(ReadFileZip::Zipped(GzipReaderPiped::from_filename(path)?))
    } else {
        R::from_readable(ReadFileZip::File(file))
    }
}

/// Creates one or two readers from one or two files using the strategy given by
/// [`reader_from_file`].
#[inline]
pub(crate) fn readers_from_files<R: RecordReader<ReadFileZip>, P: AsRef<Path>>(
    path1: P, path2: Option<P>,
) -> Result<(R, Option<R>), ReaderFromFileError> {
    let reader1 = reader_from_file(&path1).map_err(|e| ReaderFromFileError::file1(&path1, e))?;
    let reader2 = match path2 {
        Some(path2) => Some(reader_from_file(&path2).map_err(|e| ReaderFromFileError::file2(&path2, e))?),
        None => None,
    };
    Ok((reader1, reader2))
}

/// Creates a [`WriteFileZipStdout`], using `path` to determine whether a
/// regular file, zipped file, or stdout should be used.
///
/// The path is included in the error if the file cannot be created.
#[inline]
pub(crate) fn create_writer<P: AsRef<Path>>(path: Option<P>) -> std::io::Result<WriteFileZipStdout> {
    let writer = match path {
        Some(ref p) => {
            let is_gz = p.as_ref().extension().is_some_and(|ext| ext == "gz");
            let file = File::create(p).map_failed_write(p.as_ref())?;
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
/// A struct containing two writers for paired reads: one for left reads and one
/// for the right.
pub struct PairedWriters<W: Write> {
    pub(crate) writer1: W,
    pub(crate) writer2: W,
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
        let writer1 = create_writer(writer1.as_ref())?;
        let writer2 = match writer2 {
            Some(path) => Some(create_writer(Some(path.as_ref()))?),
            None => None,
        };
        Ok(Self::new(writer1, writer2))
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

/// A trait unifying [`FastQReader`] and [`FastXReader`].
pub trait RecordReader<R: Read>: Sized {
    /// Creates a reader from `read`.
    fn from_readable(read: R) -> std::io::Result<Self>;
}

impl<R: Read> RecordReader<R> for FastQReader<R> {
    #[inline]
    fn from_readable(read: R) -> std::io::Result<Self> {
        FastQReader::from_readable(read)
    }
}

impl<R: Read> RecordReader<R> for FastXReader<R> {
    #[inline]
    fn from_readable(read: R) -> std::io::Result<Self> {
        FastXReader::from_readable(read)
    }
}
