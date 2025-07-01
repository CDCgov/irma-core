use crate::io::{FastXReader, FromFilename, ReaderFromFileError, is_gz};
use flate2::read::MultiGzDecoder;
use std::{
    fs::File,
    io::{PipeReader, Read, Stdin, stdin},
    path::Path,
    thread::{self, JoinHandle},
};
use zoe::{
    define_whichever,
    prelude::{FastQReader, FastaReader},
};

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
    #[inline]
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

        // Check for EOF
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
    // TODO: Implement reading from stdin for select processes
    /// An enum for the input types [`File`] and [`Stdin`].
    ///
    /// To construct this, either use [`from_filename`] or
    /// [`from_optional_filename`].
    ///
    /// [`from_filename`]: FromFilename::from_filename
    /// [`from_optional_filename`]: FromFilename::from_optional_filename
    pub(crate) enum ReadFileStdin {
        /// A regular uncompressed file.
        File(File),
        /// The standard input stream.
        Stdin(Stdin),
    }

    impl Read for ReadFileStdin {}
}

define_whichever! {
   /// An enum for the different input types [`File`] and a gzip compressed
   /// file.
   ///
   /// For the [`Zipped`] variant, this will perform all unzipping lazily via an
   /// iterator. To instead perform unzipping eagerly on a separate thread, use
   /// [`ReadFileZipPipe`].
   ///
   /// To construct this, use [`from_filename`]. The [`Zipped`] variant is
   /// chosen if the file has extension `gz`.
   ///
   /// [`from_filename`]: FromFilename::from_filename
   /// [`Zipped`]: ReadFileZip::Zipped
    pub(crate) enum ReadFileZip {
        /// A regular uncompressed file.
        File(File),
        /// A gzip compressed file, using lazy decoding.
        Zipped(MultiGzDecoder<File>),
    }

    impl Read for ReadFileZip {}
}

define_whichever! {
    /// An enum for the different input types [`File`] and a gzip compressed
    /// file, using a separate thread for decoding.
    ///
    /// For the [`Zipped`] variant, this will perform all unzipping eagerly on a
    /// separate thread. See [`GzipReaderPiped`] for more details. To instead
    /// perform unzipping lazily on the same thread (e.g., if the entire file
    /// will not be read, or the results will be slurped into memory), use
    /// [`ReadFileZip`].
    ///
    /// To construct this, use [`from_filename`]. The [`Zipped`] variant is
    /// chosen if the file has extension `gz`.
    ///
    /// [`from_filename`]: FromFilename::from_filename
    /// [`Zipped`]: ReadFileZip::Zipped
    pub(crate) enum ReadFileZipPipe {
        /// A regular uncompressed file.
        File(File),
        /// A gzip compressed file, using eager decoding on a separate thread.
        Zipped(GzipReaderPiped),
    }

    impl Read for ReadFileZipPipe {}
}

impl FromFilename for ReadFileStdin {
    /// Creates a new [`ReadFileStdin::File`].
    #[inline]
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        File::open(path).map(Self::File)
    }
}

impl Default for ReadFileStdin {
    /// Defaults to [`ReadFileStdin::Stdin`].
    #[inline]
    fn default() -> Self {
        Self::Stdin(stdin())
    }
}

impl FromFilename for ReadFileZip {
    /// Creates a new [`ReadFileZip`] from the file.
    ///
    /// [`ReadFileZip::File`] is returned unless the file ends with extension
    /// `gz`, in which case [`ReadFileZip::Zipped`] is returned.
    #[inline]
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        let file = File::open(&path)?;
        if is_gz(path) {
            Ok(Self::Zipped(MultiGzDecoder::new(file)))
        } else {
            Ok(Self::File(file))
        }
    }
}

impl FromFilename for ReadFileZipPipe {
    /// Creates a new [`ReadFileZipPipe`] from the file.
    ///
    /// [`ReadFileZipPipe::File`] is returned unless the file ends with
    /// extension `gz`, in which case [`ReadFileZipPipe::Zipped`] is returned.
    #[inline]
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        if is_gz(&path) {
            Ok(ReadFileZipPipe::Zipped(GzipReaderPiped::from_filename(path)?))
        } else {
            Ok(ReadFileZipPipe::File(File::open(path)?))
        }
    }
}

/// Readers for a set of possibly-paired records.
///
/// This stores a single reader, and an optional second reader for paired reads.
pub struct RecordReaders<R> {
    /// The reader for the first input file.
    pub reader1: R,
    /// The optional reader for the second input file (for paired reads).
    pub reader2: Option<R>,
}

impl<R> RecordReaders<R>
where
    R: FromFilename,
{
    /// Returns [`RecordReaders`] from a path and an optional second path.
    ///
    /// It may be necessary to specify the generic on [`RecordReaders`] before
    /// calling this function.
    pub fn from_filenames(path1: impl AsRef<Path>, path2: Option<impl AsRef<Path>>) -> std::io::Result<Self> {
        let reader1 = R::from_filename(&path1).map_err(|e| ReaderFromFileError::file1(path1, e))?;
        let reader2 = match path2 {
            Some(path2) => Some(R::from_filename(&path2).map_err(|e| ReaderFromFileError::file2(path2, e))?),
            None => None,
        };
        Ok(RecordReaders { reader1, reader2 })
    }
}

impl<R> FromFilename for FastQReader<R>
where
    R: Read + FromFilename,
{
    #[inline]
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        FastQReader::from_readable(R::from_filename(path)?)
    }
}

impl<R> FromFilename for FastaReader<R>
where
    R: Read + FromFilename,
{
    #[inline]
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        FastaReader::from_readable(R::from_filename(path)?)
    }
}

impl<R> FromFilename for FastXReader<R>
where
    R: Read + FromFilename,
{
    #[inline]
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        FastXReader::from_readable(R::from_filename(path)?)
    }
}
