use crate::io::is_gz;
use flate2::read::MultiGzDecoder;
use std::{
    fs::File,
    io::{PipeReader, Read, Stdin, stdin},
    path::Path,
    thread::{self, JoinHandle},
};
use zoe::{data::err::ResultWithErrorContext, define_whichever};

/// A reader for a [gzip file](https://www.rfc-editor.org/rfc/rfc1952#page-5)
/// that may have multiple members, spawning a separate thread for the unzipping
/// process.
///
/// This uses a separate thread for unzipping to minimize IO overhead, and it
/// uses an anonymous pipe to communicate the unzipped data.
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
    ///
    /// ## Errors
    ///
    /// Any IO errors occurring when forming the pipe are propagated with
    /// context. Errors occurring during decoding appear when reading from the
    /// [`GzipReaderPiped`].
    pub fn from_readable<R>(readable: R) -> std::io::Result<Self>
    where
        R: Read + Send + 'static, {
        let (reader, mut writer) =
            std::io::pipe().with_context("Failed to intiailize the pipe for decoding the gzip data")?;

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

    /// Opens a new [`GzipReaderPiped`] from a path.
    ///
    /// The path should contain
    /// [gzip](https://www.rfc-editor.org/rfc/rfc1952#page-5) encoded data.
    ///
    /// ## Errors
    ///
    /// Any IO errors occurring when opening the file or forming the pipe are
    /// propagated.
    #[inline]
    #[allow(dead_code)]
    pub fn open<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        File::open(&path).and_then(Self::from_readable)
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
    pub(crate) enum ReadFileStdin {
        /// A regular uncompressed file.
        File(File),
        /// The standard input stream.
        Stdin(Stdin),
    }

    impl Read for ReadFileStdin {}
}

impl ReadFileStdin {
    /// Opens a [`ReadFileStdin`] from an optional path. If a path is not
    /// provided, [`ReadFileStdin::Stdin`] is used.
    ///
    /// ## Errors
    ///
    /// If a path is provided, any IO errors when opening the file are
    /// propagated. If no path is provided, this method is infallible.
    pub fn open(path: Option<impl AsRef<Path>>) -> std::io::Result<Self> {
        match path {
            Some(path) => File::open(&path).map(Self::File),
            None => Ok(ReadFileStdin::Stdin(stdin())),
        }
    }
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

impl ReadFileZip {
    /// Opens a [`ReadFileZip`] from a path.
    ///
    /// The file is determined to be zipped if it ends in `.gz`.
    ///
    /// ## Errors
    ///
    /// Any IO errors when opening the file are propagated.
    pub fn open(path: impl AsRef<Path>) -> std::io::Result<Self> {
        let file = File::open(&path)?;

        if is_gz(path) {
            Ok(Self::Zipped(MultiGzDecoder::new(file)))
        } else {
            Ok(Self::File(file))
        }
    }
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

impl ReadFileZipPipe {
    /// Opens a [`ReadFileZipPipe`] from a path.
    ///
    /// The file is determined to be zipped if it ends in `.gz`.
    ///
    /// ## Errors
    ///
    /// Any IO errors when opening the file or forming the pipe are propagated.
    pub fn open(path: impl AsRef<Path>) -> std::io::Result<Self> {
        let file = File::open(&path)?;

        if is_gz(&path) {
            Ok(ReadFileZipPipe::Zipped(GzipReaderPiped::from_readable(file)?))
        } else {
            Ok(ReadFileZipPipe::File(file))
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
