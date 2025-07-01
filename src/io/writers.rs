use crate::io::{FromFilename, MapFailedWriteExt, is_gz};
use flate2::{Compression, write::GzEncoder};
use std::{
    fs::File,
    io::{BufWriter, Stdout, Write, stdout},
    path::Path,
};
use zoe::define_whichever;

define_whichever! {
    /// An enum for the different acceptable output types. A [`BufWriter`] is used for all variants.
    #[derive(Debug)]
    pub(crate) enum  WriteFileZipStdout {
        /// A writer for a regular uncompressed file.
        File(BufWriter<File>),
        /// A writer for a gzip compressed file.
        Zipped(GzEncoder<BufWriter<File>>),
        /// A writer for uncompressed data to stdout.
        Stdout(BufWriter<Stdout>),
    }

    impl Write for WriteFileZipStdout {}
}

/// A struct containing two writers for paired reads: one for left reads and one
/// for the right.
///
/// This is compatible with the [`WriteRecord`] trait, so that two-tuples and
/// length-two arrays of records can be written to [`PairedWriters`].
///
/// [`WriteRecord`]: crate::io::WriteRecord
pub struct PairedWriters<W> {
    pub writer1: W,
    pub writer2: W,
}

impl<W> PairedWriters<W> {
    /// Creates a new [`PairedWriters`] from two writers.
    #[inline]
    pub fn new(writer1: W, writer2: W) -> Self {
        Self { writer1, writer2 }
    }
}

impl<W: Write> PairedWriters<W> {
    /// Flushes the two stored writers.
    #[inline]
    pub fn flush(&mut self) -> std::io::Result<()> {
        self.writer1.flush()?;
        self.writer2.flush()
    }
}

/// An enum for holding either a single writer (single reads) or
/// [`PairedWriters`] (for paired reads).
///
/// This is compatible with the [`WriteRecords`] trait, so that given an
/// iterator of reads, they can be written to either two writer if present, or
/// written to one and interleaved. This involves a single match statement, and
/// so does not incur significant overhead. This is *not* compatible with the
/// [`WriteRecord`] trait, since performing a match on each write would be
/// inefficient.
///
/// [`WriteRecord`]: crate::io::WriteRecord
/// [`WriteRecords`]: crate::io::WriteRecords
pub enum RecordWriters<W> {
    /// A single writer for single end reads.
    SingleEnd(W),
    /// A pair of writers for paired reads.
    PairedEnd(PairedWriters<W>),
}

impl<W: Write> RecordWriters<W> {
    /// Creates a new [`RecordWriters`] object to represent either a single
    /// writer (single reads) or two writers (paired reads).
    ///
    /// This is used for parsing clap arguments.
    #[inline]
    pub fn new(writer1: W, writer2: Option<W>) -> Self {
        match writer2 {
            Some(writer2) => Self::PairedEnd(PairedWriters::new(writer1, writer2)),
            None => Self::SingleEnd(writer1),
        }
    }

    /// Creates a new [`RecordWriters`] object to represent either a single
    /// writer (single reads) or two writer (paired reads).
    ///
    /// ## Errors
    ///
    /// Any IO errors when opening the files are propagated.
    #[allow(dead_code)]
    pub fn from_filenames(path1: impl AsRef<Path>, path2: Option<impl AsRef<Path>>) -> std::io::Result<Self>
    where
        W: FromFilename, {
        Ok(Self::new(W::from_filename(path1)?, path2.map(W::from_filename).transpose()?))
    }

    /// Creates a new [`RecordWriters`] object, where the first path is optional
    /// and defaults to [`Default::default`] if not provided (often [`stdout`]).
    ///
    /// If the second path is not included, then the [`RecordWriters`] only
    /// holds a single writer, as with [`from_filenames`].
    ///
    /// [`from_filenames`]: RecordWriters::from_filenames
    #[inline]
    pub fn from_optional_filenames(
        path1: Option<impl AsRef<Path>>, path2: Option<impl AsRef<Path>>,
    ) -> std::io::Result<Self>
    where
        W: FromFilename + Default, {
        Ok(Self::new(
            W::from_optional_filename(path1)?,
            path2.map(W::from_filename).transpose()?,
        ))
    }
}

impl FromFilename for WriteFileZipStdout {
    fn from_filename<P>(path: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        let file = File::create(&path).map_failed_write(&path)?;
        let bufwriter = BufWriter::new(file);

        let writer = if is_gz(path) {
            Self::Zipped(GzEncoder::new(bufwriter, Compression::default()))
        } else {
            Self::File(bufwriter)
        };

        Ok(writer)
    }
}

impl Default for WriteFileZipStdout {
    #[inline]
    fn default() -> Self {
        Self::Stdout(BufWriter::new(stdout()))
    }
}
