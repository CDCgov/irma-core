use crate::utils::whichever::define_whichever;
use flate2::{Compression, read::MultiGzDecoder, write::GzEncoder};
use std::{
    fs::File,
    io::{BufWriter, Error as IOError, PipeReader, Stdin, Stdout, stdout},
    path::Path,
    thread::{self, JoinHandle},
};
use zoe::prelude::FastQReader;

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

/// Open a single FASTQ file. If it ends in `gz`, return a [`FastQReader`]
/// backed by [`ReadFileZip::Zipped`], and return the thread handle for error
/// propagation. Otherwise, return a [`FastQReader`] backed by
/// [`ReadFileZip::File`].
#[inline]
pub(crate) fn open_fastq_file<P: AsRef<Path>>(path: P) -> std::io::Result<(FastQReader<ReadFileZip>, Option<IoThread>)> {
    let file = File::open(&path)?;

    let is_gz = path.as_ref().extension().is_some_and(|ext| ext == "gz");

    if is_gz {
        let (pipe, thread) = spawn_decoder(path)?;
        Ok((FastQReader::from_readable(ReadFileZip::Zipped(pipe))?, Some(thread)))
    } else {
        Ok((FastQReader::from_readable(ReadFileZip::File(file))?, None))
    }
}

/// Open one or two FASTQ files using the strategy given by [`open_fastq_file`].
/// The thread handles are grouped together in [`IoThreads`].
#[inline]
pub(crate) fn open_fastq_files<P: AsRef<Path>>(
    path1: P, path2: Option<P>,
) -> std::io::Result<(FastQReader<ReadFileZip>, Option<FastQReader<ReadFileZip>>, IoThreads)> {
    let Some(path2) = path2 else {
        let (reader, thread) = open_fastq_file(path1)?;
        let threads = IoThreads(thread, None);

        return Ok((reader, None, threads));
    };

    let (reader1, thread1) = open_fastq_file(path1)?;
    let (reader2, thread2) = open_fastq_file(path2)?;
    let threads = IoThreads(thread1, thread2);
    Ok((reader1, Some(reader2), threads))
}

/// Creates a [`WriteFileZipStdout`], using `path` to determine whether a regular
/// file, zipped file, or stdout should be used.
#[inline]
pub(crate) fn create_writer<P: AsRef<Path>>(path: Option<P>) -> Result<WriteFileZipStdout, IOError> {
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
