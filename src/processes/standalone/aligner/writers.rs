//! Traits and structs for writing the output of aligner

use crate::aligner::{AlignerConfig, AlignmentAndSeqs, Strand};
use std::io::Write;
use zoe::{
    data::{fasta::FastaSeq, sam::SamDataView},
    math::AnyInt,
    prelude::{DataOwned, NucleotidesView, QualityScores, QualityScoresView},
};

#[cfg(feature = "dev_no_rayon")]
use crate::io::WriteFileZipStdout;

#[cfg(not(feature = "dev_no_rayon"))]
use std::{error::Error, fmt::Display};
#[cfg(not(feature = "dev_no_rayon"))]
use zoe::data::err::{ErrorWithContext, GetCode};

/// The type of error yielded by a failed write. This depends on whether
/// `dev-no-rayon` is set.
#[cfg(not(feature = "dev_no_rayon"))]
pub type WriterError = ThreadedWriteError;

/// The type of error yielded by a failed write. This depends on whether
/// `dev-no-rayon` is set.
#[cfg(feature = "dev_no_rayon")]
pub type WriterError = std::io::Error;

/// A clonable writer supporting writing from multiple threads via an [`mpsc`]
/// channel.
///
/// A single dedicated thread is used for writing to the file to avoid
/// interleaved writes. The handle to this thread is stored in the first
/// constructed [`AlignmentWriterThreaded`], and all subsequently cloned copies
/// do not hold the thread handle. It is important to call [`flush`] on the
/// original writer to properly finalize the thread.
///
/// [`flush`]: AlignmentWriterThreaded::flush
/// [`mpsc`]: std::sync::mpsc
#[cfg(not(feature = "dev_no_rayon"))]
pub struct AlignmentWriterThreaded {
    sender:        std::sync::mpsc::Sender<String>,
    writer_thread: Option<std::thread::JoinHandle<std::io::Result<()>>>,
}

#[cfg(not(feature = "dev_no_rayon"))]
impl Clone for AlignmentWriterThreaded {
    #[inline]
    fn clone(&self) -> Self {
        Self {
            sender:        self.sender.clone(),
            writer_thread: None,
        }
    }
}

/// An error that could arise when forming and writing an alignment in the
/// multithreaded case.
#[derive(Debug)]
#[cfg(not(feature = "dev_no_rayon"))]
pub enum ThreadedWriteError {
    /// An IO error, which is used to hold any concrete error whose cause is
    /// known.
    IoError(std::io::Error),
    /// An error due to the writer thread being unavailable for receiving
    /// messages. This signals that an error occurred, but the cause of the
    /// error is unknown.
    ReceiverDeallocated,
}

#[cfg(not(feature = "dev_no_rayon"))]
impl From<std::io::Error> for ThreadedWriteError {
    fn from(value: std::io::Error) -> Self {
        Self::IoError(value)
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl From<ErrorWithContext> for ThreadedWriteError {
    fn from(value: ErrorWithContext) -> Self {
        Self::IoError(value.into())
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl Display for ThreadedWriteError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ThreadedWriteError::IoError(e) => write!(f, "{e}"),
            ThreadedWriteError::ReceiverDeallocated => {
                write!(f, "The receiver associated with the sender has been deallocated")
            }
        }
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl Error for ThreadedWriteError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        match self {
            ThreadedWriteError::IoError(e) => Some(e),
            ThreadedWriteError::ReceiverDeallocated => None,
        }
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl GetCode for ThreadedWriteError {
    fn get_code(&self) -> i32 {
        match self {
            ThreadedWriteError::IoError(e) => e.get_code(),
            ThreadedWriteError::ReceiverDeallocated => 1,
        }
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl AlignmentWriterThreaded {
    /// Constructs a [`AlignmentWriterThreaded`] from a regular writer by moving
    /// it into a thread and creating a channel.
    #[inline]
    #[must_use]
    pub fn from_writer<W>(mut writer: W) -> Self
    where
        W: Write + Send + 'static, {
        let (sender, receiver) = std::sync::mpsc::channel();
        let writer_thread = std::thread::spawn(move || -> std::io::Result<()> {
            while let Ok(string) = receiver.recv() {
                writeln!(writer, "{string}")?;
            }
            writer.flush()
        });

        Self {
            sender,
            writer_thread: Some(writer_thread),
        }
    }

    /// Writes a string to the [`AlignmentWriterThreaded`], properly handling
    /// errors if they occur.
    ///
    /// ## Errors
    ///
    /// If the cause of the failed write can be determined (e.g., by joining the
    /// thread), then that error is propagated as
    /// [`ThreadedWriteError::IoError`]. Otherwise,
    /// [`ThreadedWriteError::ReceiverDeallocated`] is returned.
    #[inline]
    pub fn write(&mut self, string: String) -> Result<(), ThreadedWriteError> {
        self.sender.send(string).map_err(|_| {
            if let Some(thread) = std::mem::take(&mut self.writer_thread)
                && let Err(e) = thread.join().unwrap()
            {
                ThreadedWriteError::IoError(e)
            } else {
                ThreadedWriteError::ReceiverDeallocated
            }
        })
    }

    /// Finalizes the writing by closing the thread and propagating any errors.
    ///
    /// This **must** be called on at least the originally created writer in
    /// order to properly handle all errors.
    #[inline]
    pub fn flush(self) -> std::io::Result<()> {
        if let Some(thread) = self.writer_thread {
            drop(self.sender);
            thread.join().unwrap()
        } else {
            Ok(())
        }
    }
}

/// Encapsulates the necessary logic in order for a writer to work with
/// `aligner`.
///
/// This is specifically designed to share logic between a multi-threaded
/// `AlignmentWriterThreaded` and a single-threaded `WriteFileZipStdout`.
pub trait AlignmentWriter: Sized {
    /// Given an unmapped alignment in a [`SamDataView`], write the alignment.
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> Result<(), WriterError>;

    /// Given an alignment in a [`SamDataView`] along with an alignment score,
    /// write the alignment.
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> Result<(), WriterError>;

    /// Writes an alignment in SAM format.
    ///
    /// The alignment should either correspond to:
    ///
    /// - The alignment of the query against the reference (if
    ///   [`Strand::Forward`] is passed)
    /// - The alignment of the reverse complement of the query against the
    ///   reference (if [`Strand::Reverse`]) is passed)
    ///
    /// The `MAPQ` field is not used and is set to 255. The optional `AS` tag
    /// for the score is included when the read is mapped. The query and
    /// reference name are truncated to only include the characters before the
    /// first whitespace. A trailing linebreak is not included.
    fn write_alignment<'q, 'r>(
        &mut self, alignment: AlignmentAndSeqs<'q, 'r>, config: &AlignerConfig,
    ) -> Result<(), WriterError> {
        let qname = process_header(&alignment.query.header);

        match alignment.mapping {
            Some(mapping) if mapping.inner.score > 0 => {
                let rname = process_header(&alignment.reference.name);
                let pos = mapping.inner.ref_range.start + 1;
                let mapq = 255;
                let cigar = mapping.inner.states.to_cigar_unchecked();

                match mapping.strand {
                    Strand::Forward => {
                        let flag = 0;
                        let seq = &alignment.query.sequence;
                        let qual = alignment
                            .query
                            .quality
                            .as_ref()
                            .map_or(QualityScoresView::try_from(b"*").unwrap(), DataOwned::as_view);
                        let record =
                            SamDataView::new(qname, flag, rname, pos, mapq, cigar.as_view(), seq.as_slice().into(), qual);
                        return self.write_record(record, mapping.inner.score);
                    }
                    Strand::Reverse => {
                        let flag = 16;
                        let seq = NucleotidesView::from(alignment.query.sequence.as_slice())
                            .to_reverse_complement()
                            .into_vec();
                        let qual = alignment
                            .query
                            .quality
                            .as_ref()
                            .map_or(QualityScores::try_from(b"*").unwrap(), |qual| qual.to_reverse());
                        let record = SamDataView::new(
                            qname,
                            flag,
                            rname,
                            pos,
                            mapq,
                            cigar.as_view(),
                            seq.as_slice().into(),
                            qual.as_view(),
                        );
                        return self.write_record(record, mapping.inner.score);
                    }
                };
            }
            _ => {
                if !config.exclude_unmapped {
                    return self.write_unmapped(SamDataView::unmapped(qname, "*"));
                }
            }
        };
        Ok(())
    }
}

#[cfg(feature = "dev_no_rayon")]
impl AlignmentWriter for WriteFileZipStdout {
    #[inline]
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> std::io::Result<()> {
        writeln!(self, "{record}")?;
        Ok(())
    }

    #[inline]
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> std::io::Result<()> {
        writeln!(self, "{record}\tAS:i:{score}")?;
        Ok(())
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl AlignmentWriter for AlignmentWriterThreaded {
    #[inline]
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> Result<(), ThreadedWriteError> {
        self.write(format!("{record}"))
    }

    #[inline]
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> Result<(), ThreadedWriteError> {
        self.write(format!("{record}\tAS:i:{score}"))
    }
}

/// Processes a header by removing everything after the first whitespace, or
/// using '*' if the header is unavailable.
#[inline]
fn process_header(header: &str) -> &str {
    header.split_ascii_whitespace().next().unwrap_or("*")
}

/// Writes a SAM-style header to the writer, containing the `HD` and `SQ` lines.
#[inline]
pub fn write_header<W: Write>(writer: &mut W, references: &[FastaSeq]) -> std::io::Result<()> {
    writeln!(writer, "@HD\tVN:1.4")?;
    for reference in references {
        writeln!(
            writer,
            "@SQ\tSN:{name}\tLN:{len}",
            name = process_header(&reference.name),
            len = reference.sequence.len()
        )?;
    }
    Ok(())
}
