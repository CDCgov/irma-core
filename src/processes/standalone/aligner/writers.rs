//! Traits and structs for writing the output of aligner

use crate::{
    aligner::{AlignerConfig, Strand},
    io::{FastX, WriteFileZipStdout},
};
use std::io::Write;
use zoe::{
    alignment::Alignment,
    data::{fasta::FastaSeq, sam::SamDataView},
    math::AnyInt,
    prelude::{DataOwned, NucleotidesView, QualityScores, QualityScoresView},
};

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

/// A trait alias containing any relevant bounds needed for an
/// [`AlignmentWriter`], depending on the `dev_no_rayon` feature.
#[cfg(not(feature = "dev_no_rayon"))]
pub trait AdditionalBounds: Clone + Send {}

#[cfg(not(feature = "dev_no_rayon"))]
impl<T: Clone + Send> AdditionalBounds for T {}

#[cfg(feature = "dev_no_rayon")]
pub trait AdditionalBounds {}
#[cfg(feature = "dev_no_rayon")]
impl<T> AdditionalBounds for T {}

/// Encapsulates the necessary logic in order for a writer to work with
/// `aligner`.
///
/// This is specifically designed to share logic between a multi-threaded
/// `AlignmentWriterThreaded` and a single-threaded `WriteFileZipStdout`.
pub trait AlignmentWriter: Sized {
    /// Given an unmapped alignment in a [`SamDataView`], write the alignment.
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> std::io::Result<()>;

    /// Given an alignment in a [`SamDataView`] along with an alignment score,
    /// write the alignment.
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> std::io::Result<()>;

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
    fn write_alignment<T: AnyInt>(
        &mut self, alignment: Option<(Alignment<T>, Strand)>, query: &FastX, reference: &FastaSeq, config: &AlignerConfig,
    ) -> std::io::Result<()> {
        let qname = process_header(&query.header);

        match alignment {
            Some((alignment, strand)) if alignment.score > T::ZERO => {
                let rname = process_header(&reference.name);
                let pos = alignment.ref_range.start + 1;
                let mapq = 255;
                let cigar = alignment.states.to_cigar_unchecked();

                match strand {
                    Strand::Forward => {
                        let flag = 0;
                        let seq = &query.sequence;
                        let qual = query
                            .quality
                            .as_ref()
                            .map_or(QualityScoresView::try_from(b"*").unwrap(), DataOwned::as_view);
                        let record =
                            SamDataView::new(qname, flag, rname, pos, mapq, cigar.as_view(), seq.as_slice().into(), qual);
                        return self.write_record(record, alignment.score);
                    }
                    Strand::Reverse => {
                        let flag = 16;
                        let seq = NucleotidesView::from(query.sequence.as_slice())
                            .to_reverse_complement()
                            .into_vec();
                        let qual = query
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
                        return self.write_record(record, alignment.score);
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

impl AlignmentWriter for WriteFileZipStdout {
    #[inline]
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> std::io::Result<()> {
        writeln!(self, "{record}")
    }

    #[inline]
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> std::io::Result<()> {
        writeln!(self, "{record}\tAS:i:{score}")
    }
}

#[cfg(not(feature = "dev_no_rayon"))]
impl AlignmentWriter for AlignmentWriterThreaded {
    #[inline]
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> std::io::Result<()> {
        self.sender
            .send(format!("{record}"))
            .expect("The receiver associated with `sender` has been deallocated");
        Ok(())
    }

    #[inline]
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> std::io::Result<()> {
        self.sender
            .send(format!("{record}\tAS:i:{score}"))
            .expect("The receiver associated with `sender` has been deallocated");

        Ok(())
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
