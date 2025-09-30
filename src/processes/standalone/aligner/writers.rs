//! Traits and structs for writing the output of aligner

use crate::{
    aligner::{AlignerConfig, Strand},
    io::{FastX, FromFilename, WriteFileZipStdout},
};
use std::{
    io::Write,
    path::PathBuf,
    sync::mpsc,
    thread::{self, JoinHandle},
};
use zoe::{
    alignment::Alignment,
    data::{fasta::FastaSeq, sam::SamDataView},
    math::AnyInt,
    prelude::{DataOwned, NucleotidesView, QualityScores, QualityScoresView},
};

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
/// channel ([`mpsc::Sender`]) and a single-threaded [`WriteFileZipStdout`].
pub trait AlignmentWriter: Sized {
    /// Any additional handles or data needed to finalize the writer.
    type Handle;

    /// Creates a new [`AlignmentWriter`], while also returning any additional
    /// handles or data needed to finalize the writer.
    fn new_writer(path: Option<PathBuf>) -> std::io::Result<(Self, Self::Handle)>;

    /// Finalizes the writer, so that all contents are flushed, threads are
    /// joined, and/or errors are propagated.
    fn finalize_writer(self, handle: Self::Handle) -> std::io::Result<()>;

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

impl AlignmentWriter for mpsc::Sender<String> {
    type Handle = JoinHandle<std::io::Result<()>>;

    #[inline]
    fn new_writer(path: Option<PathBuf>) -> std::io::Result<(Self, Self::Handle)> {
        let mut output = WriteFileZipStdout::from_optional_filename(path)?;
        let (sender, receiver) = mpsc::channel();
        let writer_thread = thread::spawn(move || -> std::io::Result<()> {
            while let Ok(string) = receiver.recv() {
                writeln!(output, "{string}")?;
            }
            output.flush()
        });
        Ok((sender, writer_thread))
    }

    #[inline]
    fn finalize_writer(self, handle: Self::Handle) -> std::io::Result<()> {
        // Must drop the sender before the writer thread can finish
        drop(self);
        handle.join().unwrap()
    }

    #[inline]
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> std::io::Result<()> {
        self.send(format!("{record}"))
            .expect("The receiver associated with `sender` has been deallocated");
        Ok(())
    }

    #[inline]
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> std::io::Result<()> {
        self.send(format!("{record}\tAS:i:{score}"))
            .expect("The receiver associated with `sender` has been deallocated");

        Ok(())
    }
}

impl AlignmentWriter for WriteFileZipStdout {
    type Handle = ();

    #[inline]
    fn new_writer(path: Option<PathBuf>) -> std::io::Result<(Self, Self::Handle)> {
        let writer = WriteFileZipStdout::from_optional_filename(path)?;
        Ok((writer, ()))
    }

    #[inline]
    fn finalize_writer(mut self, _handle: Self::Handle) -> std::io::Result<()> {
        self.flush()
    }

    #[inline]
    fn write_unmapped<'a>(&mut self, record: SamDataView<'a>) -> std::io::Result<()> {
        writeln!(self, "{record}")
    }

    #[inline]
    fn write_record<'a, T: AnyInt>(&mut self, record: SamDataView<'a>, score: T) -> std::io::Result<()> {
        writeln!(self, "{record}\tAS:i:{score}")
    }
}

/// Processes a header by removing everything after the first whitespace, or
/// using '*' if the header is unavailable.
#[inline]
fn process_header(header: &str) -> &str {
    header.split_ascii_whitespace().next().unwrap_or("*")
}
