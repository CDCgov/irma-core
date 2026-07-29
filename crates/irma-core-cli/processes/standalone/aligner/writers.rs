//! Traits and structs for writing the output of aligner

use crate::aligner::{AlignerConfig, AlignmentAndSeqs, Strand};
#[cfg(not(feature = "dev_no_rayon"))]
use irma_records::io::WriterThreaded;
use std::io::Write;
use zoe::{
    data::{fasta::FastaSeq, sam::SamDataView},
    math::AnyInt,
    prelude::{AsView, NucleotidesView, QualityScores, QualityScoresView},
};

#[cfg(feature = "dev_no_rayon")]
use irma_records::io::WriteFileZipStdout;

/// Encapsulates the necessary logic in order for a writer to work with
/// `aligner`.
///
/// This is specifically designed to share logic between multi-threaded
/// [`WriterThreaded`] and single-threaded [`WriteFileZipStdout`] writers.
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
    fn write_alignment<'q, 'r>(
        &mut self, alignment: AlignmentAndSeqs<'q, 'r>, config: &AlignerConfig,
    ) -> std::io::Result<()> {
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
                            .map_or(QualityScoresView::try_from(b"*").unwrap(), AsView::as_view);
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
impl AlignmentWriter for WriterThreaded {
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
