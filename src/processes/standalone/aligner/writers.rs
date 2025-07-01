//! Traits and structs for writing the output of aligner

use crate::{
    aligner::{AlignerConfig, Strand},
    io::FastX,
};
use std::sync::mpsc;
use zoe::{
    alignment::Alignment,
    data::{fasta::FastaSeq, sam::SamDataView},
    math::AnyInt,
    prelude::{DataOwned, NucleotidesView, QualityScores, QualityScoresView},
};

/// Sends an alignment to `sender` in SAM format.
///
/// The alignment should either correspond to:
///
/// - The alignment of the query against the reference (if [`Strand::Forward`]
///   is passed)
/// - The alignment of the reverse complement of the query against the reference
///   (if [`Strand::Reverse`]) is passed)
///
/// The `MAPQ` field is not used and is set to 255. The optional `AS` tag for
/// the score is included when the read is mapped. The query and reference name
/// are truncated to only include the characters before the first whitespace. A
/// trailing linebreak is not included.
pub fn send_alignment<T: AnyInt>(
    sender: &mut mpsc::Sender<String>, alignment: Option<(Alignment<T>, Strand)>, query: &FastX, reference: &FastaSeq,
    config: &AlignerConfig,
) {
    let qname = process_header(&query.header);

    match alignment {
        Some((alignment, strand)) if alignment.score > T::ZERO => {
            let rname = process_header(&reference.name);
            let pos = alignment.ref_range.start + 1;
            let mapq = 255;
            let cigar = alignment.states.to_cigar_unchecked();

            let line = match strand {
                Strand::Forward => {
                    let flag = 0;
                    let seq = &query.sequence;
                    let qual = query
                        .quality
                        .as_ref()
                        .map_or(QualityScoresView::try_from(b"*").unwrap(), DataOwned::as_view);
                    let record =
                        SamDataView::new(qname, flag, rname, pos, mapq, cigar.as_view(), seq.as_slice().into(), qual);
                    format!("{record}\tAS:i:{score}", score = alignment.score)
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
                    format!("{record}\tAS:i:{score}", score = alignment.score)
                }
            };

            sender
                .send(line)
                .expect("The receiver associated with `sender` has been deallocated");
        }
        _ => {
            if !config.exclude_unmapped {
                sender
                    .send(SamDataView::unmapped(qname, "*").to_string())
                    .expect("The receiver associated with `sender` has been deallocated");
            }
        }
    };
}

/// Processes a header by removing everything after the first whitespace, or
/// using '*' if the header is unavailable.
fn process_header(header: &str) -> &str {
    header.split_ascii_whitespace().next().unwrap_or("*")
}
