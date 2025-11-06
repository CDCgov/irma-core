//! Provides the [`WriteRecord`] and [`WriteRecords`] traits for abstracting
//! over the record types (and other related structs) which can be written to a
//! single output file or paired output files.

use crate::io::{FlushWriter, PairedWriters, RecordWriters};
use std::io::Write;
use zoe::{
    data::fasta::FastaSeq,
    prelude::{FastQ, FastQView, FastQViewMut},
};

/// A trait providing the ability for a record or record-like struct to write
/// itself, where the supported writer is given by `W`.
pub trait WriteRecord<W> {
    /// Writes the record to the writer of supported type `W`.
    fn write_record(self, writer: &mut W) -> std::io::Result<()>;
}

impl<W: Write> WriteRecord<W> for FastQ {
    /// Writes a [`FastQ`] record to a single writer.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{self}")
    }
}

impl<W: Write> WriteRecord<W> for FastQView<'_> {
    /// Writes a [`FastQView`] record to a single writer.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{self}")
    }
}

impl<W: Write> WriteRecord<W> for FastQViewMut<'_> {
    /// Writes a [`FastQViewMut`] record to a single writer.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{self}")
    }
}

impl<W: Write> WriteRecord<W> for FastaSeq {
    /// Writes a [`FastaSeq`] record to a single writer.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{self}")
    }
}

impl<W, E> WriteRecord<W> for Result<FastQ, E>
where
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes a [`FastQ`] record to a single writer, propagating an error if
    /// present.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{}", self?)
    }
}

impl<W, E> WriteRecord<W> for Result<FastQView<'_>, E>
where
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes a [`FastQView`] record to a single writer, propagating an error
    /// if present.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{}", self?)
    }
}

impl<W, E> WriteRecord<W> for Result<FastQViewMut<'_>, E>
where
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes a [`FastQViewMut`] record to a single writer, propagating an
    /// error if present.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{}", self?)
    }
}

impl<W, E> WriteRecord<W> for Result<FastaSeq, E>
where
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes a [`FastaSeq`] record to a single writer, propagating an error if
    /// present.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{}", self?)
    }
}

impl<A, B, W> WriteRecord<W> for (A, B)
where
    A: WriteRecord<W>,
    B: WriteRecord<W>,
    W: Write,
{
    /// Writes a tuple of two individual records to a single writer,
    /// interleaving them.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        let (read1, read2) = self;
        read1.write_record(writer)?;
        read2.write_record(writer)
    }
}

impl<A, B, W> WriteRecord<PairedWriters<W>> for (A, B)
where
    A: WriteRecord<W>,
    B: WriteRecord<W>,
    W: Write,
{
    /// Writes a tuple of two individual records to a pair of writers in a
    /// [`PairedWriters`] struct.
    #[inline]
    fn write_record(self, writer: &mut PairedWriters<W>) -> std::io::Result<()> {
        let (read1, read2) = self;
        read1.write_record(&mut writer.writer1)?;
        read2.write_record(&mut writer.writer2)
    }
}

impl<A, B, W, E> WriteRecord<W> for Result<(A, B), E>
where
    A: WriteRecord<W>,
    B: WriteRecord<W>,
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes a tuple of two individual records to a single writer,
    /// interleaving them and propagating an error if present.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        let (read1, read2) = self?;
        read1.write_record(writer)?;
        read2.write_record(writer)
    }
}

impl<A, B, W, E> WriteRecord<PairedWriters<W>> for Result<(A, B), E>
where
    A: WriteRecord<W>,
    B: WriteRecord<W>,
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes a tuple of two individual records to a pair of writers in a
    /// [`PairedWriters`] struct, propagating an error if present.
    #[inline]
    fn write_record(self, writer: &mut PairedWriters<W>) -> std::io::Result<()> {
        let (read1, read2) = self?;
        read1.write_record(&mut writer.writer1)?;
        read2.write_record(&mut writer.writer2)
    }
}

impl<A, W> WriteRecord<W> for [A; 2]
where
    A: WriteRecord<W>,
    W: Write,
{
    /// Writes an array of two individual records to a single writer,
    /// interleaving them.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        let [read1, read2] = self;
        read1.write_record(writer)?;
        read2.write_record(writer)
    }
}

impl<A, W> WriteRecord<PairedWriters<W>> for [A; 2]
where
    A: WriteRecord<W>,
    W: Write,
{
    /// Writes an array of two individual records to a pair of writers in a
    /// [`PairedWriters`] struct.
    #[inline]
    fn write_record(self, writer: &mut PairedWriters<W>) -> std::io::Result<()> {
        let [read1, read2] = self;
        read1.write_record(&mut writer.writer1)?;
        read2.write_record(&mut writer.writer2)
    }
}

impl<A, W, E> WriteRecord<W> for Result<[A; 2], E>
where
    A: WriteRecord<W>,
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes an array of two individual records to a single writer,
    /// interleaving them and propagating an error if present.
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        let [read1, read2] = self?;
        read1.write_record(writer)?;
        read2.write_record(writer)
    }
}

impl<A, W, E> WriteRecord<PairedWriters<W>> for Result<[A; 2], E>
where
    A: WriteRecord<W>,
    W: Write,
    std::io::Error: From<E>,
{
    /// Writes an array of two individual records to a pair of writers in a
    /// [`PairedWriters`] struct, propagating an error if present.
    #[inline]
    fn write_record(self, writer: &mut PairedWriters<W>) -> std::io::Result<()> {
        let [read1, read2] = self?;
        read1.write_record(&mut writer.writer1)?;
        read2.write_record(&mut writer.writer2)
    }
}

/// An extension trait for iterators allowing all of their records to be written
/// to a supported writer type `W`.
///
/// This relies on the [`WriteRecord`] trait.
pub trait WriteRecords<W> {
    /// Writes all the records in the iterator to the `writer` of supported type
    /// `W`.
    fn write_records(self, writer: W) -> std::io::Result<()>;
}

impl<I, W> WriteRecords<W> for I
where
    I: Iterator<Item: WriteRecord<W>>,
    W: FlushWriter,
{
    /// Writes all the records in the iterator to the `writer` of supported type
    /// `W`.
    #[inline]
    fn write_records(mut self, mut writer: W) -> std::io::Result<()> {
        self.try_for_each(|record| record.write_record(&mut writer))?;
        writer.flush_all()
    }
}

impl<I, W> WriteRecords<RecordWriters<W>> for I
where
    I: Iterator + WriteRecords<W> + WriteRecords<PairedWriters<W>>,
    W: Write,
{
    /// Writes all the records in the iterator to the one or two writes
    /// contained in the [`RecordWriters`] enum.
    ///
    /// The reads will be interleaved if there is only one writer.
    #[inline]
    fn write_records(self, writer: RecordWriters<W>) -> std::io::Result<()> {
        match writer {
            RecordWriters::SingleEnd(writer) => self.write_records(writer),
            RecordWriters::PairedEnd(writer) => self.write_records(writer),
        }
    }
}
