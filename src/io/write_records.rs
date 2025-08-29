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
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{self}")
    }
}

impl<W: Write> WriteRecord<W> for FastQView<'_> {
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{self}")
    }
}

impl<W: Write> WriteRecord<W> for FastQViewMut<'_> {
    #[inline]
    fn write_record(self, writer: &mut W) -> std::io::Result<()> {
        write!(writer, "{self}")
    }
}

impl<W: Write> WriteRecord<W> for FastaSeq {
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
    #[inline]
    fn write_record(self, writer: &mut PairedWriters<W>) -> std::io::Result<()> {
        let [read1, read2] = self?;
        read1.write_record(&mut writer.writer1)?;
        read2.write_record(&mut writer.writer2)
    }
}

/// An extension trait for iterators allowing all of their records to be written
/// to a supported writer type `W`. This relies on the [`WriteRecord`] trait.
pub trait WriteRecords<W> {
    /// Write all the records the iterator to the `writer` of supported type
    /// `W`.
    fn write_records(self, writer: W) -> std::io::Result<()>;
}

impl<I, W> WriteRecords<W> for I
where
    I: Iterator<Item: WriteRecord<W>>,
    W: FlushWriter,
{
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
    #[inline]
    fn write_records(self, writer: RecordWriters<W>) -> std::io::Result<()> {
        match writer {
            RecordWriters::SingleEnd(writer) => self.write_records(writer),
            RecordWriters::PairedEnd(writer) => self.write_records(writer),
        }
    }
}
