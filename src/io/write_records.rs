//! Provides the [`WriteRecord`] and [`WriteRecords`] traits for abstracting
//! over the record types (and other related structs) which can be written to a
//! single output file or paired output files.
//!
//! The core building block is [`WriteRecord`], which enables writing many
//! different types, such as [`FastQ`], [`FastaSeq`], the corresponding view
//! types, two-tuples of these, length-two arrays of these, and results of
//! these. The trait is generic over `W`, which can be a single end writer
//! (implementing [`Write`]) or a paired end writer ([`PairedWriters`]).
//!
//! For single output files, the reads in tuples and arrays are interleaved. For
//! paired output files, only the tuple/array values are permitted; the first
//! read is written to the first file, and the second is written to the second
//! file. Any errors in the inputs are propagated.
//!
//! This trait does *not* support [`RecordWriters`], an enum containing either
//! single or paired end writers. This is because performing a match statement
//! on every write is inefficient and not idiomatic.
//!
//! Instead, a single match statement should be performed, followed by writing
//! all the records from an iterator. This is provided via the [`WriteRecords`]
//! trait, which is an extension trait for iterators. It supports the same
//! combinations of types as [`WriteRecord`], but it also allows `W` to be
//! [`RecordWriters`] (and the implementation ensures that only a single match
//! statement is performed).
//!
//! ## Specifying Trait Bounds for [`WriteRecords`]
//!
//! The benefit of these traits is that it enables functions to be generic over
//! the types of inputs and the types of writers. If a function accepts an
//! iterator of type `I` and directly writes these records to `W` without any
//! iterator adaptors, then the trait bound `I: WriteRecords<W>` can be used.
//! However, more often, a chain of iterator adaptors is used to filter or
//! otherwise transform the input iterator. The iterator type for which
//! `write_records` is called may not even be nameable, if closures are
//! involved.
//!
//! In this case, the trait bound should be specified in terms of the item `A`
//! that is being written, rather than in terms of `I`. Two trait bounds are
//! necessary:
//!
//! 1. Ensuring `W` satisfies [`SequenceWriter`]. This trait bound may not be
//!    needed if a specific writer type is being supported (e.g.,
//!    [`RecordWriters`], [`PairedWriters`], or just a single-output writer
//!    implementing [`Write`]). However, if the function is generic over the
//!    type of writer, this is the proper trait bound.
//! 2. Ensuring `A` satisfies [`WriteRecordCompatibleItem<W>`]. This is a marker
//!    trait ensuring that the record can in fact be written to the writer and
//!    any of its variants.

use crate::io::{PairedWriters, RecordWriters};
use std::io::Write;
use zoe::{
    data::fasta::FastaSeq,
    prelude::{FastQ, FastQView, FastQViewMut},
};

/// A trait providing the ability for a record or record-like struct to write
/// itself, where the supported writer is given by `W`.
pub trait WriteRecord<W> {
    /// Writes the record to the writer of supported type `W`.
    ///
    /// The two types of writes for `W` currently supported are:
    ///
    /// - A generic `W` implementing [`Write`], which represents writing to a
    ///   single output file (single end reads)
    /// - [`PairedWriters`], which represents writing to two output files
    ///   (paired end reads)
    ///
    /// This function is available for [`FastQ`], [`FastaSeq`], the
    /// corresponding view types, two-tuples of these, length-two arrays of
    /// these, and results of these. For single end reads, the reads in tuples
    /// and arrays are interleaved. For paired end reads, the first is written
    /// to the first file, and the second is written to the second file. Any
    /// errors in the inputs are propagated.
    ///
    /// `W` being [`RecordWriters`] is not supported. [`RecordWriters`] is an
    /// enum holding either a single writer or a pair of writers, and hence each
    /// [`write_record`] call would incur a match statement. The idiomatic
    /// solution is to use the [`WriteRecords`] trait, which would do a single
    /// match statement and then write all the records in the iterator.
    ///
    /// [`write_record`]: WriteRecord::write_record
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

/// A trait marking that the writer is either a basic writer implementing
/// [`Write`], or is [`PairedWriters`]. In other words, the writer is
/// potentially compatible with [`WriteRecord`].
///
/// This trait also provides a method for flushing the writer(s), for use in
/// [`SequenceWriter`].
trait BasicSequenceWriter {
    /// Flushes all the writers contained in the struct.
    fn flush_writer(&mut self) -> std::io::Result<()>;
}

impl<W: Write> BasicSequenceWriter for W {
    #[inline]
    fn flush_writer(&mut self) -> std::io::Result<()> {
        self.flush()
    }
}

impl<W: Write> BasicSequenceWriter for PairedWriters<W> {
    #[inline]
    fn flush_writer(&mut self) -> std::io::Result<()> {
        self.writer1.flush()?;
        self.writer2.flush()
    }
}

/// A trait for unifying writers of a known type (implementing
/// [`BasicSequenceWriter`], such as single output writers and
/// [`PairedWriters`]) with writers of an unknown type ([`RecordWriters`], which
/// is an enum containing either one).
pub trait SequenceWriter: Sized {
    /// The type of the first possible writer.
    type Writer1;

    /// The type of the second possible writer, or the same as [`Writer1`] if
    /// not applicable.
    ///
    /// [`Writer1`]: SequenceWriter::Writer1
    type Writer2;

    /// Writes all records in the iterator to the writer.
    ///
    /// For most use-cases, [`WriteRecords::write_records`] will be a more
    /// useful syntax.
    fn write_records_from<I>(self, iterator: I) -> std::io::Result<()>
    where
        I: Iterator<Item: WriteRecordCompatibleItem<Self>>;
}

impl<W: BasicSequenceWriter> SequenceWriter for W {
    type Writer1 = W;
    type Writer2 = W;

    fn write_records_from<I>(mut self, mut iterator: I) -> std::io::Result<()>
    where
        I: Iterator<Item: WriteRecord<Self::Writer1> + WriteRecord<Self::Writer2>>, {
        iterator.try_for_each(|record| record.write_record(&mut self))?;
        self.flush_writer()
    }
}

impl<W> SequenceWriter for RecordWriters<W>
where
    W: BasicSequenceWriter,
    PairedWriters<W>: BasicSequenceWriter,
{
    type Writer1 = W;
    type Writer2 = PairedWriters<W>;

    fn write_records_from<I>(self, mut iterator: I) -> std::io::Result<()>
    where
        I: Iterator<Item: WriteRecord<Self::Writer1> + WriteRecord<Self::Writer2>>, {
        match self {
            RecordWriters::SingleEnd(mut writer) => {
                iterator.try_for_each(|record| record.write_record(&mut writer))?;
                writer.flush_writer()
            }
            RecordWriters::PairedEnd(mut writer) => {
                iterator.try_for_each(|record| record.write_record(&mut writer))?;
                writer.flush_writer()
            }
        }
    }
}

/// A marker trait for types which can be written to `W` with [`WriteRecords`].
///
/// This is similar to a [`WriteRecord`] trait bound, but it also allows `W` to
/// be [`RecordWriters`]. So, for example, `[FastQ; 2]` implements
/// `WriteRecordCompatibleItem<File>` (interleaved),
/// `WriteRecordCompatibleItem<PairedWriters<File>>` (paired files), and
/// `WriteRecordCompatibleItem<RecordWriters<File>>` (maybe interleaved or maybe
/// paired files).
///
/// This is useful for concisely specifying trait bounds.
pub trait WriteRecordCompatibleItem<W: SequenceWriter>: WriteRecord<W::Writer1> + WriteRecord<W::Writer2> {}

impl<T, W> WriteRecordCompatibleItem<W> for T
where
    W: SequenceWriter,
    T: WriteRecord<W::Writer1> + WriteRecord<W::Writer2>,
{
}

/// An extension trait for iterators allowing all of their records to be written
/// to a supported writer type `W`.
///
/// This relies on the [`WriteRecord`] trait.
pub trait WriteRecords<W> {
    /// Write all the records in the iterator to the `writer` of supported type
    /// `W`.
    ///
    /// See [`WriteRecord::write_record`] for more details. This function can be
    /// called with a writer implementing [`SequenceWriter`], or it can be
    /// called with a [`RecordWriters`] struct (which will entail a single match
    /// statement to determine whether the writer is single or paired end
    /// reads).
    fn write_records(self, writer: W) -> std::io::Result<()>;
}

impl<I, W> WriteRecords<W> for I
where
    I: Iterator<Item: WriteRecordCompatibleItem<W>>,
    W: SequenceWriter,
{
    #[inline]
    fn write_records(self, writer: W) -> std::io::Result<()> {
        SequenceWriter::write_records_from(writer, self)
    }
}
