use crate::io::{PairedWriters, RecordWriters, WriteFileZipStdout, WriterWithContext};
use flate2::write::GzEncoder;
use std::{
    collections::VecDeque,
    fs::File,
    io::{BufWriter, Cursor, Empty, LineWriter, PipeWriter, Sink, Stderr, Stdout, Write},
    process::ChildStdin,
};
use zoe::data::err::ResultWithErrorContext;

/// A trait allowing a writer to be finished, which may include flushing,
/// writing any footers, etc.
///
/// This does not have a blanket impl over all of [`Write`] to prevent
/// conflicts. For example, [`GzEncoder`] implements [`Write`] but also has an
/// inherent [`finish`] method which must be called. If additional unsupported
/// use-cases are encountered, implementations may need to be added within
/// IRMA-core, or a wrapper type can be used.
///
/// [`finish`]: GzEncoder::finish
pub trait Finish {
    /// Finalizes the writer, performing flushing, writing any footers, etc.
    fn finish(self) -> std::io::Result<()>;
}

impl<W: Write + Finish> Finish for GzEncoder<W> {
    fn finish(self) -> std::io::Result<()> {
        let inner = self.finish()?;
        inner.finish()
    }
}

impl<W: Finish> Finish for PairedWriters<W> {
    fn finish(self) -> std::io::Result<()> {
        let res1 = self.writer1.finish();
        let res2 = self.writer2.finish();
        res1.and(res2)
    }
}

impl<W: Finish> Finish for RecordWriters<W> {
    fn finish(self) -> std::io::Result<()> {
        match self {
            RecordWriters::SingleEnd(inner) => inner.finish(),
            RecordWriters::PairedEnd(inner) => inner.finish(),
        }
    }
}

impl Finish for WriteFileZipStdout {
    fn finish(self) -> std::io::Result<()> {
        match self {
            WriteFileZipStdout::File(writer) => writer.finish(),
            WriteFileZipStdout::Zipped(writer) => writer.finish(),
            WriteFileZipStdout::Stdout(writer) => writer.finish(),
        }
    }
}

// Stdlib impls for implementors of Write

macro_rules! finish_for_write_no_ref {
     ($($ty:ty),*) => {
        $(
            impl Finish for $ty {
                fn finish(mut self) -> std::io::Result<()> {
                    self.flush()
                }
            }

            impl Finish for &mut $ty {
                fn finish(self) -> std::io::Result<()> {
                    self.flush()
                }
            }
        )*
    };
}

macro_rules! finish_for_write {
    ($($ty:ty),*) => {
        finish_for_write_no_ref!{$($ty),*}

        $(
            impl Finish for &$ty {
                fn finish(mut self) -> std::io::Result<()> {
                    self.flush()
                }
            }
        )*
    };
}

finish_for_write! {
    File,
    Stdout,
    Stderr,
    ChildStdin,
    PipeWriter,
    Empty,
    Sink
}

finish_for_write_no_ref! {
    Cursor<&mut [u8]>,
    Vec<u8>,
    VecDeque<u8>
}

impl Finish for &mut [u8] {
    fn finish(mut self) -> std::io::Result<()> {
        self.flush()
    }
}

impl<const N: usize> Finish for Cursor<[u8; N]> {
    fn finish(mut self) -> std::io::Result<()> {
        self.flush()
    }
}

impl<const N: usize> Finish for &mut Cursor<[u8; N]> {
    fn finish(self) -> std::io::Result<()> {
        self.flush()
    }
}

// Blanket impls

impl<W: Finish + Write> Finish for BufWriter<W> {
    fn finish(self) -> std::io::Result<()> {
        // Flush the buffer into the underlying writer
        let inner = self.into_inner()?;
        // Finish the inner writer
        inner.finish()
    }
}

impl<W: Finish + Write> Finish for LineWriter<W> {
    fn finish(self) -> std::io::Result<()> {
        // Flush the buffer into the underlying writer
        let inner = self.into_inner()?;
        // Finish the inner writer
        inner.finish()
    }
}

impl<W: Finish> Finish for WriterWithContext<W> {
    fn finish(self) -> std::io::Result<()> {
        self.writer.finish().with_context(self.description)?;
        Ok(())
    }
}

impl<W: Finish> Finish for Box<W> {
    fn finish(self) -> std::io::Result<()> {
        (*self).finish()
    }
}
