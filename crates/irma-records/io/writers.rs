use crate::io::{Finish, WithContext, WriterWithContext, is_gz};
use flate2::{Compression, write::GzEncoder};
use std::{
    error::Error,
    fmt::Display,
    fs::File,
    io::{BufWriter, Stdout, Write, stdout},
    path::Path,
    sync::{
        Arc, OnceLock,
        mpsc::{RecvError, SendError},
    },
    thread::JoinHandle,
};
use zoe::{define_whichever, impl_traits};

/// A writer enabling writing to either a file or stdout.
///
/// This includes buffering via [`BufWriter`] and context with in
/// [`WithContext`] to add context to write errors.
pub struct WriteFileStdout(WithContext<BufWriter<WriteFileStdoutInner>>);

define_whichever! {
    /// An enum unifying writing to a regular file or stdout, without buffering
    /// or context.
    ///
    /// Since all variants require buffering, we do not include [`BufWriter`]
    /// inside the variants, and instead wrap the enum in a [`BufWriter`]. This
    /// allows dispatch to _only_ occur when writing a full buffer, as compared
    /// to on every write.
    #[derive(Debug)]
    enum WriteFileStdoutInner {
        /// A writer for a regular uncompressed file.
        File(File),
        /// A writer for uncompressed data to stdout.
        Stdout(Stdout),
    }

    impl Write for WriteFileStdoutInner {}
}

impl_traits! {
    impl Write for WriteFileStdout {}
}

define_whichever! {
    /// An enum unifying writing to a regular file, a gzip compressed file, or
    /// stdout.
    ///
    /// A [`BufWriter`] is used for all variants, and all variants are wrapped
    /// in [`WriterWithContext`] to add context to write errors.
    ///
    /// ## Limitations
    ///
    /// Currently only the default level of compression is exposed.
    #[derive(Debug)]
    pub enum WriteFileZipStdout {
        /// A writer for a regular uncompressed file.
        File(WithContext<BufWriter<File>>),
        /// A writer for a gzip compressed file.
        Zipped(WithContext<GzEncoder<BufWriter<File>>>),
        /// A writer for uncompressed data to stdout.
        Stdout(WithContext<BufWriter<Stdout>>),
    }

    impl Write for WriteFileZipStdout {}
}

impl WriteFileStdout {
    /// Creates a new [`WriteFileStdout`] from an optional filename. If a path
    /// is not provided, `stdout` is used.
    ///
    /// ## Errors
    ///
    /// If a path is provided, any IO errors when creating the file are
    /// propagated. If no path is provided, this method is infallible. Any
    /// failed writes will have context added including the path if available.
    pub fn create(path: Option<impl AsRef<Path>>) -> std::io::Result<Self> {
        match path {
            Some(path) => {
                let file = File::create(&path)?;
                Ok(Self(
                    BufWriter::new(WriteFileStdoutInner::File(file)).with_path_context("Failed to write to file", path),
                ))
            }
            None => Ok(Self(
                BufWriter::new(WriteFileStdoutInner::Stdout(stdout())).with_context("Failed to write to stdout"),
            )),
        }
    }

    /// Similar to [`WriteFileStdout::create`], but uses a specified `capacity`
    /// for the underlying [`BufWriter`].
    pub fn with_capacity(capacity: usize, path: Option<impl AsRef<Path>>) -> std::io::Result<Self> {
        match path {
            Some(path) => {
                let file = File::create(&path)?;
                Ok(Self(
                    BufWriter::with_capacity(capacity, WriteFileStdoutInner::File(file))
                        .with_path_context("Failed to write to file", path),
                ))
            }
            None => Ok(Self(
                BufWriter::with_capacity(capacity, WriteFileStdoutInner::Stdout(stdout()))
                    .with_context("Failed to write to stdout"),
            )),
        }
    }
}

impl WriteFileZipStdout {
    /// Creates a new [`WriteFileZipStdout`] from an optional filename. If a
    /// path is not provided, [`WriteFileZipStdout::Stdout`] is used.
    ///
    /// ## Errors
    ///
    /// If a path is provided, any IO errors when creating the file are
    /// propagated. If no path is provided, this method is infallible. Any
    /// failed writes will have context added including the path if available.
    pub fn create(path: Option<impl AsRef<Path>>) -> std::io::Result<Self> {
        match path {
            Some(path) => {
                let file = File::create(&path)?;
                let bufwriter = BufWriter::new(file);

                let writer = if is_gz(&path) {
                    Self::Zipped(
                        GzEncoder::new(bufwriter, Compression::default())
                            .with_path_context("Failed to write to zipped file", path),
                    )
                } else {
                    Self::File(bufwriter.with_path_context("Failed to write to file", path))
                };

                Ok(writer)
            }
            None => Ok(WriteFileZipStdout::Stdout(
                BufWriter::new(stdout()).with_context("Failed to write to stdout"),
            )),
        }
    }

    /// Similar to [`WriteFileZipStdout::create`], but uses a specified
    /// `capacity` for the underlying [`BufWriter`].
    pub fn with_capacity(capacity: usize, path: Option<impl AsRef<Path>>) -> std::io::Result<Self> {
        match path {
            Some(path) => {
                let file = File::create(&path)?;
                let bufwriter = BufWriter::with_capacity(capacity, file);

                let writer = if is_gz(&path) {
                    Self::Zipped(
                        GzEncoder::new(bufwriter, Compression::default())
                            .with_path_context("Failed to write to zipped file", path),
                    )
                } else {
                    Self::File(bufwriter.with_path_context("Failed to write to file", path))
                };

                Ok(writer)
            }
            None => Ok(WriteFileZipStdout::Stdout(
                BufWriter::with_capacity(capacity, stdout()).with_context("Failed to write to stdout"),
            )),
        }
    }
}

impl Finish for WriteFileStdoutInner {
    fn finish(self) -> std::io::Result<()> {
        match self {
            WriteFileStdoutInner::File(writer) => writer.finish(),
            WriteFileStdoutInner::Stdout(writer) => writer.finish(),
        }
    }
}

impl Finish for WriteFileStdout {
    fn finish(self) -> std::io::Result<()> {
        self.0.finish()
    }
}

/// A struct containing two writers for paired reads: one for left reads and one
/// for the right.
///
/// This is compatible with the [`WriteRecord`] trait, so that two-tuples and
/// length-two arrays of records can be written to [`PairedWriters`].
///
/// [`WriteRecord`]: crate::io::WriteRecord
pub struct PairedWriters<W> {
    pub writer1: W,
    pub writer2: W,
}

impl<W> PairedWriters<W> {
    /// Creates a new [`PairedWriters`] from two writers.
    #[inline]
    pub fn new(writer1: W, writer2: W) -> Self {
        Self { writer1, writer2 }
    }
}

/// An enum for holding either a single writer (single reads) or
/// [`PairedWriters`] (for paired reads).
///
/// This is compatible with the [`WriteRecords`] trait, so that given an
/// iterator of reads, they can be written to either two writer if present, or
/// written to one and interleaved. This involves a single match statement, and
/// so does not incur significant overhead. This is *not* compatible with the
/// [`WriteRecord`] trait, since performing a match on each write would be
/// inefficient.
///
/// [`WriteRecord`]: crate::io::WriteRecord
/// [`WriteRecords`]: crate::io::WriteRecords
pub enum RecordWriters<W> {
    /// A single writer for single end reads.
    SingleEnd(W),
    /// A pair of writers for paired reads.
    PairedEnd(PairedWriters<W>),
}

impl<W> RecordWriters<W> {
    /// Creates a new [`RecordWriters`] object to represent either a single
    /// writer (single reads) or two writers (paired reads).
    ///
    /// This is used for parsing clap arguments.
    #[inline]
    pub fn new(writer1: W, writer2: Option<W>) -> Self {
        match writer2 {
            Some(writer2) => Self::PairedEnd(PairedWriters::new(writer1, writer2)),
            None => Self::SingleEnd(writer1),
        }
    }
}

/// A clonable writer supporting writing from multiple threads via an [`mpsc`]
/// channel.
///
/// A single dedicated thread is used for writing to the underlying writer to
/// avoid interleaved writes. Calling [`write`] is non-blocking, instead
/// queueing the bytes in the channel, or in a per-clone buffer when batching is
/// enabled. Flush before cloning if earlier writes must precede clone writes.
/// Calling [`flush`] blocks until all previously queued data has been written
/// to the underlying writer.
///
/// To ensure that errors are not silently ignored, for each linked
/// [`WriterThreaded`], ensure that either [`flush`] is called on it, or it is
/// dropped before a [`flush`] call on a different clone of the
/// [`WriterThreaded`].
///
/// [`flush`]: WriterThreaded::flush
/// [`write`]: WriterThreaded::write
/// [`mpsc`]: std::sync::mpsc
#[derive(Debug)]
pub struct WriterThreaded {
    /// Declared before `shared` so pending bytes are sent before the last clone
    /// finishes the writer thread.
    batch: Option<PendingBatch>,

    /// The shared state between [`WriterThreaded`] clones.
    ///
    /// This enables a proper [`Finish`] implementation, since we can check
    /// whether a single writer clone remains and we can store the
    /// [`JoinHandle`].
    shared: Arc<WriterThreadedShared>,

    /// The first IO error produced by the underlying writer in a
    /// commonly-accessible location. The outer [`Arc`] ensures that the memory
    /// is shared between all the threads. The [`OnceLock`] acts as an `Option`
    /// that is written to once. The [`SharedIoError`] contains another [`Arc`]
    /// to allow the inner error to be clone-able.
    ///
    /// This is not in [`WriterThreadedShared`] since the thread cannot be
    /// formed until the error is wrapped in [`Arc`] (so that a clone can be
    /// moved into the thread), but the [`Arc`] around [`WriterThreadedShared`]
    /// cannot be formed until the thread handle is available.
    writer_error: Arc<OnceLock<SharedIoError>>,
}

#[derive(Debug)]
struct PendingBatch {
    bytes:    Vec<u8>,
    capacity: usize,
    // Kept here because WriterThreaded cannot implement Drop while finish moves shared out.
    sender:   std::sync::mpsc::Sender<Msg>,
}

const BATCH_HEADROOM: usize = 1024;

impl PendingBatch {
    fn buffer(capacity: usize) -> Vec<u8> {
        let reserve = if capacity == 0 {
            0
        } else {
            capacity.checked_add(BATCH_HEADROOM).expect("batch capacity overflow")
        };
        Vec::with_capacity(reserve)
    }
}

impl Drop for PendingBatch {
    fn drop(&mut self) {
        if !self.bytes.is_empty() {
            let _ = self.sender.send(Msg::Data(std::mem::take(&mut self.bytes)));
        }
    }
}

impl Clone for WriterThreaded {
    fn clone(&self) -> Self {
        Self {
            shared:       self.shared.clone(),
            writer_error: self.writer_error.clone(),
            batch:        self.batch.as_ref().map(|batch| PendingBatch {
                bytes:    PendingBatch::buffer(batch.capacity),
                capacity: batch.capacity,
                sender:   batch.sender.clone(),
            }),
        }
    }
}

/// The inner fields of a [`WriterThreaded`].
#[derive(Debug)]
struct WriterThreadedShared {
    /// The sending portion of the channel.
    sender: std::sync::mpsc::Sender<Msg>,

    /// A handle to the thread, used to enable having both a [`Finish`] and
    /// [`Drop`] implementation without using unsafe code like [`ManuallyDrop`].
    ///
    /// This only becomes `None` within [`drop`] or [`finish`]. Otherwise, this
    /// is guaranteed to be `Some`.
    ///
    /// [`finish`]: Finish::finish
    /// [`ManuallyDrop`]: std::mem::ManuallyDrop
    handle: Option<JoinHandle<std::io::Result<()>>>,
}

enum Msg {
    /// A line of data to be written to the underlying writer.
    Data(Vec<u8>),
    /// A request to flush the underlying writer and all previously received
    /// messages in the channel. The sender is used to communicate when the
    /// flush has been successfully completed.
    Flush(std::sync::mpsc::Sender<()>),
    /// A request to finish the underlying writer and terminate the thread.
    Finish,
}

impl WriterThreaded {
    /// Enables per-clone byte batching. Each batch reserves at least `capacity`
    /// bytes and is sent once it reaches `capacity` bytes; remaining bytes are
    /// sent on flush, finish, or drop. A zero `capacity` reserves nothing and
    /// sends every write immediately. Unlike [`new`](Self::new), writes copy into
    /// the batch.
    #[must_use]
    pub fn with_batch_capacity<W>(writer: W, capacity: usize) -> Self
    where
        W: Write + Finish + Send + 'static, {
        let mut threaded = Self::new(writer);
        threaded.batch = Some(PendingBatch {
            bytes: PendingBatch::buffer(capacity),
            capacity,
            sender: threaded.shared.sender.clone(),
        });
        threaded
    }

    fn writer_failure(&self) -> std::io::Error {
        self.writer_error.get().map_or_else(
            || std::io::Error::new(std::io::ErrorKind::BrokenPipe, "writer thread is no longer running"),
            Into::into,
        )
    }

    fn send_batch(&mut self) -> std::io::Result<()> {
        let Some(batch) = self.batch.as_mut() else { return Ok(()) };
        if batch.bytes.is_empty() {
            return Ok(());
        }
        let bytes = std::mem::replace(&mut batch.bytes, PendingBatch::buffer(batch.capacity));
        self.shared
            .sender
            .send(Msg::Data(bytes))
            .map_err(|SendError(_)| self.writer_failure())
    }

    fn append_batch(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if self
            .batch
            .as_ref()
            .is_some_and(|batch| !batch.bytes.is_empty() && batch.bytes.len().saturating_add(buf.len()) > batch.capacity)
        {
            self.send_batch()?;
        }
        let batch = self.batch.as_mut().expect("batching is enabled");
        batch.bytes.extend_from_slice(buf);
        if batch.bytes.len() >= batch.capacity {
            self.send_batch()?;
        }
        Ok(buf.len())
    }

    /// Writes an owned buffer without copying when batching is disabled.
    pub fn write_vec(&mut self, buf: Vec<u8>) -> std::io::Result<usize> {
        let len = buf.len();

        if self.batch.is_some() {
            return self.append_batch(&buf);
        }

        match self.shared.sender.send(Msg::Data(buf)) {
            Ok(()) => Ok(len),
            Err(SendError(_)) => Err(self.writer_failure()),
        }
    }
}

impl Write for WriterThreaded {
    /// Writes a buffer into the writer, returning how many bytes were queued
    /// for writing.
    ///
    /// Without batching, this allocates a vector; use [`WriterThreaded::write_vec`]
    /// to avoid copying an existing vector. With batching, bytes are copied into
    /// the batch instead.
    fn write(&mut self, buf: &[u8]) -> std::io::Result<usize> {
        if self.batch.is_some() {
            return self.append_batch(buf);
        }
        self.write_vec(buf.to_vec())
    }

    /// Flushes the contents of the queue, which includes all writes made by
    /// this clone before this call was made. The relative ordering between this
    /// flush call and other clones of [`WriterThreaded`] cannot be guaranteed.
    /// This function will block until the flush is successful.
    ///
    /// To ensure that errors are not silently ignored, for each linked
    /// [`WriterThreaded`], ensure that either [`flush`] is called on it, or it
    /// is dropped before a [`flush`] call on a different clone of the
    /// [`WriterThreaded`]. This is similar to a [`BufWriter`], where failing to
    /// call [`flush`] may result in errors getting silently ignored.
    ///
    /// [`BufWriter`]: std::io::BufWriter
    /// [`flush`]: WriterThreaded::flush
    fn flush(&mut self) -> std::io::Result<()> {
        self.send_batch()?;
        let (sender, receiver) = std::sync::mpsc::channel();

        match self.shared.sender.send(Msg::Flush(sender)) {
            Ok(()) => match receiver.recv() {
                Ok(()) => Ok(()),
                Err(RecvError) => Err(self.writer_failure()),
            },
            Err(SendError(_)) => Err(self.writer_failure()),
        }
    }

    /// Writes a buffer into the writer.
    ///
    /// Without batching, this allocates a vector; use [`WriterThreaded::write_vec`]
    /// to avoid copying an existing vector. With batching, bytes are copied into
    /// the batch instead.
    fn write_all(&mut self, buf: &[u8]) -> std::io::Result<()> {
        self.write(buf)?;
        Ok(())
    }

    fn write_fmt(&mut self, args: std::fmt::Arguments<'_>) -> std::io::Result<()> {
        if let Some(batch) = self.batch.as_mut() {
            Write::write_fmt(&mut batch.bytes, args)?;
            if batch.bytes.len() >= batch.capacity {
                self.send_batch()?;
            }
            return Ok(());
        }
        // We use write_vec to ensure only a single allocation is performed
        if let Some(s) = args.as_str() {
            self.write_vec(s.as_bytes().to_vec()).map(|_| ())
        } else {
            self.write_vec(std::fmt::format(args).into_bytes()).map(|_| ())
        }
    }
}

impl WriterThreaded {
    /// Constructs a [`WriterThreaded`] from a regular writer by moving it into
    /// a thread and creating a channel.
    #[inline]
    #[must_use]
    pub fn new<W>(mut writer: W) -> Self
    where
        W: Write + Finish + Send + 'static, {
        let (sender, receiver) = std::sync::mpsc::channel();

        let writer_error = Arc::new(OnceLock::new());

        let thread_writer_error = writer_error.clone();

        // We do not bind the thread handle, but the thread will still run to
        // completion (it will run until an error, or until all senders are
        // dropped)
        let handle = std::thread::spawn(move || -> std::io::Result<()> {
            while let Ok(msg) = receiver.recv() {
                match msg {
                    Msg::Data(bytes) => {
                        let res = writer.write_all(&bytes);
                        if let Err(err) = res {
                            let err: SharedIoError = err.into();
                            thread_writer_error.get_or_init(|| err.clone());
                            return Err(err.into());
                        }
                    }
                    Msg::Flush(done) => {
                        let res = writer.flush();
                        if let Err(err) = res {
                            let err: SharedIoError = err.into();
                            thread_writer_error.get_or_init(|| err.clone());
                            return Err(err.into());
                        }

                        // Signal to the thread that called flush that the flush
                        // has now been completed. This should never fail since
                        // the thread calls `recv` immediately after sending the
                        // Msg, but if it does for some pathological case, we
                        // may as well poison as much as possible.
                        let res = done.send(());
                        if let Err(SendError(())) = res {
                            let err: SharedIoError = std::io::Error::other("Failed to confirm successful flush").into();
                            thread_writer_error.get_or_init(|| err.clone());
                            return Err(err.into());
                        }
                    }
                    Msg::Finish => return writer.finish(),
                }
            }

            // This will be reached only after all senders have hung up. The
            // final WriterThreaded clone should call finish when it gets
            // dropped, which will cause the thread to return before reaching
            // this point. Hence, this is likely unreachable but to be safe we
            // include it.
            writer.finish()
        });

        Self {
            shared: Arc::new(WriterThreadedShared {
                sender,
                handle: Some(handle),
            }),
            writer_error,
            batch: None,
        }
    }
}

impl Drop for WriterThreadedShared {
    /// Drops the [`WriterThreadedShared`].
    ///
    /// [`flush`]: WriterThreaded::flush
    fn drop(&mut self) {
        let Some(handle) = self.handle.take() else {
            return;
        };

        let _ = self.sender.send(Msg::Finish);
        let _ = handle.join();
    }
}

impl Finish for WriterThreaded {
    fn finish(mut self) -> std::io::Result<()> {
        self.send_batch()?;
        // Ensure this is the final WriterThreaded clone, and bind the shared
        // state
        let mut shared = Arc::into_inner(self.shared)
            .ok_or(std::io::Error::other("Cannot finish WriterThreaded while clones still exist"))?;

        // Extract the join handle
        let Some(handle) = shared.handle.take() else {
            // This should not be reachable, but it means the thread was already
            // joined (and an error propagated), so we can just end
            return Ok(());
        };

        // Send a finish message. An Ok value indicates that either the message
        // was successfully received (which will cause finish to be called on
        // the inner writer), or that the receiver hung up immediately after
        // send returned Ok. The receiver could hang up due to a different error
        // (propagated when we join the thread handle), a panic (propagated via
        // unwrap on the join), or the process getting killed (in which case
        // this parent process will also be killed). So, we can safely ignore
        // the result.
        let _ = shared.sender.send(Msg::Finish);

        // Drop the final shared writer so that the thread can terminate
        drop(shared);

        // Join the handle, waiting for the thread to finish and propagating any
        // errors
        handle.join().unwrap()
    }
}

/// A [`std::io::Error`] for which [`Clone`] can be called (allowing many
/// threads to yield the same underlying error while preserving any source
/// errors).
#[derive(Clone, Debug)]
struct SharedIoError(Arc<std::io::Error>);

impl Display for SharedIoError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.0.fmt(f)
    }
}

impl Error for SharedIoError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        self.0.source()
    }
}

impl From<&SharedIoError> for std::io::Error {
    fn from(value: &SharedIoError) -> Self {
        std::io::Error::new(value.0.kind(), value.clone())
    }
}

impl From<SharedIoError> for std::io::Error {
    fn from(value: SharedIoError) -> Self {
        std::io::Error::new(value.0.kind(), value)
    }
}

impl From<std::io::Error> for SharedIoError {
    fn from(value: std::io::Error) -> Self {
        Self(Arc::new(value))
    }
}