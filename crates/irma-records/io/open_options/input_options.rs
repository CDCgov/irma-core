use crate::io::{
    FastXReader, GzipReaderInThread, InputContext, IterWithContext, IterWithErrorContext, OptionalPaths, PairedErrors,
    ReadFileStdin, ReadFileZip, ReadFileZipInThread, ReaderType, ReaderWithContext, RecordReaders,
    open_options::PairedStruct,
};
use std::{
    fs::File,
    io::{BufReader, Read, Stdin, stdin},
    path::Path,
};
use zoe::{
    data::sam::SAMReader,
    prelude::{FastQReader, FastaReader},
};

/// A builder pattern for opening input files in IRMA-core.
///
/// This supports many features, such as:
///
/// - Handling unpaired or paired inputs
/// - Interpreting the paths in multiple ways, such as [`File`],
///   [`ReadFileZip`], [`ReadFileZipInThread`], and [`ReadFileStdin`]
/// - Parsing the inputs into [`FastQReader`], [`FastaReader`], [`FastXReader`],
///   or [`SAMReader`]
/// - Automatically adding context including the path and record type (if
///   applicable) to any errors while opening the file or reading from the file
/// - Automatically adding context including the path and record type to each
///   item of the fallible iterator (for [`FastQReader`], [`FastaReader`], etc.)
///
/// To use this, use the following steps:
///
/// 1. Select the appropriate constructor:
///    - [`InputOptions::new_from_path`]: This is used when a single input file
///      is being opened from a path, compatible with [`File`], [`ReadFileZip`],
///      and [`ReadFileZipInThread`].
///    - [`InputOptions::new_from_opt_path`]: This is used when a single input
///      file is being opened from an optional path, compatible with
///      [`ReadFileStdin`]. If the path is not provided, stdin is used.
///    - [`InputOptions::new_from_paths`]: This is used when potentially paired
///      input files are being opened. This is compatible with [`File`],
///      [`ReadFileZip`], and [`ReadFileZipInThread`].
///    - [`InputOptions::new_from_opt_paths`]: This is used when potentially
///      paired input files are being opened, and the first path is optional.
///      This is compatible with [`ReadFileStdin`].
/// 2. Call a method to interpret the path as something readable. The options
///    may differ depending on the constructor used.
///    - `use_file`: Interpret the path as a regular file ([`File`])
///    - `use_file_or_zip`: Interpret the path as a regular or zipped file
///      ([`ReadFileZip`])
///    - `use_file_or_stdin`: Interprets the optional path as a regular file or
///      stdin if no path is provided ([`ReadFileStdin`])
/// 3. If using `use_file_or_zip`, optionally specify that the decoding should
///    happen eagerly on a separate thread with `decode_in_thread`. This is
///    useful when the contents are processed in a streamed manner (i.e., not
///    collected).
/// 4. Optionally parse the input as a particular record format:
///    - `parse_fastq`: Parses the input as a FASTQ file ([`FastQReader`])
///    - `parse_fasta`: Parses the input as a FASTA file ([`FastaReader`])
///    - `parse_fastx`: Parses the input as either a FASTQ or FASTA file
///      ([`FastXReader`])
///    - `parse_sam`: Parses the input as a SAM file ([`SAMReader`])
/// 5. Call the `open` method to retrieve the inputs, with context automatically
///    added to any errors.
///
/// See the uses of this in IRMA-core for examples.
///
/// Some tips for using this effectively:
///
/// - When using [`FastXReader`], this builder will return an
///   `IterWithContext<FastXReader<...>>`, which works as an iterator but does
///   not facilitate matching on the FASTQ and FASTA variants. To dispatch based
///   on these variants, call [`dispatch`].
/// - The final reader or iterator returned has context fully added, for errors
///   when opening, iterating, or reading. Hence, it is good to avoid adding
///   redundant context. For example, when using this to open a file and then
///   perform some type of custom parsing on it, it is best to encapsulate the
///   parsing logic in a separate function (e.g., with traits [`FromStr`] or
///   [`TryFrom`]). The errors within the parsing function should include
///   descriptions of the parsing error, and path context can be added when the
///   parsing function is called.
///
/// [`FromStr`]: std::str::FromStr
/// [`dispatch`]: crate::io::IterWithContext::dispatch
pub struct InputOptions<'a, R> {
    /// Any context needed to properly display error context.
    context: InputContext<'a>,

    /// The input constructed so far, as a result. This can be a [`Path`],
    /// `Option<Path>`, [`File`], an enum such as [`ReadFileZip`], an iterator
    /// such as [`FastQReader`], or a pair of these with [`RecordReaders`] or
    /// [`OptionalPaths`].
    ///
    /// The [`Result`] is sometimes guaranteed to be [`Ok`], such as when `R` is
    /// [`Path`] or [`OptionalPaths`]. Similarly, the [`PairedErrors`] is
    /// guaranteed to be [`PairedErrors::Err1`] except when [`RecordReaders`] is
    /// used.
    ///
    /// By having a simpler user-facing type for `R` (that hides the underlying
    /// [`Result`]), we lose these type-state guarantees. However, this was
    /// deemed worthwhile.
    input: Result<R, PairedErrors>,
}

impl<'a> InputOptions<'a, &'a Path> {
    /// Creates a new [`InputOptions`] from a specified path.
    ///
    /// This can then be interpreted as [`File`], [`ReadFileZip`], or
    /// [`ReadFileZipInThread`].
    pub fn new_from_path<P>(path: &'a P) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        let path = path.as_ref();
        Self {
            context: InputContext::new(Some(path), None),
            input:   Ok(path),
        }
    }

    /// Interprets the path using [`File`] for reading.
    pub fn use_file(self) -> InputOptions<'a, File> {
        InputOptions {
            context: self.context,
            input:   self.input.and_then(|path| File::open(path).map_err(PairedErrors::Err1)),
        }
    }

    /// Interprets the path using [`ReadFileZip`], which supports regular files
    /// and [gzip files](https://www.rfc-editor.org/rfc/rfc1952#page-5).
    ///
    /// The file is determined to be zipped if the path ends in `.gz`.
    pub fn use_file_or_zip(self) -> InputOptions<'a, ReadFileZip> {
        InputOptions {
            context: self.context,
            input:   self
                .input
                .and_then(|path| ReadFileZip::open(path).map_err(PairedErrors::Err1)),
        }
    }
}

impl<'a> InputOptions<'a, Stdin> {
    /// Creates a new [`InputOptions`] for reading from [`Stdin`].
    pub fn new_stdin() -> Self {
        Self {
            context: InputContext::default(),
            input:   Ok(stdin()),
        }
    }
}

impl<'a> InputOptions<'a, Option<&'a Path>> {
    /// Creates a new [`InputOptions`] from an optional path.
    ///
    /// This can then be interpreted as a [`ReadFileStdin`].
    pub fn new_from_opt_path<P>(path: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        let path = path.map(AsRef::as_ref);
        Self {
            context: InputContext::new(path, None),
            input:   Ok(path),
        }
    }

    /// Interprets the optional path using [`ReadFileStdin`], which supports
    /// regular files and stdin (in the case that no path was provided).
    #[allow(dead_code)]
    pub fn use_file_or_stdin(self) -> InputOptions<'a, ReadFileStdin> {
        InputOptions {
            context: self.context,
            input:   self
                .input
                .and_then(|path| ReadFileStdin::open(path).map_err(PairedErrors::Err1)),
        }
    }
}

impl<'a> InputOptions<'a, RecordReaders<&'a Path>> {
    /// Creates a new [`InputOptions`] from potentially paired paths.
    ///
    /// This is designed to dynamically support both unpaired and paired
    /// behavior, with the decision made at runtime. For paired inputs, `path2`
    /// should be `Some`, and for an unpaired input, `path2` should be `None`.
    ///
    /// The paths can then be interpreted as [`File`], [`ReadFileZip`], or
    /// [`ReadFileZipInThread`].
    pub fn new_from_paths<P>(path1: &'a P, path2: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        let path1 = path1.as_ref();
        let path2 = path2.map(AsRef::as_ref);
        Self {
            context: InputContext::new(Some(path1), path2),
            input:   Ok(RecordReaders {
                reader1: path1,
                reader2: path2,
            }),
        }
    }

    /// Interprets the path(s) using [`File`] for reading.
    pub fn use_file(self) -> InputOptions<'a, RecordReaders<File>> {
        InputOptions {
            context: self.context,
            input:   self.input.and_then(|readers| readers.try_map(File::open)),
        }
    }

    /// Interprets the path(s) using [`ReadFileZip`], which supports regular
    /// files and [gzip files](https://www.rfc-editor.org/rfc/rfc1952#page-5).
    ///
    /// Each file is determined to be zipped if the path end in `.gz`.
    pub fn use_file_or_zip(self) -> InputOptions<'a, RecordReaders<ReadFileZip>> {
        InputOptions {
            context: self.context,
            input:   self.input.and_then(|readers| readers.try_map(ReadFileZip::open)),
        }
    }
}

impl<'a> InputOptions<'a, OptionalPaths<'a>> {
    /// Creates a new [`InputOptions`] from two optional, potentially paired
    /// paths.
    ///
    /// This is designed to dynamically support both unpaired and paired
    /// behavior, with the decision made at runtime. For paired inputs, `path2`
    /// should be `Some`, and for an unpaired input, `path2` should be `None`.
    ///
    /// The paths can then be interpreted as [`ReadFileStdin`]. Only `path1` has
    /// the potential of being [`ReadFileStdin::Stdin`], since if `path2` is
    /// `None`, this corresponds to unpaired input.
    pub fn new_from_opt_paths<P1, P2>(path1: Option<&'a P1>, path2: Option<&'a P2>) -> Self
    where
        P1: AsRef<Path> + ?Sized,
        P2: AsRef<Path> + ?Sized, {
        let path1 = path1.map(AsRef::as_ref);
        let path2 = path2.map(AsRef::as_ref);
        Self {
            context: InputContext::new(path1, path2),
            input:   Ok(OptionalPaths { path1, path2 }),
        }
    }

    /// Interprets the optional path(s) using [`ReadFileStdin`].
    ///
    /// Only `path1` has the potential of being [`ReadFileStdin::Stdin`], since
    /// if `path2` is `None`, this corresponds to unpaired input.
    pub fn use_file_or_stdin(self) -> InputOptions<'a, RecordReaders<ReadFileStdin>> {
        InputOptions {
            context: self.context,
            input:   self.input.and_then(|paths| paths.try_map_readers(ReadFileStdin::open)),
        }
    }
}

impl<'a> InputOptions<'a, ReadFileZip> {
    /// Performs decoding of a GZIP file in a separate thread and sends them
    /// back using an anonymous pipe, which can speed up streamed processing of
    /// the contents.
    ///
    /// If the contents are collected before the program continues, this is
    /// unlikely to add benefit.
    pub fn decode_in_thread(self) -> InputOptions<'a, ReadFileZipInThread> {
        InputOptions {
            context: self.context,
            input:   self.input.and_then(|input| match input {
                ReadFileZip::File(file) => Ok(ReadFileZipInThread::File(file)),
                ReadFileZip::Zipped(decoder) => Ok(ReadFileZipInThread::Zipped(
                    GzipReaderInThread::from_decoder(decoder).map_err(PairedErrors::Err1)?,
                )),
            }),
        }
    }
}

impl<'a> InputOptions<'a, RecordReaders<ReadFileZip>> {
    /// Performs decoding of a GZIP file in a separate thread and sends them
    /// back using an anonymous pipe, which can speed up streamed processing of
    /// the contents.
    ///
    /// If the contents are collected before the program continues, this is
    /// unlikely to add benefit.
    pub fn decode_in_thread(self) -> InputOptions<'a, RecordReaders<ReadFileZipInThread>> {
        let srcs = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        let input = srcs.try_map(|input| match input {
            ReadFileZip::File(file) => Ok(ReadFileZipInThread::File(file)),
            ReadFileZip::Zipped(decoder) => Ok(ReadFileZipInThread::Zipped(GzipReaderInThread::from_decoder(decoder)?)),
        });

        InputOptions {
            context: self.context,
            input,
        }
    }
}

impl<'a, R> InputOptions<'a, R>
where
    R: Read,
{
    /// Parses the input as a FASTQ file, via the iterator [`FastQReader`].
    pub fn parse_fastq(self) -> InputOptions<'a, FastQReader<R>> {
        let src = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        InputOptions {
            context: self.context.with_reader1(ReaderType::FastQ),
            input:   FastQReader::from_readable(src).map_err(PairedErrors::Err1),
        }
    }

    /// Parses the input as a FASTA file, via the iterator [`FastaReader`].
    pub fn parse_fasta(self) -> InputOptions<'a, FastaReader<R>> {
        let src = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        InputOptions {
            context: self.context.with_reader1(ReaderType::Fasta),
            input:   FastaReader::from_readable(src).map_err(PairedErrors::Err1),
        }
    }

    /// Parses the input as either a FASTQ file or a FASTA file, via the
    /// iterator [`FastXReader`].
    pub fn parse_fastx(mut self) -> InputOptions<'a, FastXReader<R>> {
        let src = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        let src = FastXReader::from_readable(src);

        if let Ok(reader) = &src {
            match reader {
                FastXReader::Fastq(_) => self.context = self.context.with_reader1(ReaderType::FastQ),
                FastXReader::Fasta(_) => self.context = self.context.with_reader1(ReaderType::Fasta),
            }
        } else {
            self.context = self.context.with_reader1(ReaderType::FastX);
        }

        InputOptions {
            context: self.context,
            input:   src.map_err(PairedErrors::Err1),
        }
    }

    /// Parses the input as a SAM file, via the iterator [`SAMReader`].
    pub fn parse_sam(self) -> InputOptions<'a, SAMReader<R, true>> {
        let src = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        InputOptions {
            context: self.context.with_reader1(ReaderType::Sam),
            input:   SAMReader::from_readable(src).map_err(PairedErrors::Err1),
        }
    }
}

impl<'a, R> InputOptions<'a, RecordReaders<R>>
where
    R: Read,
{
    /// Parses the input(s) as FASTQ files, via the iterator [`FastQReader`].
    pub fn parse_fastq(self) -> InputOptions<'a, RecordReaders<FastQReader<R>>> {
        let srcs = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        let input = srcs.try_map(FastQReader::from_readable);

        InputOptions {
            context: self.context.with_reader1(ReaderType::FastQ).with_reader2(ReaderType::FastQ),
            input,
        }
    }

    /// Parses the input(s) as FASTA file(s), via the iterator [`FastaReader`].
    #[allow(dead_code)]
    pub fn parse_fasta(self) -> InputOptions<'a, RecordReaders<FastaReader<R>>> {
        let srcs = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        let input = srcs.try_map(FastaReader::from_readable);

        InputOptions {
            context: self.context.with_reader1(ReaderType::Fasta).with_reader2(ReaderType::Fasta),
            input,
        }
    }

    /// Parses the input(s) as either FASTQ files or FASTA files, via the
    /// iterator [`FastXReader`].
    pub fn parse_fastx(mut self) -> InputOptions<'a, RecordReaders<FastXReader<R>>> {
        let srcs = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        let input = srcs.try_map(FastXReader::from_readable);

        if let Ok(readers) = &input {
            match readers.reader1 {
                FastXReader::Fastq(_) => self.context = self.context.with_reader1(ReaderType::FastQ),
                FastXReader::Fasta(_) => self.context = self.context.with_reader1(ReaderType::Fasta),
            }

            if let Some(reader2) = &readers.reader2 {
                match reader2 {
                    FastXReader::Fastq(_) => self.context = self.context.with_reader2(ReaderType::FastQ),
                    FastXReader::Fasta(_) => self.context = self.context.with_reader2(ReaderType::Fasta),
                }
            }
        } else {
            self.context = self.context.with_reader1(ReaderType::FastX).with_reader2(ReaderType::FastX);
        }

        InputOptions {
            context: self.context,
            input,
        }
    }

    /// Parses the input(s) as SAM files, via the iterator [`SAMReader`].
    pub fn parse_sam(self) -> InputOptions<'a, RecordReaders<SAMReader<R, true>>> {
        let srcs = match self.input {
            Ok(src) => src,
            Err(e) => {
                return InputOptions {
                    context: self.context,
                    input:   Err(e),
                };
            }
        };

        let input = srcs.try_map(SAMReader::from_readable);

        InputOptions {
            context: self.context.with_reader1(ReaderType::Sam).with_reader2(ReaderType::Sam),
            input,
        }
    }
}

impl<'a, R> InputOptions<'a, R>
where
    R: Read,
{
    /// A helper function for opening something implementing [`Read`].
    ///
    /// For opening an iterator ([`FastQReader`], [`FastaReader`], etc.), use
    /// `open_iter`.
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type if a `parse` method was
    /// called and the error originated during parsing. Any failed reads from
    /// the reader will also have similar context due to the
    /// [`ReaderWithContext`] wrapper.
    fn open_readable(self) -> std::io::Result<ReaderWithContext<BufReader<R>>> {
        match self.input {
            Ok(reader) => Ok(InputContext::add_reader_context(BufReader::new(reader), self.context.input1)),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}

impl<'a, R> InputOptions<'a, RecordReaders<R>>
where
    R: Read,
{
    /// A helper function for opening [`RecordReaders`] of something
    /// implementing [`Read`].
    ///
    /// For opening an iterator ([`FastQReader`], [`FastaReader`], etc.), use
    /// `open_iter`.
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type if a `parse` method was
    /// called and the error originated during parsing. Any failed reads from
    /// the readers will also have similar context due to the
    /// [`ReaderWithContext`] wrapper.
    fn open_readable(self) -> std::io::Result<RecordReaders<ReaderWithContext<R>>> {
        match self.input {
            Ok(readers) => Ok(self.context.add_paired_reader_context(readers)),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}

impl<'a, I> InputOptions<'a, I>
where
    I: IterWithErrorContext,
{
    /// A helper function for opening an iterator ([`FastQReader`],
    /// [`FastaReader`], [`FastXReader`], or [`SAMReader`]).
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available and the record type if the error originated during
    /// parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    fn open_iter(self) -> std::io::Result<IterWithContext<I>> {
        match self.input {
            Ok(reader) => Ok(InputContext::add_iter_context(
                reader,
                self.context.reader1,
                self.context.input1,
            )),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}

impl<'a, I> InputOptions<'a, RecordReaders<I>>
where
    I: IterWithErrorContext,
{
    /// A helper function for opening a [`RecordReaders`] of iterators
    /// ([`FastQReader`], [`FastaReader`], [`FastXReader`], or [`SAMReader`]).
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available and the record type if the error originated during
    /// parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    fn open_iter(self) -> std::io::Result<RecordReaders<IterWithContext<I>>> {
        match self.input {
            Ok(readers) => Ok(self.context.add_paired_iter_context(readers)),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}

impl InputOptions<'_, File> {
    /// Opens the [`File`] for reading, wrapping it in a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// IO errors when opening the file are propagated. Context is added that
    /// includes the path. Any failed reads from the reader will also have
    /// similar context due to the [`ReaderWithContext`] wrapper.
    pub fn open(self) -> std::io::Result<ReaderWithContext<BufReader<File>>> {
        self.open_readable()
    }
}

impl InputOptions<'_, ReadFileZip> {
    /// Opens the [`ReadFileZip`], wrapping it in a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// IO errors when opening the file are propagated. Context is added that
    /// includes the path. Any failed reads from the reader will also have
    /// similar context due to the [`ReaderWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<ReaderWithContext<BufReader<ReadFileZip>>> {
        self.open_readable()
    }
}

impl InputOptions<'_, ReadFileZipInThread> {
    /// Opens the [`ReadFileZipInThread`], wrapping it in a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// IO errors when opening the file or forming the pipe are propagated.
    /// Context is added that includes the path. Any failed reads from the
    /// reader will also have similar context due to the [`ReaderWithContext`]
    /// wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<ReaderWithContext<BufReader<ReadFileZipInThread>>> {
        self.open_readable()
    }
}

impl InputOptions<'_, ReadFileStdin> {
    /// Opens the [`ReadFileStdin`], wrapping it in a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// If a path was provided, IO errors when opening the file are propagated.
    /// Context is added that includes the path. Any failed reads from the
    /// reader will also have similar context due to the [`ReaderWithContext`]
    /// wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<ReaderWithContext<BufReader<ReadFileStdin>>> {
        self.open_readable()
    }
}

impl<R> InputOptions<'_, FastQReader<R>>
where
    R: Read,
{
    /// Opens the [`FastQReader`].
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type FASTQ if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<IterWithContext<FastQReader<R>>> {
        self.open_iter()
    }
}

impl<R> InputOptions<'_, FastaReader<R>>
where
    R: Read,
{
    /// Opens the [`FastaReader`].
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type FASTA if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    pub fn open(self) -> std::io::Result<IterWithContext<FastaReader<R>>> {
        self.open_iter()
    }
}

impl<R> InputOptions<'_, FastXReader<R>>
where
    R: Read,
{
    /// Opens the [`FastXReader`].
    ///
    /// The output (after handling the error case) can be directly used as an
    /// iterator, and context will be added since it is wrapped in
    /// [`IterWithContext`]. To instead match on whether the reader is a
    /// [`FastaReader`] or a [`FastQReader`], call [`IterWithContext::dispatch`]
    /// first to move the [`IterWithContext`] wrapper separately into each
    /// variant.
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    pub fn open(self) -> std::io::Result<IterWithContext<FastXReader<R>>> {
        self.open_iter()
    }
}

impl<R> InputOptions<'_, SAMReader<R, true>>
where
    R: Read,
{
    /// Opens the [`SAMReader`].
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type SAM if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<IterWithContext<SAMReader<R, true>>> {
        self.open_iter()
    }
}

impl InputOptions<'_, RecordReaders<File>> {
    /// Opens the potentially paired [`File`] inputs, wrapping each in a
    /// [`BufReader`].
    ///
    /// ## Errors
    ///
    /// IO errors when opening the files are propagated. Context is added that
    /// includes the path. Any failed reads from the readers will also have
    /// similar context due to the [`ReaderWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<RecordReaders<BufReader<ReaderWithContext<File>>>> {
        self.open_readable().map(|readers| readers.map(BufReader::new))
    }
}

impl InputOptions<'_, RecordReaders<ReadFileZip>> {
    /// Opens the potentially paired [`ReadFileZip`] inputs, wrapping each in a
    /// [`BufReader`].
    ///
    /// ## Errors
    ///
    /// IO errors when opening the files are propagated. Context is added that
    /// includes the path. Any failed reads from the readers will also have
    /// similar context due to the [`ReaderWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<RecordReaders<BufReader<ReaderWithContext<ReadFileZip>>>> {
        self.open_readable().map(|readers| readers.map(BufReader::new))
    }
}

impl InputOptions<'_, RecordReaders<ReadFileZipInThread>> {
    /// Opens the potentially paired [`ReadFileZipInThread`] inputs, wrapping
    /// each in a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// IO errors when opening the files or forming the pipe are propagated.
    /// Context is added that includes the path. Any failed reads from the
    /// readers will also have similar context due to the [`ReaderWithContext`]
    /// wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<RecordReaders<BufReader<ReaderWithContext<ReadFileZipInThread>>>> {
        self.open_readable().map(|readers| readers.map(BufReader::new))
    }
}

impl InputOptions<'_, RecordReaders<ReadFileStdin>> {
    /// Opens the potentially paired [`ReadFileStdin`] inputs, wrapping each in
    /// a [`BufReader`].
    ///
    /// ## Errors
    ///
    /// If a path was provided for the first input, IO errors when opening the
    /// file are propagated. Context is added that includes the path. Any failed
    /// reads from the readers will also have similar context due to the
    /// [`ReaderWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<RecordReaders<BufReader<ReaderWithContext<ReadFileStdin>>>> {
        self.open_readable().map(|readers| readers.map(BufReader::new))
    }
}

impl<R> InputOptions<'_, RecordReaders<FastQReader<R>>>
where
    R: Read,
{
    /// Opens the potentially paired [`FastQReader`] inputs.
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type FASTQ if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    pub fn open(self) -> std::io::Result<RecordReaders<IterWithContext<FastQReader<R>>>> {
        self.open_iter()
    }
}

impl<R> InputOptions<'_, RecordReaders<FastaReader<R>>>
where
    R: Read,
{
    /// Opens the potentially paired [`FastaReader`] inputs.
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type FASTA if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<RecordReaders<IterWithContext<FastaReader<R>>>> {
        self.open_iter()
    }
}

impl<R> InputOptions<'_, RecordReaders<FastXReader<R>>>
where
    R: Read,
{
    /// Opens the potentially paired [`FastXReader`] inputs.
    ///
    /// The output readers (after handling the error case) can be directly used
    /// as iterators, and context will be added since they are wrapped in
    /// [`IterWithContext`]. To instead match on whether each reader is a
    /// [`FastaReader`] or a [`FastQReader`], call [`IterWithContext::dispatch`]
    /// first on each reader to move the [`IterWithContext`] wrapper separately
    /// into each variant.
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    pub fn open(self) -> std::io::Result<RecordReaders<IterWithContext<FastXReader<R>>>> {
        self.open_iter()
    }
}

impl<R> InputOptions<'_, RecordReaders<SAMReader<R, true>>>
where
    R: Read,
{
    /// Opens the potentially paired [`SAMReader`] inputs.
    ///
    /// ## Errors
    ///
    /// Any IO or parsing error is propagated. Context is added that includes
    /// the path if available, and the record type SAM if the error originated
    /// during parsing. Any items that are errors in the iterator will also have
    /// similar context due to the [`IterWithContext`] wrapper.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<RecordReaders<IterWithContext<SAMReader<R, true>>>> {
        self.open_iter()
    }
}
