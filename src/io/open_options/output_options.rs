use crate::io::{OptionalPaths, OutputContext, PairedErrors, RecordWriters, WriteFileZipStdout, open_options::PairedStruct};
use std::{fs::File, io::BufWriter, path::Path};

/// A builder pattern for creating output files in IRMA-core.
///
/// This supports many features, such as:
///
/// - Handling unpaired or paired inputs
/// - Interpreting the paths in multiple ways, such as [`File`] and
///   [`WriteFileZipStdout`]
/// - Automatically adding context including the path and record type (if
///   applicable) to any errors
/// - Altering the capacity of the [`BufWriter`]
///
/// To use this, use the following steps:
///
/// 1. Select the appropriate constructor:
///    - [`OutputOptions::new_from_path`]: This is used when a single output
///      file is being created from a path, compatible with [`File`].
///    - [`OutputOptions::new_from_opt_path`]: This is used when a single output
///      file is being created from an optional path, compatible with
///      [`WriteFileZipStdout`]. If the path is not provided, stdout is used.
///    - [`OutputOptions::new_from_paths`]: This is used when potentially paired
///      output files are being created. This is compatible with [`File`].
///    - [`OutputOptions::new_from_opt_paths`]: This is used when potentially
///      paired output files are being created, and the first path is optional.
///      This is compatible with [`WriteFileZipStdout`].
/// 2. Optionally set the capacity for the [`BufWriter`] to create using
///    `with_capacity`.
/// 3. Call a method to interpret the path as something readable. The options
///    may differ depending on the constructor used.
///    - `use_file`: Interpret the path as a regular file ([`File`])
///    - `use_file_zip_or_stdout`: Interpret the path as a regular file, zipped
///      file, or stdout if no path is provided ([`WriteFileZipStdout`])
/// 4. Call the `open` method to create the outputs, with context automatically
///    added to any errors.
///
/// See the uses of this in IRMA-core for examples.
///
/// When using [`FastXReader`], this builder will return an
/// `IterWithContext<FastXReader<...>>`, which works as an iterator but does not
/// facilitate matching on the FASTQ and FASTA variants. To dispatch based on
/// these variants, call [`dispatch`].
///
/// [`FastXReader`]: crate::io::fastx::FastXReader
/// [`dispatch`]: crate::io::IterWithContext::dispatch
pub struct OutputOptions<'a, W> {
    /// Any context needed to properly display error context.
    context: OutputContext<'a>,

    /// The output constructed so far, as a result. This can be a [`Path`],
    /// `Option<Path>`, [`File`], an enum such as [`WriteFileZipStdout`], or a
    /// pair of these with [`RecordWriters`] or [`OptionalPaths`].
    ///
    /// The [`Result`] is sometimes guaranteed to be [`Ok`], such as when `W` is
    /// [`Path`] or [`OptionalPaths`]. Similarly, the [`PairedErrors`] is
    /// guaranteed to be [`PairedErrors::Err1`] except when [`RecordWriters`] is
    /// used.
    ///
    /// By having a simpler user-facing type for `W` (that hides the underlying
    /// [`Result`]), we lose these type-state guarantees. However, this was
    /// deemed worthwhile.
    output: Result<W, PairedErrors>,

    /// The capacity for the [`BufWriter`] to create.
    capacity: Option<usize>,
}

impl<'a> OutputOptions<'a, &'a Path> {
    /// Creates a new [`OutputOptions`] from a specified path.
    ///
    /// This can then be interpreted as [`File`].
    #[allow(dead_code)]
    pub fn new_from_path<P>(path: &'a P) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        let path = path.as_ref();
        Self {
            context:  OutputContext::new(Some(path), None),
            output:   Ok(path),
            capacity: None,
        }
    }

    /// Sets the capacity for the [`BufWriter`] to use when opening the writer.
    pub fn with_capacity(mut self, capacity: usize) -> Self {
        self.capacity = Some(capacity);
        self
    }

    /// Interprets the path using [`File`] for writing.
    #[allow(dead_code)]
    pub fn use_file(self) -> OutputOptions<'a, BufWriter<File>> {
        let output = self.output.and_then(|path| {
            if let Some(capacity) = self.capacity {
                File::create(path).map(|file| BufWriter::with_capacity(capacity, file))
            } else {
                File::create(path).map(BufWriter::new)
            }
            .map_err(PairedErrors::Err1)
        });

        OutputOptions {
            context: self.context,
            output,
            capacity: self.capacity,
        }
    }
}

impl<'a> OutputOptions<'a, Option<&'a Path>> {
    /// Creates a new [`OutputOptions`] from an optional path.
    ///
    /// This can then be interpreted as a [`WriteFileZipStdout`].
    pub fn new_from_opt_path<P>(path: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        let path = path.map(AsRef::as_ref);
        Self {
            context:  OutputContext::new(path, None),
            output:   Ok(path),
            capacity: None,
        }
    }

    /// Sets the capacity for the [`BufWriter`] to use when opening the writer.
    #[allow(dead_code)]
    pub fn with_capacity(mut self, capacity: usize) -> Self {
        self.capacity = Some(capacity);
        self
    }

    /// Interprets the optional path using [`WriteFileZipStdout`], which
    /// supports regular files, [gzip
    /// files](https://www.rfc-editor.org/rfc/rfc1952#page-5), and stdout (in
    /// the case that no path was provided).
    ///
    /// The output is zipped if the path ends in `.gz`.
    pub fn use_file_zip_or_stdout(self) -> OutputOptions<'a, WriteFileZipStdout> {
        let output = self.output.and_then(|path| {
            if let Some(capacity) = self.capacity {
                WriteFileZipStdout::with_capacity(capacity, path)
            } else {
                WriteFileZipStdout::create(path)
            }
            .map_err(PairedErrors::Err1)
        });

        OutputOptions {
            context: self.context,
            output,
            capacity: self.capacity,
        }
    }
}

impl<'a> OutputOptions<'a, RecordWriters<&'a Path>> {
    /// Creates a new [`OutputOptions`] from potentially paired paths.
    ///
    /// This is designed to dynamically support both unpaired and paired
    /// behavior, with the decision made at runtime. For paired outputs, `path2`
    /// should be `Some`, and for an unpaired output, `path2` should be `None`.
    ///
    /// The paths can then be interpreted as [`File`].
    #[allow(dead_code)]
    pub fn new_from_paths<P>(path1: &'a P, path2: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        let path1 = path1.as_ref();
        let path2 = path2.map(AsRef::as_ref);
        Self {
            context:  OutputContext::new(Some(path1), path2),
            output:   Ok(RecordWriters::new(path1, path2)),
            capacity: None,
        }
    }

    /// Sets the capacity for the [`BufWriter`] to use when opening the writer.
    #[allow(dead_code)]
    pub fn with_capacity(mut self, capacity: usize) -> Self {
        self.capacity = Some(capacity);
        self
    }

    /// Interprets the path(s) using [`File`] for writing.
    #[allow(dead_code)]
    pub fn use_file(self) -> OutputOptions<'a, RecordWriters<BufWriter<File>>> {
        OutputOptions {
            context:  self.context,
            output:   self.output.and_then(|writers| {
                writers.try_map(|path| {
                    if let Some(capacity) = self.capacity {
                        File::create(path).map(|file| BufWriter::with_capacity(capacity, file))
                    } else {
                        File::create(path).map(BufWriter::new)
                    }
                })
            }),
            capacity: self.capacity,
        }
    }
}

impl<'a> OutputOptions<'a, OptionalPaths<'a>> {
    /// Creates a new [`OutputOptions`] from two optional, potentially paired
    /// paths.
    ///
    /// This is designed to dynamically support both unpaired and paired
    /// behavior, with the decision made at runtime. For paired outputs, `path2`
    /// should be `Some`, and for an unpaired output, `path2` should be `None`.
    ///
    /// The paths can then be interpreted as [`WriteFileZipStdout`]. Only
    /// `path1` has the potential of being [`WriteFileZipStdout::Stdout`], since
    /// if `path2` is `None`, this corresponds to unpaired output.
    pub fn new_from_opt_paths<P>(path1: Option<&'a P>, path2: Option<&'a P>) -> Self
    where
        P: AsRef<Path> + ?Sized, {
        let path1 = path1.map(AsRef::as_ref);
        let path2 = path2.map(AsRef::as_ref);
        Self {
            context:  OutputContext::new(path1, path2),
            output:   Ok(OptionalPaths { path1, path2 }),
            capacity: None,
        }
    }

    /// Sets the capacity for the [`BufWriter`] to use when opening the writer.
    #[allow(dead_code)]
    pub fn with_capacity(mut self, capacity: usize) -> Self {
        self.capacity = Some(capacity);
        self
    }

    /// Interprets the optional path(s) using [`WriteFileZipStdout`].
    ///
    /// Only `path1` has the potential of being [`WriteFileZipStdout::Stdout`],
    /// since if `path2` is `None`, this corresponds to unpaired output.
    pub fn use_file_zip_or_stdout(self) -> OutputOptions<'a, RecordWriters<WriteFileZipStdout>> {
        OutputOptions {
            context:  self.context,
            output:   self.output.and_then(|paths| {
                paths.try_map_writers(|path| {
                    if let Some(capacity) = self.capacity {
                        WriteFileZipStdout::with_capacity(capacity, path)
                    } else {
                        WriteFileZipStdout::create(path)
                    }
                })
            }),
            capacity: self.capacity,
        }
    }
}

impl<'a> OutputOptions<'a, BufWriter<File>> {
    /// Opens the [`File`] for writing, wrapping it in a [`BufWriter`].
    ///
    /// ## Errors
    ///
    /// IO errors when opening the file are propagated. Context is added that
    /// includes the path.
    pub fn open(self) -> std::io::Result<BufWriter<File>> {
        match self.output {
            Ok(writer) => Ok(writer),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}

impl<'a> OutputOptions<'a, WriteFileZipStdout> {
    /// Opens the [`WriteFileZipStdout`].
    ///
    /// ## Errors
    ///
    /// If a path was provided, IO errors when creating the file are propagated.
    /// Context is added that includes the path.
    pub fn open(self) -> std::io::Result<WriteFileZipStdout> {
        match self.output {
            Ok(writer) => Ok(writer),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}

impl<'a> OutputOptions<'a, RecordWriters<BufWriter<File>>> {
    /// Opens the potentially paired [`File`]s for writing, wrapping each in a
    /// [`BufWriter`].
    ///
    /// ## Errors
    ///
    /// IO errors when creating the files are propagated. Context is added that
    /// includes the path.
    #[allow(dead_code)]
    pub fn open(self) -> std::io::Result<RecordWriters<BufWriter<File>>> {
        match self.output {
            Ok(writer) => Ok(writer),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}

impl<'a> OutputOptions<'a, RecordWriters<WriteFileZipStdout>> {
    /// Opens the potentially paired [`WriteFileZipStdout`] outputs.
    ///
    /// ## Errors
    ///
    /// If a path was provided for the first input, IO errors when creating the
    /// file are propagated. Context is added that includes the path.
    pub fn open(self) -> std::io::Result<RecordWriters<WriteFileZipStdout>> {
        match self.output {
            Ok(writer) => Ok(writer),
            Err(e) => Err(self.context.add_context(e).into()),
        }
    }
}
