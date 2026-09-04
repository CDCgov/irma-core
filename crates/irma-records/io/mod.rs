use std::path::{Path, PathBuf};
use zoe::data::err::ResultWithErrorContext;

mod error_context;
mod fastx;
mod finish;
mod open_options;
mod readers;
mod write_records;
mod writers;

pub use error_context::*;
pub use fastx::*;
pub use finish::*;
pub use open_options::*;
pub use readers::*;
pub use write_records::*;
pub use writers::*;

/// Checks whether a path represents a device file, such as `/dev/null`.
fn is_linux_device(path: &Path) -> bool {
    path.starts_with("/dev/")
}

/// A trait for validating that input and output paths do not have conflicts.
pub trait ValidatePaths {
    /// Returns the paths that will be read from by the process.
    fn inputs(&self) -> impl IntoIterator<Item = &PathBuf>;

    /// Returns the paths that will be written to by the process.
    fn outputs(&self) -> impl IntoIterator<Item = &PathBuf>;

    /// Validates that no path is both an input and an output, and that all
    /// output paths are distinct.
    ///
    /// Device files (paths beginning with `/dev`/) are ignored.
    ///
    /// ## Errors
    ///
    /// All input paths must exist, and the parent directories of the output
    /// paths must exist. The paths must be successfully canonicalized. All
    /// output paths must be distinct and cannot also be input paths.
    fn validate_paths(&self) -> std::io::Result<()> {
        let inputs = self
            .inputs()
            .into_iter()
            .filter(|path| !is_linux_device(path))
            .map(|path| std::fs::canonicalize(path).with_path_context("Failed to canonicalize path", path));

        let outputs = self
            .outputs()
            .into_iter()
            .filter(|path| !is_linux_device(path))
            .map(|path| {
                // If the output path already exists (including as a symlink),
                // canonicalize it directly so that aliases such as symlinks
                // pointing at an input file are resolved to their real path.
                if path.exists() {
                    return Ok(std::fs::canonicalize(path).with_path_context("Failed to canonicalize path", path)?);
                }

                let filename = path.file_name().ok_or_else(|| {
                    std::io::Error::other(format!("Failed to find filename of path: {path}", path = path.display()))
                })?;
                let parent = match path.parent() {
                    Some(parent) if !parent.as_os_str().is_empty() => parent,
                    _ => Path::new("."),
                };
                let canonical_parent =
                    std::fs::canonicalize(parent).with_path_context("Failed to canonicalize parent path", parent)?;
                Ok(canonical_parent.join(filename))
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        for input1 in inputs {
            let input1 = input1?;

            for output in &outputs {
                if &input1 == output {
                    return Err(std::io::Error::other(format!(
                        "Found a file as both an input and an output: {input1}",
                        input1 = input1.display()
                    )));
                }
            }
        }

        for (i, output1) in outputs.iter().enumerate() {
            let rest = &outputs[i + 1..];
            for output2 in rest {
                if output1 == output2 {
                    return Err(std::io::Error::other(format!(
                        "Two output files were the same: {output1}",
                        output1 = output1.display()
                    )));
                }
            }
        }

        Ok(())
    }
}

/// Checks whether a file is a [gzip
/// file](https://www.rfc-editor.org/rfc/rfc1952#page-5).
///
/// This is currently done naively by seeing if it ends with a `gz` extension.
#[inline]
pub fn is_gz<P: AsRef<Path>>(path: P) -> bool {
    path.as_ref().extension().is_some_and(|ext| ext == "gz")
}
