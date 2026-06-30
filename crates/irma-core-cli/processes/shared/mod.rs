use jiff::Zoned;
use std::error::Error;
use zoe::{data::err::DisplayErrStack, search::ByteSubstringMut};

pub mod trimming;

/// Replaces tabs with spaces in a String.
///
/// The legacy xflate format (XFL) uses tabs as a delimiter and these cannot be
/// present in the FASTQ header.
pub(crate) fn replace_tabs_with_spaces(header: &mut String) {
    // Safety: tab to space replacements are both ASCII and therefore do not
    // disrupt UTF-8 boundaries
    unsafe {
        header.as_mut_vec().replace_all_bytes(b'\t', b' ');
    }
}

/// An extension trait for an error enabling it to be printed alongside a
/// message.
pub(crate) trait PrintWarning {
    /// Prints the error and its backtrace alongside a message to either
    /// `stdout` or `stderr`.
    ///
    /// The message includes a timestamp, indentation based on the shell level
    /// (`SHLVL` environmental variable), and the word `WARNING`.
    fn warn(&self, program: &str, message: &str, use_stderr: bool);
}

impl<E> PrintWarning for E
where
    E: Error + 'static,
{
    fn warn(&self, program: &str, message: &str, use_stderr: bool) {
        // Pad the timestamp based on the shell level using environment variable
        // SHLVL
        let shlvl = std::env::var("SHLVL").ok().and_then(|v| v.parse::<usize>().ok()).unwrap_or(0);
        let pad = "  ".repeat(shlvl.saturating_sub(1));

        if use_stderr {
            eprintln!(
                "[{now}] {pad}{program} WARNING :: {message}",
                now = Zoned::now().strftime("%Y-%m-%d %k:%M:%S")
            );
            eprint!("{}", self.display_stack())
        } else {
            println!(
                "[{now}] {pad}{program} WARNING :: {message}",
                now = Zoned::now().strftime("%Y-%m-%d %k:%M:%S")
            );
            print!("{}", self.display_stack())
        }
    }
}
