use crate::Cli;
use clap::{CommandFactory, error::ErrorKind};

pub(crate) mod clipping;

/// Aborts clap with a given error `message` due to a custom parsing error.
///
/// The subcommand should be specified as a lowercase string with `subcommand`
/// if available. This ensures the help message is as informative as possible.
/// If an invalid subcommand is passed, it will be ignored.
pub(crate) fn abort_clap(kind: ErrorKind, message: impl std::fmt::Display, subcommand: Option<&str>) -> ! {
    let mut command = Cli::command();

    if let Some(subcommand) = subcommand
        && let Some(c) = command.get_subcommands_mut().find(|c| c.get_name() == subcommand)
    {
        c.error(kind, message).exit();
    } else {
        command.error(kind, message).exit()
    }
}
