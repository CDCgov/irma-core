use std::{
    error::Error,
    fmt::{Debug, Display},
};
use zoe::data::err::GetCode;

/// An error arising from paired headers which were incorrect (either poorly
/// formatted or mismatching).
#[derive(Copy, Clone, Eq, PartialEq, Hash, Debug)]
pub enum PairedHeaderError {
    /// The first header failed to be parsed.
    ParsingErrorFirst,
    /// The second header failed to be parsed.
    ParsingErrorSecond,
    /// There was a mismatch between the headers.
    Mismatch,
}

impl Display for PairedHeaderError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PairedHeaderError::ParsingErrorFirst => {
                write!(f, "Failed to parse the first header")
            }
            PairedHeaderError::ParsingErrorSecond => {
                write!(f, "Failed to parse the second header")
            }
            PairedHeaderError::Mismatch => write!(f, "Mismatching IDs found!"),
        }
    }
}

impl Error for PairedHeaderError {}
impl GetCode for PairedHeaderError {}
