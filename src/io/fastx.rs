use std::{
    fs::File,
    io::{BufRead, ErrorKind, Read},
    path::Path,
};
use zoe::prelude::{FastQReader, FastaReader};

/// A reader for either FASTA or FASTQ.
pub enum FastXReader<R: Read> {
    Fastq(FastQReader<R>),
    Fasta(FastaReader<R>),
}

impl<R: std::io::Read> FastXReader<R> {
    /// Creates an iterator over either FASTA or FASTQ data from a type
    /// implementing [`Read`], automatically wrapping it with a buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if an IO error occurs, or if the file does not contain
    /// `@` or `>` as its first non-whitespace character.
    pub fn from_readable(inner: R) -> std::io::Result<Self> {
        let mut buffer = std::io::BufReader::new(inner);
        let start = loop {
            let start_buffer = buffer.fill_buf()?;
            if start_buffer.is_empty() {
                break None;
            }
            let start = start_buffer.iter().find(|x| !x.is_ascii_whitespace());
            if start.is_some() {
                break start;
            }
        };

        match start {
            Some(b'>') => Ok(FastXReader::Fasta(FastaReader::from_bufreader(buffer)?)),
            Some(b'@') => Ok(FastXReader::Fastq(FastQReader::from_bufreader(buffer)?)),
            _ => Err(std::io::Error::new(
                ErrorKind::InvalidData,
                "Unable to determine whether the file is FASTA or FASTQ!",
            )),
        }
    }
}

impl FastXReader<std::fs::File> {
    /// Creates an iterator over either FASTA or FASTQ data from a file, using a
    /// buffered reader.
    ///
    /// ## Errors
    ///
    /// Will return `Err` if file or permissions do not exist, or if the file
    /// does not contain `@` or `>` as its first non-whitespace character.
    #[allow(dead_code)]
    pub fn from_filename<P>(filename: P) -> std::io::Result<Self>
    where
        P: AsRef<Path>, {
        FastXReader::from_readable(File::open(filename)?)
    }
}
