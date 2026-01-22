use std::{
    fmt::{Display, Write},
    io::{BufRead, ErrorKind, Read},
};
use zoe::{
    data::{
        CheckSequence,
        err::WithErrorContext,
        fasta::FastaSeq,
        records::{HeaderReadable, SequenceReadable},
    },
    define_whichever,
    prelude::{FastQ, FastQReader, FastaReader, QualityScores},
};

define_whichever! {
    /// A reader over either FASTA or FASTQ data, determined automatically based
    /// on the first non-whitespace character.
    ///
    /// This can be used as an [`Iterator`] directly, in which case it yields
    /// [`FastX`] records and performs a match statement on every [`next`] call.
    /// Using [`for_each`] or other iterator methods can eliminate extra match
    /// statements.
    ///
    /// Or, this can be matched, with each type of reader handled manually.
    ///
    /// [`next`]: Iterator::next
    /// [`for_each`]: Iterator::for_each
    pub enum FastXReader<R: Read> {
        /// A FASTQ reader (records include quality scores).
        Fastq(FastQReader<R>),
        /// A FASTA reader (records do not include quality scores).
        Fasta(FastaReader<R>),
    }

    #[map(FastX::map_result)]
    impl<U: Read> Iterator for FastXReader<U> {
        type Item = std::io::Result<FastX>;

        // The rest of the methods are added by define_whichever
    }
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
            Some(other) => Err(std::io::Error::new(
                ErrorKind::InvalidData,
                format!("File started with invalid character: {}", *other as char),
            )
            .with_context("Unable to determine whether the file is FASTA or FASTQ!"))?,
            None => Err(std::io::Error::new(
                ErrorKind::InvalidData,
                "No FASTQ or FASTA data was found!",
            )),
        }
    }
}

/// A record type for either FASTQ or FASTA data.
pub struct FastX {
    pub header:   String,
    pub sequence: Vec<u8>,
    pub quality:  Option<QualityScores>,
}

impl FastX {
    /// A helper function for mapping a result of [`FastQ`] or [`FastaSeq`] to a
    /// result of [`FastX`].
    #[inline]
    fn map_result<T>(result: std::io::Result<T>) -> std::io::Result<Self>
    where
        FastX: From<T>, {
        result.map(|record| record.into())
    }
}

impl From<FastQ> for FastX {
    #[inline]
    fn from(value: FastQ) -> Self {
        Self {
            header:   value.header,
            sequence: value.sequence.into_vec(),
            quality:  Some(value.quality),
        }
    }
}

impl From<FastaSeq> for FastX {
    #[inline]
    fn from(value: FastaSeq) -> Self {
        Self {
            header:   value.name,
            sequence: value.sequence,
            quality:  None,
        }
    }
}

impl Display for FastX {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if let Some(quality) = &self.quality {
            writeln!(f, "@{}", self.header)?;

            if self.sequence.is_ascii_simd::<16>() {
                // SAFETY: we just checked it is ASCII using our fast SIMD function.
                // ASCII is valid UTF8.
                f.write_str(unsafe { std::str::from_utf8_unchecked(&self.sequence) })?;
            } else {
                f.write_str(&String::from_utf8_lossy(&self.sequence))?;
            }

            writeln!(f, "\n+\n{quality}")
        } else {
            writeln!(f, ">{}", self.header)?;

            if self.sequence.is_ascii_simd::<16>() {
                // SAFETY: we just checked it is ASCII using our fast SIMD function.
                // ASCII is valid UTF8.
                f.write_str(unsafe { std::str::from_utf8_unchecked(&self.sequence) })?;
            } else {
                f.write_str(&String::from_utf8_lossy(&self.sequence))?;
            }

            f.write_char('\n')
        }
    }
}

impl HeaderReadable for FastX {
    #[inline]
    fn header(&self) -> &str {
        &self.header
    }
}

impl SequenceReadable for FastX {
    #[inline]
    fn sequence_bytes(&self) -> &[u8] {
        &self.sequence
    }
}
