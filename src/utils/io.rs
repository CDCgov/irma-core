use crate::utils::whichever::define_whichever;
use flate2::{Compression, bufread::MultiGzDecoder, write::GzEncoder};
use std::{
    fs::File,
    io::{BufReader, BufWriter, Error as IOError, Stdin, Stdout, stdout},
    path::Path,
};

define_whichever! {
    @match_compressed_reader

    #[doc="An enum for the different acceptable input types"]
    pub(crate) enum ReadFileZip {
        File(BufReader<File>),
        Zipped(MultiGzDecoder<BufReader<File>>),
    }

    impl Read for ReadFileZip {}
}

define_whichever! {
    @match_compressed_writer

    #[doc="An enum for the different acceptable output types"]
    pub(crate) enum  WriteFileZipStdout {
        File(BufWriter<File>),
        Zipped(GzEncoder<BufWriter<File>>),
        Stdout(BufWriter<Stdout>),
    }

    impl Write for WriteFileZipStdout {}
}

define_whichever! {
    @match_file_stdin

    pub(crate) enum ReadFileStdin {
        File(File),
        Stdin(Stdin),
    }

    impl Read for ReadFileStdin {}
}

pub(crate) fn open_fastq_file<P: AsRef<Path>>(path: P) -> Result<ReadFileZip, IOError> {
    let file = File::open(&path)?;
    let buf_reader = BufReader::new(file);

    let is_gz = path.as_ref().extension().is_some_and(|ext| ext == "gz");

    let reader = if is_gz {
        ReadFileZip::Zipped(MultiGzDecoder::new(buf_reader))
    } else {
        ReadFileZip::File(buf_reader)
    };

    Ok(reader)
}

pub(crate) fn create_writer<P: AsRef<Path>>(path: Option<P>) -> Result<WriteFileZipStdout, IOError> {
    let writer = match path {
        Some(ref p) => {
            let is_gz = p.as_ref().extension().is_some_and(|ext| ext == "gz");
            let file = File::create(p)?;
            let buf_writer = BufWriter::new(file);

            if is_gz {
                WriteFileZipStdout::Zipped(GzEncoder::new(buf_writer, Compression::default()))
            } else {
                WriteFileZipStdout::File(buf_writer)
            }
        }
        None => WriteFileZipStdout::Stdout(BufWriter::new(stdout())),
    };

    Ok(writer)
}
