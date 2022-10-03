#![allow(unused_variables, irrefutable_let_patterns, dead_code, unused_imports)]

use clap::{Args, ValueHint};
use either::*;
use std::fs::File;
use std::fs::OpenOptions;
use std::io::prelude::*;
use std::io::{stdin, BufReader, BufWriter, Stdin};
use std::path::PathBuf;

/*

         . "\t\t-C|--complement-add\t\t\t\n"
         . "\t\t-O|--ordinal-headers\t\t\t\n"
         . "\t\t-F|--file-id <STR>\t\t\t\n"
         . "\t\t-S|--save-quality <STR>\t\t\t\n"
         . "\t\t-A|--save-stats <STR>\t\t\t\n"
         . "\t\t-K|--skip-remaining\t\t\t\n"
         . "\t\t-H|--keep-header\t\t\t\n"
         . "\t\t-c|--clip-adapter <STR>\t\t\tClip adapter.\n"
         . "\t\t-m|--mask-adapter <STR>\t\t\tMask adapter.\n"
         . "\t\t-Z|--fuzzy-adapter\t\t\tAllow one mismatch.\n"
         . "\t\t-U|--uracil-to-thymine\t\t\t\n"
         . "\t\t-E|--enforce-clipped-length\t\t\n"

*/

#[derive(Args)]
pub struct FastqConverterArgs {
    fastq_input_file: Option<PathBuf>,

    #[arg(short = 'Q', long)]
    /// Outputs fastQ instead of fastA format.
    fastq_output: bool,
    #[arg(short = 'C', long)]
    /// Take the reverse complement and add to data.
    complement_add: bool,
    #[arg(short = 'H', long)]
    /// Keep the fastq header as usual.
    keep_header: bool,

    #[arg(short = 'T', long, default_value_t = 0)]
    /// Specify the read quality threshold (geometric mean, median).
    read_quality: u8,
    #[arg(short = 'M', long)]
    /// Interprets the threshold (-T) as the median, not the geometric mean.
    use_median: bool,

    #[arg(short = 'L', long, default_value_t = 0)]
    /// Minimum length of sequence read data, filtered otherwise.
    min_length: usize,
    #[arg(short = 'E', long)]
    /// The minimum length threshold (-L) is enforced when adapter clipped (-c).
    enforce_clipped_length: bool,

    #[arg(short = 'm', long)]
    /// Specify adapter sequence and mask when found in reads.
    mask_adapter: Option<String>,
    #[arg(short = 'c', long)]
    /// Specify adapter sequence and clip appropriate ends when found in reads.
    clip_adapter: Option<String>,
    #[arg(short = 'Z', long)]
    /// Allow up to one mismatch for adapter clipping (-c) or masking (-m).
    fuzzy_adapter: bool,

    #[arg(short = 'U', long)]
    /// Covert uracil (U) to thymine (T) in read sequences.
    uracil_to_thymine: bool,

    #[arg(short = 'S', long, value_hint = ValueHint::FilePath)]
    /// Save quality file for back-mapping.
    save_quality_file: Option<PathBuf>,

    #[arg(short = 'G', long, value_hint = ValueHint::FilePath)]
    /// Quality control log path and filename.
    log_file: Option<PathBuf>,
    #[arg(short = 'g', long)]
    /// Specify log ID tag (integer) for output collation.
    log_id: Option<usize>,

    #[arg(short = 'O', long)]
    /// Replace header with monotonically increasing ordinal headers.
    ordinal_headers: bool,
    #[arg(short = 'F', long)]
    /// File id tag for ordinal header mode (-O).
    ordinal_file_id: Option<String>,

    #[arg(short = 'A', long, value_hint = ValueHint::FilePath)]
    /// Save quality vs. length statistics file for analysis.
    save_stats: Option<PathBuf>,
    #[arg(short = 'K', long)]
    /// Do not output data FASTA/FASTQ data (assumes -A).
    skip_remaining: bool,
}

pub fn fastq_process(args: FastqConverterArgs) {
    // Improve with Zoe
    let mut fastq_file_reader = if let Some(file_path) = args.fastq_input_file {
        BufReader::new(Either::Left(
            OpenOptions::new()
                .read(true)
                .create(true)
                .open(&file_path)
                .unwrap_or_else(|_| panic!("Cannot open log file! See: {}", file_path.display())),
        ))
    } else {
        BufReader::new(Either::Right(stdin()))
    };

    let log_file_writer = if let Some(file_path) = args.log_file {
        Some(BufWriter::new(
            OpenOptions::new()
                .append(true)
                .create(true)
                .open(&file_path)
                .unwrap_or_else(|_| panic!("Cannot open log file! See: {}", file_path.display())),
        ))
    } else {
        None
    };

    let save_stats_writer = if let Some(file_path) = args.save_stats {
        Some(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create(true)
                .open(&file_path)
                .unwrap_or_else(|_| {
                    panic!("Cannot open statistics file! See: {}", file_path.display())
                }),
        ))
    } else {
        None
    };

    let quality_file_writer = if let Some(file_path) = args.save_quality_file {
        Some(BufWriter::new(
            OpenOptions::new()
                .write(true)
                .create(true)
                .open(&file_path)
                .unwrap_or_else(|_| {
                    panic!(
                        "Cannot open quality file for writing! See: {}",
                        file_path.display()
                    )
                }),
        ))
    } else {
        None
    };

    let file_id = if let Some(id) = args.ordinal_file_id {
        if args.complement_add {
            "|".to_owned() + &id + "|"
        } else {
            id + "|"
        }
    } else {
        "".to_owned()
    };

    let (forward_adapter, reverse_adapter, adapter_mask) = if let Some(ref a) = args.mask_adapter {
        let forward = a.as_bytes().to_ascii_uppercase();
        let length = forward.len();
        let reverse = forward.iter().copied().rev().collect(); // add reverse complement using Zoe

        (forward, reverse, [b'N'].repeat(length))
    } else if let Some(ref a) = args.clip_adapter {
        let forward = a.as_bytes().to_ascii_uppercase();
        let reverse = forward.iter().copied().rev().collect(); // add reverse complement using Zoe

        (forward, reverse, Vec::new())
    } else {
        (Vec::new(), Vec::new(), Vec::new())
    };

    //if args.fuzzy_adapter && (args.mask_adapter.is_some() || args.clip_adapter.is_some()) {}

    //let (pp,nn,dnp)

    let mut reads_passing_qc: u32 = 0;
    let mut ordinal_id: u32 = 0;
    let mut fastq_buffer = Vec::new();

    // improve with Zoe
    let msg = "FASTQ reader failure!";
    loop {
        let bytes = fastq_file_reader
            .read_until(b'\n', &mut fastq_buffer)
            .expect(msg);
        if bytes == 0 {
            break;
        }

        let header = if fastq_buffer.len() > 1 {
            if fastq_buffer.ends_with(b"\n") {
                fastq_buffer.pop();
            }

            String::from_utf8(fastq_buffer.clone()).expect(msg)
        } else {
            String::new()
        };
        fastq_buffer.clear();

        let bytes = fastq_file_reader
            .read_until(b'\n', &mut fastq_buffer)
            .expect(msg);
        if bytes == 0 {
            break;
        }

        let sequence = if fastq_buffer.len() > 1 {
            if fastq_buffer.ends_with(b"\n") {
                fastq_buffer.pop();
            }

            fastq_buffer
                .iter()
                .copied()
                .map(|c| c.to_ascii_uppercase())
                .collect()
        } else {
            Vec::new()
        };
        fastq_buffer.clear();

        // Junk line
        let bytes = fastq_file_reader
            .read_until(b'\n', &mut fastq_buffer)
            .expect(msg);
        if bytes == 0 {
            break;
        }
        fastq_buffer.clear();

        let bytes = fastq_file_reader
            .read_until(b'\n', &mut fastq_buffer)
            .expect(msg);
        if bytes == 0 {
            break;
        }

        let quality = if fastq_buffer.len() > 1 {
            if fastq_buffer.ends_with(b"\n") {
                fastq_buffer.pop();
            }

            fastq_buffer
                .iter()
                .copied()
                .map(|c| c.to_ascii_uppercase())
                .collect()
        } else {
            Vec::new()
        };
        fastq_buffer.clear();

        print!(
            "(1){}\n(2){}\n(3){}\n",
            header,
            String::from_utf8_lossy(&sequence),
            String::from_utf8_lossy(&quality)
        );
    }
}
/*

    while ( my $hdr = <> ) {
        chomp($hdr);

        my $seq = <>;
        chomp($seq);

        my $junk = <>;
        chomp($junk);

        my $quality = <>;
        chomp($quality);



my @fuzzyAdaptersFWD = ();
my @fuzzyAdaptersREV = ();
if ( $fuzzyAdapter && ( $maskAdapter || $clipAdapter ) ) {
    my $L = length $forwardAdapter;
    if ( $L > 0 ) {
        my $N = '[ATCGN]';
        $L--;
        for my $i ( 0 .. $L ) {
            my $tmp = $forwardAdapter;
            substr( $tmp, $i, 1, $N );
            push( @fuzzyAdaptersFWD, $tmp );
            $tmp = $reverseAdapter;
            substr( $tmp, $i, 1, $N );
            push( @fuzzyAdaptersREV, $tmp );
        }
    }
} else {
    $fuzzyAdapter = 0;
}

if ( $minLength < 0 ) {
    die("ERROR: minimum length must be a non-negative integer.\n");
} else {
    $minLength = int($minLength);
}

local $RS = "\n";

my $reads_passing_qc = 0;
my $id               = 0;
my $pp               = q{};
my $nn               = q{};
my $dnp              = q{};

if ($complementAndAdd) {
    $pp  = 'P';    # position strand
    $nn  = 'N';    # negative strand
    $dnp = '_';    # delimiter
}

while ( my $hdr = <> ) {
    chomp($hdr);

    my $seq = <>;
    chomp($seq);

    my $junk = <>;
    chomp($junk);

    my $quality = <>;
    chomp($quality);

    if ( length $seq < $minLength ) {
        next;
    }

    if ($uracilToThymine) {
        $seq =~ tr/uU/tT/;
    }

    if ($clipAdapter) {
        if ( $seq =~ /$reverseAdapter/ismx ) {
            $seq     = substr( $seq,     0, $LAST_MATCH_START[0] );
            $quality = substr( $quality, 0, $LAST_MATCH_START[0] );
        } elsif ( $seq =~ /$forwardAdapter/ismx ) {
            $seq     = substr( $seq,     $LAST_MATCH_END[0] );
            $quality = substr( $quality, $LAST_MATCH_END[0] );
        } elsif ($fuzzyAdapter) {
            my $findFuzzyFWD = 1;
            foreach my $tmp (@fuzzyAdaptersREV) {
                if ( $seq =~ /$tmp/ismx ) {
                    $seq          = substr( $seq,     0, $LAST_MATCH_START[0] );
                    $quality      = substr( $quality, 0, $LAST_MATCH_START[0] );
                    $findFuzzyFWD = 0;
                    last;
                }
            }

            if ($findFuzzyFWD) {
                foreach my $tmp (@fuzzyAdaptersFWD) {
                    if ( $seq =~ /$tmp/ismx ) {
                        $seq     = substr( $seq,     $LAST_MATCH_END[0] );
                        $quality = substr( $quality, $LAST_MATCH_END[0] );
                        last;
                    }
                }
            }
        }

        if ( $clippedMinLength && length $seq < $minLength ) {
            next;
        }
    } elsif ($maskAdapter) {
        if ( $seq =~ /$reverseAdapter/ismx ) {
            $seq =~ s/$reverseAdapter/$adapterMask/ismx;
        } elsif ( $seq =~ /$forwardAdapter/ismx ) {
            $seq =~ s/$forwardAdapter/$adapterMask/ismx;
        } elsif ($fuzzyAdapter) {
            my $findFuzzyFWD = 1;
            foreach my $tmp (@fuzzyAdaptersREV) {
                if ( $seq =~ /$tmp/ismx ) {
                    $seq =~ s/$tmp/$adapterMask/ismx;
                    $findFuzzyFWD = 0;
                    last;
                }
            }

            if ($findFuzzyFWD) {
                foreach my $tmp (@fuzzyAdaptersFWD) {
                    if ( $seq =~ /$tmp/ismx ) {
                        $seq =~ s/$tmp/$adapterMask/ismx;
                        last;
                    }
                }
            }
        }
    }

    my @a = unpack( "c* i*", $quality );
    my $q = 0;
    my $n = scalar(@a);
    $id++;

    if ($useMedian) {
        my @sorted = sort(@a);
        if ( $n % 2 == 0 ) {
            $q = ( $sorted[$n / 2] + $sorted[( $n / 2 ) - 1] ) / 2 - 33;
        } else {
            $q = $sorted[( $n - 1 ) / 2] - 33;
        }
    } else {
        ## use geometric mean otherwise
        foreach my $x (@a) {
            $q += $x;
        }
        $q = ( $q - $n * 33 ) / $n;
    }

    if ($saveStats) {
        print $STAT $id, "\t", $q, "\t", $n, "\n";
        if ($skipRemaining) {
            next;
        }
    }

    if ( $q >= $qualityThreshold ) {
        $reads_passing_qc++;
        $hdr = substr( $hdr, 1 );

        if ( !$keepHeader ) { $hdr =~ tr/ /_/; }

        if ($ordinal) {
            if ($fastQformat) {
                print STDOUT $pp, $fileID, $id, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
            } else {
                print STDOUT '>', $pp, $fileID, $id, "\n", $seq, "\n";
                if ($saveQualityFile) {
                    print $QUA $pp, $fileID, $id, "\t", $hdr, "\t", $quality, "\n";
                }
            }
        } else {
            if ($fastQformat) {
                print STDOUT '@', $hdr, $dnp, $pp, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
            } else {
                print STDOUT '>', $hdr, $dnp, $pp, '|', sprintf( "%.1f", $q ), '|', $n, "\n", $seq, "\n";
                if ($saveQualityFile) {
                    print $QUA $hdr, $dnp, $pp, "\t", $quality, "\n";
                }
            }
        }

        # Take the reverse complement and add it
        # courtesy http://reverse-complement.com/ambiguity.html
        if ($complementAndAdd) {
            $seq = reverse($seq);
            $seq =~ tr/gcatrykmbvdhuGCATRYKMBVDHU/cgtayrmkvbhdaCGTAYRMKVBHDA/;
            $quality = reverse($quality);
            if ($ordinal) {
                if ($fastQformat) {
                    print STDOUT $nn, $fileID, $id, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
                } else {
                    print STDOUT '>', $nn, $fileID, $id, "\n", $seq, "\n";
                    if ($saveQualityFile) {
                        print $QUA $nn, $fileID, $id, "\t", $hdr, "\t", $quality, "\n";
                    }
                }
            } else {
                if ($fastQformat) {
                    print STDOUT $hdr, $dnp, $nn, "\n", $seq, "\n", $junk, "\n", $quality, "\n";
                } else {
                    print STDOUT '>', $hdr, $dnp, $nn, '|', sprintf( "%.1f", $q ), '|', $n, "\n", $seq, "\n";
                    if ($saveQualityFile) {
                        print $QUA $hdr, $dnp, $nn, "\t", $quality, "\n";
                    }
                }
            }
        }
    }
}

if ( $logFile ne q{} ) {
        if ( $logID ne q{} ) {
        $logID = ':' . $logID;
    }
    my $reads_passing_length = $id;
    print $LOG_OUT $logFile, "$logID\t", $reads_passing_length, "\t", $reads_passing_qc, "\t", $qualityThreshold, "\t",
      $minLength,
      "\t", $useMedian, "\n";
    close $LOG_OUT or croak("Cannot close file: $OS_ERROR\n");
}

*/
