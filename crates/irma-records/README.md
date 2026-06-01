# irma-records

`irma-records` is a small library for working with sequencing records in the
formats that IRMA handles most often. It is not tied to influenza analysis or to
the full `irma-core` workflow: it provides reusable readers, writers, record
transforms, paired-read handling, and SAM pair merging for tools that need to
process FASTA, FASTQ, or SAM data.

FASTA stores sequence records. FASTQ stores sequence records together with per-base
quality scores. SAM stores sequence alignments. These formats are common in NGS
(next generation sequencing), where many short reads are produced from a sample
and later cleaned, paired, aligned, or summarized.

## What It Provides

The [`io`] module contains convenience types for opening plain or compressed
inputs, selecting FASTA or FASTQ parsing from the data itself, writing records,
and attaching useful path context to errors. These helpers are useful for command
line tools that accept either files or standard streams and want consistent
error messages.

The [`fastq`] module adds operations for editing reads before analysis, including
hard clipping, masking, adapter or barcode removal, primer handling, poly-G
cleanup, canonical base recoding, and read-quality summaries. These operations
are exposed through [`fastq::ReadTransforms`] for FASTQ records.

The [`paired`] module helps keep paired-end reads together. It can zip two read
streams into checked pairs, split interleaved paired reads back into two streams,
and inspect common read-name conventions to determine whether a record belongs
to the first or second side of a pair.

The [`sam`] module provides pair merging for aligned SAM records and summary
statistics describing agreement and disagreement between mates.

Finally, the [`hashing`] module provides the seeded hash builder used by IRMA
tools when reproducible hashing is requested.

## How IRMA Uses It

In `irma-core`, this crate is the record-handling layer beneath the command line
workflows. Read quality control uses [`fastq`] transforms to remove or mask parts
of reads and to evaluate read quality before keeping, filtering, or compacting
records. Commands that sample, interleave, deinterleave, align, or write sequence
data use [`io`] and [`paired`] so that FASTA and FASTQ files can be handled with a
shared set of checks and output paths.

IRMA also uses the [`sam`] functionality after alignment to combine paired-end
evidence against a reference sequence and collect high-level merge statistics.
Those details may evolve, but the role of this crate remains the same: keep
format-specific record work in one place so higher-level workflows can focus on
analysis steps.

## Where To Start

For format-agnostic FASTA or FASTQ input, start with [`io::FastXReader`]. For
file-opening helpers, see [`io::InputOptions`] and [`io::OutputOptions`]. For
FASTQ read cleanup, see [`fastq::ReadTransforms`]. For paired-end streams, see
[`paired::ZipPairedReadsExt`] and [`paired::DeinterleavedPairedReadsExt`]. For
aligned read-pair merging, see [`sam::SamMergeablePairs`].

## License and Notices

This crate inherits all licenses, disclaimers and notices described in the [project
README](https://github.com/CDCgov/irma-core#notices).
