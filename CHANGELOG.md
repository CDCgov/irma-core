# IRMA-core Changelog

All notable changes to this project will be documented in this file. The format
is roughly based on [Keep a Changelog], and this project tries to adheres to
[Semantic Versioning]. For IRMA vs IRMA-core compatibility, please see the
[version matrix](VERSION_MATRIX.md).

## [0.9.0] - TBD

### Added

- Adds `verbose` flag to `sampler` for printing total number of reads in
  original input and number of downsampled records to `stderr`
- Adds more detailed error messages with more context, including file paths for many IO errors

### Fixes

- Fixes a bug where leading whitespace in a record file can cause the program to stall (in `sampler` and `aligner`)
- Allows repeated input/output files from processes if they are device files (e.g., `/dev/null`)

## [0.8.1] - 2025-12-02

### Fixes

- `aligner` with `--rev-comp` and `--best-match` now produces the correct positions and CIGAR strings

## [0.8.0] - 2025-11-21

### Added

- Adds standalone `aligner` process for performing pairwise local sequence
  alignment with the Striped Smith Waterman algorithm

### Changed

- The `trimmer` process now can perform deinterleaving, as well as widow filtering with interleaved inputs and outputs

## [0.7.0] - 2025-09-17

### Added

- Adds standalone `xleave` process for interleaving and de-interleaving
  paired-end FastQ and FASTA files

## [0.6.2] - 2025-09-25

### Fixes

- `xflate` recently produced invalid FASTQ during inflation. This has been
  fixed. This bug was not surfaced in IRMA.

## [0.6.1] - 2025-09-09

### Changed

- The `fastq-converter` process has been removed.
- More context is added for IO errors in `sampler`

### Fixes

- Fixes bug where SRA paired-end headers could be could be incorrectly flagged
  as unmated when `--filter-widows` was selected in `trimmer`, or with paired
  reads in `sampler`
- Fixes a reproducibility issue in `sampler` where the seed wasn't being used in
  one case

## [0.6.0] - 2025-08-29

### Added

- Adds a new standalone process called `sampler` that provides support for
randomly downsampling FastQ or FASTA files. Deinterleaving is also supported.

### Changed

- Use threads and anonymous pipes (in Rust) when decompressing zipped inputs
- Improve error messages to include file names in `trimmer` and `preprocess`
- Other dependencies updated, such as a faster foldhash.

### Fixes

- `mergeSAMpairs`: updated Zoe dependency (in v0.0.19, using v0.0.20) fixes hard clipping in merged reads. Reader
  error messages are also improved generally.

## [0.5.1] - 2025-06-03

### Fixes

- Fixes a bug where `--polyg-trim` would trim the sequence but not the quality scores.

## [0.5.0] - 2025-05-30

### Added

- `preprocess` (formerly `qc-trim-deflate`) now supports zipped inputs (`fastq.gz` files).

### Changed

- Re-factor and specify dependencies.
- The process `qc-trim-deflate` is now called `preprocess`, and it includes most
  of the same functionality as `trimmer` (note that its arguments have been
  conformed to match `trimmer` rather than the legacy `fastq_converter`
  process).
- The deprecated `fastq_converter` will be removed in the next version.
- For `trimmer`/`preprocess`, `--filter-widows` can no longer be specified when
  there is only a single input file.
- For `fastq_converter`/`preprocess`, `--keep-header` has been removed since
  IRMA has made this an unmodifiable default since v1.1.2.

## [0.4.3] - 2025-05-01

### Fixes

- Fixed bug in `trimmer` to allow handling of multi-member `fastq.gz` files

## [0.4.2] - 2025-04-30

### Changed

- The `Dockerfile` has been updated for nextflow compatibility (it currently
  needed: awk, grep, ps, sed)
- The `num-procs` process now allows `LOCAL_PROCS_OVERRIDE` to be greater than
  the available cores. Use `IFX_LOCAL_PROCS` if you would like it capped at the
  available cores.

## [0.4.1] - 2025-04-28

### Added

- The `trimmer` process now supports compressed `fastq.gz` files for both input and output.
- Adds `flate2` crate as a dependency for handling `.gz` files

### Changed

- The release process has been improved to produce better and more artifacts.

## [0.4.0] - 2025-04-16

### Added

- `trimmer` now supports `-f, --filter-widows` for filtering widowed/orphaned reads in a paired-end context, with support for two input and two output `fastq` files.
- Github actions added to generate downloadable releases

### Changed

- Updates dependencies, particularly with `foldhash` and `zoe` improvements
- For `trimmer`, long argument `--fastq-output-file` is now `--fastq-output`

### Fixes

- For `num-procs`, the environemntal var `LOCAL_PROCS_OVERRIDE` now caps at the
  max available cores, consistent with the warning message.

## [0.3.1] - 2025-03-21

- **Added**: The `num-procs` process now has the option `--cap-cores-using-env`
  so environmental variables like `NSLOTS` and `IFX_LOCAL_PROCS` can further cap
  resources. The var `LOCAL_PROCS_OVERRIDE` trumps other settings up to the
  available cores.

## [0.3.0] - 2025-03-18

### Added

- Adds `trimmer` standalone process. Read the [docs](docs/TRIMMER.md).
- Adds integrated `num-procs` process for portable core counts.
- Adds documents for open-sourcing.

### Changed

- Refactors into integrated / standalone processes.
- Updates dependencies.

### Fixes

- Incorporates fix for Zoe that affected IUPAC base recoding.

## [0.2.0] - 2024-10-30

- **Added**: Updates the new `qc-trim-deflate` (name subject to change) process with various fixes and improvements.
- **Changed**: Slims down and corrects raw read counts in the the fastq-converter process.

## [0.1.6] - 2024-10-11

- **Added**: xflate and qc-trim-deflate processes (initial sketch) and rewritten the flow for the fastq-converter process.

## [0.1.5] - 2024-09-04

- **Added**: custom inexact matching algorithm from [Zoe]

<!-- Versions -->
[0.9.0]: https://github.com/CDCgov/irma-core/compare/v0.8.1...v0.9.0
[0.8.1]: https://github.com/CDCgov/irma-core/compare/v0.8.0...v0.8.1
[0.8.0]: https://github.com/CDCgov/irma-core/compare/v0.7.0...v0.8.0
[0.7.0]: https://github.com/CDCgov/irma-core/compare/v0.6.1...v0.7.0
[0.6.2]: https://github.com/CDCgov/irma-core/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/CDCgov/irma-core/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/CDCgov/irma-core/compare/v0.5.1...v0.6.0
[0.5.1]: https://github.com/CDCgov/irma-core/compare/v0.5.0...v0.5.1
[0.5.0]: https://github.com/CDCgov/irma-core/compare/v0.4.3...v0.5.0
[0.4.3]: https://github.com/CDCgov/irma-core/compare/v0.4.2...v0.4.3
[0.4.2]: https://github.com/CDCgov/irma-core/compare/v0.4.1...v0.4.2
[0.4.1]: https://github.com/CDCgov/irma-core/compare/v0.4.0...v0.4.1
[0.4.0]: https://github.com/CDCgov/irma-core/compare/v0.3.1...v0.4.0
[0.3.1]: https://github.com/CDCgov/irma-core/compare/v0.3.0...v0.3.1
[0.3.0]: https://github.com/CDCgov/irma-core/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/CDCgov/irma-core/compare/v0.1.6...v0.2.0
[0.1.6]: https://github.com/CDCgov/irma-core/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/CDCgov/irma-core/compare/IRMA@v1.1.5...v0.1.5

<!-- Links -->
[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html
[Zoe]: https://github.com/CDCgov/zoe
