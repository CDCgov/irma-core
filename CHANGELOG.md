# IRMA-core Changelog

All notable changes to this project will be documented in this file. The format is roughly based on [Keep a Changelog], and this project tries to adheres to [Semantic Versioning]. For IRMA vs IRMA-core compatibility, please see the [version matrix](VERSION_MATRIX.md).

## [0.3.0] - 2025-03-18

## Added

- Adds `trimmer` standalone process. Read the [docs](docs/TRIMMER.md).
- Adds integrated `num-procs` process for portable core counts.
- Adds documents for open-sourcing.

### Changed

- Refactors into integrated / standalone processes.
- Updates dependencies.

### Fixes

- Incorporates fix for Zoe that affected IUPAC base recoding.

## [0.2.0] - 2024-10-30

-- **Added**: Updates the new `qc-trim-deflate` (name subject to change) process with various fixes and improvements.
-- **Changed**: Slims down and corrects raw read counts in the the fastq-converter process.

## [0.1.6] - 2024-10-11

- **Added**: xflate and qc-trim-deflate processes (initial sketch) and rewritten the flow for the fastq-converter process.

## [0.1.5] - 2024-09-04

- **Added**: custom inexact matching algorithm from [Zoe]

<!-- Versions -->
[0.3.0]: https://github.com/CDCgov/irma-core/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/CDCgov/irma-core/compare/v0.1.6...v0.2.0
[0.1.6]: https://github.com/CDCgov/irma-core/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/CDCgov/irma-core/compare/IRMA@v1.1.5...v0.1.5

<!-- Links -->
[keep a changelog]: https://keepachangelog.com/en/1.0.0/
[semantic versioning]: https://semver.org/spec/v2.0.0.html
[Zoe]: https://github.com/CDCgov/zoe
