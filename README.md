# IRMA-core

**IRMA-core** is suite of subcommands for both stand-alone and integrated use with [IRMA](https://cdcgov.github.io/irma) for *select portions* of the virus sequence assembly process. Processes integrated with IRMA are subject to more rapid internal revision and change while stand-alone processes will attempt to avoid breaking CLI changes where possible. More primitive bioinformatic functionality in IRMA-core can sometimes be pushed down to [The Zoe Anthology](https://github.com/CDCgov/zoe) library for general use cases (PRs welcomed).

**This binary compiles using Rust *nightly*.**

## Processes and their usage

### Standalone subcommands

| Process   | Description                                                                                             | Usage                      |
| --------- | ------------------------------------------------------------------------------------------------------- | -------------------------- |
| `trimmer` | Used for removing adapters, barcodes, and primers among other things. [Read the docs](docs/TRIMMER.md). | `irma-core trimmer --help` |

### Integrated with IRMA

| Process            | Description                                                                                                        | Usage                              |
| ------------------ | ------------------------------------------------------------------------------------------------------------------ | ---------------------------------- |
| `qc-trim-deflate`  | Performs all-in-one FastQ quality control, trimming, and deflation to XFL and FASTA formats. Similar to `trimmer`. | `irma-core qc-trim-deflate --help` |
| `fastq-converter`ø | Performs FastQ quality control, file conversions, and adapter trimming.                                            | `irma-core fastq-converter --help` |
| `merge-sam`†       | Merges Illumina paired-end reads with parsimonious error correction and detection                                  | `irma-core merge-sam --help`       |
| `xflate`†          | Deflates FastQ files to deduplicated Fasta files, or reinflates deduplicated Fasta files to FastQ files.           | `irma-core xflate --help`          |
| `num-procs`        | Provides the physical or logical cores of a CPU portably.                                                          | `irma-core num-procs --help`       |

For compatibility notes between IRMA-core and IRMA, see the [version matrix](docs/VERSION_MATRIX.md).

  *† May be combined into a future process, deprecated and removed.*\
  *ø Deprecated, will be removed.*

## FAQ

### How do I install IRMA-core?

The *correct* version of IRMA-core will ship with IRMA and [is packaged in all releases since v1.1.0](https://github.com/CDCgov/irma/releases). IRMA-core became non-optional in IRMA v1.3 and later. You can also download IRMA-core for standalone use cases from the [IRMA-core release page](https://github.com/CDCgov/irma-core/releases).

To compile and install IRMA-core yourself, first install [rustup](https://forge.rust-lang.org/infra/other-installation-methods.html), and then:

```bash
rustup toolchain install nightly
git clone https://github.com/CDCgov/irma-core
cd irma-core
cargo +nightly b -r

# Install here:
ls -l target/release/irma-core
```

For RHEL 8 compatible Linux distributions, you can also re-build IRMA-core using the [latest builder image](docs/BUILDER-README.md):

```bash
docker pull ghcr.io/your-org/irma-core/builder:latest
```

### How should this application be cited?

One can cite the IRMA [manuscript](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-3030-6) in general, but a software citation of IRMA-core specifically is also possible for stand-alone use cases. We include some citations for convenience in [`CITATION.bib`](https://github.com/CDCgov/irma-core/blob/main/BIBLIOGRAPHY.bib)).

### Why does IRMA-core have a difference license than IRMA?

IRMA [packages](https://wonder.cdc.gov/amd/flu/irma/disclaimer.html) other binaries and that affects its license. IRMA-core is an indendent project from IRMA and so can use more flexible licensing.

### Will IRMA-core replace IRMA?

No, IRMA contains data + code/binary + configuration that are appropriate for use in that repository. However, in the future it is possible IRMA-core could eliminate dependence on components and make IRMA easier to use and maintain.

### What are some goals of the project?

In no particular order:

- Retire technical debt and make the IRMA project more maintainable going forward.
- Increase IRMA's portability to systems like Windows.
- Eliminate unneeded dependencies and components.
- Relax licensing for IRMA itself (see above).
- Accelerate execution and/or reduce memory pressure within IRMA.
- Ease the integration of modern features.
- Provide stand-alone functionality outside of IRMA.

## Notices

### Contact Info

For direct correspondence on the project, feel free to contact: [Samuel S. Shepard](mailto:sshepard@cdc.gov), Centers for Disease Control and Prevention or reach out to other [contributors](CONTRIBUTORS.md).

### Public Domain Standard Notice

This repository constitutes a work of the United States Government and is not subject to domestic copyright protection under 17 USC § 105. This repository is in the public domain within the United States, and copyright and related rights in the work worldwide are waived through the [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).  All contributions to this repository will be released under the CC0 dedication.  By submitting a pull request you are agreeing to comply with this waiver of copyright interest.

### License Standard Notice

The repository utilizes code licensed under the terms of the Apache Software License and therefore is licensed under ASL v2 or later. This source code in this repository is free: you can redistribute it and/or modify it under the terms of the Apache Software License version 2, or (at your option) any later version. This source code in this repository is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the Apache Software License for more details. You should have received a copy of the Apache Software License along with this program. If not, see: <http://www.apache.org/licenses/LICENSE-2.0.html>. The source code forked from other open source projects will inherit its license.

### Privacy Standard Notice

This repository contains only non-sensitive, publicly available data and information. All material and community participation is covered by the [Disclaimer](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md). For more information about CDC's privacy policy, please visit <http://www.cdc.gov/other/privacy.html>.

### Contributing Standard Notice

Anyone is encouraged to contribute to the repository by [forking](https://help.github.com/articles/fork-a-repo) and submitting a pull request. (If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git).) By contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users under the terms of the [Apache Software License v2](http://www.apache.org/licenses/LICENSE-2.0.html) or later.

All comments, messages, pull requests, and other submissions received through CDC including this GitHub page may be subject to applicable federal law, including but not limited to the Federal Records Act, and may be archived. Learn more at [http://www.cdc.gov/other/privacy.html](http://www.cdc.gov/other/privacy.html).

### Records Management Standard Notice

This repository is not a source of government records, but is a copy to increase collaboration and collaborative potential. All government records will be published through the [CDC web site](http://www.cdc.gov).

## Additional Standard Notices

Please refer to [CDC's Template Repository](https://github.com/CDCgov/template) for more information about [contributing to this repository](https://github.com/CDCgov/template/blob/main/CONTRIBUTING.md), [public domain notices and disclaimers](https://github.com/CDCgov/template/blob/main/DISCLAIMER.md), and [code of conduct](https://github.com/CDCgov/template/blob/main/code-of-conduct.md).
