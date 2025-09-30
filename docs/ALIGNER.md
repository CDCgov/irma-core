# IRMA-core Aligner README

## Motivation and Goals

IRMA-core's `aligner` provides an efficient and exact local sequence alignment routine using Striped Smith-Waterman. This is particularly useful for aligning Next Generation Sequencing (NGS) reads and performing full alignment between short genomes.

## Multithreading

`aligner` uses `rayon` to perform multithreading to enable higher throughput. To specify the number of threads, set the `RAYON_NUM_THREADS` environmental variable as described [here](https://docs.rs/rayon/latest/rayon/fn.max_num_threads.html). Or, to limit to a single thread, pass `--single-thread` to `aligner`.

For benchmarking or scenarios where a single thread is always used, the `dev_no_rayon` feature can be enabled in IRMA-core to remove the use of channels. This will also remove the `--single-thread` flag since it is no longer relevant. This feature may be removed in future releases, and so should not be relied upon except for testing.

## Inputs and Outputs

The first positional argument is a FASTA file containing the reference sequence(s), and the second argument is a FASTA or FASTQ file containing the queries. Either file may be gzip-compressed, in which case it is assumed to end in `.gz`. The output is in the [SAM alignment format](https://samtools.github.io/hts-specs/SAMv1.pdf). The score is reported with the `AS` tag for mapped reads, and the `MAPQ` field is not used (it is set to 255). The output file is specified with `--out` or `--output` flags. If not specified, output is directed to `STDOUT`. If the provided file ends in `.gz`, the output will be zipped.

As an example, consider the following inputs:

- `reference.fasta`

  ```text
  >reference
  GACTCAGTAAGACACGGTCTAGCTGACTGT
  ```

- `query.fastq`

  ```text
  @query1
  AGTAACACGGTC
  +
  IIIIIIIIIIII
  @query2
  GTCTAGGTGACTA
  +
  IIIIIIIIIIIII
  ```

Aligner could be run with a command such as:

```bash
irma-core aligner \
    reference.fasta query.fastq \
    --out alignments.sam \
    -m 1 -x 1 -o 2 -e 1
```

This produces:

```text
@query1 0       reference       6       255     5M2D7M  *       0       0       AGTAACACGGTC    IIIIIIIIIIII    AS:i:9
@query2 0       reference       17      255     12M1S   *       0       0       GTCTAGGTGACTA   IIIIIIIIIIIII   AS:i:10
```

## Scoring

IRMA-core `aligner` supports both DNA and amino acid alignments. For DNA, alignment is case-insensitive over the alphabet `ACGTN`, with all other symbols being treated as `N`. For amino acid alignment, the case-insensitive alphabet is `ACDEFGHIKLMNPQRSTVWY*BJZX`, with all other symbols being treated as `X`. The alphabet is specified with `--alphabet dna` (default) or `--alphabet aa`.

Scoring is done with either:

- A simple weight matrix, where a fixed score is used for matching residues (`--matching` or `-m`) and mismatching residues (`--mismatch` or `-x`). For DNA alignment, `--ignore-n` can be passed to make all comparisons involving an `N` have a score of 0.
- A protein substitution matrix can be specified using `--matrix <MAT>`, where `<MAT>` may be:
  - BLOSUM: `blosum30`, `blosum35`, `blosum40`, `blosum45`, `blosum50`, `blosum55`, `blosum60`, `blosum62`, `blosum65`, `blosum70`, `blosum75`, `blosum80`, `blosum85`, `blosum90`, `blosum95`, `blosum100`
  - PAM: `pam30`, `pam40`, `pam70`, `pam120`, `pam200`, `pam250`

All other combinations are invalid and will result in an error. Some more examples:

| Example                   | Scoring Scheme                                                                      |
| ------------------------- | ----------------------------------------------------------------------------------- |
| Nothing                   | A simple DNA weight matrix with default weights (`--matching 2 --mismatch -5`)      |
| `-m 1 -x 1 -o 2 -e 1`     | A simple DNA weight matrix with some weights and penalties specified                |
| `--alphabet aa`           | The default protein substitution matrix `BLOSUM_62`                                 |
| `--matrix pam250`         | A named protein weight matrix from *Zoe*, in this case `PAM250`                     |
| `--alphabet aa -m 5 -x 2` | A simple protein weight matrix with user-specified match score and mismatch penalty |

If only one of `--matching` or `--mismatch` are specified, then the other weights are set to the default.

An affine gap penalty is used with `--gap-open` and `--gap-extend`, which default to `10` and `1` respectively. Note that the gap open penalty must be at least the gap extend penalty. We use the affine gap formula, $W(k) = u(k-1) + v$, where $k$ is the gap length, $u$ is the gap extend penalty, $v$ is the gap open penalty, and $W(k)$ is the total penalty for the gap. In order to use $W(k) = uk + v$, simply pass gap open as the gap open plus gap extension.

| Parameter             | Default    | Kind                                                 | Description                                                      |
| --------------------- | ---------- | ---------------------------------------------------- | ---------------------------------------------------------------- |
| `--matching` (`-m`)   | 2          | integer $0\leq x\leq 127$                            | The score for matching residues in a simple weight matrix        |
| `--mismatch` (`-x`)   | 5          | integer $0\leq x\leq 127$                            | The penalty for mismatching residues in a simple weight matrix   |
| `--gap-open` (`-o`)   | 10         | integer $0\leq x\leq 127$                            | The penalty for opening a gap                                    |
| `--gap-extend` (`-e`) | 1          | integer $0\leq x\leq 127$                            | The penalty for extending a gap                                  |
| `--ignore-n`          | False      |                                                      | Use a score of 0 when `N` is being compared for the DNA alphabet |
| `--matrix`            | `blosum62` | [`blosum30`, `blosum35`, ..., `pam30`, `pam40`, ...] | The protein substitution matrix to use for scoring               |
| `--alphabet`          | `dna`      | [`dna`, `aa`]                                        | The alphabet to interpret the inputs as                          |

## Other Options

For DNA alignments, passing `--rev-comp` or `-r` will also check the alignment against the reverse complement and return whichever is better. SAM uses the 5th bit (16 or 0b0001 0000) to indicate that the best alignment was against the reverse complement of the reference. To exclude unmapped (zero-scoring) alignments from the output, use `--exclude-unmapped`.

By default, `aligner` will align all references against all queries and output each result. To instead only output the best match for each query, use `--best-match`.

`aligner` by default builds the striped profile from the query sequences, using them once and then discarding them. For certain data, it may be more efficienct to build larger profiles from the reference(s) and then reuse them. To do this, use `--profile-from-ref`.

| Parameter            | Description                                                                                       |
| -------------------- | ------------------------------------------------------------------------------------------------- |
| `--rev-comp` (`-r`)  | Also checks alignments against the reverse complement, outputting whichever has the highest score |
| `--exclude-unmapped` | Excludes unmapped alignments from the output file                                                 |
| `--best-match`       | The best matching alignment for each query is output, instead of all of them                      |
| `--profile-from-ref` | Which sequence to build the striped profiles from                                                 |
| `--single-thread`    | Sets the number of `rayon` threads to 1. See [here](#features) for more details                   |
