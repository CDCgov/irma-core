# IRMA-core "Sampler" README

## Motivation and Goals

Next generation sequencers often produce an excess of sequences, for reasons such as to create redundancy for enabling consensus generation, due to requiring a minimum amount of materials per run, or other reasons. Running the large amount of sequences through a pipeline can be computationally expensive, so it may be desirable to transform the data into a smaller subset before further manipulation to increase efficiency.

IRMA-core's sampler provides efficient, random, and fully-representative downsampling (also referred to as subsampling), as well as some other useful functionality.

## Downsampling Targets

IRMA-core's sampler requires a target for downsampling. This can either be provided as `--subsample-target`, which is the exact number of reads to be in the subsampled output, or `--percent-target`, which is the percentage of the original amount of reads to be in the subsampled output.

Percent targets must be provided as an integer [0-100] and may not provide an exact percentage downsampled in the cases of streamed or compressed input.

If a `--percent-target` of 100 is provided, no downsampling will occur. This could be useful for de-interleaving without downsampling.

If a `--subsample-target` is provided that is *greater* than the amount of sequences in the input, the process will succeed and give an output that is identical to the input, but provide a warning for the user.

## Inputs and Outputs

Sampler can downsample `FASTQ` and `FASTA` formats. Inputs are provided as positional arguments, with sampler accepting either a single file, or as a pair of paired-read files. The files may also be a stream (e.g., from a process substitution) or a `.gz` compressed file.

For outputs, provide one output with `-o` or two outputs with `-1` and `-2`. If no output is provided, IRMA-core will output the subsampled data to `stdout`. *Note* for paired-read input, if no output is provided when there are two inputs, the `stdout` output will be interleaved.

### Streamed and Piped inputs

Inputs can be provided via pipes or with command substitution.
The following will stream data from the `input.fasta` file, randomly downsample it to 10000 sequences, and print the output to `stdout`.

```bash
irma-core sampler \
    <(cat input.fasta) \
    --subsample-target 10000
```

### Zipped  Input/Output

The following will take a zipped `.fastq.gz` input, perform downsampling to approximately 10% of the file size, and output a standard `.fastq` file.

```bash
irma-core sampler \
    intput.fastq.gz \
    --output_file1 sampled.fastq
    --percent-target 10
```

## Paired Reads

Some sequencers (including Illumina sequencers) generate reads from both ends of the DNA fragments, resulting in two FASTQ files of paired reads. To handle these, you can provide an additional paired input file as a positional argument after the first, and/or provide an additional output file with `--output_file2` or `-2`. The paired reads will be downsampled simultaneously, ensuring that pairs get either kept or removed together.

The following will take paired input, downsample them to 10000 paired sequences, and output two zipped files (each file will contain 10000 sequences).

```bash
irma-core sampler \
    input.fastq \
    input2.fastq \
    --output-file1 sampled.fastq.gz \
    --output-file2 sampled2.fastq.gz\
    --subsample-target 10000
```

### Interleaving and De-interleaving

Paired reads can sometimes be in interleaved formats, where alternating sequences in the file will belong to the "left" and "right" paired ends.

For example:

```fastq
@Read1 1:N:0
AGCAGATAAAGACAAATCAGCTGGTTTTCCATNNN
+
FFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFF###
@Read1 2:N:0
ATGGAAAACCAGCTGATTTGTCTTTATGTGCTNNN
+
FFFFFFFFFFFFFFFF:FFFFFFFFFF,FFFF###
@Read2 1:N:0
TCTGTTTCTTCACCTGATGCTGTTGTTTACCGAGGT
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
@Read2 2:N:0
ACCTCGGTAAACAACAGCATCAGGTGAAGAAACAGA
+
FFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFF:
```

Here two sequences are paired.

Sampler can handle two paired inputs, and output a single interleaved file:

```bash
irma-core sampler \
    input1.fastq \
    input2.fastq \
    --subsample-target 10000 \
    -o interleaved.fastq
```

Sampler can also take an interleaved input and output two *de-interleaved* files:

```bash
irma-core sampler \
    interleaved.fastq.gz \
    --percent-target 10 \
    -1 output1.fastq \
    -2 output2.fastq
```

Currently, if a single input and single output are provided, each read will be treated individually; interleaved input with interleaved output is *not* supported.