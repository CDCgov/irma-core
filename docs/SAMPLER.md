# IRMA-core "Sampler" README

## Motivation and Goals

Next generation sequencers often produce an excess of reads, for reasons such as to create redundancy for enabling consensus generation, due to requiring a minimum amount of reagents per run, or other reasons. Running the full read set through a pipeline can be computationally expensive, so it may be desirable to transform the data into a smaller subset to save time.

IRMA-core's sampler process provides efficient, random, and fully-representative downsampling (also referred to as subsampling), as well as some other useful functionality.

## Downsampling Targets

IRMA-core's sampler requires a target for downsampling. This can either be provided as `--subsample-target`, which is the exact number of reads to be in the subsampled output, or `--percent-target`, which is the percentage of the original amount of reads to be in the subsampled output.

Percent targets must be provided as an integer [0-100] and may not provide an exact percentage downsampled in the cases of streamed or compressed input.

- If a `--percent-target` of 100 is provided, no downsampling will occur. This could be useful for de-interleaving without downsampling.
- If a `--subsample-target` is provided that is *greater* than the amount of sequences in the input, the process will succeed and give an output that is identical to the input, but provide a warning for the user.

## Inputs and Outputs

Sampler can downsample `FASTQ` and `FASTA` formats. Inputs are provided as positional arguments, with sampler accepting either a single file, or as a pair of paired-read files. The files may also be a stream (e.g., from a process substitution) or a `.gz` compressed file.

For outputs, you can select one output file with `-o` or two output files with `-1` and `-2`. If no output is provided, IRMA-core will output the subsampled data to `stdout`.

*Note:* if only one output is selected (`stdout` or a file) for paired-end inputs, that output will be interleaved.

### Streamed Inputs

Inputs can be provided via a file (or command substitution for multiple files). The following will stream data from the `input.fasta` file, randomly downsample it to 10,000 sequences, and print the output to `stdout`.

*Note:* in cases of multiple files being provided, `sampler` will treat these as one continuous stream of the files concatenated together.

```bash
# simple file
irma-core sampler input.fasta --subsample-target 10000

# lots of files together
irma-core sampler <(cat *.fasta *.fa *.fas) --subsample-target 10000
```

### Zipped Input/Output

The following will take a zipped `.fastq.gz` input, perform downsampling to approximately 10% of the file size, and output a standard `.fastq` file.

```bash
irma-core sampler intput.fastq.gz \
    --output_file1 sampled.fastq \
    --percent-target 10
```

## Paired Reads

Some sequencers (including Illumina sequencers) generate reads from both ends of the DNA fragments, resulting in two paired-end read files. To handle these, you can provide an additional paired input file as a positional argument after the first, and/or provide an additional output file with `--output_file2` or `-2`. Paired-end reads will be **downsampled together**, ensuring that pairs get either kept or removed together.

The following will take paired input, downsample them to 10000 paired sequences, and output two zipped files (each file will contain 10000 sequences).

```bash
irma-core sampler \
    input_R1.fastq input_R2.fastq \
    --output-file1 sampled_R1.fastq.gz \
    --output-file2 sampled_R2.fastq.gz \
    --subsample-target 10000
```

### Interleaving and De-interleaving

Paired-end reads can sometimes be in interleaved formats (for example, downloading from the SRA web UI), where alternating sequences in the file will belong to the "left" and "right" paired ends.

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

**Two** paired-end inputs / **one** interleaved output:

```bash
irma-core sampler \
    input_R1.fastq input_R2.fastq \
    --subsample-target 10000 \
    -o interleaved.fastq
```

**One** interleaved input / **two** deinterleaved outputs:

```bash
irma-core sampler \
    interleaved.fastq.gz \
    --percent-target 10 \
    -1 output_R1.fastq -2 output_R2.fastq
```

If a single input and output are provided, no interleaving is assumed.

## Verbose Output

An optional flag of `--verbose` or `-v` can be used to print diagnostics to `stderr`. The output is of the form:
`Downsampled 177564 total records to 35512 (20.00 %).`