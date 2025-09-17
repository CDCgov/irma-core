# IRMA-core "XLeave" README

## Motivation and Goals

This is a small, standalone process for interleaving two paired-end FastQ/FASTA files, or de-interleaving one interleaved paired-end FastQ/FASTA into two separate files. These tasks can also be performed in conjunction with primer/adapter trimming or random downsampling using `trimmer` or `sampler`, respectively.

## Inputs and Outputs

`xleave` will infer whether to interleave or de-interleave based on the number of inputs and outputs. If one input file is provided with two outputs, it will attempt to de-interleave the input file. If two inputs are provided with one output, it will attempt to interleave the two input files. If only one input and output, or two inputs and two outputs are provided, `xleave` will throw an error, as no interleaving or de-interleaving can occur in these cases.

`xleave` can handle FastQ and FASTA formats. Inputs are provided as positional arguments, with `xleave` accepting either a single file, or a pair of paired-end files. The files may also be a stream (e.g., from a process substitution) or a `.gz` compressed file.

For outputs, you can select one output file with `-o` or output files with `-1` and `-2`. If no output is provided, IRMA-core will output the interleaved data to `stdout`. Note that `stdout` can only be used for single, interleaved output.

The following will take a zipped `.fastq.gz` input and de-interleave it into two files.

```bash
irma-core xleave \
    input.fastq.gz \
    -1 out_R1.fastq \
    -2 out_R2.fastq
```

The following will take two paired `.fastq` files and interleave them into one zipped file.

```bash
irma-core xleave \
    input_R1.fastq input_R2.fastq \
    -o interleaved_out.fastq.gz
```

## Paired Headers

In cases of both interleaving and de-interleaving, validation of headers is performed. As `xleave` reads the input(s), it checks each pair of paired headers to ensure that they match. If a mismatch is found, the process will exit early.
