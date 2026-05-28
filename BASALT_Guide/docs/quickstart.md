# Quick Start

This page gets you running BASALT in a few minutes.

## Prerequisites

- Conda environment created (see [Installation](installation.md))
- Input files in the current working directory

## Minimal Command

```bash
conda activate basalt_env
BASALT -a assembly.fasta -s sample_R1.fq,sample_R2.fq -t 32 -m 128
```

This runs the full pipeline (autobinning + refinement + reassembly) using default settings.

## Common Workflows

=== "Short-read only"

    ```bash
    BASALT -a as1.fa,as2.fa \
        -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq \
        -t 60 -m 250
    ```

=== "Short + Long reads"

    ```bash
    BASALT -a as1.fa \
        -s sample_R1.fq,sample_R2.fq \
        -l nanopore.fastq \
        -t 60 -m 250
    ```

=== "Short + HiFi reads"

    ```bash
    BASALT -a as1.fa \
        -s sample_R1.fq,sample_R2.fq \
        -hf hifi_reads.fastq \
        -t 32 -m 128
    ```

=== "Fast mode"

    ```bash
    BASALT -a as1.fa \
        -s sample_R1.fq,sample_R2.fq \
        -t 32 -m 128 \
        --sensitive quick --refinepara quick
    ```

## Input File Requirements

1. **Place all input files in the working directory** — BASALT does not support absolute paths.
2. **Assembly files** must be in FASTA format (`.fa`, `.fna`, `.fasta`).
3. **Read files** must be in FASTQ format (`.fq`, `.fastq`). Compressed files (`.gz`, `.tar.gz`, `.zip`) are automatically decompressed.
4. **Each paired-end sample** requires read pairs separated by `,`. Multiple samples are separated by `/`.

## Resuming a Run

BASALT uses checkpoints. If your run is interrupted:

```bash
BASALT --mode continue
```

For a fresh start (discarding previous progress):

```bash
BASALT -a assembly.fasta -s reads_R1.fq,reads_R2.fq -t 32 -m 128 --mode new
```

## Expected Output

See the [Output Files](output.md) page for details on BASALT's output structure.

## What's Next?

- [Usage](usage.md) — full command-line reference
- [Tutorial](tutorial.md) — step-by-step guide with demo data
- [FAQ](faq.md) — common questions and answers
