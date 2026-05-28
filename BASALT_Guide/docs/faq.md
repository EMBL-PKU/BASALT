# FAQ

## General

### What is BASALT?

BASALT (Binning Across a Series of Assemblies Toolkit) is a tool for binning and post-binning refinement of metagenomic assemblies. It generates high-quality metagenome-assembled genomes (MAGs) from short-read, long-read, and hybrid datasets.

### What makes BASALT different from metaWRAP or DASTool?

Three main differences:

1. **Multi-assembly support**: BASALT accepts multiple single assemblies and co-assemblies in a single run with automatic dereplication.
2. **Deep learning refinement**: Neural network-based contamination detection and removal at the individual contig level.
3. **Long-read integration**: Efficient utilisation of long reads for contig retrieval and polishing.

### What data types does BASALT support?

| Mode | Short Reads | Long Reads (ONT/PacBio) | PacBio HiFi |
|---|---|---|---|
| SRS only | :material-check: | — | — |
| SRS + LRS | :material-check: | :material-check: | — |
| SRS + HiFi | :material-check: | — | :material-check: |
| HiFi only | — | — | :material-check: |

Long-read-only (ONT/CLR) without short reads is not yet supported.

---

## Performance

### How long does BASALT take to run?

This depends on dataset size, complexity, and chosen parameters. A typical run with `--sensitive sensitive` on a 32-core workstation with 256 GB RAM may take 12–24 hours for a moderate-complexity metagenome (e.g., human gut). The demo dataset completes in ~6 hours on a 32-core machine.

### How can I speed up BASALT?

- Use `--sensitive quick --refinepara quick` for the fastest preset
- Run only the module you need: `--module autobinning`, `--module refinement`, or `--module reassembly`
- Increase thread count (`-t`)
- Ensure sufficient RAM to avoid swap overhead
- For multiple assemblies, BASALT is already more efficient than running multiple single-assembly jobs

### Does BASALT support GPU?

Yes. BASALT v1.2.0 supports GPU acceleration for Semibin2 and deep learning model inference. Ensure CUDA-compatible PyTorch is installed.

---

## Workflow

### Can I run only part of the pipeline?

Yes. Use `--module`:

- `--module autobinning`: Run only binning and bin selection
- `--module refinement`: Run only contamination removal and contig retrieval (on existing bins)
- `--module reassembly`: Run only reassembly (on existing refined bins)
- `--module all`: Run the full pipeline (default)

### Can I use my own bins with BASALT refinement?

Yes. There are two approaches:

1. **Standalone refinement** (`-r`): Point BASALT to a single binset folder.
2. **Data Feeding** (`-d`): Import multiple external binsets with automatic coverage recalculation.

### Can the same contig appear in multiple bins?

Under SA + CA mode, redundant bins can be generated. BASALT's Bin Selection module identifies and removes redundancies automatically through dereplication. The final best binset is non-redundant.

### What if my run is interrupted?

BASALT records progress in `Basalt_checkpoint.txt`. Resume with:

```bash
BASALT --mode continue
```

For a completely fresh restart:

```bash
BASALT --mode new
```

---

## Input & Output

### Does BASALT support absolute file paths?

No. Place all input files in the current working directory, or create soft links using `ln -s`.

### Can I specify an output directory?

BASALT generates all output in the current working directory. The `-o` flag sets the output folder name prefix, not the path.

### Where are the final MAGs?

In `<output_prefix>_final_binset/` (default: `Final_binset_final_binset/`). Each bin is a single multi-FASTA file.

---

## Troubleshooting

### CheckM2 database not found

```bash
checkm2 database --download
```

Refer to [CheckM2 documentation](https://github.com/chklovski/CheckM2) for details.

### BASALT: command not found

Ensure the conda environment is activated and scripts have correct permissions:

```bash
conda activate basalt_env
chmod -R 755 $(dirname $(which BASALT))/*
```

### IndexError: list index out of range

This typically means BASALT cannot parse file paths. Ensure all input files are in the working directory (no absolute paths) and the `-s` argument uses the correct delimiter format (`/` between samples, `,` between read pairs).

### quality_report.tsv not found

This occurs when too few bins pass quality thresholds, so CheckM2 produces no output. Try:

- Lowering `--min-cpn` (e.g., `--min-cpn 20`)
- Using `--sensitive more-sensitive` for more thorough binning
- Ensuring sufficient sequencing coverage

---

## Installation

### Can I install BASALT without Conda?

It is strongly recommended to use Conda. BASALT depends on many bioinformatics tools that are easier to manage through Bioconda. Alternatively, use the Singularity image.

### How do I install on macOS?

BASALT is designed for Linux x64 systems. macOS is not officially supported. Use the Singularity image in a Linux VM or Docker container on macOS.
