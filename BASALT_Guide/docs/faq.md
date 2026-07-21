# Troubleshooting and FAQ

Start from evidence. Preserve the complete command, stdout, stderr, `Basalt_log.txt`, `Basalt_checkpoint.txt`, software versions, database configuration, and input checksums before changing files or rerunning a stage.

## Fast diagnostic sequence

```bash
pwd
conda info --envs
command -v BASALT checkm2 metabat2 SemiBin2 bowtie2 samtools spades.py
test -n "$BASALT_WEIGHT" && find "$BASALT_WEIGHT" -maxdepth 1 -name '*_ensemble.csv'
tail -n 30 Basalt_checkpoint.txt
tail -n 100 Basalt_log.txt
grep -Ei 'error|failed|warning|traceback' *.log 2>/dev/null || true
```

## Installation and environment

### `BASALT: command not found`

Activate the environment used during installation and inspect the launcher:

```bash
conda activate basalt_env
command -v BASALT
ls -l "$CONDA_PREFIX/bin/BASALT" "$CONDA_PREFIX/bin/BASALT.py"
```

If the launcher is absent, return to the cloned repository and reinstall:

```bash
conda activate basalt_env
cd /path/to/BASALT
bash install.sh
```

### BASALT models are missing

The pipeline expects `BASALT_WEIGHT` to name a directory containing five ensemble descriptor files and the associated checkpoints.

```bash
export BASALT_WEIGHT=/absolute/path/to/BASALT_WEIGHT
find "$BASALT_WEIGHT" -maxdepth 1 -name '*_ensemble.csv' | wc -l
```

If the count is not `5`, rerun the model downloader and inspect extraction errors:

```bash
python /path/to/BASALT/BASALT_models_download.py --path "$BASALT_WEIGHT"
```

### `CheckM2 database not found`

```bash
checkm2 database --download
export CHECKM2DB=/absolute/path/to/CheckM2_database/uniref100.KO.1.dmnd
```

Test the installed CheckM2 version with a small independent input before rerunning a long BASALT analysis.

### `libcrypto.so` or another shared library is missing

Do not create an ABI-changing symbolic link as a first response. That can make a program start while remaining binary-incompatible.

Inspect the executable and solve the environment consistently:

```bash
command -v samtools
ldd "$(command -v samtools)" | grep 'not found' || true
conda list samtools openssl
```

Reinstall compatible packages in a clean environment or use the project container.

## Inputs and paths

### `IndexError: list index out of range` near input parsing

Check the paired-end grammar:

```text
one sample:       R1.fastq,R2.fastq
two samples:      S1_R1.fastq,S1_R2.fastq/S2_R1.fastq,S2_R2.fastq
```

Also confirm that filenames contain no commas or `/` characters beyond the intended delimiters.

### Do absolute paths work?

Absolute paths are not handled consistently across all legacy shell-command paths. Use a dedicated working directory and simple local symlinks. BASALT-Air is the separate implementation intended for first-class absolute-path and separate work/output-directory support.

### Can compressed reads be used?

The pipeline contains handlers for `.gz`, `.tar.gz`, and `.zip`. Plain FASTA and FASTQ inputs remain easier to audit. Test the exact archive layout on a small run because automated extraction can introduce unexpected filenames or overwrite existing files.

### Can I combine single assemblies and co-assemblies?

Yes, provided the supplied read coverage is biologically compatible with each assembly. BASALT compares selected candidates across assemblies to reduce redundancy. Document which samples contributed to every assembly and coverage matrix.

## Runtime and checkpoints

### A scheduler reports success, but the output is incomplete

Some external tools are invoked through shell calls, so a failed subcommand may not always produce a non-zero final BASALT status. Audit:

1. the last checkpoint;
2. errors in captured stderr;
3. optional-binner warnings;
4. the final FASTA count and sizes;
5. the final quality report.

### How do I resume?

From the same working directory and unchanged environment:

```bash
BASALT --mode continue
```

Do not resume after changing inputs, code, model files, databases, or major dependency versions. Start a new run directory instead.

### The run is too slow or uses too much storage

- use `--sensitive quick --refinepara quick`;
- remove optional binners;
- run a pilot on one assembly;
- request adequate local scratch storage;
- monitor memory and I/O rather than assuming CPU is the bottleneck.

Changing the preset changes the candidate space. Report the change as a method setting, not only as a performance optimization.

## Outputs and quality

### Where are final MAGs?

For the default CheckM2 route, they are in the directory named by `-o`. Without `-o`, use `Final_binset/`.

### Why is `quality_report.tsv` missing?

Possible causes include:

- CheckM2 database or executable failure;
- no bins reaching a stage;
- an external stage failing before quality assessment;
- output being written under a stage-specific directory;
- cleanup or manual movement.

Search logs before lowering scientific quality thresholds. Lowering `--min-cpn` changes which bins enter refinement and is not a general fix for a missing database or failed program.

### Does a folder named `refined` or `MAGs` guarantee acceptable genomes?

No. Directory names encode workflow state. Apply explicit completeness, contamination, taxonomic, size, and study-specific criteria to the final report.

### Can the same contig appear in several candidate bins?

Yes, especially before within- and cross-assembly selection. BASALT reduces candidate redundancy, but the final biological interpretation still requires a stated dereplication and taxonomic framework.

## Data types and acceleration

### Which read designs are supported?

BASALT accepts assembly FASTA files with paired-end short reads, ONT or PacBio CLR reads, PacBio HiFi reads, or compatible combinations. Reassembly and polishing stages depend on which read types are present and can be skipped when their required inputs are absent.

### Does BASALT use a GPU?

Some dependencies, including SemiBin 2, PyTorch, and optional VAMB configurations, may use compatible accelerators. Availability depends on the installed packages and drivers. Do not assume GPU use from the BASALT command alone. Record device utilization and software versions if acceleration affects runtime or reproducibility.

## Optional binners

### BASALT finished after an optional binner failed

This is expected failure handling in the current CheckM2 orchestration. Optional-adapter exceptions are logged and the pipeline continues. Inspect the log and verify which candidates actually entered selection.

See [Extra binners](extra-binners.md).

## BASALT versus BASALT-Air

BASALT-Air is a separate Pixi-based implementation with different path and output controls. The Conda BASALT commands on this site do not support `--workdir`, `--outdir`, `--version`, or `--check-deps`.

Do not assume two implementations produce byte-identical outputs without a version-matched regression test. Record which repository, release, dependency lock, models, and databases were used.

## Reporting a bug

Open a [GitHub issue](https://github.com/PKU-EMBL/BASALT/issues) with:

- BASALT commit or release;
- operating system and resource allocation;
- Conda environment export;
- complete command with sensitive paths redacted consistently;
- `Basalt_log.txt`, `Basalt_checkpoint.txt`, stdout, and stderr;
- model and database configuration;
- relevant external-tool versions;
- a minimal shareable reproducer when possible.

Do not upload human, clinical, or otherwise restricted sequencing data to a public issue.
