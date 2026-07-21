# Quick start

This example runs the full Conda-based BASALT workflow on one assembly and one paired-end short-read sample. It is designed to expose path, environment, database, and checkpoint problems before a larger analysis.

## Before you run

Confirm that:

- [BASALT is installed](installation.md) in an active environment;
- `BASALT_WEIGHT` points to the extracted model directory;
- CheckM2 can locate its database;
- the assembly and both read files belong to the intended sample or study design;
- the working directory is new or contains only this run.

## 1. Create an isolated working directory

```bash
mkdir -p /project/basalt_runs/study_01
cd /project/basalt_runs/study_01
```

BASALT writes checkpoints and intermediates to the current working directory. One directory per run prevents checkpoint collisions and makes provenance easier to audit.

## 2. Link the input files

```bash
ln -s /project/data/assembly.fasta .
ln -s /project/data/sample_R1.fastq .
ln -s /project/data/sample_R2.fastq .
```

Check the links and record input checksums:

```bash
test -s assembly.fasta
test -s sample_R1.fastq
test -s sample_R2.fastq

sha256sum \
  assembly.fasta \
  sample_R1.fastq \
  sample_R2.fastq \
  > input.sha256
```

## 3. Run BASALT

```bash
conda activate basalt

BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -t 32 \
  -m 128 \
  --sensitive sensitive \
  --refinepara quick \
  --min-cpn 35 \
  --max-ctn 20 \
  -q checkm2 \
  --mode new \
  -o study_01_basalt \
  > basalt.stdout.log 2> basalt.stderr.log
```

`-t` controls requested threads. `-m` reports available RAM in GB to BASALT and influences some internal parallelism; it is not a hard operating-system memory limit.

## 4. Monitor without changing the run

```bash
tail -f Basalt_log.txt
```

In another shell, inspect the last checkpoint:

```bash
tail -n 5 Basalt_checkpoint.txt
```

Do not launch a second BASALT process in the same directory.

## 5. Inspect completion and outputs

A completed full run should contain the final output directory named by `-o`:

```bash
test -d study_01_basalt
find study_01_basalt -maxdepth 1 -type f -name '*.fa' | wc -l
```

Retain these provenance files with the final MAGs:

```text
BASALT_command.txt
Basalt_checkpoint.txt
Basalt_log.txt
basalt.stdout.log
basalt.stderr.log
input.sha256
study_01_basalt/
```

The presence of FASTA files establishes that the workflow produced bins. It does not establish that every bin meets a study-specific quality or taxonomic criterion. Inspect the final CheckM2 report and document all downstream filters.

## Resume an interrupted run

Resume from the same directory, environment, code version, models, and databases:

```bash
BASALT --mode continue \
  > basalt.resume.stdout.log 2> basalt.resume.stderr.log
```

The current CLI defaults omitted arguments to empty values and relies on checkpointed intermediate files during continuation. Preserve the original command in `BASALT_command.txt`. If inputs, parameters, code, models, or databases changed, start a new run directory instead of resuming.

## Common input variants

::::{tab-set}
:::{tab-item} Short and long reads
```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -l nanopore.fastq \
  -t 64 -m 256 --mode new -o study_01_basalt
```
:::

:::{tab-item} PacBio HiFi
```bash
BASALT \
  -a assembly.fasta \
  -hf hifi.fastq \
  -t 32 -m 128 --mode new -o study_01_basalt
```
:::

:::{tab-item} Two short-read samples
```bash
BASALT \
  -a assembly.fasta \
  -s sample_1_R1.fastq,sample_1_R2.fastq/sample_2_R1.fastq,sample_2_R2.fastq \
  -t 64 -m 256 --mode new -o study_01_basalt
```
:::
::::

Read [Input formats and paths](inputs.md) before using several samples or assemblies. The delimiter describes file grouping; it does not encode biological pairing between assemblies and samples.

## Next steps

- [Tutorial](tutorial.md) for the public demo dataset and an HPC submission template.
- [Command-line reference](usage.md) for every option and its actual default.
- [Outputs and quality control](output.md) for stage-specific artifacts.
- [Reproducibility and reporting](reproducibility.md) before interpreting or publishing a MAG set.
