# Tutorial

This tutorial uses the public BASALT demo dataset as a regression test. It is not a benchmark for new biological conclusions, and exact outputs can change with BASALT, database, model, and dependency versions.

## Demo dataset

The archived dataset is available from Figshare:

> Qiu, Z. BASALT demo files. Figshare (2023). [https://doi.org/10.6084/m9.figshare.22323424](https://doi.org/10.6084/m9.figshare.22323424)

The record contains:

| File | Purpose |
|---|---|
| `Data.tar.gz` | Short reads, long reads, and an OPERA-MS assembly |
| `Final_bestbinset.tar.gz` | Historical expected final bins |
| `basalt.sh` | Historical demonstration command |

The complete download is approximately 967 MB. Preserve the Figshare version and file checksums in the run record.

## 1. Prepare the environment

Follow [Installation](installation.md), then verify:

```bash
eval "$(micromamba shell hook --shell bash)"
micromamba activate basalt
# Conda alternative: conda activate basalt
BASALT --help
test -n "$BASALT_WEIGHT"
test -d "$BASALT_WEIGHT"
command -v checkm2 metabat2 SemiBin2 bowtie2 samtools spades.py
```

For activation-free execution, prefix environment commands with `micromamba run -n basalt`. The remaining commands in this tutorial assume an active `basalt` environment.

## 2. Download and inspect the demo

Download all files from the Figshare record into a new directory. Do not overwrite an existing BASALT run.

```bash
mkdir -p /project/basalt_demo
cd /project/basalt_demo
```

After downloading, record checksums before extraction:

```bash
sha256sum Data.tar.gz Final_bestbinset.tar.gz basalt.sh \
  > figshare-input.sha256
```

Inspect archive paths and the historical command before executing anything:

```bash
tar -tzf Data.tar.gz | sed -n '1,40p'
sed -n '1,200p' basalt.sh
```

The archived `basalt.sh` may reflect the 2023 CLI and output names. Use the current command grammar documented below unless reproducing the historical software environment.

## 3. Extract the data

```bash
mkdir -p data expected
tar -xzf Data.tar.gz -C data
tar -xzf Final_bestbinset.tar.gz -C expected
```

Inventory the extracted filenames:

```bash
find data -maxdepth 2 -type f -print | sort
find expected -maxdepth 2 -type f -print | sort
```

Identify the assembly, paired-end files, and long-read file from that inventory. The variables below are placeholders; replace them with the extracted basenames.

## 4. Create a run directory

```bash
mkdir -p run
cd run

ln -s ../data/<assembly.fasta> assembly.fasta
ln -s ../data/<short_R1.fastq> short_R1.fastq
ln -s ../data/<short_R2.fastq> short_R2.fastq
ln -s ../data/<long_reads.fastq> long_reads.fastq

sha256sum \
  assembly.fasta short_R1.fastq short_R2.fastq long_reads.fastq \
  > input.sha256
```

Use simple link names so the command is independent of archive-specific paths.

## 5. Run BASALT

```bash
BASALT \
  -a assembly.fasta \
  -s short_R1.fastq,short_R2.fastq \
  -l long_reads.fastq \
  -t 32 \
  -m 128 \
  --sensitive sensitive \
  --refinepara quick \
  --min-cpn 35 \
  --max-ctn 20 \
  -q checkm2 \
  --mode new \
  -o demo_basalt \
  > basalt.stdout.log 2> basalt.stderr.log
```

The Figshare record reports a runtime within 6 hours on a 32-core Intel Xeon Gold 5218 workstation for the historical demo environment. Treat that value as context, not a service-level expectation. Current dependency versions, storage, database placement, and accelerator availability can change runtime.

## 6. Audit completion

```bash
tail -n 30 Basalt_checkpoint.txt
tail -n 80 Basalt_log.txt
grep -Ei 'error|failed|warning|traceback' \
  basalt.stderr.log Basalt_log.txt || true

test -d demo_basalt
find demo_basalt -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l
```

Do not rely on process exit status alone. Some external commands are invoked through shell calls, and optional-binner failures can be handled as warnings.

## 7. Compare with the historical expected set

The archived expected set is useful for regression, but filenames and selected bins can differ across releases and quality databases. Begin with structural comparisons:

```bash
find demo_basalt -maxdepth 1 -type f -name '*.fa' | wc -l
find ../expected -type f \( -name '*.fa' -o -name '*.fasta' -o -name '*.fna' \) | wc -l

find demo_basalt -maxdepth 1 -type f -name '*.fa' -exec sha256sum {} + \
  | sort > observed-bins.sha256
```

Interpret a difference only after checking:

- BASALT commit or release;
- model files;
- CheckM2 version and database;
- external-binner and assembler versions;
- sensitivity and refinement settings;
- warnings and skipped stages.

A checksum mismatch does not by itself establish a regression. FASTA ordering, identifier normalization, and dependency changes can alter bytes without changing the underlying sequence set.

## HPC example

```bash
#!/bin/bash
#SBATCH --job-name=basalt_demo
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -u

export MAMBA_ROOT_PREFIX=/absolute/path/to/micromamba-root
export BASALT_WEIGHT=/absolute/path/to/BASALT_WEIGHT
export CHECKM2DB=/absolute/path/to/uniref100.KO.1.dmnd
MICROMAMBA_BIN=/absolute/path/to/micromamba
cd /project/basalt_demo/run

"$MICROMAMBA_BIN" run -n basalt BASALT \
  -a assembly.fasta \
  -s short_R1.fastq,short_R2.fastq \
  -l long_reads.fastq \
  -t "$SLURM_CPUS_PER_TASK" \
  -m 128 \
  --sensitive sensitive \
  --refinepara quick \
  -q checkm2 \
  --mode new \
  -o demo_basalt
```

Replace every placeholder with a validated site path. Conda users can source the site `conda.sh`, activate `basalt`, and run `BASALT` directly instead. `set -u` detects unset shell variables but does not change BASALT's handling of external-command failures. Perform the completion audit after the scheduler job ends.

## Archive the tutorial result

```bash
micromamba env export -n basalt > basalt-environment.yml
# Conda alternative: conda env export -n basalt --no-builds > basalt-environment.yml
git -C /path/to/BASALT rev-parse HEAD > basalt-git-commit.txt

(
  cd "$BASALT_WEIGHT"
  find . -type f \
    \( -name '*.pth' -o -name '*_ensemble.csv' \) -print0 \
    | sort -z \
    | xargs -0 sha256sum
) > basalt-models.sha256

tar -czf demo-provenance.tar.gz \
  BASALT_command.txt \
  Basalt_checkpoint.txt \
  Basalt_log.txt \
  basalt.stdout.log \
  basalt.stderr.log \
  input.sha256 \
  basalt-environment.yml \
  basalt-git-commit.txt \
  basalt-models.sha256
```

Proceed to [Reproducibility and reporting](reproducibility.md) before adapting the command to study data.
