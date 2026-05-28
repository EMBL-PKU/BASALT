# Tutorial

A step-by-step guide for running BASALT with demo data.

---

## Prerequisites

- BASALT installed and environment activated (see [Installation](installation.md))
- Demo data downloaded from the repository
- At least 32 cores and 128 GB RAM recommended

---

## Demo Files

The `BASALT demo files/` directory contains:

| File | Description |
|---|---|
| `Data.tar.gz` | Short-read, long-read FASTQ files and an OPERA-MS assembled contig file |
| `Final_bestbinset.tar.gz` | Expected final output bins for validation |
| `basalt.sh` | Shell script to run the demo |

Extract the data:

```bash
tar -zxf BASALT\ demo\ files/Data.tar.gz
```

A 32-core workstation with 128 GB RAM is expected to complete the demo in ~6 hours.

---

## Example 1: Short-Read Only, Single Assembly

### Step 1: Prepare Working Directory

```bash
mkdir basalt_test && cd basalt_test
ln -s /path/to/assembly.fasta .
ln -s /path/to/sample_R1.fq .
ln -s /path/to/sample_R2.fq .
```

### Step 2: Activate Environment

```bash
conda activate basalt_env
```

### Step 3: Run BASALT

```bash
BASALT -a assembly.fasta \
    -s sample_R1.fq,sample_R2.fq \
    -t 60 -m 250
```

### Step 4: Check Results

```bash
ls Final_binset_final_binset/
# List of MAG FASTA files
```

---

## Example 2: Multi-Assembly with Short + Long Reads

### Step 1: Run

```bash
BASALT -a as1.fa,as2.fa,as3.fa \
    -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq/s3_R1.fq,s3_R2.fq \
    -l lr1.fastq,lr2.fastq,lr3.fastq \
    -t 64 -m 256
```

### Step 2: What happens

1. **Autobinning**: MetaBAT2, CONCOCT, Semibin2 (sensitive preset) run on all 3 assemblies.
2. **Bin Selection**: Within each assembly, the best bins are selected. Then cross-assembly dereplication removes redundancies.
3. **Refinement**: DL model removes contamination. Paired-end tracking and long-read connections retrieve missed contigs. Long reads also trigger polishing (Racon).
4. **Reassembly**: SPAdes reassembles each bin to improve contiguity.

---

## Example 3: Quick Mode

For large or complex datasets where runtime is a concern:

```bash
BASALT -a as1.fa,as2.fa \
    -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq \
    -t 64 -m 256 \
    --sensitive quick \
    --refinepara quick
```

This uses fewer binners and parameter combinations, and skips deep contig retrieval.

---

## Example 4: Only Run Refinement on Existing Bins

If you already have bins from another tool and want to use BASALT's deep-learning refinement:

```bash
BASALT -a assembly.fa \
    -s sample_R1.fq,sample_R2.fq \
    -r My_Bins_Folder \
    -c Coverage_matrix_for_binning_assembly.fa.txt \
    -t 32 -m 128
```

---

## Example 5: Data Feeding (Multiple External Binsets)

Import bins from VAMB, manual curation, or other tools for refinement:

```bash
BASALT -s sample_R1.fq,sample_R2.fq \
    -d vamb_output,my_curated_bins \
    --binset-index 1 \
    -t 32 -m 128
```

This will:
1. Reindex bins from both folders
2. Recompute coverage matrices
3. Generate paired-end connections
4. Run CheckM2 quality assessment
5. Proceed to Refinement

---

## Running on HPC Clusters

### SLURM Example

```bash
#!/bin/bash
#SBATCH --job-name=basalt
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G
#SBATCH --time=48:00:00

source ~/.bashrc
conda activate basalt_env

BASALT -a as1.fa,as2.fa \
    -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq \
    -t 64 -m 250 \
    --mode new
```

### Singularity with SLURM

```bash
#!/bin/bash
#SBATCH --job-name=basalt_sif
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G

singularity run -B $PWD basalt.sif BASALT \
    -a as1.fa -s s1_R1.fq,s1_R2.fq \
    -t 32 -m 128
```

---

## Verifying Results

After a successful run, check the quality of your MAGs:

```bash
# Check completeness and contamination
head Final_binset_final_binset/quality_report.tsv

# Count MAGs meeting quality thresholds
awk -F'\t' 'NR>1 && $3>=50 && $4<10 {count++} END {print count " MAGs ≥50% complete, <10% contamination"}' \
    Final_binset_final_binset/quality_report.tsv
```

---

## Troubleshooting Common Issues

### Run hangs or crashes

1. Check `Basalt_log.txt` for the last completed step.
2. Ensure sufficient RAM: BASALT is memory-intensive, especially during binning and reassembly.
3. Resume with `BASALT --mode continue`.

### Low MAG count

1. Try `--sensitive more-sensitive` for more thorough binning.
2. Lower `--min-cpn` threshold (e.g., `--min-cpn 20`).
3. Ensure sufficient sequencing depth — low-coverage samples produce fewer recoverable MAGs.
4. Consider using `--refinepara deep` for more aggressive contig retrieval.

### Error: CheckM2 database not found

```bash
checkm2 database --download
```

Refer to the [FAQ](faq.md) for more troubleshooting guidance.
