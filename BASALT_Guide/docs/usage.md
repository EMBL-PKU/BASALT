# Usage

## Command-Line Interface

```
BASALT -a <assemblies> -s <short_reads> -t <threads> -m <RAM_GB> [options]
```

---

## Required Arguments

| Argument | Type | Description | Example |
|---|---|---|---|
| `-a`, `--assemblies` | `str` | Comma-separated assembly FASTA files | `-a as1.fa,as2.fa,as3.fa` |
| `-s`, `--shortreads` | `str` | Paired-end reads. Read pairs separated by `,`; datasets by `/` | `-s s1_r1.fq,s1_r2.fq/s2_r1.fq,s2_r2.fq` |
| `-t`, `--threads` | `int` | Number of CPU threads | `-t 64` |
| `-m`, `--ram` | `int` | Maximum RAM in GB (minimum: 32) | `-m 250` |

**Supported file formats:**

- Assemblies: `.fa`, `.fna`, `.fasta`
- Reads: `.fq`, `.fastq`
- Compression: `.gz`, `.tar.gz`, `.zip`

---

## Optional Arguments

| Argument | Type | Default | Choices | Description |
|---|---|---|---|---|
| `-l`, `--longreads` | `str` | *none* | — | Comma-separated long-read files (ONT/PacBio, excludes HiFi) |
| `-hf`, `--hifi` | `str` | *none* | — | Comma-separated PacBio HiFi read files |
| `-e`, `--extra_binner` | `str` | *none* | `m`, `v`, `l` | Additional binners: `m`=MetaBinner, `v`=VAMB, `l`=LorBin |
| `-o`, `--out` | `str` | `Final_binset` | — | Output folder name prefix |
| `-q`, `--quality-check` | `str` | `checkm2` | `checkm`, `checkm2` | Quality assessment software |
| `--min-cpn` | `int` | `35` | — | Minimum completeness (%) for refinement |
| `--max-ctn` | `int` | `20` | — | Maximum contamination (%) for refinement |
| `--mode` | `str` | `continue` | `new`, `continue` | Start fresh or resume from checkpoint |
| `--module` | `str` | `all` | `autobinning`, `refinement`, `reassembly`, `all` | Pipeline stage to run |
| `--sensitive` | `str` | `sensitive` | `quick`, `sensitive`, `more-sensitive` | Binning sensitivity preset |
| `--refinepara` | `str` | `quick` | `quick`, `deep` | Refinement depth |
| `-r`, `--refinement-binset` | `str` | *none* | — | Binset folder for standalone refinement |
| `-c`, `--coverage-list` | `str` | *none* | — | Coverage files for standalone refinement |
| `-b`, `--binsets-list` | `str` | *none* | — | Binset folders for dereplication |
| `-d`, `--data-feeding-folder` | `str` | *none* | — | External binset folders for Data Feeding |
| `--binset-index` | `int` | `500` | — | Start index for extra binsets in Data Feeding |

---

## Sensitivity Presets

The `--sensitive` flag controls which binners are used and their parameter ranges:

| Preset | Binners | Parameters | Runtime |
|---|---|---|---|
| `quick` | MetaBAT2, Semibin2 | MetaBAT2: 200/300/400/500; Semibin2: 100 | :material-rabbit: Fastest |
| `sensitive` | MetaBAT2, CONCOCT, Semibin2 | As above + CONCOCT (1–2 settings) | :material-tortoise: Balanced |
| `more-sensitive` | Maxbin2, MetaBAT2, CONCOCT, Semibin2 | Maxbin2: 0.3/0.5/0.7/0.9 + full settings | :material-snail: Most thorough |

---

## Refinement Depth

| Option | Description |
|---|---|
| `quick` | Standard contig retrieval. Faster, suitable for most datasets. |
| `deep` | Extended contig retrieval at the sequence retrieval step. Recovers more contigs but takes longer. |

---

## Extra Binners

Binners enabled via `-e` run **in addition to** the default set (MetaBAT2, Maxbin2, CONCOCT, Semibin2):

| Flag | Tool | Method | Reference |
|---|---|---|---|
| `m` | MetaBinner | k-mer composition + coverage clustering | — |
| `v` | VAMB | Variational autoencoder | *Nat. Biotechnol.* (2021) |
| `l` | LorBin | Long-read multi-scale adaptive clustering | *Nat. Commun.* (2025) |

Combine multiple extra binners: `-e m,l` or `-e m,v,l`.

---

## Examples

### Basic Short-Read Run

```bash
BASALT -a assembly.fasta \
    -s sample_R1.fq,sample_R2.fq \
    -t 32 -m 128
```

### Multi-Assembly with Long Reads

```bash
BASALT -a as1.fa,as2.fa,as3.fa \
    -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq \
    -l long.fastq \
    -t 60 -m 250
```

### Fast Mode

```bash
BASALT -a assembly.fasta \
    -s sample_R1.fq,sample_R2.fq \
    -t 32 -m 128 \
    --sensitive quick --refinepara quick
```

### High-Quality Filtering with Custom Thresholds

```bash
BASALT -a as1.fa \
    -s sample_R1.fq,sample_R2.fq \
    -t 60 -m 250 \
    --min-cpn 50 --max-ctn 10
```

### Short + Long Reads, Sensitive Mode, CheckM2

```bash
BASALT -a as1.fa,as2.fa \
    -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq \
    -l lr1.fastq,lr2.fastq \
    -t 64 -m 256 \
    --sensitive more-sensitive --refinepara deep \
    -q checkm2
```

### Autobinning Only (Skip Refinement + Reassembly)

```bash
BASALT -a as1.fa,as2.fa \
    -s s1_R1.fq,s1_R2.fq/s2_R1.fq,s2_R2.fq \
    -t 60 -m 250 \
    --module autobinning
```

### Refinement Only on Existing Bins

```bash
BASALT -a assembly.fa \
    -s sample_R1.fq,sample_R2.fq \
    -r My_MAGs_Folder \
    -c Coverage_matrix.txt \
    -t 32 -m 128
```

### Data Feeding: Import External Bins

```bash
BASALT -s sample_R1.fq,sample_R2.fq \
    -d vamb_bins,metabat_bins \
    --binset-index 1 \
    -t 32 -m 128
```

### With Extra Binners

```bash
BASALT -a as1.fa \
    -s sample_R1.fq,sample_R2.fq \
    -t 32 -m 128 \
    -e m,l \
    --sensitive more-sensitive
```
