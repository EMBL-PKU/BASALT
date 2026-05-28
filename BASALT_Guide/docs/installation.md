# Installation

## Conda Installation (Recommended)

### 1. Clone the Repository

```bash
git clone https://github.com/EMBL-PKU/BASALT.git
cd BASALT
```

### 2. Create the Conda Environment

```bash
conda create -n basalt_env -c conda-forge -c bioconda \
    python=3.12 \
    megahit metabat2 maxbin2 concoct prodigal semibin \
    bedtools blast bowtie2 diamond checkm2 \
    unicycler spades samtools racon pplacer pilon \
    ncbi-vdb minimap2 miniasm idba hmmer entrez-direct \
    biopython uv --yes

conda activate basalt_env
```

### 3. Install Python Packages

```bash
uv pip install tensorflow torch torchvision tensorboard tensorboardx \
    lightgbm scikit-learn numpy==1.26.4 python-igr \
    scipy pandas matplotlib cython biolib joblib tqdm requests checkm-genome
```

### 4. Download Model Weights

```bash
python BASALT_models_download.py --path "/path/to/model/folder"
```

### 5. Install BASALT Scripts

```bash
chmod +x install.sh
bash install.sh
chmod +x /path/to/basalt/bin/*
```

### 6. Set Up Environment Variables

Add the following to your `~/.bashrc`:

```bash
export CHECKM2DB=/path/to/checkm2db/CheckM2_database/uniref100.KO.1.dmnd
export CHECKM_DATA_PATH=/path/to/checkmdb
export BASALT_WEIGHT=/path/to/BASALT
```

Then reload:

```bash
source ~/.bashrc
```

The CheckM and CheckM2 databases, along with the latest Singularity image, are available from [Google Drive](https://drive.google.com/drive/folders/1d0e_2FpYRBAZLwKXl8fA-yDK4b5PBA_E?usp=sharing).

---

## Singularity Installation

For users who prefer containers, or for users in China mainland with network limitations, a prebuilt Singularity image is available.

### Singularity Image

The image (`basalt.sif`) includes all dependencies: CheckM, CheckM2, Semibin2, Bowtie2, BWA, etc.

**Run BASALT directly:**

```bash
# When basalt.sif is in your home directory
singularity run basalt.sif BASALT -a as1.fa \
    -s S1_R1.fq,S1_R2.fq/S2_R1.fq,S2_R2.fq \
    -t 32 -m 128
```

**With bind mount for custom paths:**

```bash
singularity run -B /media/emma basalt.sif BASALT -h
```

**Run in background with screen:**

```bash
screen -dmS basalt_job bash -c 'singularity run basalt.sif BASALT -a as1.fa -s reads_R1.fq,reads_R2.fq -t 32 -m 128 > log_basalt'
```

**Invoke individual tools inside the image:**

```bash
singularity run basalt.sif bowtie2 -h
singularity run basalt.sif checkm2 predict -h
singularity run basalt.sif samtools --help
```

---

## Installing from China Mainland

Chinese users experiencing slow network speeds can use mirror sources:

```bash
site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
conda config --add channels ${site}/pkgs/free/
conda config --add channels ${site}/pkgs/main/
conda config --add channels ${site}/cloud/conda-forge/
conda config --add channels ${site}/cloud/bioconda/
```

Alternatively, use the [Singularity image](#singularity-installation) or download model weights from alternative sources (see [Release Notes](Release_notes.md)).

---

## System Requirements

| Resource | Minimum | Recommended |
|---|---|---|
| Operating System | Linux x64 | Linux x64 |
| CPU Cores | 8 | 32+ |
| RAM | 128 GB | 256 GB+ |
| Storage | 100 GB free | 500 GB+ |
| Python | 3.12 | 3.12 |

---

## Verifying the Installation

To verify that BASALT is correctly installed:

```bash
conda activate basalt_env
BASALT -h
```

This should print the help message with all available options.

To test with the demo dataset, refer to the [Tutorial](tutorial.md).
