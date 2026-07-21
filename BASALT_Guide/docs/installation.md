# Installation

This page installs the Conda-based BASALT repository. BASALT depends on several compiled bioinformatics programs, two quality-control databases, and trained model files. A successful Python import alone is therefore not a complete installation test.

## Platform and resources

| Requirement | Supported or practical value |
|---|---|
| Operating system | Linux x86-64 |
| Python | 3.12 |
| CPU | 8 threads minimum; 32 or more recommended for routine datasets |
| Memory | 128 GB practical starting point; 256 GB or more for large or complex runs |
| Storage | Dataset-dependent; allow space for mappings, candidate bins, reassemblies, and archives |

macOS and Windows are not supported as native execution platforms. Use a Linux host, cluster, virtual machine, or compatible container runtime.

## Conda installation

### 1. Clone the repository

```bash
git clone https://github.com/PKU-EMBL/BASALT.git
cd BASALT
```

Record the version you installed:

```bash
git rev-parse HEAD
```

### 2. Create the environment

```bash
conda create -n basalt_env -c conda-forge -c bioconda \
  python=3.12 \
  megahit metabat2 maxbin2 concoct prodigal semibin \
  bedtools blast bowtie2 bwa diamond checkm2 \
  unicycler spades samtools racon pilon ncbi-vdb \
  minimap2 miniasm idba hmmer entrez-direct biopython uv \
  --yes

conda activate basalt_env
```

Install the Python dependencies into the active environment:

```bash
uv pip install \
  torch tensorboardx pillow scikit-learn \
  numpy==1.26.4 scipy pandas matplotlib tqdm requests
```

`basalt_environment.yml` is the sole repository environment definition and is provided as an alternative to the explicit command above. It currently uses Tsinghua mirror channels. Review and replace the channel URLs if those mirrors are unsuitable for your site.

### 3. Install the BASALT command

Run the installer while `basalt_env` is active:

```bash
bash install.sh
```

The script copies the BASALT programs into `$CONDA_PREFIX/bin` and creates the `BASALT` launcher. Re-run the installer after updating the repository because this is a copied installation, not an editable package.

### 4. Download the trained models

Choose a persistent model directory and keep the path free of shell metacharacters:

```bash
BASALT_models_download.py --path "$PWD/BASALT_WEIGHT"
export BASALT_WEIGHT="$PWD/BASALT_WEIGHT"
```

Add the export to your shell initialization file only after verifying the path. The extracted directory should contain five `*_ensemble.csv` files and their associated model checkpoints.

```bash
find "$BASALT_WEIGHT" -maxdepth 1 -name '*_ensemble.csv' | wc -l
```

Expected output: `5`.

### 5. Configure CheckM2

CheckM2 is the default quality-control backend:

```bash
checkm2 database --download
```

If CheckM2 does not discover the database automatically, set the path to its DIAMOND database:

```bash
export CHECKM2DB=/absolute/path/to/CheckM2_database/uniref100.KO.1.dmnd
```

The legacy CheckM backend is selected with `-q checkm`. CheckM can impose a different Python and dependency stack from the maintained CheckM2 route. Validate it in a compatible environment or container, install its database, and set the location required by that installation, commonly `CHECKM_DATA_PATH`.

### 6. Verify the executable stack

```bash
BASALT --help

command -v \
  BASALT checkm2 metabat2 SemiBin2 bowtie2 bwa samtools \
  minimap2 spades.py unicycler blastn makeblastdb
```

Then record the environment for provenance:

```bash
conda env export --no-builds > basalt-environment.yml
```

:::{warning}
The Conda edition documented here does not implement `BASALT --version` or `BASALT --check-deps`. Those options belong to BASALT-Air. Use `BASALT --help`, `command -v`, explicit tool version commands, and the tutorial smoke test.
:::

## Singularity execution

The upstream repository has historically described a prebuilt `basalt.sif` image, but its former Google Drive endpoint returned HTTP 404 during this documentation update. Use the Conda installation above unless you already possess a trusted image or the maintainers publish a new verified artifact. Follow the [upstream repository](https://github.com/PKU-EMBL/BASALT) for distribution updates.

If you already have an image, record where it came from and verify its checksum before production use.

Verify the image before production use:

```bash
singularity run basalt.sif BASALT --help
singularity run basalt.sif checkm2 --help
singularity inspect basalt.sif > singularity-inspect.txt
```

Run from a dedicated directory and bind all required storage explicitly:

```bash
cd /project/basalt_run

singularity run \
  --bind /project:/project \
  basalt.sif \
  BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -t 32 -m 128 --mode new -o study_basalt
```

Preserve an image checksum with each analysis:

```bash
sha256sum basalt.sif > basalt.sif.sha256
```

On a cluster, make the bind explicit even if the runtime normally auto-binds the current directory. Confirm that inputs, the working directory, databases, and `BASALT_WEIGHT` are visible inside the container before starting a long job. Run under the site scheduler or a persistent session, capture stdout and stderr, and do not rely on terminal persistence as a substitute for checkpoint and output auditing.

## China mainland

The repository includes `basalt_environment.yml` with Tsinghua mirror channels. Model files are also mirrored on [Baidu Netdisk](https://pan.baidu.com/s/1ouKqabxHYr1GmvpquQCzqw?pwd=embl) with extraction code `embl`.

Do not mix packages from several mirror and upstream channel stacks within the same solved environment unless necessary. Export the final environment and verify each external executable after installation.

## Updating BASALT

The installer copies scripts into the environment. Update and reinstall them explicitly:

```bash
cd /path/to/BASALT
git pull --ff-only
git rev-parse HEAD
conda activate basalt_env
bash install.sh
```

Start new analyses in a new working directory after an update. Resuming a checkpoint with different code, model weights, databases, or external-tool versions weakens reproducibility and may produce incompatible intermediates.

## Next step

Run the [quick start](quickstart.md), then the [tutorial](tutorial.md) with the public demo dataset.
