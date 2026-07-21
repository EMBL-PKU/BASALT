# Installation

This page installs the Conda-compatible BASALT repository. BASALT depends on several compiled bioinformatics programs, quality-control databases, and trained model files. A successful Python import alone is therefore not a complete installation test.

## Platform and resources

| Requirement | Supported or practical value |
|---|---|
| Operating system | Linux x86-64 |
| Python | 3.12 |
| CPU | 8 threads minimum; 32 or more recommended for routine datasets |
| Memory | 128 GB practical starting point; 256 GB or more for large or complex runs |
| Storage | Dataset-dependent; allow space for mappings, candidate bins, reassemblies, and archives |

macOS and Windows are not supported as native execution platforms. Use a Linux host, cluster, virtual machine, or compatible container runtime.

## Conda-compatible installation

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

Micromamba is the recommended package manager because it uses the libmamba solver, downloads packages in parallel, and does not require a Python-equipped base environment:

```bash
micromamba create -n basalt -f basalt_environment.yml \
  --strict-channel-priority --yes
```

Current Conda is an equivalent fallback:

```bash
CONDA_CHANNEL_PRIORITY=strict \
  conda env create -n basalt -f basalt_environment.yml \
  --solver libmamba --yes
```

`basalt_environment.yml` is the sole repository environment definition. It uses only `conda-forge` and `bioconda`; the commands above enforce strict channel priority without editing global configuration. The file pins NumPy for the trained-model stack and omits redundant direct dependencies. The legacy CheckM package is intentionally excluded because its dependency stack conflicts with the maintained Python 3.12/CheckM2 route; install it in a separate compatible environment only when `-q checkm` is required.

If package access is slow or unreliable in mainland China, use the [network-aware installer](#china-mainland) rather than changing this portable YAML or the global Conda configuration.

### 3. Install the BASALT command

Install directly through the selected manager without depending on shell activation:

```bash
micromamba run -n basalt bash install.sh
# or: conda run -n basalt bash install.sh
```

The script copies the BASALT programs into `$CONDA_PREFIX/bin` and creates the `BASALT` launcher. Re-run the installer after updating the repository because this is a copied installation, not an editable package.

Optional extra binners are not part of the base installation. In particular, `-e l` uses the [PKU-EMBL/LorBin-BASALT-Extrabinner](https://github.com/PKU-EMBL/LorBin-BASALT-Extrabinner) in-house fork, which must be installed and validated separately. Do not add it to the portable BASALT environment YAML without re-solving and testing the complete stack; follow the pinned-source procedure in [Extra binners](extra-binners.md#lorbin).

### 4. Download the trained models

Choose a persistent model directory and keep the path free of shell metacharacters. Automatic mode first uses the official [PKU-EMBL/BASALT_WEIGHT repository](https://huggingface.co/PKU-EMBL/BASALT_WEIGHT), whose client supports concurrent and resumable downloads, and then falls back to the legacy Figshare ZIP:

```bash
micromamba run -n basalt \
  BASALT_models_download.py \
  --source auto \
  --path "$PWD/BASALT_WEIGHT"

export BASALT_WEIGHT="$PWD/BASALT_WEIGHT"
```

Conda users can replace `micromamba run` with `conda run`. A public Hugging Face repository does not require login. To force one source or use a file obtained through another route:

```bash
# Official Hugging Face only
BASALT_models_download.py --source huggingface --path "$BASALT_WEIGHT"

# Existing ZIP, including a manually downloaded Baidu Netdisk copy
BASALT_models_download.py --source archive \
  --archive /path/to/BASALT.zip \
  --path "$BASALT_WEIGHT"

# Site-operated object storage or another trusted URL
BASALT_models_download.py --source url \
  --url https://example.org/BASALT.zip \
  --path "$BASALT_WEIGHT"
```

Add the export to your shell initialization file only after verifying the path. The downloader checks both the five top-level `*_ensemble.csv` descriptors and their corresponding checkpoint directories.

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
micromamba env export -n basalt > basalt-environment.yml
# or: conda env export -n basalt --no-builds > basalt-environment.yml
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

`BASALT_setup_China_mainland.py` handles the environment, script installation, and model weights as one recoverable workflow. It prefers an installed `micromamba`, then `mamba`, then Conda. The following command downloads the official standalone micromamba only when micromamba is unavailable, stores its shared package cache under `~/micromamba`, probes several mirror endpoints, and downloads weights from Hugging Face first:

```bash
python3 BASALT_setup_China_mainland.py \
  --manager micromamba \
  --bootstrap-micromamba \
  --mirror auto \
  --env-name basalt
```

The bootstrapped binary is placed at `~/.local/bin/micromamba`; the installer does not initialize or edit the shell. It also does not run `conda config`, edit `~/.condarc`, edit pip configuration, or change environment files to mode `777`. Activate later with the path printed at completion, or continue using `micromamba run -n basalt COMMAND`.

Inspect the exact commands without changing the system:

```bash
python3 BASALT_setup_China_mainland.py \
  --manager micromamba \
  --bootstrap-micromamba \
  --mirror tuna \
  --model-source none \
  --dry-run
```

### Mirror selection

| Option | Conda channels | Matching temporary PyPI index | Intended use |
|---|---|---|---|
| `--mirror auto` | probes TUNA, BFSU, USTC, then upstream | follows the selected preset | default for variable networks |
| `--mirror tuna` | TUNA `conda-forge` and `bioconda` | TUNA PyPI | stable manual selection |
| `--mirror bfsu` | BFSU `conda-forge` and `bioconda` | BFSU PyPI | alternative education-network route |
| `--mirror ustc` | USTC dynamic caches for `conda-forge` and `bioconda` | USTC PyPI | community-channel-only environment |
| `--mirror upstream` | official `conda.anaconda.org` | inherited/default PyPI | VPN, proxy, or international route |

The presets follow the current [TUNA](https://mirrors.tuna.tsinghua.edu.cn/help/anaconda/), [BFSU](https://mirrors.bfsu.edu.cn/help/anaconda/), and [USTC](https://mirrors.ustc.edu.cn/help/anaconda.html) Anaconda mirror guidance. USTC currently documents dynamic caches for `conda-forge` and `bioconda`, which is sufficient because the maintained BASALT environment no longer depends on `defaults` or the historical `esteinig` channel.

For an institutional mirror, repeat `--channel` in decreasing-priority order and optionally set the PyPI endpoint:

```bash
python3 BASALT_setup_China_mainland.py \
  --channel https://mirror.example.edu/anaconda/cloud/conda-forge \
  --channel https://mirror.example.edu/anaconda/cloud/bioconda \
  --pip-index-url https://mirror.example.edu/pypi/simple
```

Do not combine packages from different channel stacks in one solved environment. If a mirror is reachable but lacks a required artifact, rerun with `--update --mirror OTHER` after confirming the selected route, or create a new environment under a different name. Existing package caches can be used with `--offline`.

### Weight fallbacks

Automatic mode uses [Hugging Face](https://huggingface.co/PKU-EMBL/BASALT_WEIGHT) first and Figshare second. To create only the software environment, pass `--model-source none`. Model files are additionally available from [Baidu Netdisk](https://pan.baidu.com/s/1ouKqabxHYr1GmvpquQCzqw?pwd=embl) with extraction code `embl`; after manual download, keep the original ZIP and run:

```bash
python3 BASALT_setup_China_mainland.py \
  --model-source archive \
  --model-archive /absolute/path/BASALT.zip \
  --model-dir /persistent/path/BASALT_WEIGHT
```

The same local archive works independently of the installer:

```bash
BASALT_models_download.py \
  --source archive \
  --archive /absolute/path/BASALT.zip \
  --path /persistent/path/BASALT_WEIGHT
```

Use `--hf-endpoint` only for a trusted Hugging Face-compatible service operated or approved by your institution. Proxy variables such as `HTTPS_PROXY` are inherited automatically.

After installation, export the final environment, record the selected mirror, and verify every external executable. Database downloads are separate from package and model installation; run `checkm2 database --download` through the network route approved for the compute site.

## Updating BASALT

The installer copies scripts into the environment. Update and reinstall them explicitly:

```bash
cd /path/to/BASALT
git pull --ff-only
git rev-parse HEAD
python3 BASALT_setup_China_mainland.py \
  --source-dir . \
  --env-name basalt \
  --update \
  --model-source none
```

Start new analyses in a new working directory after an update. Resuming a checkpoint with different code, model weights, databases, or external-tool versions weakens reproducibility and may produce incompatible intermediates.

## Next step

Run the [quick start](quickstart.md), then the [tutorial](tutorial.md) with the public demo dataset.
