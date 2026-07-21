# BASALT: Binning Across a Series of Assemblies Toolkit

[![License: MIT](https://img.shields.io/badge/license-MIT-2f855a.svg)](LICENSE)
[![Python 3.12](https://img.shields.io/badge/Python-3.12-3776ab.svg)](https://www.python.org/)
[![Nature Communications](https://img.shields.io/badge/Nature%20Communications-2024-b31b1b.svg)](https://doi.org/10.1038/s41467-024-46539-7)
[![Documentation](https://img.shields.io/badge/docs-MkDocs-526cfe.svg)](BASALT_Guide/docs/index.md)

BASALT is a metagenomic workflow for recovering and refining metagenome-assembled genomes (MAGs) across one or more assemblies. It combines multi-binner candidate generation, within- and between-assembly bin selection, contig-level contamination screening, contig retrieval, and reassembly in one checkpointed pipeline.

The published evaluation reports improved MAG recovery and downstream genome-resolved analyses relative to the tested workflows and datasets. These results define the evidence boundary: performance depends on community complexity, sequencing depth, assembly quality, read type, and the selected parameters. See the [Nature Communications article](https://doi.org/10.1038/s41467-024-46539-7) for the complete experimental design.

<p align="center">
  <img src="fig/workflow.png" alt="BASALT workflow from assemblies and reads to refined MAGs" width="86%">
</p>

## Documentation

The [BASALT guide](BASALT_Guide/docs/index.md) is the source of truth for installation, inputs, commands, outputs, troubleshooting, and reproducible reporting.

- [Installation](BASALT_Guide/docs/installation.md)
- [Quick start](BASALT_Guide/docs/quickstart.md)
- [Inputs and study design](BASALT_Guide/docs/inputs.md)
- [Command-line reference](BASALT_Guide/docs/usage.md)
- [Tutorial](BASALT_Guide/docs/tutorial.md)
- [Reproducibility and reporting](BASALT_Guide/docs/reproducibility.md)

## What BASALT does

1. **Autobinning.** BASALT runs MetaBAT 2 and SemiBin 2 across parameter settings. The `sensitive` and `more-sensitive` presets add CONCOCT, and `more-sensitive` also adds MaxBin 2.0.
2. **Bin selection.** Candidate bins are compared within each assembly and then across assemblies to reduce redundancy.
3. **Refinement.** A five-model multilayer-perceptron ensemble screens contigs using sequence-composition and coverage-derived features. Read connectivity is then used to retrieve candidate contigs.
4. **Reassembly.** BASALT can reassemble and polish selected bins when the required read types are available.

Optional adapters support MetaBinner (`-e m`), VAMB (`-e v`), and LorBin (`-e l`). Their installation and input requirements are documented separately in the [extra-binner guide](BASALT_Guide/docs/extra-binners.md).

## Supported input designs

Every run requires at least one assembly and at least one source of read coverage. BASALT accepts the following designs:

| Design | Assembly | Short reads | ONT or PacBio CLR | PacBio HiFi |
|---|---:|---:|---:|---:|
| Short-read only | required | required | no | no |
| Long-read only | required | no | required | no |
| HiFi only | required | no | no | required |
| Short reads + ONT or CLR | required | required | required | no |
| Short reads + HiFi | required | required | no | required |

Input paths are resolved from the working directory. Use filenames or relative paths without shell metacharacters; for the most reliable execution, place or symlink inputs into a dedicated run directory. Paired-end files within a sample are separated by a comma, and multiple short-read samples are separated by `/`.

## Requirements

- Linux x86-64
- Python 3.12
- 8 CPU threads and 128 GB RAM as a practical starting point
- 32 or more threads and 256 GB RAM for larger or more complex datasets
- Sufficient temporary storage for mapping, binning, and reassembly intermediates

Resource use is dataset-dependent. Run a pilot or request more memory for large assemblies, many samples, or `more-sensitive` mode.

## Installation

The Conda installation below is the maintained route for this repository. A prebuilt Singularity image is also available; see the [installation guide](BASALT_Guide/docs/installation.md) for databases, model weights, container use, and verification.

```bash
git clone https://github.com/PKU-EMBL/BASALT.git
cd BASALT

conda create -n basalt_env -c conda-forge -c bioconda \
  python=3.12 megahit metabat2 maxbin2 concoct prodigal semibin \
  bedtools blast bowtie2 bwa diamond checkm2 \
  unicycler spades samtools racon pilon ncbi-vdb \
  minimap2 miniasm idba hmmer entrez-direct biopython uv --yes

conda activate basalt_env
uv pip install torch tensorboardx pillow scikit-learn \
  numpy==1.26.4 scipy pandas matplotlib tqdm requests

bash install.sh
```

Download the trained BASALT models and expose their directory:

```bash
python BASALT_models_download.py --path "$PWD/BASALT_WEIGHT"
export BASALT_WEIGHT="$PWD/BASALT_WEIGHT"
```

Set the CheckM2 database path if it is not already configured by your environment:

```bash
checkm2 database --download
export CHECKM2DB=/absolute/path/to/CheckM2_database/uniref100.KO.1.dmnd
```

Verify the command and external tools before a production run:

```bash
BASALT --help
command -v checkm2 metabat2 SemiBin2 bowtie2 samtools spades.py
```

## Quick start

Create a clean working directory and link the inputs into it:

```bash
mkdir -p basalt_run
cd basalt_run
ln -s /absolute/path/assembly.fasta .
ln -s /absolute/path/sample_R1.fastq .
ln -s /absolute/path/sample_R2.fastq .
```

Run the full short-read workflow:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -t 32 \
  -m 128 \
  --mode new \
  -o study_basalt
```

The final MAGs and final CheckM2 report are written under `study_basalt/`. Intermediate files, `BASALT_command.txt`, `Basalt_log.txt`, and `Basalt_checkpoint.txt` remain in the working directory.

To resume an interrupted run from that same directory:

```bash
BASALT --mode continue
```

Do not combine files from different runs in one working directory. `--mode new` initializes checkpoint state, but it does not provide complete cleanup or provenance isolation.

## Common commands

Short reads across two assemblies and two samples:

```bash
BASALT \
  -a assembly_1.fasta,assembly_2.fasta \
  -s sample_1_R1.fastq,sample_1_R2.fastq/sample_2_R1.fastq,sample_2_R2.fastq \
  -t 64 -m 256 --mode new -o study_basalt
```

Short reads plus ONT or PacBio CLR reads:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -l nanopore.fastq \
  -t 64 -m 256 --mode new -o study_basalt
```

PacBio HiFi reads:

```bash
BASALT \
  -a assembly.fasta \
  -hf hifi.fastq \
  -t 32 -m 128 --mode new -o study_basalt
```

Faster candidate generation and refinement:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  --sensitive quick --refinepara quick \
  -t 32 -m 128 --mode new -o study_basalt
```

See the [command-line reference](BASALT_Guide/docs/usage.md) before using standalone refinement, data feeding, dereplication, or module-specific execution. These routes require intermediate files from a compatible BASALT run.

## Reproducible use

For each analysis, archive:

- the exact BASALT command and BASALT Git commit or release;
- the Conda environment export and versions of external executables;
- CheckM or CheckM2 database identity and quality thresholds;
- input assembly and read checksums;
- `BASALT_command.txt`, `Basalt_log.txt`, and `Basalt_checkpoint.txt`;
- the final quality report and all exclusions applied after BASALT.

The [reporting guide](BASALT_Guide/docs/reproducibility.md) provides a Methods template and a minimum reporting checklist.

## Citation

If BASALT contributes to an analysis, cite:

> Qiu, Z., Yuan, L., Lian, C.-A. *et al.* BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. *Nature Communications* **15**, 2179 (2024). [https://doi.org/10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7)

```bibtex
@article{qiu2024basalt,
  author  = {Qiu, Zhiguang and Yuan, Li and Lian, Chun-Ang and Lin, Bin and
             Chen, Jie and Mu, Rong and Qiao, Xuejiao and Zhang, Liyu and
             Xu, Zheng and Fan, Lu and Zhang, Yunzeng and Wang, Shanquan and
             Li, Junyi and Cao, Huiluo and Li, Bing and Chen, Baowei and
             Song, Chi and Liu, Yongxin and Shi, Lili and Tian, Yonghong and
             Ni, Jinren and Zhang, Tong and Zhou, Jizhong and Zhuang, Wei-Qin and
             Yu, Ke},
  title   = {BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis},
  journal = {Nature Communications},
  volume  = {15},
  pages   = {2179},
  year    = {2024},
  doi     = {10.1038/s41467-024-46539-7}
}
```

If an optional binner materially contributes to the reported analysis, cite that software and record its version and parameters.

## Support and licence

Report reproducible bugs through [GitHub Issues](https://github.com/PKU-EMBL/BASALT/issues). Include the command, log, checkpoint, software versions, database configuration, and the smallest shareable input that reproduces the problem.

Scientific or collaboration enquiries:

- Ke Yu: [yuke.sz@pku.edu.cn](mailto:yuke.sz@pku.edu.cn)
- Zhaorui (Elijah) Jiang: [zrjiang25@stu.pku.edu.cn](mailto:zrjiang25@stu.pku.edu.cn)

BASALT is released under the [MIT License](LICENSE).
