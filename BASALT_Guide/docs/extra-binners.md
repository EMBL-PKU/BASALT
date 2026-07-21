# Extra binners

BASALT can add candidate bins from MetaBinner, VAMB, or LorBin to the default candidate pool. These adapters are optional and have independent dependencies, interfaces, and failure modes.

| Flag | Adapter | Intended evidence source |
|---|---|---|
| `-e m` | MetaBinner | Sequence composition and abundance profiles |
| `-e v` | VAMB | Variational autoencoder representation and abundance |
| `-e l` | [LorBin–BASALT Extra-binner](https://github.com/PKU-EMBL/LorBin-BASALT-Extrabinner) | Long-read-oriented binning through the BASALT-compatible in-house fork |

Combine flags with commas, for example `-e m,v`. An adapter runs in addition to the binners selected by `--sensitive`.

## Scientific use

Adding a binner expands the candidate set presented to BASALT selection. It does not guarantee more non-redundant, high-quality, or biologically correct MAGs. Choose an adapter because its assumptions fit the data, not solely to maximize candidate count.

For every optional binner, record:

- software version and installation source;
- complete command and parameters;
- input assembly and coverage files;
- whether CPU or GPU execution was used;
- adapter warnings or failures;
- candidate count before BASALT selection;
- the software citation.

## Failure semantics

The current CheckM2 orchestration catches failures from optional adapters, writes a warning, and continues with other candidates. A BASALT run can therefore finish even when an optional binner failed.

Search the logs explicitly:

```bash
grep -E '\[WARNING\]|\[ERROR DETAILS\]|Extra binner' \
  Basalt_log.txt basalt.stderr.log
```

Confirm that the expected adapter output directory exists and contains non-empty FASTA files. Do not report an optional binner as part of the method if it failed before contributing candidates.

## MetaBinner

MetaBinner requires its own executable workflow and a compatible `run_metabinner.sh` installation. BASALT derives a coverage profile and k-mer table before invoking the adapter.

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -e m \
  -t 32 -m 128 --mode new -o study_metabinner
```

Verify the MetaBinner installation outside BASALT before using this route. The adapter currently discovers installation paths through the active Conda environments, which may be site-dependent.

## VAMB

The repository contains an adapter for the VAMB 5 command form. VAMB uses mapped BAM files produced during BASALT preprocessing. VAMB is not part of the validated base environment: VAMB 5.0.4 pins PyTorch 2.6.0, whereas the current BASALT model stack uses the CPU-generic PyTorch 2.4 release. A plain `pip install vamb` inside `basalt` can therefore replace the tested PyTorch and scientific-Python packages.

Create an isolated VAMB environment and verify it independently:

```bash
micromamba create -n basalt-vamb python=3.12 pip --yes
micromamba run -n basalt-vamb \
  python -m pip install 'vamb==5.0.4'
micromamba run -n basalt-vamb vamb --help
```

The safest route is to run VAMB independently, retain its version and environment export, and import its FASTA bins through [data feeding](#manual-recovery-through-data-feeding). The direct `-e v` adapter requires a validated `vamb` executable on the BASALT process `PATH`; a separately activated environment is not automatically visible. Do not expose it by installing or upgrading packages inside the validated `basalt` environment. If a site administrator provides a wrapper to the isolated executable, test that wrapper on a small dataset and record it as part of the workflow.

Run the adapter:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -e v \
  -t 32 -m 128 --mode new -o study_vamb
```

VAMB GPU behaviour depends on its isolated PyTorch installation and local accelerator stack. Record `vamb`, Python, PyTorch, CUDA, and driver versions when GPU execution is used.

## LorBin

LorBin is intended for long-read metagenomic binning. For `-e l`, the maintained integration target is [PKU-EMBL/LorBin-BASALT-Extrabinner](https://github.com/PKU-EMBL/LorBin-BASALT-Extrabinner), an in-house fork of the upstream `LorMeBioAI/LorBin` repository. It retains the `LorBin` executable and algorithmic workflow while adapting the package to the BASALT Python 3.12 environment and carrying BASALT-oriented fixes for argument parsing, marker evaluation, runtime name errors, and related integration defects.

:::{important}
Do not report this route only as “LorBin.” Record that the BASALT-compatible in-house fork was used, together with its Git commit. The LorBin paper describes the method; the fork repository and commit identify the implementation that actually produced the candidate bins.
:::

### Install into the BASALT environment

The fork declares Python 3.12 and exposes `LorBin=lorbin.lorbin:main`. Its repository documents installation on top of the BASALT environment. Clone and pin the source before installing it:

```bash
git clone https://github.com/PKU-EMBL/LorBin-BASALT-Extrabinner.git
cd LorBin-BASALT-Extrabinner
git rev-parse HEAD | tee lorbin-basalt-commit.txt

micromamba run -n basalt \
  python -m pip install --no-deps .

micromamba run -n basalt LorBin --help
```

Conda users can replace `micromamba run` with `conda run`. `--no-deps` prevents pip from silently replacing the solved BASALT packages; the repository's `setup.py` does not currently declare Python-package dependencies. The BASALT environment already supplies the required sequence-processing tools, but its transitive PyTorch and scientific-Python versions can differ from the fork's standalone reference environment. Verify the exact resolved stack rather than assuming compatibility:

```bash
micromamba run -n basalt python -c \
  'import Bio, joblib, numpy, pandas, scipy, sklearn, torch; \
print("Biopython", Bio.__version__); print("NumPy", numpy.__version__); \
print("pandas", pandas.__version__); print("SciPy", scipy.__version__); \
print("scikit-learn", sklearn.__version__); print("joblib", joblib.__version__); \
print("PyTorch", torch.__version__)'

micromamba run -n basalt \
  sh -c 'command -v LorBin minimap2 samtools hmmsearch prodigal bedtools'
```

Run the fork's public test dataset before a production analysis. The repository links the test data from [Zenodo record 13883404](https://zenodo.org/records/13883404).

If the fork cannot coexist with the validated BASALT environment, use its standalone `lorbin_env.yaml`, run LorBin independently, and import the resulting FASTA bins through [data feeding](#manual-recovery-through-data-feeding). A separately activated LorBin environment is not automatically visible to a BASALT process running in another environment.

### Run through BASALT

The adapter depends on the installed `LorBin` executable and compatible BAM files generated during BASALT preprocessing.

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -l nanopore.fastq \
  -e l \
  -t 32 -m 128 --mode new -o study_lorbin
```

The adapter invokes the fork's `LorBin bin` command, then normalizes returned FASTA files into the BASALT candidate-bin naming scheme. Confirm that `Basalt_log.txt` records successful completion and that the corresponding `*_LorBin_genomes/` directory contains non-empty FASTA files. A final BASALT result is not evidence that LorBin succeeded because optional-binner failures are warning-only.

### Report the implementation

A minimal Methods statement is:

> Long-read-oriented candidate bins were generated with the PKU-EMBL LorBin–BASALT Extra-binner in-house fork (package version 0.1.0; Git commit [SHA]) through the BASALT `-e l` adapter. Candidate FASTA files were subsequently evaluated within the BASALT selection and refinement workflow.

Replace the version and commit with the installed values. Do not imply that all fork-generated candidates were retained in the final bin set.

## Manual recovery through data feeding

If an adapter fails but the underlying tool can be run independently, preserve its bin FASTA files and import them through BASALT data feeding. Contig identifiers must remain compatible with the BASALT assembly.

```bash
BASALT \
  -s sample_R1.fastq,sample_R2.fastq \
  -d external_binner_bins \
  --binset-index 1 \
  -t 32 -m 128 -q checkm2 -o imported
```

The CheckM2 data-feeding route writes `imported_data_feeded/` for this example.

## Citation

Cite the exact optional binner used. For the LorBin adapter, cite the algorithm paper and identify the forked code version:

> Xue, W. *et al.* LorBin: efficient binning of long-read metagenomes by multiscale adaptive clustering and evaluation. *Nature Communications* **16**, 9353 (2025).

> LorBin source-code archive. [https://doi.org/10.5281/zenodo.13864645](https://doi.org/10.5281/zenodo.13864645). Record the corresponding [LorBin–BASALT Extra-binner](https://github.com/PKU-EMBL/LorBin-BASALT-Extrabinner) Git commit used for the BASALT run.

Use each software project's recommended citation and version-specific documentation for MetaBinner and VAMB.
