# Release notes

Release notes summarize intended changes. For reproducible analysis, also record the exact Git commit, model files, quality-control database, and dependency environment.

## Unreleased

This development snapshot hardens the code paths audited after the `V1.2.2` tag. It must receive a new release tag before these changes are described as part of a published version.

### Fixed

- Corrected MetaBinner, VAMB 5, and LorBin command construction, executable discovery, sorted-BAM inputs, output-table parsing, and normalized candidate-folder names.
- Added strict `BASALT_WEIGHT` validation down to every descriptor-referenced checkpoint, rejected unsafe descriptor/archive paths, and removed the unpinned legacy model fallback.
- Added CLI validation for positive resource values, percentage ranges, paired-read grammar, optional-binner identifiers, and route-specific input requirements.
- Fixed undefined variables and imports identified by Python 3.12 compilation and Ruff critical-error checks in comparison, retrieval, mapping, dataset, and dereplication paths.
- Made the installed launcher independent of the caller's active `$CONDA_PREFIX` and changed installation to fail clearly when MetaBAT 2 does not supply `jgi_summarize_bam_contig_depths`.
- Removed an ineffective Pilon probe and redundant post-reassembly moves that could emit errors or reference an undefined final loop item.

### Repository maintenance

- Removed the obsolete general-purpose `BASALT_setup.py`, duplicate `BASALT/BASALT` launcher, tracked bytecode and Finder metadata, and the stale bundled x86-64 `jgi_summarize_bam_contig_depths` binary.
- Added `unzip` to the canonical environment because compressed read inputs use it at runtime, and removed the unused direct Pilon dependency.
- Standardized installed executable permissions on mode `0755`; world-writable `0777` permissions are neither required nor recommended.
- Added Linux CI for Python, Bash and Perl syntax, critical Ruff diagnostics, release-smoke tests, launcher isolation, and permission behavior.

## BASALT 1.2.2

Released 21 July 2026.

BASALT 1.2.2 improves installation reliability, model reproducibility, and documentation. It retains the Python 3.12 workflow introduced in 1.2.1 while preventing optional packages from silently replacing the validated BASALT runtime.

### Highlights

- Resolves the Python 3.12 base runtime in one Conda transaction using the CPU-generic PyTorch stack.
- Downloads models from `PKU-EMBL/BASALT_WEIGHT` with a pinned revision, resumable transfers, connection timeouts, and documented fallback routes.
- Improves mainland-China installation with micromamba and selectable TUNA, BFSU, USTC, upstream, or institutional mirrors without modifying global Conda or pip configuration.
- Keeps VAMB and other conflicting optional adapters outside the validated base environment.
- Expands the README and Read the Docs guide, including LorBin integration, multi-assembly study design, troubleshooting, reproducibility, and BASALT-Air guidance.

### Fixed

- Removed the duplicate uv/pip scientific-stack transaction that could download a second PyTorch build and large CUDA dependencies.
- Added `--hf-timeout` so blocked Hugging Face routes fail promptly instead of hanging indefinitely.
- Pinned the default model revision to `bc98b102522d1c80dd8c2594df4ab3155438320e`.
- Isolated VAMB because VAMB 5.0.4 requires PyTorch 2.6, whereas the validated BASALT model runtime uses CPU-generic PyTorch 2.4.
- Standardized model discovery around one validated absolute `BASALT_WEIGHT` path, preferably configured once in `~/.bashrc`.
- Removed obsolete channel requirements and separated legacy CheckM from the maintained CheckM2/Python 3.12 route.

### Validation

The release was installed and tested on Ubuntu 22.04 x86-64 with micromamba. The BASALT launcher, Python stack, CheckM2 database, required executables, 81 model checkpoints, five ensemble descriptors, and CPU model loading all passed validation. The complete documentation also passed strict Sphinx and external-link checks. These checks do not constitute a new biological-performance benchmark.

### Compatibility and upgrade notes

- Recreate or update the environment from `basalt_environment.yml`, then rerun `install.sh` because BASALT scripts are copied into the environment.
- Keep VAMB, LorBin, MetaBinner, legacy CheckM, and GPU-specific stacks isolated unless independently validated.
- Record the BASALT commit, model revision and checksums, CheckM2 database, environment export, commands, and warnings for production analyses.
- Start updated analyses in a new working directory rather than resuming checkpoints created with a different software stack.

## BASALT 1.2.1

Tagged 15 May 2026.

### Changes

- Updated the VAMB adapter for the VAMB 5 command family.
- Added VAMB 5 to the Python 3.12 environment definition.
- Added warning-based failure handling for optional binners so other candidate generators can continue.

### Compatibility boundaries

- Optional-adapter success must be verified from logs and output FASTA files.
- VAMB, PyTorch, CUDA, and driver versions remain part of the analysis provenance.
- The Conda edition does not implement `BASALT --version`; record the Git tag and commit directly.

## BASALT 1.2.0

Released 16 December 2025.

### Changes

- Updated the documented environment to Python 3.12.
- Added a LorBin adapter under `-e l`.
- Selected CheckM2 as the default quality-control backend.
- Added accelerator-compatible dependency paths for SemiBin 2 and model inference.
- Made the model directory configurable through `BASALT_WEIGHT`.
- Renamed the historical `--autopara` option to `--sensitive`.

### Compatibility boundaries

- The Conda edition does not implement `BASALT --version` or `BASALT --check-deps`.
- Optional binner adapters require independent installation and validation.
- Resuming checkpoints created with earlier code or database versions is not recommended.

## BASALT 1.1.0

Released 12 June 2024.

- Revised installation workflows.
- Added China-mainland mirror support.
- Added a Singularity distribution route.
- Included maintenance and stability changes.

## BASALT 1.0.0

Released 18 August 2023.

- Introduced multi-assembly candidate generation and bin selection.
- Introduced model-based contig screening and read-supported retrieval.
- Added bin-level reassembly and OLC stages.
- Added external-bin data feeding.
- Added CheckM-based quality assessment and textual checkpoints.

## Related project: BASALT-Air 1.0.0

BASALT-Air v1.0.0 is now available from the separate [BASALT-Air repository](https://github.com/PKU-EMBL/BASALT-Air). See [BASALT-Air: the new implementation](basalt-air.md) for the implementation comparison and migration boundary.

BASALT-Air uses Pixi and adds first-class absolute paths, separate work and output directories, and implementation-specific dependency checks. Its CLI and dependency model differ from the Conda-based BASALT documented on this site.

Do not infer byte-identical output between implementations without a version-matched regression test.

## Publication

> Qiu, Z., Yuan, L., Lian, C.-A. *et al.* BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. *Nature Communications* **15**, 2179 (2024). [https://doi.org/10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7)
