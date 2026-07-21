# Release notes

Release notes summarize intended changes. For reproducible analysis, also record the exact Git commit, model files, quality-control database, and dependency environment.

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

BASALT-Air is maintained in the separate [BASALT-Air repository](https://github.com/PKU-EMBL/BASALT-Air).

BASALT-Air uses Pixi and adds first-class absolute paths, separate work and output directories, and implementation-specific dependency checks. Its CLI and dependency model differ from the Conda-based BASALT documented on this site.

Do not infer byte-identical output between implementations without a version-matched regression test.

## Publication

> Qiu, Z., Yuan, L., Lian, C.-A. *et al.* BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. *Nature Communications* **15**, 2179 (2024). [https://doi.org/10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7)
