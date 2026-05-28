# Release Notes

## BASALT-Air v1.0.0 (2026-04-16)

BASALT-Air is a lightweight, modular reimplementation of BASALT using [Pixi](https://pixi.sh) for dependency management.

### Key Improvements

- **Absolute path support**: No need to `cd` into data directories. Run BASALT from anywhere.
- **Separated work/output directories**: `--workdir` for intermediate files, `--outdir` for final results.
- **Pixi-based installation**: Faster, more reproducible than Conda. No environment activation headaches.
- **Hugging Face model hosting**: Download weights via `huggingface-cli` or `pixi run download-weights`.
- **Built-in dependency check**: `BASALT --check-deps` and `BASALT --version`.
- **Configurable CUDA version**: Set in `pixi.toml` instead of manual CUDA toolkit management.
- **Semicolon separator**: Multiple datasets separated by `;` instead of `/` for cleaner CLI.

### Breaking Changes from BASALT (Conda)

| | BASALT (Conda) | BASALT-Air |
|---|---|---|
| Package manager | Conda | Pixi |
| Activation | `conda activate basalt_env` | `pixi shell` |
| Dataset separator | `/` (slash) | `;` (semicolon) |

> BASALT-Air shares the same core pipeline and produces identical output. The change is in packaging and CLI ergonomics only.

### Links

- GitHub: [https://github.com/PKU-EMBL/BASALT-Air](https://github.com/PKU-EMBL/BASALT-Air)
- Hugging Face: [https://huggingface.co/PKU-EMBL/BASALT_WEIGHT](https://huggingface.co/PKU-EMBL/BASALT_WEIGHT)

---

## v1.2.0 (2025-12-16)

### Major Changes

- **Python 3.12 support**: Upgraded from Python 3.8 to 3.12 with a streamlined installation script.
- **LorBin integration**: Added LorBin (*Nat. Commun.*, 2025) as an extra binner (`-e l`), the state-of-the-art long-read binning tool.
- **CheckM2 as default**: Changed the default quality assessment from CheckM to CheckM2 v1.1.0.
- **GPU acceleration**: Added GPU support for Semibin2 and deep learning model inference.
- **Flexible model weights**: Model weights path is now configurable via `BASALT_WEIGHT` environment variable instead of being fixed to `~/.cache/BASALT`.
- **VAMB 5.0 support**: Updated VAMB integration for compatibility with VAMB 5.0 (Python 3.12).

### Breaking Changes

- Python 3.8 is no longer supported. Python 3.12 is required.
- The default QC backend is now CheckM2 (previously CheckM). Use `-q checkm` to revert to CheckM.
- `--autopara` argument has been renamed to `--sensitive`.

---

## v1.1.0 (2024-06-12)

### Changes

- Improved installation workflow with `BASALT_setup.py`.
- Added China mainland mirror support (`BASALT_setup_China_mainland.py`).
- Added Singularity image support.
- Various bug fixes and stability improvements.

---

## v1.0.0 (2023-08-18)

Initial public release.

### Features

- Multi-assembly binning with four binners (MetaBAT2, Maxbin2, CONCOCT, Semibin2).
- Deep-learning-based contamination detection and removal (S5 outlier remover).
- Contig retrieval via paired-end tracking (S6) and long-read connections (S7).
- Reassembly with SPAdes (S9) and Unicycler hybrid (S9p).
- Data Feeding module for importing external bins.
- CheckM-based quality assessment.
- Checkpoint-based resumption.

---

## Publication

> Qiu, Z., Yuan, L., Lian, CA. et al. **BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis.** *Nature Communications* **15**, 2179 (2024).
