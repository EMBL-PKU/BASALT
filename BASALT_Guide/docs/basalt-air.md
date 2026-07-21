# BASALT-Air: the new implementation

:::{admonition} BASALT-Air v1.0.0 is available
:class: important

[PKU-EMBL/BASALT-Air](https://github.com/PKU-EMBL/BASALT-Air) is a newer, separately maintained implementation of the BASALT workflow. Its official repository currently identifies the software as **v1.0.0**.
:::

## Why BASALT-Air?

BASALT-Air retains the multi-assembly binning and refinement model while modernizing deployment, paths, and run provenance. Its official repository describes:

- a Pixi environment and lock file for Python and external bioinformatics tools;
- first-class absolute input paths;
- separate `--workdir` and `--outdir` locations;
- short-read, ONT or PacBio CLR, PacBio HiFi, and hybrid workflows;
- modular autobinning, refinement, reassembly, and data-feeding execution;
- dependency checks, timestamped logs, and a structured run manifest.

The installed entry points are lowercase `basalt` and `basalt_models_download`.

## Which implementation should I use?

| Situation | Recommended implementation |
|---|---|
| New deployment that benefits from absolute paths and isolated work/output locations | BASALT-Air |
| New reproducible environment managed from a lock file | BASALT-Air |
| Existing analysis created with this repository's checkpoints | Continue with the same pinned BASALT version |
| Reproduction of a published or archived Conda BASALT run | Use the recorded BASALT commit, dependencies, models, and databases |
| Cross-implementation comparison | Run both independently and compare with prespecified metrics |

Choosing BASALT-Air is an operational recommendation, not a guarantee of more or better MAGs. The scientific evidence boundary remains defined by the BASALT publication and the datasets and parameters tested.

## Important compatibility differences

| Interface | Conda BASALT documented here | BASALT-Air v1.0.0 |
|---|---|---|
| Executable | `BASALT` | `basalt` |
| Environment | Conda plus installation script | Pixi environment and lock file |
| Absolute paths | Not handled consistently in legacy command paths | Supported directly |
| Work/output control | Current working directory plus `-o` name | `--workdir`, `--outdir`, and `-o` |
| Multiple paired-end datasets | `/` separates pairs | BASALT-Air documentation uses `;` for multiple absolute-path pairs |
| Dependency audit | Manual tool checks | `basalt --check-deps` |
| Version check | Record Git tag and commit | `basalt --version` |
| Run metadata | Logs, command file, and checkpoint | Adds `BASALT_run_manifest.json` and timestamped logs |

Do not paste a BASALT-Air command into this repository's `BASALT` executable. Do not resume a checkpoint produced by one implementation with the other.

## Start with the authoritative repository

Use the BASALT-Air repository for current installation, dependency, and CLI details:

- [BASALT-Air repository and README](https://github.com/PKU-EMBL/BASALT-Air)
- [BASALT-Air issues](https://github.com/PKU-EMBL/BASALT-Air/issues)

A minimal orientation sequence is:

```bash
git clone https://github.com/PKU-EMBL/BASALT-Air.git
cd BASALT-Air
pixi install
pixi run version
pixi run check-deps
```

Configure model-weight and CheckM2 database paths exactly as described in the BASALT-Air repository before starting an analysis.

## Reporting BASALT-Air

State explicitly that BASALT-Air was used and report:

- BASALT-Air version and Git commit;
- `pixi.lock` checksum or archived lock file;
- complete `basalt` command;
- model-weight and database identifiers;
- input and output checksums;
- `BASALT_run_manifest.json`, logs, and any resumed or skipped stages.

The BASALT-Air repository asks users to cite the original BASALT publication:

> Qiu, Z., Yuan, L., Lian, C.-A. *et al.* BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. *Nature Communications* **15**, 2179 (2024). [https://doi.org/10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7)
