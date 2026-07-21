# Developer reference

BASALT is a command-line workflow rather than a stable importable Python API. Internal functions, filenames, argument lists, and intermediate formats may change without deprecation. External software should invoke a pinned `BASALT` command and treat the documented files as workflow outputs.

## Entry points

| File | Responsibility |
|---|---|
| `BASALT.py` | Parse the CLI and dispatch normal, data-feeding, refinement, and dereplication routes |
| `BASALT_main_d.py` | Orchestrate the default CheckM2 workflow |
| `BASALT_main_c.py` and `BASALT_main_c_*` | Legacy CheckM orchestration |
| `Data_feeding.py` | Normalize external binsets and generate compatible mapping and quality inputs |

## Pipeline modules

| Stage | CheckM2 implementation | Legacy CheckM implementation |
|---|---|---|
| Candidate generation | `S1_Autobinners_2qc_11152023.py` | shared |
| Optional binners | `S1e_extra_binners.py` | shared |
| Coverage and connectivity | `S2_BinsAbundance_PE_connections_multiple_processes_pool_10032023.py` | `S2_BinsAbundance_PE_connections_multiple_processes_pool_checkm.py` |
| Within-assembly selection | `S3_Bins_comparator_within_group_10042023.py` | `S3_Bins_comparator_within_group_checkm.py` |
| Cross-assembly selection | `S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_10242023.py` | `S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm.py` |
| Contig screening | `S5_Outlier_remover_DL_11012023.py` | `S5_Outlier_remover_DL_checkm.py` |
| PE-supported retrieval | `S6_retrieve_contigs_from_PE_contigs_10302023.py` | `S6_retrieve_contigs_from_PE_contigs_checkm.py` |
| Additional retrieval | `S7_Contigs_retrieve_within_group_10262023.py` | `S7_Contigs_retrieve_within_group_checkm.py` |
| Long-read retrieval and polishing | `S7lr_finding_sr_contigs_basing_lr_and_polishing_11022023.py` | `S7lr_finding_sr_contigs_basing_lr_and_polishing_checkm.py` |
| OLC | `S8_OLC_new_10262023.py` | `S8_OLC_new_checkm.py` |
| Reassembly | `S9_Reassembly_10262023.py`, `S9p_Hybrid_Reassembly_10262023.py` | matching `_checkm.py` files |
| Final OLC | `S10_OLC_new_10262023.py` | `S10_OLC_new_checkm.py` |

## Model components

| File | Responsibility |
|---|---|
| `model.py` | Multilayer-perceptron and residual-block definitions |
| `ensemble.py` | Load model descriptors and combine predictions |
| `my_dataset.py` | PyTorch dataset wrappers for contig features |
| `utils.py` | Training and evaluation utilities |
| `BASALT_models_download.py` | Download model descriptors and checkpoints |

The trained model files are external assets. Code version alone is insufficient to reproduce model-based refinement.

## Development constraints

- Many modules execute work at the filesystem and process level.
- Several commands are assembled as shell strings, so filenames require conservative handling.
- Intermediate filenames are part of checkpoint coordination.
- Broad exception handlers and shell exit-code handling can obscure failures.
- Dated module filenames are historical identifiers, not semantic versions.

Changes to intermediate naming, checkpoint text, or contig identifiers can break continuation and downstream stages. Test full and resumed workflows when modifying them.

## Documentation build

Documentation sources are under `BASALT_Guide/docs/` and use MkDocs Material.

```bash
python -m venv .venv-docs
. .venv-docs/bin/activate
python -m pip install -r BASALT_Guide/requirements.txt
mkdocs build --strict --config-file BASALT_Guide/mkdocs.yml
```

Read the Docs uses the repository-root `.readthedocs.yaml`. A local strict build should pass before documentation changes are merged.

## Contribution checklist

1. Keep CLI defaults synchronized between `BASALT.py`, the README, and `usage.md`.
2. Add a bounded test or demo case for changed workflow behavior.
3. Preserve or explicitly migrate checkpoint semantics.
4. Document required external-tool and database versions.
5. Run a strict documentation build.
6. Report limitations and failure modes with the feature description.
