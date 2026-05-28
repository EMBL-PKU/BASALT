# Pipeline Architecture

BASALT is organised into three main functional modules operating in sequence:

```mermaid
flowchart LR
    A[Raw Reads + Assemblies] --> B[Autobinning]
    B --> C[Bin Selection]
    C --> D[Refinement<br/>DL + Contig Retrieval]
    D --> E[Reassembly]
    E --> F[Final MAGs]
```

## Module 1: Autobinning + Bin Selection

Runs multiple binning algorithms in parallel across all assemblies, evaluates quality, and selects the optimal non-redundant binset.

### Steps

| Step | Script | Description |
|---|---|---|
| S1 | `S1_Autobinners_2qc_11152023.py` | Run MetaBAT2, Maxbin2, CONCOCT, Semibin2 with multiple parameter sets |
| S1e | `S1e_extra_binners.py` | Run optional extra binners (MetaBinner, VAMB, LorBin) |
| S2 | `S2_BinsAbundance_PE_connections_*.py` | Compute abundance profiles and paired-end connectivity |
| S3 | `S3_Bins_comparator_within_group_*.py` | Compare bins within each assembly group, select best |
| S4 | `S4_Multiple_Assembly_Comparitor_*.py` | Cross-assembly dereplication to remove redundant bins |

### Binners and Parameters

| Binner | Parameter Values | Sensitivity |
|---|---|---|
| MetaBAT2 | 200, 300, 400, 500 | All presets |
| Maxbin2 | 0.3, 0.5, 0.7, 0.9 | `more-sensitive` only |
| CONCOCT | 1–3 flexible thresholds | `sensitive`, `more-sensitive` |
| Semibin2 | 100 (default) | All presets |
| MetaBinner | 100 | Extra (`-e m`) |
| VAMB | 100 | Extra (`-e v`) |
| LorBin | 100 | Extra (`-e l`) |

## Module 2: Refinement

Removes contamination from individual bins using a deep learning ensemble, then retrieves missed contigs.

### Steps

| Step | Script | Description |
|---|---|---|
| S5 | `S5_Outlier_remover_DL_*.py` | DL ensemble classifies each contig as Real or Contaminated; contaminated contigs are removed |
| S6 | `S6_retrieve_contigs_from_PE_contigs_*.py` | Retrieve missing contigs using paired-end connectivity patterns |
| S7 | `S7_Contigs_retrieve_within_group_*.py` | Intra-group contig retrieval and OLC refinement |
| S7lr | `S7lr_finding_sr_contigs_basing_lr_and_polishing_*.py` | Long-read-based contig retrieval + Racon polishing |
| S8 | `S8_OLC_new_*.py` | Overlap-Layout-Consensus contig merging |

### Deep Learning Architecture

The S5 outlier remover uses an **ensemble of 5 MLP models** with the following architecture:

```
Input (TNF + Coverage + Coverage Change Ratio)
    → Linear(hidden_size)
    → LBR Block × 2 (Linear → BatchNorm → ReLU → residual skip)
    → Linear(num_classes=2)
    → Softmax → Real / Contaminated
```

Each LBR block expands to 4× hidden size and contracts back with a residual connection. The 5 models are trained independently and their predictions are combined via majority voting.

## Module 3: Reassembly

Reassembles refined bins to improve genome quality and contiguity.

### Steps

| Step | Script | Description |
|---|---|---|
| S9 | `S9_Reassembly_*.py` | Short-read reassembly using SPAdes |
| S9p | `S9p_Hybrid_Reassembly_*.py` | Hybrid reassembly using Unicycler (requires long reads) |
| S10 | `S10_OLC_new_*.py` | Final Overlap-Layout-Consensus refinement |

## Checkpoint System

BASALT records progress after each step in `Basalt_checkpoint.txt`. The checkpoint file tracks:

```
1st autobinner done!
2nd bin selection within group done!
3rd bin selection within multiple groups done!
4th outlier removal done!
5th contig retrieval done!
6th long-read contig retrieval done!
7th intra-group contig retrieval done!
8th reassembly done!
9th final OLC done!
```

Use `--mode continue` to resume from the last completed step. Use `--mode new` to start fresh.

## Data Feeding Workflow

The Data Feeding module (`Data_feeding.py`) allows users to bypass the Autobinning module entirely:

1. Rename and reindex external bins to be compatible with BASALT
2. Recompute coverage matrices using Bowtie2 mapping
3. Generate paired-end connection files
4. Run CheckM/CheckM2 quality assessment
5. Feed the prepared data directly into the Refinement module
