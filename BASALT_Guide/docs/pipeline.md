# Pipeline and methods

BASALT is organized as a sequence of candidate generation, selection, refinement, and reassembly stages. The description below distinguishes executable operations from biological interpretation.

```mermaid
flowchart LR
    A["Assemblies and reads"] --> B["Autobinning"]
    B --> C["Within-assembly selection"]
    C --> D["Cross-assembly dereplication"]
    D --> E["Contig-level refinement"]
    E --> F["Read-supported retrieval"]
    F --> G["Reassembly and OLC"]
    G --> H["Final MAG set and quality report"]
```

## Stage 0: preprocessing and coverage

BASALT filters assembly sequences shorter than 1,000 bp, maps the supplied reads, and derives depth and connectivity files. Paired-end reads are mapped with Bowtie 2. Long reads use long-read-aware mapping paths, including minimap2 where applicable.

**Interpretive boundary.** Coverage is informative only when the supplied reads and assemblies belong to a coherent study design. BASALT does not test whether files were paired correctly at the biological level.

## Stage 1: candidate generation

The selected `--sensitive` preset controls candidate diversity:

| Preset | Candidate generators |
|---|---|
| `quick` | MetaBAT 2 and SemiBin 2 |
| `sensitive` | MetaBAT 2, CONCOCT, and SemiBin 2 |
| `more-sensitive` | MaxBin 2.0, MetaBAT 2, CONCOCT, and SemiBin 2 |

MetaBAT 2 is run across four internal threshold settings. CONCOCT cluster counts are derived from preliminary bin counts. MaxBin 2.0 uses four probability thresholds in `more-sensitive` mode. Optional adapters can add MetaBinner, VAMB, or LorBin candidates.

**Interpretive boundary.** Candidate multiplicity increases the space evaluated by later selection stages. It does not make every candidate independent or biologically valid.

## Stages 2–4: bin characterization and selection

BASALT computes abundance and paired-end connectivity summaries, compares candidates within an assembly, and then compares selected bins across assemblies. The selection logic uses combinations of contig overlap, coverage, sequence composition, and CheckM or CheckM2 quality estimates.

These stages reduce redundant candidates and retain representatives for downstream refinement.

**Interpretive boundary.** Computational dereplication is not a taxonomic or strain-definition procedure. Use an explicit downstream species or strain criterion if the scientific question requires one.

## Stage 5: contig-level contamination screening

The outlier-removal stage derives tetranucleotide-frequency and coverage-related features for contigs. An ensemble of five multilayer perceptrons classifies contigs into retained and predicted-contaminant classes. The ensemble combines model predictions by voting.

The model architecture contains linear, batch-normalization, and rectified-linear layers with residual blocks. Model weights are distributed separately and are located through `BASALT_WEIGHT`.

**Interpretive boundary.** The output is a model prediction conditioned on the training distribution and available features. It is not direct experimental proof that a contig is foreign to a population genome. Report the model version or checksum when refinement affects the main result.

## Stages 6–8: contig retrieval and OLC

BASALT evaluates unbinned or excluded contigs using paired-end connectivity, long-read connectivity, coverage consistency, and sequence-composition criteria. Compatible contigs may be added to a bin. Overlap–layout–consensus (OLC) operations can merge or extend sequences in applicable routes.

`--refinepara deep` extends retrieval in compatible code paths. The effect depends on the available read types and checkpoint state.

**Interpretive boundary.** Retrieved contigs should be treated as read- and feature-supported candidates. Evaluate quality changes, taxonomic consistency, and unexpected sequence composition before interpreting an expanded bin.

## Stages 9–10: reassembly and final OLC

When paired-end short reads are available, BASALT performs bin-level reassembly using SPAdes-related paths. Hybrid routes use long reads where supported. Reassembled and original candidates are compared, followed by another OLC step in compatible workflows.

Reassembly is skipped when the required short-read inputs are absent. Individual bins may also lack sufficient reads for successful reassembly.

**Interpretive boundary.** Higher contiguity is not sufficient evidence of higher genome accuracy. Examine completeness, contamination, assembly size, N50, read support, and taxonomic coherence together.

## Finalization

The final selected directory is moved to the name supplied by `-o`. FASTA files are normalized to `.fa` names where possible, and a final CheckM2 report is generated when one is not already present. Cleanup then archives or removes several intermediate groups.

Retain logs and checksums before applying additional filtering or renaming.

## Checkpoint model

`Basalt_checkpoint.txt` records textual completion markers for major stages. `--mode continue` reads these markers and attempts to resume from the last recognized state.

Checkpointing provides operational recovery, not environment reproducibility. A valid resume requires compatible inputs, code, models, databases, software versions, and intermediate filenames.

## Implementation map

| Stage | Main implementation |
|---|---|
| CLI and dispatch | `BASALT.py` |
| CheckM2 orchestration | `BASALT_main_d.py` |
| Autobinning | `S1_Autobinners_2qc_11152023.py` |
| Optional binners | `S1e_extra_binners.py` |
| Abundance and PE connectivity | `S2_BinsAbundance_PE_connections_multiple_processes_pool_10032023.py` |
| Within-assembly comparison | `S3_Bins_comparator_within_group_10042023.py` |
| Cross-assembly comparison | `S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_10242023.py` |
| Outlier screening | `S5_Outlier_remover_DL_11012023.py` |
| PE-supported retrieval | `S6_retrieve_contigs_from_PE_contigs_10302023.py` |
| Additional retrieval and polishing | `S7_Contigs_retrieve_within_group_10262023.py`, `S7lr_finding_sr_contigs_basing_lr_and_polishing_11022023.py` |
| OLC | `S8_OLC_new_10262023.py` |
| Reassembly | `S9_Reassembly_10262023.py`, `S9p_Hybrid_Reassembly_10262023.py` |
| Final OLC | `S10_OLC_new_10262023.py` |

The dated filenames are implementation identifiers, not software release numbers. Cite the BASALT release or commit used for analysis.

## Methodological reporting

A Methods section should state:

- BASALT release or Git commit;
- input assemblies and read datasets, including preprocessing;
- selected quality-control backend and database;
- sensitivity and refinement settings;
- completeness and contamination thresholds;
- optional binners and their versions;
- compute environment and external-tool versions;
- post-BASALT filtering, dereplication, taxonomy, and quality criteria.

See [Reproducibility and reporting](reproducibility.md) for a reusable template.
