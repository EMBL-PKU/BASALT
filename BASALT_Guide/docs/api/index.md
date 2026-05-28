# API Reference

BASALT is primarily a command-line tool. This section provides reference documentation for the core Python modules.

## Core Modules

| Module | Description |
|---|---|
| `BASALT.py` | CLI argument parsing and workflow dispatch |
| `BASALT_main_d.py` | Main pipeline orchestration (CheckM2 branch) |
| `model.py` | MLP neural network architecture for contamination detection |
| `ensemble.py` | Model ensemble inference for S5 outlier removal |
| `utils.py` | Utility functions for training and evaluation |
| `my_dataset.py` | PyTorch Dataset for loading contig features |

## Pipeline Steps

| Module | Pipeline Stage |
|---|---|
| `S1_Autobinners_2qc_11152023.py` | S1: Multi-binner execution and initial QC |
| `S1e_extra_binners.py` | S1e: Extra binner integration (MetaBinner, VAMB, LorBin) |
| `S2_BinsAbundance_PE_connections_*.py` | S2: Abundance profiling and PE connectivity |
| `S3_Bins_comparator_within_group_*.py` | S3: Within-assembly bin comparison |
| `S4_Multiple_Assembly_Comparitor_*.py` | S4: Cross-assembly dereplication |
| `S5_Outlier_remover_DL_*.py` | S5: DL-based contamination removal |
| `S6_retrieve_contigs_from_PE_contigs_*.py` | S6: PE-based contig retrieval |
| `S7_Contigs_retrieve_within_group_*.py` | S7: Within-group contig retrieval |
| `S7lr_finding_sr_contigs_basing_lr_and_polishing_*.py` | S7lr: Long-read contig retrieval + polishing |
| `S8_OLC_new_*.py` | S8: Overlap-Layout-Consensus |
| `S9_Reassembly_*.py` | S9: Short-read reassembly (SPAdes) |
| `S9p_Hybrid_Reassembly_*.py` | S9p: Hybrid reassembly (Unicycler) |
| `S10_OLC_new_*.py` | S10: Final OLC refinement |

## Helper Modules

| Module | Description |
|---|---|
| `Data_feeding.py` | External binset import and coverage recalculation |
| `gen_kmer.py` | k-mer frequency computation for MetaBinner |
| `Cleanup.py` | Intermediate file cleanup |
| `Final_drep.py` | Final dereplication of binsets |
| `Cytoscapeviz.pl` | Cytoscape connection visualization (Perl) |
| `calc.kmerfreq.pl` | k-mer frequency calculation (Perl) |
| `jgi_summarize_bam_contig_depths` | BAM depth summary (JGI tool) |
