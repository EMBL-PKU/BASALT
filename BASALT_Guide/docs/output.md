# Output Files

All BASALT outputs are generated under the current working directory.

## Final Output

The final curated MAGs are located in:

```
<output_prefix>_final_binset/
```

Default: `Final_binset_final_binset/`

## Output Directory Structure

| Path | Description |
|---|---|
| `Final_binset_final_binset/` | Final curated MAGs (one FASTA per bin) |
| `Basalt_checkpoint.txt` | Checkpoint file used for `--mode continue` |
| `Basalt_log.txt` | Detailed runtime log |
| `BASALT_command.txt` | Record of the BASALT command and parameters used |

## Intermediate Outputs

BASALT generates intermediate outputs at each stage. These can be useful for debugging or if you want to use intermediate results:

### Autobinning Stage

| Path | Description |
|---|---|
| `1_<assembly>_<param>_<binner>_genomes/` | Raw bins from each binner/parameter combination |
| `*_checkm2/` or `*_checkm/` | Quality assessment for each binset |
| `Bins_folder.txt` | Mapping of assembly names to binset folders |
| `Depth_total.txt` | Coverage depth data for all assemblies |
| `Connections_total_dict.txt` | Paired-end connection data for all assemblies |
| `Assembly_MoDict.txt` | Modified assembly file mappings |

### BestBins Selection

| Path | Description |
|---|---|
| `BestBinsSet/` | Non-redundant bins selected after within-assembly and cross-assembly dereplication |
| `BestBinsSet_comparison_files/` | Comparison data used for dereplication |
| `Coverage_matrix_list.txt` | List of coverage matrices for each assembly |
| `Bestbinset_list.txt` | List of best binsets selected |

### Refinement Stage

| Path | Description |
|---|---|
| `BestBinsSet_outlier_refined/` | Bins after DL-based contamination removal |
| `Predicted_potential_outlier.txt` | Per-contig contamination predictions |
| `BestBinsSet_outlier_refined_filtrated/` | Bins after filtration by completeness/contamination |
| `BestBinsSet_outlier_refined_filtrated_retrieved/` | Bins after PE-based and LR-based contig retrieval |
| `BestBinsSet_outlier_refined_filtrated_retrieved_MAGs/` | Final MAGs after intra-group contig retrieval |

### Reassembly Stage

| Path | Description |
|---|---|
| `*_reassembly/` | Bins reassembled with SPAdes (short-read only) |
| `*_hybrid_reassembly/` | Bins reassembled with Unicycler (hybrid) |
| `*_reassembly_OLC/` | Final bins after Overlap-Layout-Consensus refinement |

## Quality Reports

Each binset folder contains a quality report:

| File | Backend | Contents |
|---|---|---|
| `quality_report.tsv` | CheckM2 | Bin name, completeness, contamination, N50, genome size |
| `bin_stats_ext.tsv` | CheckM | Bin name, marker lineage, completeness, contamination, genome size, mean scaffold length |

## Log Files

| File | Description |
|---|---|
| `Basalt_log.txt` | Runtime event log with timestamps |
| `Basalt_checkpoint.txt` | Last completed step (1st–9th) |

## Cleanup

BASALT performs automatic cleanup of intermediate files (temporary BAM files, SAM files, etc.) between steps. Intermediate binsets from earlier stages are preserved on disk.
