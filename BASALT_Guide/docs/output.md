# Outputs and quality control

BASALT writes final and intermediate outputs into the current working directory. The exact intermediate set depends on input read types, selected modules, quality-control backend, optional binners, and whether stages were skipped or resumed.

## Final output

For a normal CheckM2 run, `-o` names the final directory:

```text
<value supplied to -o>/
```

With the default, this is:

```text
Final_binset/
```

The directory contains final bin FASTA files, generally with a `.fa` suffix, and a final quality report when CheckM2 completes successfully.

!!! warning "Historical names"
    Earlier documentation referred to `Final_binset_final_binset/` or `Final_bestbinset/`. Those names do not describe the current `-o` behaviour in the default CheckM2 orchestration.

## Run-level provenance

| File | Role | Retain? |
|---|---|---:|
| `BASALT_command.txt` | Recorded input and parameter summary | yes |
| `Basalt_log.txt` | Runtime events and warnings | yes |
| `Basalt_checkpoint.txt` | Textual stage completion markers | yes |
| `Autobinner_checkpoint.txt` | Candidate-generation progress | recommended |
| captured stdout and stderr | External-tool messages not always duplicated in the BASALT log | yes |

Archive these files with the final MAGs. They are required to distinguish a completed result from a partially successful run.

## Principal intermediate groups

Names can vary, but the following patterns reflect the current workflow.

### Candidate generation

| Pattern | Contents |
|---|---|
| `*_metabat_genomes/` | MetaBAT 2 candidate bins |
| `*_maxbin2_genomes/` | MaxBin 2.0 candidates in `more-sensitive` mode |
| `*_concoct_genomes/` | CONCOCT candidates in non-quick modes |
| `*_semibin_genomes/` | SemiBin 2 candidates |
| `*_checkm2/` or `*_checkm/` | Stage-specific quality estimates |
| `Coverage_matrix_*.txt` | Coverage matrices used by downstream stages |
| `condense_connections_*.txt` | Summarized paired-end connectivity |

### Selection and refinement

| Pattern | Contents |
|---|---|
| `*BestBinsSet*` or `BestBinset*` | Selected or dereplicated binsets at successive stages |
| `*outlier_refined*` | Bins after model-based contig screening |
| `*filtrated_retrieved*` | Bins after thresholding and contig retrieval |
| `Predicted_potential_outlier.txt` | Contig-level model predictions where produced |

### Reassembly and polishing

| Pattern | Contents |
|---|---|
| `*_re-assembly*` | Reassembled bin candidates |
| `*_OLC*` | Bins processed by OLC comparison or extension |
| `Remained_seq.tar.gz` | Unassigned or remaining reads from compatible polishing paths |

Cleanup can archive or remove intermediate groups after a full run. Copy required audit artifacts before manual cleanup.

## Quality reports

CheckM2 commonly writes `quality_report.tsv`. BASALT may copy the final report to `Final_bestbinset_quality_report.tsv` inside the output directory when a final report is generated after selection.

CheckM2 columns and their precise meaning are defined by the installed CheckM2 version. Do not assume a fixed column order in downstream scripts without inspecting the header.

For the legacy CheckM path, stage outputs include CheckM lineage-workflow files such as `bin_stats_ext.tsv`.

## Interpreting quality estimates

Completeness and contamination are model- or marker-based estimates. They are not direct measurements of genome truth. Interpret them together with:

- genome size and sequence count;
- N50 or other contiguity summaries;
- taxonomic assignment and marker consistency;
- read support and coverage distribution;
- unexpected composition or duplicated regions;
- the biological and environmental context.

State the software version, database, and threshold logic whenever quality classes are reported.

## Completion audit

Before accepting the result:

```bash
test -d study_01_basalt
test -s Basalt_log.txt
test -s Basalt_checkpoint.txt

find study_01_basalt -maxdepth 1 -type f -name '*.fa' -size +0c | wc -l
tail -n 20 Basalt_checkpoint.txt
tail -n 50 Basalt_log.txt
```

Then verify that:

1. the expected final stage completed or was intentionally skipped;
2. optional-binner failures were understood;
3. each retained FASTA is non-empty and parseable;
4. a matching final quality report exists;
5. all post-processing filters are scripted and recorded;
6. final files have checksums.

```bash
find study_01_basalt -maxdepth 1 -type f -print0 \
  | sort -z \
  | xargs -0 sha256sum \
  > study_01_basalt.sha256
```

## Do not infer from directory names

Intermediate names encode program stages, not validated biological categories. A folder containing `MAGs`, `refined`, `polished`, `retrieved`, or `Best` is not evidence that every contained genome passes a publication-specific standard. Apply and report explicit acceptance criteria.
