# Reproducibility and reporting

A reproducible BASALT analysis requires more than the final FASTA files. The result depends on input preprocessing, assembly provenance, read-to-assembly relationships, BASALT code, trained models, external executables, quality-control databases, parameters, checkpoint history, and post-processing decisions.

## Minimum provenance record

### Inputs

- biological sample identifiers and inclusion or exclusion criteria;
- sequencing platform, library preparation, and read preprocessing;
- assembly software, version, parameters, and source reads;
- exact assembly and read filenames with checksums;
- mapping between each assembly and the reads used as coverage evidence.

### Software and models

- BASALT release and Git commit;
- installation route, operating system, and hardware allocation;
- Conda environment export;
- versions of binners, mappers, assemblers, polishers, and quality-control tools;
- source repository and Git commit for every forked optional binner, including the LorBin–BASALT Extra-binner used by `-e l`;
- checksums of the BASALT ensemble descriptors and checkpoints;
- accelerator, driver, CUDA, and PyTorch versions when a GPU is used.

### Databases

- CheckM or CheckM2 version;
- database name, release or acquisition date, path, and checksum when available;
- any taxonomic database used after BASALT.

### Parameters and workflow state

- complete BASALT command;
- `--sensitive`, `--refinepara`, `--min-cpn`, `--max-ctn`, and `-q`;
- optional binners and their parameters;
- whether the run started new or resumed;
- skipped, failed, manually rerun, or replaced stages;
- `BASALT_command.txt`, `Basalt_log.txt`, and `Basalt_checkpoint.txt`.

### Outcome definition

- criteria used to call a sequence set a MAG;
- completeness, contamination, size, contiguity, and taxonomy filters;
- downstream dereplication threshold and software;
- treatment of strain-level redundancy;
- all manual curation and exclusions;
- final quality report and output checksums.

## Capture commands

Run these from the analysis directory and adapt tool names to the installed environment:

```bash
git -C /path/to/BASALT rev-parse HEAD > basalt-git-commit.txt
# Required when the LorBin adapter (`-e l`) contributed candidates:
git -C /path/to/LorBin-BASALT-Extrabinner rev-parse HEAD \
  > lorbin-basalt-git-commit.txt
micromamba env export -n basalt > basalt-environment.yml
# Conda alternative: conda env export -n basalt --no-builds > basalt-environment.yml

{
  python --version
  checkm2 --version
  metabat2 --help 2>&1 | head -n 2
  SemiBin2 --version
  bowtie2 --version | head -n 1
  samtools --version | head -n 1
  minimap2 --version
  spades.py --version
} > software-versions.txt 2>&1

sha256sum input/* > input.sha256
(
  cd "$BASALT_WEIGHT"
  find . -type f \
    \( -name '*.pth' -o -name '*_ensemble.csv' \) -print0 \
    | sort -z \
    | xargs -0 sha256sum
) > basalt-models.sha256
```

Some tools print versions only through help output or stderr. Inspect `software-versions.txt` rather than assuming every command succeeded.

## Claim–evidence–boundary audit

Before reporting a result, map each major claim to direct evidence and a boundary.

| Claim type | Minimum evidence | Boundary to state |
|---|---|---|
| More MAGs recovered | Defined MAG criteria, matched comparator, same input data, counts and quality distributions | Tested datasets and parameter settings |
| Lower contamination | Same quality estimator and database, before–after paired comparison | Estimator uncertainty and model dependence |
| Improved contiguity | N50 or sequence-count changes with genome size and contamination retained | Reassembly and read-support conditions |
| Novel taxa or functions | Versioned taxonomy or annotation workflow with appropriate controls | Database coverage and classification threshold |
| Faster workflow | Wall time, CPU time, memory, storage, hardware, and comparable stages | Hardware and dataset scale |

Avoid phrases such as “high quality,” “improved,” “faster,” or “non-redundant” unless the metric, comparator, threshold, and evaluation scope are explicit.

## Methods template

Adapt the bracketed fields without adding unsupported claims:

> Metagenomic assemblies were processed with BASALT [release; Git commit] on [operating system and hardware]. Assemblies were generated with [assembler, version, parameters] from [samples and read types]. BASALT received [number] assemblies, [number] paired-end datasets, [number] ONT or CLR datasets, and [number] HiFi datasets. Candidate generation used the `[quick|sensitive|more-sensitive]` preset with [optional binners, versions, source repository or fork commit, and parameters]. Contig refinement used the `[quick|deep]` setting and retained bins entering refinement at a minimum estimated completeness of [value]% and maximum estimated contamination of [value]%. Genome quality was estimated with [CheckM or CheckM2, version, database identifier]. The run used [threads] threads and [RAM] GB and was [started new|resumed from stated checkpoint]. Final bins were filtered using [complete criteria] and dereplicated with [software, version, threshold]. Commands, environment files, logs, database identifiers, model checksums, and output checksums are available at [repository or archive].

## Results reporting template

Report observations before interpretation:

> BASALT produced [number] final bins. Under the prespecified criteria of [criteria], [number] were retained for downstream analysis. Median estimated completeness was [value] and median estimated contamination was [value], calculated with [software and database]. Relative to [matched comparator], BASALT yielded [effect with units and uncertainty] under the same input and quality-control protocol. This comparison is limited to [datasets, read types, parameter settings, and exclusions].

Do not claim a mechanism from a quality-score change alone. If a contig-level explanation is proposed, support it with read mappings, composition, taxonomy, and an appropriate comparison.

## Data and code availability

A durable release package should include:

- input accession identifiers or a controlled-access procedure;
- scripts that generate the BASALT command and post-processing results;
- BASALT commit and environment lock or export;
- model and database identifiers;
- logs and checkpoints;
- final MAG FASTA files and quality tables when data governance permits;
- checksums and a manifest describing every file;
- the licence and citation requirements.

Large or restricted read files need not be duplicated when stable accessions and exact preprocessing steps are available.

## Reproducibility boundary

Byte-identical output may not be feasible across hardware accelerators, multithreaded external tools, and database releases. State the intended reproducibility level:

- **computational replay**, using the same images, files, and checksums;
- **workflow reproduction**, using the same versions and parameters with tolerable ordering differences;
- **scientific replication**, recovering the same bounded conclusion under an independently executed protocol.

Do not describe a workflow as reproducible without specifying which level was tested.
