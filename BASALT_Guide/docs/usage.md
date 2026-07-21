# Command-line reference

This page documents the Conda-based BASALT CLI in this repository. BASALT-Air has a different path model and sample delimiter.

## Command structure

```text
BASALT \
  -a <assembly[,assembly...]> \
  [-s <R1,R2[/R1,R2...]>] \
  [-l <long_read[,long_read...]>] \
  [-hf <hifi_read[,hifi_read...]>] \
  [options]
```

A normal new run requires one or more assemblies and at least one read source, although the argument parser does not enforce every valid combination.

## Input and resource options

| Option | Parsed default | Meaning |
|---|---:|---|
| `-a`, `--assemblies` | none | Comma-separated assembly FASTA files |
| `-s`, `--shortreads` | none | Paired-end files separated by `,`; samples separated by `/` |
| `-l`, `--longreads` | none | Comma-separated ONT or PacBio CLR read files |
| `-hf`, `--hifi` | none | Comma-separated PacBio HiFi read files |
| `-t`, `--threads` | `4` | Requested CPU threads |
| `-m`, `--ram` | `32` | Available RAM in GB used by internal scheduling decisions |

`-m` is not a hard memory cap. The operating system or scheduler must enforce a limit if one is required.

Read [Inputs and study design](inputs.md) for filename, compression, and sample-matching constraints.

## Workflow options

| Option | Default | Accepted values | Meaning |
|---|---:|---|---|
| `-o`, `--out` | `Final_binset` | directory name | Name of the final output directory |
| `-q`, `--quality-check` | `checkm2` | `checkm2`, `checkm` | Quality-control backend |
| `--min-cpn` | `35` | integer percentage | Minimum estimated completeness for bins entering refinement |
| `--max-ctn` | `20` | integer percentage | Maximum estimated contamination for bins entering refinement |
| `--mode` | `continue` | `new`, `continue` | Initialize or resume checkpoint state |
| `--module` | `all` | `autobinning`, `refinement`, `reassembly`, `all` | Requested pipeline section |
| `--sensitive` | `sensitive` | `quick`, `sensitive`, `more-sensitive` | Candidate-generation preset |
| `--refinepara` | `quick` | `quick`, `deep` | Contig-retrieval depth |
| `-e`, `--extra_binner` | none | `m`, `v`, `l` | Optional MetaBinner, VAMB, or LorBin adapter; combine with commas |

!!! note "Output naming"
    In the current CheckM2 workflow, the final directory is the value supplied to `-o`. For example, `-o study_01_basalt` produces `study_01_basalt/`. The option is a name, not a portable absolute output path.

## Advanced input routes

| Option | Default | Meaning |
|---|---:|---|
| `-d`, `--data-feeding-folder` | none | Comma-separated external binset directories for data feeding |
| `--binset-index` | `500` | Starting index assigned to imported binsets |
| `-r`, `--refinement-binset` | empty | Existing binset for standalone outlier screening |
| `-c`, `--coverage-list` | none | Comma-separated compatible coverage matrices |
| `-b`, `--binsets-list` | none | Comma-separated existing binsets for dereplication |

These routes require BASALT-compatible identifiers and intermediate files. They are not interchangeable shortcuts for a normal new run. Preserve the source assembly, reads, coverage-generation method, and original binning provenance.

## Sensitivity presets

| Preset | Candidate generators | Consequence |
|---|---|---|
| `quick` | MetaBAT 2 and SemiBin 2 | Fewer candidate sets and lower runtime |
| `sensitive` | MetaBAT 2, dynamically parameterized CONCOCT, and SemiBin 2 | Default balance of candidate diversity and runtime |
| `more-sensitive` | MaxBin 2.0, MetaBAT 2, dynamically parameterized CONCOCT, and SemiBin 2 | More candidate sets and higher compute and storage demand |

The presets change candidate generation, not a universal accuracy threshold. `more-sensitive` can increase opportunities for recovery but does not guarantee more acceptable MAGs.

## Refinement depth

| Value | Behaviour |
|---|---|
| `quick` | Standard retrieval path |
| `deep` | Extends retrieval in compatible code paths and generally requires more time |

Report this option together with completeness and contamination thresholds. A deeper search changes the candidate space and should not be treated as a purely computational setting.

## Recommended invocation pattern

Specify scientific and computational defaults explicitly so the run remains interpretable if software defaults change:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -t 32 -m 128 \
  --sensitive sensitive \
  --refinepara quick \
  --min-cpn 35 \
  --max-ctn 20 \
  -q checkm2 \
  --mode new \
  -o study_01_basalt
```

## Examples

### Several assemblies and short-read samples

```bash
BASALT \
  -a assembly_1.fasta,assembly_2.fasta \
  -s sample_1_R1.fastq,sample_1_R2.fastq/sample_2_R1.fastq,sample_2_R2.fastq \
  -t 64 -m 256 \
  --mode new -o study_multi_assembly
```

### Short reads plus ONT or PacBio CLR reads

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -l nanopore.fastq \
  -t 64 -m 256 \
  --mode new -o study_hybrid
```

### PacBio HiFi

```bash
BASALT \
  -a assembly.fasta \
  -hf hifi.fastq \
  -t 32 -m 128 \
  --mode new -o study_hifi
```

### Autobinning only

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  --module autobinning \
  -t 32 -m 128 --mode new
```

Module-specific execution relies on checkpoint and intermediate state. Start with a full run unless you understand the stage dependencies described on the [pipeline page](pipeline.md).

### Data feeding

```bash
BASALT \
  -s sample_R1.fastq,sample_R2.fastq \
  -d external_bins_1,external_bins_2 \
  --binset-index 1 \
  -t 32 -m 128 \
  -q checkm2 \
  -o study_external
```

With CheckM2, a non-default `-o study_external` is transformed to `study_external_data_feeded/` for this data-feeding route.

### Standalone outlier screening

```bash
BASALT \
  -s sample_R1.fastq,sample_R2.fastq \
  -l nanopore.fastq \
  -a source_assembly.fasta \
  -r existing_bins \
  -c Coverage_matrix_for_binning_source_assembly.fasta.txt \
  -t 32 -m 128 \
  -q checkm2
```

The supplied assembly, reads, bin contig identifiers, and coverage matrix must be mutually compatible.

## Resume semantics

`--mode continue` reads `Basalt_checkpoint.txt` and reuses intermediates in the current directory. It does not prove that those files were generated by the same software version or command.

Before resuming, confirm:

- the same working directory and inputs;
- the same BASALT commit or release;
- the same model files and database;
- compatible external-tool versions;
- no manual renaming or partial cleanup of intermediates.

If any condition changed, begin a new run in a new directory.

## Exit status and logs

Some external programs are invoked through shell commands, and not every failure is propagated as a non-zero BASALT exit status. A scheduler state of `COMPLETED` is therefore insufficient evidence of a successful scientific result.

Check all of the following:

1. the last entries in `Basalt_log.txt` and `Basalt_checkpoint.txt`;
2. the presence and size of final FASTA files;
3. the final quality report;
4. stderr for failed external commands;
5. warnings about optional binners or skipped stages.
