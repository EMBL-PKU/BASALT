# Inputs and study design

BASALT can combine evidence across assemblies and read datasets, but it cannot repair a mismatched study design. Confirm sample identity, read provenance, assembly provenance, and intended coverage relationships before launching the workflow.

## Required logical inputs

A normal run requires:

1. one or more metagenomic assemblies in FASTA format; and
2. at least one coverage source supplied as paired-end short reads, ONT or PacBio CLR reads, or PacBio HiFi reads.

The CLI does not validate all biologically invalid combinations before starting. Treat the checks below as part of the analysis protocol.

## Assemblies

Supply assemblies with `-a` as a comma-separated list:

```text
-a assembly_1.fasta,assembly_2.fasta
```

Accepted filename suffixes are `.fa`, `.fna`, and `.fasta`. BASALT filters assembly sequences shorter than 1,000 bp during preprocessing.

Assemblies may represent single-sample assemblies, co-assemblies, or a combination of both. Each assembly should be derived from reads that are represented in the supplied coverage data. Adding an unrelated assembly can create uninformative mappings and invalid comparisons.

## Paired-end short reads

Within each sample, separate R1 and R2 with a comma. Separate samples with `/`:

```text
-s sample_1_R1.fastq,sample_1_R2.fastq/sample_2_R1.fastq,sample_2_R2.fastq
```

The order of samples is retained internally. Use unambiguous filenames and verify that every pair has consistent read identifiers and orientation.

!!! warning "Shell interpretation"
    Quote the complete `-s` value if a wrapper, scheduler, or shell treats `/` or other characters specially. BASALT itself uses `/` as the sample delimiter for this Conda edition.

## ONT and PacBio CLR reads

Supply non-HiFi long-read files with `-l` as a comma-separated list:

```text
-l nanopore_1.fastq,nanopore_2.fastq
```

BASALT maps these reads to assemblies and can use their connectivity for retrieval and polishing. Long-read support does not imply that every bin will be polished. A bin may lack sufficient mapped reads or may fail an external reassembly or polishing step.

## PacBio HiFi reads

Supply HiFi reads separately with `-hf`:

```text
-hf hifi_1.fastq,hifi_2.fastq
```

Do not place HiFi files under `-l`. BASALT uses distinct mapping and processing logic for `-hf`.

## Compressed files

The pipeline contains decompression paths for `.gz`, `.tar.gz`, and `.zip` inputs. Plain FASTA and FASTQ files are the most reproducible choice because archive layout and automatic extraction can vary. If compressed files are used, test the exact format on a small run and archive the original checksums.

## Paths and working directories

The Conda edition writes intermediates into the current working directory and constructs several shell commands from filenames. Use a dedicated run directory with simple filenames or symlinks:

```bash
mkdir -p /project/basalt_runs/study_01
cd /project/basalt_runs/study_01
ln -s /project/data/assembly.fasta .
ln -s /project/data/sample_R1.fastq .
ln -s /project/data/sample_R2.fastq .
```

Avoid whitespace, quotes, wildcard characters, leading hyphens, and shell metacharacters in input names and output prefixes. Absolute paths may work in some code paths but are not consistently handled across the pipeline. Local names or relative symlinks are the supported practice for this guide.

## External binsets

Data feeding (`-d`), standalone refinement (`-r`), and dereplication (`-b`) are advanced routes.

- A binset directory should contain one FASTA file per bin.
- Contig identifiers must match the relevant assembly and coverage files.
- Coverage matrices must originate from compatible reads and assemblies.
- Bin identifiers should be unique across imported directories.
- Preserve the original binning method, version, parameters, and quality report.

Do not use a coverage matrix from an unrelated assembly merely because the dimensions are compatible.

## Preflight checklist

Run these checks before a production job:

```bash
test -s assembly.fasta
test -s sample_R1.fastq
test -s sample_R2.fastq

awk '/^>/{n++} END{print n " assembly sequences"}' assembly.fasta
awk 'END{print NR/4 " R1 reads"}' sample_R1.fastq
awk 'END{print NR/4 " R2 reads"}' sample_R2.fastq

sha256sum assembly.fasta sample_R1.fastq sample_R2.fastq > input.sha256
```

For large files, use a validated FASTA/FASTQ checker available at your site. The line-count check assumes unwrapped four-line FASTQ records and is only a quick diagnostic.

## Study-design record

Before running, record:

| Field | Minimum description |
|---|---|
| Biological samples | Sample identifiers, source, collection design, and exclusions |
| Sequencing | Platform, library type, read processing, and retained read counts |
| Assembly | Assembler, version, parameters, and which reads contributed |
| BASALT coverage sources | Which read files were mapped to which assemblies |
| Optional binners | Software, version, parameters, and reason for inclusion |
| Quality control | CheckM or CheckM2 version, database, and thresholds |

This record is required to interpret whether a recovered MAG represents an expected population and whether comparisons across assemblies are valid.
