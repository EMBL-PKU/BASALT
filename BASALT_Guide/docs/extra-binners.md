# Extra binners

BASALT can add candidate bins from MetaBinner, VAMB, or LorBin to the default candidate pool. These adapters are optional and have independent dependencies, interfaces, and failure modes.

| Flag | Adapter | Intended evidence source |
|---|---|---|
| `-e m` | MetaBinner | Sequence composition and abundance profiles |
| `-e v` | VAMB | Variational autoencoder representation and abundance |
| `-e l` | LorBin | Long-read-oriented binning |

Combine flags with commas, for example `-e m,v`. An adapter runs in addition to the binners selected by `--sensitive`.

## Scientific use

Adding a binner expands the candidate set presented to BASALT selection. It does not guarantee more non-redundant, high-quality, or biologically correct MAGs. Choose an adapter because its assumptions fit the data, not solely to maximize candidate count.

For every optional binner, record:

- software version and installation source;
- complete command and parameters;
- input assembly and coverage files;
- whether CPU or GPU execution was used;
- adapter warnings or failures;
- candidate count before BASALT selection;
- the software citation.

## Failure semantics

The current CheckM2 orchestration catches failures from optional adapters, writes a warning, and continues with other candidates. A BASALT run can therefore finish even when an optional binner failed.

Search the logs explicitly:

```bash
grep -E '\[WARNING\]|\[ERROR DETAILS\]|Extra binner' \
  Basalt_log.txt basalt.stderr.log
```

Confirm that the expected adapter output directory exists and contains non-empty FASTA files. Do not report an optional binner as part of the method if it failed before contributing candidates.

## MetaBinner

MetaBinner requires its own executable workflow and a compatible `run_metabinner.sh` installation. BASALT derives a coverage profile and k-mer table before invoking the adapter.

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -e m \
  -t 32 -m 128 --mode new -o study_metabinner
```

Verify the MetaBinner installation outside BASALT before using this route. The adapter currently discovers installation paths through the active Conda environments, which may be site-dependent.

## VAMB

The repository contains an adapter for the VAMB 5 command form. VAMB uses mapped BAM files produced during BASALT preprocessing.

Install and verify a compatible VAMB version in the active environment:

```bash
python -m pip install 'vamb>=5,<6'
vamb --help
```

Run the adapter:

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -e v \
  -t 32 -m 128 --mode new -o study_vamb
```

VAMB GPU behaviour depends on its PyTorch installation and local accelerator stack. Record `vamb`, Python, PyTorch, CUDA, and driver versions when GPU execution is used.

## LorBin

LorBin is intended for long-read metagenomic binning. The adapter depends on a separately installed `LorBin` executable and compatible BAM files.

```bash
BASALT \
  -a assembly.fasta \
  -s sample_R1.fastq,sample_R2.fastq \
  -l nanopore.fastq \
  -e l \
  -t 32 -m 128 --mode new -o study_lorbin
```

Treat the current LorBin integration as an optional adapter that requires local validation. Run LorBin on a small test, inspect its command in stderr or logs, and verify the returned FASTA set before production use.

## Manual recovery through data feeding

If an adapter fails but the underlying tool can be run independently, preserve its bin FASTA files and import them through BASALT data feeding. Contig identifiers must remain compatible with the BASALT assembly.

```bash
BASALT \
  -s sample_R1.fastq,sample_R2.fastq \
  -d external_binner_bins \
  --binset-index 1 \
  -t 32 -m 128 -q checkm2 -o imported
```

The CheckM2 data-feeding route writes `imported_data_feeded/` for this example.

## Citation

Cite the exact optional binner used. For LorBin, see:

> Xue, W. *et al.* LorBin: efficient binning of long-read metagenomes by multiscale adaptive clustering and evaluation. *Nature Communications* **16**, 9353 (2025).

Use each software project's recommended citation and version-specific documentation for MetaBinner and VAMB.
