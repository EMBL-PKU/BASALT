<section class="basalt-hero" markdown>

<div class="basalt-hero__eyebrow">Genome-resolved metagenomics · PKU-EMBL Laboratory</div>

# BASALT documentation

<p class="basalt-hero__tagline">A reproducible workflow for recovering and refining metagenome-assembled genomes across one or more assemblies.</p>

<figure class="basalt-hero__figure" markdown="span">
  ![BASALT workflow from assemblies and reads to refined MAGs](img/workflow.png)
  <figcaption>Conceptual BASALT workflow. The pipeline page defines the executable stages and their interpretation boundaries.</figcaption>
</figure>

<p class="basalt-hero__summary">BASALT connects multi-binner candidate generation, within- and between-assembly bin selection, contig-level contamination screening, contig retrieval, and reassembly in one checkpointed analysis.</p>

<div class="basalt-hero__actions">
  <a class="basalt-button basalt-button--primary" href="quickstart/">Get started</a>
  <a class="basalt-button basalt-button--secondary" href="pipeline/">Explore the pipeline</a>
</div>

<div class="basalt-hero__badges">
  <a href="https://doi.org/10.1038/s41467-024-46539-7">Nature Communications 2024</a>
  <span>Multi-assembly</span>
  <span>Checkpointed</span>
  <span>MIT licensed</span>
</div>

</section>

!!! abstract "Evidence boundary"

    The published study reports improved MAG recovery and downstream genome-resolved analyses relative to the workflows and datasets tested. This is not a guarantee for every community or sequencing design. Dataset complexity, sequencing depth, assembly quality, read type, database choice, and parameter settings can change the outcome.

## Start here

Choose a route based on the decision you need to make. Commands and defaults on this site are checked against the current repository.

<div class="grid cards basalt-start-grid" markdown>

- :material-download-circle-outline:{ .lg .middle } **Install and verify**

    Configure BASALT, model weights, quality-control databases, and external executables.

    [Installation :octicons-arrow-right-24:](installation.md)

- :material-rocket-launch-outline:{ .lg .middle } **Run a minimal analysis**

    Launch an isolated short-read run with explicit resources, mode, and output name.

    [Quick start :octicons-arrow-right-24:](quickstart.md)

- :material-dna:{ .lg .middle } **Design the inputs**

    Match assemblies, samples, coverage sources, and optional long-read evidence.

    [Inputs and study design :octicons-arrow-right-24:](inputs.md)

- :material-transit-connection-variant:{ .lg .middle } **Understand the method**

    Follow each executable stage and distinguish software behaviour from interpretation.

    [Pipeline and methods :octicons-arrow-right-24:](pipeline.md)

- :material-chart-box-outline:{ .lg .middle } **Audit the outputs**

    Locate final MAGs, inspect quality reports, and diagnose incomplete stages.

    [Outputs and quality control :octicons-arrow-right-24:](output.md)

- :material-notebook-check-outline:{ .lg .middle } **Report the analysis**

    Archive provenance and use the minimum Methods and reporting checklist.

    [Reproducibility and reporting :octicons-arrow-right-24:](reproducibility.md)

</div>

## Method at a glance

| Stage | Function | Main evidence used | Principal boundary |
|---|---|---|---|
| Autobinning | Generate candidate bins under several algorithms or parameter settings | Sequence composition and read depth | Candidate diversity depends on the selected preset and optional binners |
| Bin selection | Compare bins within and across assemblies | Bin quality estimates, overlap, composition, and coverage | Dereplication does not establish biological strain identity |
| Refinement | Screen candidate contaminant contigs and retrieve supported contigs | Composition, coverage-derived features, and read connectivity | Predictions depend on the trained models and available coverage signals |
| Reassembly | Reassemble selected bins when compatible reads are present | Reads assigned to individual bins | Reassembly can be skipped or fail for bins with insufficient read support |

The [pipeline page](pipeline.md) separates what the software executes from how results should be interpreted.

## Supported analysis designs

Every run requires one or more assembly FASTA files and at least one read source that can provide coverage information.

| Design | Short reads | ONT or PacBio CLR | PacBio HiFi |
|---|---:|---:|---:|
| Short-read only | required | no | no |
| Long-read only | no | required | no |
| HiFi only | no | no | required |
| Short reads + ONT or CLR | required | required | no |
| Short reads + HiFi | required | no | required |

Read the [input page](inputs.md) before combining assemblies or samples. BASALT does not infer which reads belong to which biological sample.

## BASALT and BASALT-Air

This site documents the Conda-based BASALT repository. [BASALT-Air](https://github.com/PKU-EMBL/BASALT-Air) is a separate Pixi-based implementation with absolute-path handling and separate work and output directories. Its command grammar is not interchangeable with the commands on this site.

Use the documentation that belongs to the implementation and version you ran. Do not mix BASALT-Air separators or path options into a BASALT command.

## Evidence and citation

If BASALT contributes to an analysis, cite:

> Qiu, Z., Yuan, L., Lian, C.-A. *et al.* BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. *Nature Communications* **15**, 2179 (2024). [https://doi.org/10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7)

```bibtex
@article{qiu2024basalt,
  author  = {Qiu, Zhiguang and Yuan, Li and Lian, Chun-Ang and Lin, Bin and
             Chen, Jie and Mu, Rong and Qiao, Xuejiao and Zhang, Liyu and
             Xu, Zheng and Fan, Lu and Zhang, Yunzeng and Wang, Shanquan and
             Li, Junyi and Cao, Huiluo and Li, Bing and Chen, Baowei and
             Song, Chi and Liu, Yongxin and Shi, Lili and Tian, Yonghong and
             Ni, Jinren and Zhang, Tong and Zhou, Jizhong and Zhuang, Wei-Qin and
             Yu, Ke},
  title   = {BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis},
  journal = {Nature Communications},
  volume  = {15},
  pages   = {2179},
  year    = {2024},
  doi     = {10.1038/s41467-024-46539-7}
}
```

Also cite optional binners and quality-control software when they materially affect the reported result.

## Scope of this documentation

The guide distinguishes three kinds of statement:

- **Executable behaviour**, checked against the current repository.
- **Scientific interpretation**, limited to what the workflow and published evidence support.
- **Operational guidance**, identified as a recommendation when it is not a program requirement.

If documentation and observed behaviour differ, retain the full command, `Basalt_log.txt`, `Basalt_checkpoint.txt`, software versions, and a minimal reproducer. Report the discrepancy through [GitHub Issues](https://github.com/PKU-EMBL/BASALT/issues).
