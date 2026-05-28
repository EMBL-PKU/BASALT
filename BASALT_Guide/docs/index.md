---
hide:
  - navigation
  - toc
---

# Welcome to BASALT

**Binning Across a Series of Assemblies Toolkit**

<div class="grid cards" markdown>

- :material-dna: **Multi-Assembly Binning**

    ---

    Process single assemblies and co-assemblies in one run with automatic dereplication. No need for separate tools like dRep.

- :material-brain: **Deep Learning Refinement**

    ---

    Neural network ensemble (5 MLP models) identifies and removes contamination from individual bins using TNF, coverage, and k-mer features.

- :material-swap-horizontal: **Hybrid Read Support**

    ---

    Maximise long-read utilisation for contig retrieval and polishing. ~90% faster than assembly-level polishing.

- :material-rocket-launch: **GPU Acceleration**

    ---

    GPU support for Semibin2 binning and deep learning inference in BASALT v1.2.0+.

</div>

---

## What is BASALT?

BASALT is a comprehensive toolkit for **binning and post-binning refinement** of metagenomic assemblies. It produces high-quality **metagenome-assembled genomes (MAGs)** from short-read, long-read (ONT/PacBio), and HiFi metagenomic data.

The pipeline integrates multiple binning algorithms, a deep-learning-based contamination removal module, and a reassembly engine to maximise MAG recovery while minimising contamination.

<figure markdown="span">
  ![BASALT Workflow](img/workflow.png)
  <figcaption>Overview of the BASALT pipeline, from raw reads to final MAGs.</figcaption>
</figure>

---

## BASALT vs BASALT-Air

BASALT comes in two editions:

| | **BASALT** | **BASALT-Air** |
|---|---|---|
| **Package manager** | Conda / Singularity | Pixi |
| **Path support** | Working directory only | Absolute paths anywhere |
| **Dataset separator** | `/` (slash) | `;` (semicolon) |
| **Output control** | Writes to CWD | `--workdir` + `--outdir` |
| **CUDA config** | Manual | Built-in via `pixi.toml` |
| **Model download** | `BASALT_models_download.py` | Hugging Face / pixi tasks |
| **Status** | Mature (v1.2.0) | New (v1.0.0) |

> **Recommendation:** New users should start with **BASALT-Air** for easier installation and flexible path handling. Existing BASALT users can continue using the Conda edition — both versions share the same core pipeline.

---

## Key Features

| Feature | Description |
|---|---|
| **Multiple binners** | MetaBAT2, Maxbin2, CONCOCT, Semibin2 + optional MetaBinner, VAMB, LorBin |
| **Multi-assembly** | Process any number of single assemblies and co-assemblies in one run |
| **DL contamination removal** | 5-model MLP ensemble classifies real vs. contaminated contigs |
| **Contig retrieval** | PE-tracking and long-read connectivity recover missed contigs |
| **Reassembly** | SPAdes / Unicycler hybrid reassembly for improved genome quality |
| **Data feeding** | Import externally generated bins for refinement and reassembly |
| **Checkpoint resumption** | Interrupted runs resume cleanly from the last completed step |
| **Dual QC backends** | CheckM and CheckM2 quality assessment |

---

## Citation

If you use BASALT in your research, please cite:

> Qiu, Z., Yuan, L., Lian, CA. et al. **BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis.** *Nature Communications* **15**, 2179 (2024). [DOI: 10.1038/s41467-024-46539-7](https://doi.org/10.1038/s41467-024-46539-7)

```bibtex
@article{qiu2024basalt,
  title   = {BASALT refines binning from metagenomic data and increases
             resolution of genome-resolved metagenomic analysis},
  author  = {Qiu, Zhiguang and Yuan, Li and Lian, Chun-Ang and Lin, Bin
             and Chen, Jie and Mu, Rong and Qiao, Xuejiao and Zhang, Liyu
             and Xu, Zheng and Fan, Lu and others},
  journal = {Nature Communications},
  volume  = {15},
  number  = {1},
  pages   = {2179},
  year    = {2024},
  doi     = {10.1038/s41467-024-46539-7}
}
```

---

## Contact

- **Ke Yu** — [yuke.sz@pku.edu.cn](mailto:yuke.sz@pku.edu.cn)
- **Zhaorui (Elijah) Jiang** — [zrjiang25@stu.pku.edu.cn](mailto:zrjiang25@stu.pku.edu.cn)

---

## References

1. Qiu, Z. et al. BASALT refines binning from metagenomic data... *Nat. Commun.* **15**, 2179 (2024).
2. Uritskiy, G.V. et al. MetaWRAP—a flexible pipeline... *Microbiome* **6**, 1-13 (2018).
3. Sieber, C.M. et al. Recovery of genomes from metagenomes via DAS Tool... *Nat. Microbiol.* **3**, 836-843 (2018).
4. Olm, M.R. et al. dRep: a tool for fast and accurate genomic comparisons... *ISME J.* **11**, 2864-2868 (2017).
5. Xue, W. et al. LorBin: efficient binning of long-read metagenomes... *Nat. Commun.* **16**, 9353 (2025).
