# Study design patterns

BASALT can compare candidate bins recovered from several assemblies. This makes the assembly plan part of the scientific method rather than a file-management detail. The patterns below adapt the useful study-design material from the former Chinese manual to the current BASALT command line.

:::{admonition} Design principle
:class: important

Combine assemblies only when their read coverage is biologically interpretable together. More assemblies can expose complementary genome fragments, but unrelated samples can also create misleading mappings and comparisons. Test a representative pilot before scaling the design.
:::

## Decide what belongs in one run

Use one BASALT run when all of the following are true:

- the assemblies address the same biological comparison or coherent sample group;
- the supplied reads are valid coverage evidence for those assemblies;
- sample identity, time point, treatment, and assembly provenance are known;
- a cross-assembly representative set is scientifically meaningful.

Split the analysis when combining samples would obscure a prespecified contrast, when read provenance is uncertain, or when one group is likely to dominate coverage solely because it was sequenced more deeply.

| Question | If yes | If no |
|---|---|---|
| Do the assemblies represent one coherent biological stratum? | They may be compared in one run | Separate runs by stratum |
| Can every coverage file be traced to the reads and assembly that produced it? | Record the mapping explicitly | Repair provenance before running |
| Is a pooled assembly biologically justified? | Include it alongside individual assemblies in a pilot | Prefer sample-specific assemblies |
| Would pooling erase a time, treatment, host, or habitat contrast? | Keep strata separate | Pooling may be evaluated |

## Pattern 1: individual and pooled assemblies

For replicate samples from one habitat or treatment, a pilot can compare sample-specific assemblies with one pooled assembly made from the same biological group. BASALT then evaluates candidates from both assembly strategies against the supplied read evidence.

```text
replicate reads ──> sample-specific assemblies ──┐
       │                                        ├─> one BASALT run
       └───────> biologically matched pool ─────┘
```

Example for three paired-end samples and optional matched long reads:

```bash
BASALT \
  -a sample_1.fa,sample_2.fa,sample_3.fa,matched_pool.fa \
  -s sample_1_R1.fastq,sample_1_R2.fastq/sample_2_R1.fastq,sample_2_R2.fastq/sample_3_R1.fastq,sample_3_R2.fastq \
  -l sample_1_ont.fastq,sample_2_ont.fastq,sample_3_ont.fastq \
  -t 64 -m 256 \
  --mode new \
  -o habitat_A_basalt
```

Omit `-l` when long reads are not available. The pooled assembly must be produced from the stated group; it should not silently include reads from another habitat, treatment, or host.

### Why test both assembly strategies?

Sample-specific assemblies preserve individual provenance and can retain sample-restricted populations. A pooled assembly can gain coverage for genomes shared across replicates, but it can also merge closely related populations or increase assembly complexity. BASALT provides a framework for comparing the resulting candidates; it does not prove that pooling improved the assembly.

Report the assembler, version, parameters, and exact reads used for every individual and pooled assembly. Compare the pilot using prespecified MAG, redundancy, taxonomy, and read-support criteria.

## Pattern 2: longitudinal or staged sampling

For time courses or developmental stages, preserve the temporal contrast by processing each stage as a separate stratum. Within a stage, sample-specific assemblies and a matched pooled assembly may be evaluated together.

```text
stage 1 replicates ─> individual + stage-1 pool ─> BASALT run 1
stage 2 replicates ─> individual + stage-2 pool ─> BASALT run 2
stage 3 replicates ─> individual + stage-3 pool ─> BASALT run 3
```

Stage 1 example:

```bash
BASALT \
  -a stage1_sample1.fa,stage1_sample2.fa,stage1_sample3.fa,stage1_pool.fa \
  -s stage1_sample1_R1.fastq,stage1_sample1_R2.fastq/stage1_sample2_R1.fastq,stage1_sample2_R2.fastq/stage1_sample3_R1.fastq,stage1_sample3_R2.fastq \
  -t 64 -m 256 \
  --mode new \
  -o stage1_basalt
```

Repeat in new working directories for later stages. Perform downstream cross-stage taxonomy, abundance analysis, or dereplication with a declared similarity threshold. Do not infer temporal persistence merely because BASALT bins share names or broadly similar marker profiles.

## Pattern 3: several habitats or treatments

For habitats such as host-associated, water, and sediment samples, first decide whether the analysis asks for a shared catalogue or habitat-specific populations.

- For habitat-specific recovery, run each habitat separately and compare final MAGs downstream.
- For a shared catalogue, retain habitat labels, assess mapping specificity, and declare the downstream dereplication rule.
- Do not pool reads only to increase depth if the pooled coverage no longer represents a meaningful population mixture.

The former Chinese manual illustrated this pattern with aquaculture replicates. The transferable principle is biological stratification, not the specific sample names or a claim that one assembly strategy is always superior.

## When not to pool

Avoid pooled assembly or pooled coverage by default when:

- samples come from different treatments, hosts, sites, or time points that define the main contrast;
- a high-abundance sample would dominate a low-abundance group;
- strain variation is central to the research question;
- extraction, library, or sequencing protocols are incompatible;
- contamination or sample swaps have not been resolved;
- read-to-assembly provenance cannot be reconstructed.

## Minimum design manifest

Create a machine-readable table before launching BASALT:

| Column | Example | Purpose |
|---|---|---|
| `sample_id` | `stage1_sample1` | Stable biological identifier |
| `stratum` | `stage1` | Time, treatment, habitat, or host group |
| `assembly` | `stage1_sample1.fa` | Assembly supplied to BASALT |
| `assembly_type` | `individual` or `pool` | Distinguish assembly strategy |
| `assembly_reads` | accession or checksum list | Identify reads used for assembly |
| `coverage_reads` | R1/R2, ONT, CLR, or HiFi files | Identify evidence supplied to BASALT |
| `assembler` | name, version, parameters | Make assembly reproducible |
| `run_directory` | absolute archival location | Prevent checkpoint mixing |

Validate the manifest against [Input formats](inputs.md), then archive it with the materials listed in [Reproducibility and reporting](reproducibility.md).

## Evaluate the pilot

Compare alternative study designs using the same acceptance criteria and quality-control database. At minimum, inspect:

- number of retained genomes under explicit completeness and contamination thresholds;
- redundancy at a declared average nucleotide identity criterion;
- taxonomic consistency and strain resolution;
- assembly size, contiguity, and read support;
- genomes unique to pooled or individual assemblies;
- compute time, peak memory, and temporary storage.

Choose a design from those results and the biological question. Do not select a workflow solely because it produces the largest bin count.
