### TUTORIAL
#### A guide to analyzing metagenomic data with BASALT
To activate this environment, use  
```
conda activate BASALT
```
To deactivate an active environment, use  
```
conda deactivate  
```
##### Simple command:
```
BASALT -a assembly.fasta -s short.R1.fq,short.R2.fq -l nanopore.fastq -t 60 -m 230
```

##### Required parameters:
**Input data**:

`-a` &nbsp;&nbsp; ASSEMBLIES, --assemblies ASSEMBLIES &nbsp;&nbsp; List of assemblies, e.g.: as1.fa,as2.fa

`-s` &nbsp;&nbsp; --shortreads SR_DATASETS&nbsp;&nbsp;List of paired-end reads, e.g.:r1_1.fq,r1_2.fq/r2_1.fq,r2_2.fq (paried_ends reads need '/' to seperate)

**Optional parameters**:

`-h` &nbsp;&nbsp; show this help message and exit

`-l` &nbsp;&nbsp; LONG_DATASETS, --longreads LONG_DATASETS&nbsp;&nbsp; List of long reads, e.g.: lr1.fq,lr2.fq

`-e` &nbsp;&nbsp; EXTRA_BINNER, --extra_binner EXTRA_BINNER. Extra binner for binning: m: metabinner, v: vamb; for instance: -e m, means BASALT will use metabinner for binning besides metabat2, maxbin2, and concoct

`-c` &nbsp;&nbsp; HI_C_DATASET, --HIC HI_C_DATASET. List of Hi-C dataset(s), e.g.: hc1.fq,hc2.fq

`-t` &nbsp;&nbsp; Number of threads. For example, -t 120. Default thread is set at 2. 

`-m` &nbsp;&nbsp; Maximum memory limit at gigabytes (GB). For example, -r 750. Default memory limit is set at 64GB. The maximum memory limit is critical for reassembly.  

`--mode` &nbsp;&nbsp; RUNNING_MODE&nbsp;&nbsp;   Start a new project (new) or continue to run (continue). e.g. --mode continue / --mode new

`--module`  &nbsp;&nbsp; FUNCTIONAL_MODULE&nbsp;&nbsp;  Three modules: 1. autobinning; 2. refinement; 3. reassembly. Default will run all modules. But you could set the only perform modle. e.g. --module reassembly  

`--autopara` &nbsp;&nbsp; AUTOBINING_PARAMETERS&nbsp;&nbsp; Three parameters to chose: 1. more-sensitive; 2.sensitive; 3. quick. Default: more-sensitive. e.g. --autopara sensitive  

`--refinepara` &nbsp;&nbsp; REFINEMENT_PARAMTER&nbsp;&nbsp; Two refinement parameters to chose: 1. deep; 2. quick. Default: deep. e.g. --refinepara quick

`--max-ctn` &nbsp;&nbsp; Contamination cutoff in the refinement step. Default cutoff is set at 20, which means BASALT will only refine those bins with contamination at 20 or below. 

`--min-cpn` &nbsp;&nbsp; Completeness cutoff in the refinement step. Default cutoff is set at 35, which means BASALT will only refine those bins with completeness at 35 or above.  

**Detailed case**:
```
BASALT -a assembly.fasta -s short.R1.fq,short.R2.fq -l nanopore.fastq -t 60 -m 230
```

**Note**:

If you encounter an error during use the BASALT, try to resubmit the task.
