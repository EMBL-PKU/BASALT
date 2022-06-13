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

`-l` &nbsp;&nbsp; LONG_DATASETS, --longreads LONG_DATASETS&nbsp;&nbsp; List of long reads, e.g.: lr1.fq,lr2.fq

**Optional parameters**: 
`-h` &nbsp;&nbsp; show this help message and exit
`-t` &nbsp;&nbsp; Number of threads. For example, -t 120. Default thread is set at 2.  
`-m` &nbsp;&nbsp; Maximum memory limit at gigabytes (GB). For example, -r 750. Default memory limit is set at 64GB. The maximum memory limit is critical for reassembly.  
`--continue` &nbsp;&nbsp; Continue run from the last available check-point. BASALT will check the checkpoints and start from the last checkpoint. The default of this parameter is on.  
`--new` &nbsp;&nbsp; Restart binning step.  
`--autobinning` &nbsp;&nbsp; Autobinning options. Available options are “quick” (default), “sensitive” and “more-sensitive”. The computational time will increase if later options were selected, but more high-quality bins are expected to be recovered.  
`--reassembly` &nbsp;&nbsp; Reassembly options. Available options are “quick-refinement” (default) and “reassembly”. If reassembly was needed, please use option --reassembly. BASALT will perform short reads assembly if only short read sequences were provided. Alternatively, BASALT will perform hybrid assembly if both short read and long read sequences were provided.  
`--max-ctn` &nbsp;&nbsp; Contamination cutoff in the refinement step. Default cutoff is set at 20, which means BASALT will only refine those bins with contamination at 20 or below.  
`--min-cpn` &nbsp;&nbsp; Completeness cutoff in the refinement step. Default cutoff is set at 35, which means BASALT will only refine those bins with completeness at 35 or above.  

**Detailed case**:
```
BASALT -al assembly1.fa,assembly2.fa,assembly3.fa -ds dataset1_read1.fq,dataset1_read2.fq\;dataset2_read1.fq,dataset2_read2.fq --long-reads ont1.fq,ont2.fq -t 120 -r 750 --autobinning more-sensitive --reassembly --max-ctn 20 --mix-cpn 35
```

**Note**:

If you encounter an error during use the BASALT, try to resubmit the task.
