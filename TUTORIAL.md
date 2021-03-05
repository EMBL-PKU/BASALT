### TUTORIAL
#### A guide to analyzing metagenomic data with BASALT

##### Simple command:
```
BASALT -al assembly1.fa,assembly2.fa,assembly3.fa -ds dataset1_read1.fq, dataset1_read2.fq\;dataset2_read1.fq,dataset2_read2.fq
```

##### Required parameters:
**Input data**:

`-al` &nbsp;&nbsp; Multiple assembly files for binning. Assembly files can be short reads assembled fasta, long reads assembled fasta or hybrid assembled fasta. Fasta files are separated by comma.  
`-ds` &nbsp;&nbsp; Sequence files for binning. Only short read sequences are valid. Default sequence files are pair-end library. Sequence files within pair-end library are separated by comma, and different pair-end library are separated by “\;”.  

**Optional parameters**:  
`-t` &nbsp;&nbsp; Number of threads. For example, -t 120. Default thread is set at 2.  
`-r` &nbsp;&nbsp; Maximum memory limit at gigabytes (GB). For example, -r 750. Default memory limit is set at 64GB. The maximum memory limit is critical for reassembly.  
`--long-reads` &nbsp;&nbsp; Long read sequences. Sequences are separated by “,”.  
`--continue` &nbsp;&nbsp; Continue run from the last available check-point. BASALT will check the checkpoints and start from the last checkpoint. This default of this parameter is on.  
`--new` &nbsp;&nbsp; Restart binning step.  
`--autobinning` &nbsp;&nbsp; Autobinning options. Available options are “quick” (default), “sensitive” and “more-sensitive”. The computational time will increase if later options were selected, but more high-quality bins are expected to be recovered.  
`--reassembly` &nbsp;&nbsp; Reassembly options. Available options are “quick-refinement” (default) and “reassembly”. If reassembly was needed, please use option --reassembly. BASALT will perform short reads assembly if only short read sequences were provided. Alternatively, BASALT will perform hybrid assembly if both short read and long read sequences were provided.  
`--max-ctn` &nbsp;&nbsp; Contamination cutoff in the refinement step. Default cutoff is set at 20, which means BASALT will only refine those bins with contamination at 20 or below.  
`--min-cpn` &nbsp;&nbsp; Completeness cutoff in the refinement step. Default cutoff is set at 35, which means BASALT will only refine those bins with completeness at 35 or above.  

**Detailed case**:
```
BASALT -al assembly1.fa,assembly2.fa,assembly3.fa -ds dataset1_read1.fq, dataset1_read2.fq\;dataset2_read1.fq,dataset2_read2.fq --long-reads ont1.fq,ont2.fq -t 120 -r 750 --autobinning more-sensitive --reassembly --max-ctn 20 --mix-cpn 35