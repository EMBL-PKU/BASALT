## BASALT - Binning Across a Series of AssembLies Toolkit

### Type of input data for BASALT
BASALT is a versatile tool with high efficiency for binning and post-binning refinement. BASALT can generate high quality metagenome-assembled genomes (MAGs) from various input data types including: 1) assembly from short-read sequences (SRS); 2) assembly from long-read sequences (LRS); [Note: only PacBio-HiFi data is supported in v1.0.1 for long-read only assemblies, other types of LRS data will be available in v1.0.2.] 3) hybrid assembly from SRS + LRS. Specific features of BASALT are listed below:

1.	Multiple assemblies as input with dereplication function
BASALT developed a comprehensive method incorporating multiple assembly files, including single assemblies (SAs) and co-assemblies (CAs) in one run. Additionally, a dereplication step is applied after initial binning to efficiently remove redundant bins. Comparatively, prominent binning tools such as metaWRAP 1 and DASTool 2 only support single assembly file as input, where multiple binning processes are required if there are multiple assembly files in a dataset. Moreover, redundant bins generated under SA + CA mode need to be removed using dereplication tools such as dRep 3. Although BASALT takes longer time than metaWRAP and DASTool in one single run, considerable amount of time will be saved when processing multiple assembly files, and a significantly more and better quality MAGs can be generated by BASALT than other tools, based on our assessment 4.
2.	Standalone Refinement module
BASALT can effectively identify and remove potential contamination sequences in Refinement module based on neural networks. Specifically, users can import their bins along with raw sequences to BASALT to run the Refinement module independently without initial binning using BASALT.
3.	High read use efficiency of LRS
BASALT maximized the utilization of LRS in the post-binning refinement steps. Firstly, LRS will be used at Sequence retrieval step by recruiting unused sequences to target bins via pair-end tracking. After processing the Sequence retrieval function, an extra polishing step will be performed in the existence of LRS. Polishing at this step will save ~90% of the computation time than conducting at assembly step with same iterations, as the data size is largely reduced. Furthermore, LRS will be exploited again at reassembly step using the SPAdes Hybrid function (default). [Note: reassembly function is not applicable on LRS-alone dataset in v1.0.1, but will be available in later version.] Although reassembly may take a considerable amount of time, large augmentation of genome quality can be observed after reassembly.
For any issue compiling and running BASALT, as well as bug report, please do not hesitate to contact us (yuke.sz@pku.edu.cn). Thanks for using BASALT!



### SYSTEM REQUIREMENTS
1.	Required dependencies

  	Linux x64 systems, 8+ cores, and 128GB+ RAM

  	Python (>=3.0) modules: biopython, pandas, numpy, scikit-learn, copy, multiprocessing, collections，pytorch, tensorboardx

  	Perl

  	Java (>=1.7)

  	Binning tools: MetaBAT2, Maxbin2, CONCOCT, Semibin2, Metabinner

  	Note: VAMB was used to be implemented in BASALT, but due to the conflict of environment and unsatisfactory performance on environmental datasets, we temporarily removed VAMB from BASALT environment. However, bins generated using VAMB can still be imported to BASALT directly for post-binning refinements.

  	Sequences processing tools: Bowtie2, BWA, SAMtools, Prodigal, BLAST+, HMMER, Minimap2

  	Sequences assembly and polishing tools: SPAdes, IDBA-UD, Pilon, Racon, Unicycler

  	Genome quality assessment tools: CheckM, CheckM2, pplacer

  	Note: CheckM2 database is not compiled along with BASALT installation in v1.0.1. To setup CheckM2 database, please refer to CheckM2 user guide (https://github.com/chklovski/CheckM2).



### INSTALLATION
1.	Quick installation
   
  	Download BASALT_setup.py and run:
   ```
   python BASALT_setup.py
   ```
   Please remain patient, as the installation process may take an extended period.

2. Quick installation from China mainland 从中国内地快速安装BASALT
   
   For users in China mainland who may experience a network issue, please download the alternative script ‘BASALT_setup_China_mainland.py’ and run:

   中国内地且无法翻墙的用户推荐使用‘BASALT_setup_China_mainland.py’安装
   ```
   python BASALT_setup_China_mainland.py
   ```
   Then, download the trained models for neural networks BASALT.zip from Tencent iCloud (https://share.weiyun.com/r33c2gqa) and run:
   ```
   mv BASALT.zip ~/.cache
   cd ~/.cache
   unzip BASALT.zip
   ```

3. Manual installation (recommended)
   
   Install Miniconda (https://docs.anaconda.com/free/miniconda/miniconda-install/) or Anaconda (https://docs.anaconda.com/free/anaconda/install/index.html)

   Add mirrors to increase download speed of BASALT dependent software (optional):
   ```
   site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
   conda config --add channels ${site}/pkgs/free/
   conda config --add channels ${site}/pkgs/main/
   conda config --add channels ${site}/cloud/conda-forge/
   conda config --add channels ${site}/cloud/bioconda/
   ```

   Download the BASALT installation file and create a conda environment:
   ```
   git clone https://github.com/EMBL-PKU/BASALT.git
   cd BASALT
   conda env create -n BASALT --file basalt_env.yml
   ```

   Please remain patient, as the installation process may take an extended period.

   If you have encountered an error, please download 'basalt_env.yml' from Tencent iCloud (https://share.weiyun.com/xXdRiDkl) and create a conda environment:
   ```
   conda env create -n BASALT --file basalt_env.yml
   ```

   After successfully creating the conda environment, change file permissions for BASALT script files:
   ```  
   chmod -R 777 <PATH_TO_CONDA>/envs/BASALT/bin/*
   ```
   Example: To easily find your path to conda environments, simply use:
   ```
   conda info --envs
   ```
   and you can find your path to BASALT environment, such as:
   ```
   # conda environments:
   #
   base     /home/emma/miniconda3
   BASALT   /home/emma/miniconda3/envs/BASALT
   ```
   Then, change permission to BASALT script folders:
   ```
   chmod -R 777 /home/emma/miniconda/envs/BASALT/bin/*
   ```

   Download the trained models for neural networks 'BASALT.zip' from FigShare:
   ```
   wget https://figshare.com/ndownloader/files/41093033
   mv 41093033 BASALT.zip
   mv BASALT.zip ~/.cache
   cd ~/.cache
   unzip BASALT.zip
   ```
   For users from China mainland, please download the models BASALT.zip from Tencent iCloud (https://share.weiyun.com/r33c2gqa) and run:
   ```
   mv BASALT.zip ~/.cache
   cd ~/.cache
   unzip BASALT.zip
   ```

4. Another way to install BASALT in China mainland
   以singularity的方式加载BASALT的sif镜像使用BASALT,可通过微云的以下网址获得BASALT.sif镜像文件
   ```
   https://share.weiyun.com/xKmoBmrF
   ```
      
   将BASALT的singularity镜像（BASALT.sif）放置在服务器的home目录下。以执行singularity的命令运行，如
   ```
   singularity run BASALT.sif BASALT -a as1.fa -s S1_R1.fq,S1_R2.fq/S2_R1.fq,S2_R2.fq -t 32 -m 128
   ```

   如BASALT.sif不在home目录下运行需要添加 -B挂载，如
   ```
   singularity run -B /media/emma BASALT.sif BASALT -h
   ```

   需要后台挂载运行，nohup可能会出现意外，但是集群一般sbatch等提交命令的方式可以正常运行。实验室的服务器则考虑使用screen命令。
请严格参考screen命令的执行方式（除非你很熟悉screen，切勿擅自修改命令执行方式）。如
   ```
    screen -dmS session_name bash -c 'bash basalt.sh >log_basalt'
   ```
   请注意session_name要起跟自己有辨识度唯一的名字，避免发生意外情况

   BASALT.sif含有checkm1 checkm2 semibin  bowtie2 bwa等很多软件，均可以通过以下方式调用：
   ```
   singularity run BASALT.sif bowtie2 -h
   ```
   
6. Test files
   Sample demo files (see BASALT demo files) are prepared for testing whether the BASALT script can be successfully performed, and the bins can be generated. The demo files contain Data.tar.gz, Final_bestbinset.tar.gz and basalt.sh.
   ```
   Data.tar.gz -> short read and long read raw sequence files and an OPERA-MS assembled contig file.
   Final_bestbinset.tar.gz -> expected output of final bins.
   basalt.sh -> script running this demo
   ```
   A workstation with a configuration of Intel(R) Xeon(R) Gold 5218 CPU @ 2.30GHz with 32 cores is expected to complete processing of this demo dataset within 6 hours.



### USAGE
1.	General usage

  	To run BASALT, use BASALT under conda environment, or use BASALT.py for standalone users:
   ```
   BASALT [-h] [-a ASSEMBLIES] [-s SR_DATASETS] [-l LONG_DATASETS] [-hf HIFI_DATASET] [-c HI_C_DATASET] [-t THREADS] [-m RAM] [-e EXTRA_BINNER] [-qc QC_SOFTWARE] [--min-cpn MIN_COMPLETENESS] [--max-ctn MAX_CONTAMINATION] [--mode RUNNING_MODE] [--module FUNCTIONAL_MODULE] [--autopara AUTOBINING_PARAMETERS] [--refinepara REFINEMENT_PARAMTER]![image](https://github.com/EMBL-PKU/BASALT/assets/62051720/61fb5b05-2844-4867-9598-f91e0709fa9a)
   ```
   Required arguments
   
   ```
   -a	list of assemblies, e.g., -a assembly1.fa,assembly2.fa
   Files ending with .fa, .fna, and .fasta are all supported. Additionally, compressed files ending with .gz, .tar.gz, and .zip are also supported.
   
   -s	short-read datasets, e.g., -s d1_r1.fq,d1_r2.fq/d2_r1.fq,d2_r2.fq
   Please note, read files within each pair are separated with ‘,’, and read pairs are separated with ‘/’. Reads files ending with .gz, .tar.gz, and .zip are also supported.
   
   -l		long-read datasets, e.g., -l lr1.fq,lr2.fq
   
   -hf	PacBio-HiFi datasets, e.g., -hf hifi1.fq,hifi2.fq
   
   -c	Hi-C datasets, e.g., -c hc1.fq,hc2.fq
   Read files within each pair are separated with ‘,’. Reads files ending with .gz, .tar.gz, and .zip are also supported.
   
   -t	number of threads, e.g., -t 32
   
   -m	RAM, e.g., -m 128
   Suggested minimum RAM is 32G.
   ```
   Optional arguments
   ```
   --min-cpn		Minimum completeness cutoff, e.g., --min-cpn 30 (default: 35)
   
   --max-ctn		Maximum contamination cutoff, e.g., --max-ctn 25 (default: 20)
   
   --mode		Running mode. Start a new project from the beginning –-mode new or continue the previous run –-mode continue. (default: continue)
   
   --module		Running mode. Run only Autobinning + Bin Selection modules –-module autobinning, Refinement module –-module refinement, Gap filling module –-module reassembly, or running all modules –-module all. (default: all)
   
   --autopara		Autobinning parameters. 
   –-autopara more-sensitive Choose recommended binners with full parameters: Maxbin2 [0.3, 0.5, 0.7, 0.9], MetaBAT2 [200, 300, 400, 500], CONCOCT [2-3 flexible parameters based on result of MetaBAT2], and Semibin2 [100]
   –-autopara sensitive Partial binners with partial parameters: MetaBAT2 [200, 300, 400, 500], CONCOCT [1-2 flexible parameters based on result of MetaBAT2], and Semibin2 [100]
   –-autopara quick Limited binners: MetaBAT2 [200, 300, 400, 500] and Semibin2 [100]
(default: more-sensitive)

   --refinepara	Refinement parameters. 
   --refinepara deep will enable deep refinement at sequence retrieval step. Disable this function by setting the parameter with 
   –-refinepara quick. (default: deep)
   
   --hybrid_reassembly	Setting hybrid reassembly parameters. In reassembly function, BASALT uses SPAdes Hybrid function as default parameter
   –-hybrid_assembly n to process hybrid reassembly in the existence of SRS and LRS. Use –-hybrid_assembly y to use Unicycler for hybrid reassembly. Please note that it will take a considerable amount of time when using Unicycler for hybrid reassembly.

   -qc			Selection of quality check software. Use CheckM by setting this parameter at –qc checkm , or use CheckM2 by setting with –qc checkm2.
   
   -e			Enable extra binners. We temporarily disabled VAMB in BASALT v1.0.1. To enable Metabinner, use –e m in addition to other binners
   
   -h 			Help documents.
   ```

2.	Example
   Run BASALT based on SRS datasets:
   ```
   BASALT\
   -a as1.fa,as2.fa,as3.fa\
   -s srs1_r1.fq,srs1_r2.fq/srs2_r1.fq,srs2_r2.fq\
   -t 60 -m 250
   ```

   Run BASALT based on SRS + LRS datasets:
   ```
   BASALT\
   -a as1.fa,as2.fa,as3.fa\
   -s srs1_r1.fq,srs1_r2.fq/srs2_r1.fq,srs2_r2.fq\
   -l lrs1.fq,lrs2.fq -t 60 -m 250
   ```

   Run BASALT based on customized parameters:
   ```
   BASALT\
   -a as1.fa,as2.fa,as3.fa\
   -s srs1_r1.fq,srs1_r2.fq/srs2_r1.fq,srs2_r2.fq\
   -l lr1.fq,lr2.fq -hf hifi1.fq\
   -t 60 -m 250\
   --autopara sensitive --refinepara quick --min-cpn 40 --max-ctn 15 -qc checkm2
   ```



## Troubleshooting
1.	Error from SAMtools when installing BASALT:
   ```
   samtools: error while loading shared libraries: libcrypto.so1.0.0: cannot open shared object file: No such file or directory
   ```
   Troubleshooting: Check the file libcrypto.so1.0.0:
   ```
   ls <PATH_TO_CONDA>/envs/BASALT/lib/libcry*
   ```
   If libcrypto.so.1.1 was found instead of libcrypto.so.1.0.0, create a soft link of libcrypto.so.1.1 to libcrypto.so.1.0.0:
   ```
   cd <PATH_TO_CONDA>/envs/BASALT/lib
   ln -s libcrypto.so.1.1 libcrypto.so.1.0.0
   ```
   Then check if SAMtools is available:
   ```
   samtools -help
   ```
2.	Error encountered when running BASALT:
   ```
   Traceback (most recent call last):
   File "/users/.conda/envs/BASALT/bin/BASALT", line 57, in
   datasets[str(n)].append(pr[1].strip())
   IndexError: list index out of range
   ```
   This is because BASALT does not support reading files with absolute path. To address this issue, simply move/copy corresponding files to the working directory or establish soft links under the working directory.
   
3.	Error encountered when running BASALT:
   ```
   Traceback (most recent call last):
   File "/users/.conda/envs/BASALT/bin/BASALT", line 53, in 
   datasets_list=sr_datasets.split('/')
   ```
   This is because BASALT does not support LRS only mode except PacBio-HiFi reads in v1.0.1. Please use SRS + LRS instead of LRS only. LRS only mode for Nanopore/PacBio long reads will be available in v1.0.2.

4.	Error encountered when running BASALT:
   ```
   INFO: Running CheckM2 version 1.0.1
   [03/13/2024 12:56:34 PM] INFO: Running quality prediction workflow with 30 threads.
   [03/13/2024 12:56:34 PM] ERROR: DIAMOND database not found. Please download database using <checkm2 database --download>
   ```
   This is because CheckM2 database is not installed. Users can simply download CheckM2 database from their official website (https://github.com/chklovski/CheckM2) to address this issue.

5.	Error encountered when running BASALT:
   ```
   BASALT: command not found!
   ```
   This issue occurred because BASALT scripts cannot be found or administrated. Please firstly check if BASALT has been successfully installed, by checking script files under the bin folder of BASALT environment, e.g.:
   ```
   <PATH_TO_CONDA>/envs/BASALT/bin/
   ```
   If no script file was found, please download BASALT script again and copy to the bin folder of BASALT environment, by using the following commands:
   ```
   unzip BASALT_script.zip
   chmod -R 777 BASALT_script
   mv BASALT_script/* <PATH_TO_CONDA>/envs/BASALT/bin/
   ```
   If BASALT scripts are found in the bin folder of BASALT environment, try to change permissions by using the following commands:
   ```
   chmod -R 777 <PATH_TO_CONDA>/envs/BASALT/bin/*
   ```
6.	Error encountered when running BASALT:
   ```
   Traceback (most recent call last):
  File "/user/miniconda3/envs/BASALT/bin/BASALT", line 137, in <module>
    BASALT_main_d(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, autobining_parameters, refinement_paramter, max_ctn, min_cpn, pwd, QC_software)
  File "/user/miniconda3/envs/BASALT/bin/BASALT_main_d.py", line 494, in BASALT_main_d
    Contig_recruiter_main(best_binset_from_multi_assemblies, outlier_remover_folder, num_threads, continue_mode, min_cpn, max_ctn, assembly_mo_list, connections_list, lr_connection_list, coverage_matrix_list, refinement_paramter, pwd)
  File "/user/miniconda3/envs/BASALT/bin/S6_retrieve_contigs_from_PE_contigs_10302023.py", line 1819, in Contig_recruiter_main
    parse_bin_in_bestbinset(assemblies_list, binset+'_filtrated', outlier_remover_folder, PE_connections_list, lr_connection_list, num_threads, last_step, coverage_matrix_list, refinement_mode)
  File "/user/miniconda3/envs/BASALT/bin/S6_retrieve_contigs_from_PE_contigs_10302023.py", line 1695, in parse_bin_in_bestbinset
    bin_comparison(str(binset), bins_checkm, str(binset)+'_retrieved', refinement_mode, num_threads)
  File "/user/miniconda3/envs/BASALT/bin/S6_retrieve_contigs_from_PE_contigs_10302023.py", line 731, in bin_comparison
    for line in open('quality_report.tsv','r'):
FileNotFoundError: [Errno 2] No such file or directory: 'quality_report.tsv'
   ```
	
 This is possibly due to the insufficient number of bins generated due to the low coverage of datasets at the current step, which CheckM2 cannot generate quality file.


   
### FAQ
1.	Q: Same contigs may be used multiple times in the single assembly + co-assembly mode. Does this strategy affect the final output?
   A: Redundant bins can be generated under single assembly + co-assembly mode when raw reads are used multiple times in the assembly step. For example, bin1 is clustered from single assembly A1, bin2 is clustered from co-assembly A1+A2+A3, where bin1 and bin2 are the same genome clustered with same contigs. This redundancy can be identified in the bin selection module, and redundant bins will be removed at this step. The final best binset is a non-redundant binset.

2.	Q: BASALT takes longer time to finish compared to metaWRAP. Is there a way to reduce the computation time?

  	A: Indeed, BASALT will take approximately doubled amount of time compared with metaWRAP, which computation time may even prolong with the increase of sample complexity. However, BASALT can obtain more MAGs with better quality than metaWRAP at similar computation time (i.e., running BASALT without gap filling module). Moreover, BASALT may save computation time on multiple assemblies because it only takes a single run to finish. If users wish to accelerate the entire procedure, we suggest using quick mode (--autopara quick) and skip deep refinement (--refinepara quick). This will significantly reduce the computation time.

4.	Q: How do I process refinement with existing bins?

  	A: To use refinement module only, users can feed their bin sequences with raw reads into BASALT using the script data_feeding.py.
Note: A new version of BASALT (v1.0.2) will be released around late May 2024. We will aim to optimize the above functions to make it more user-friendly.

6.	Q: Is there an output parameter?

  	A: No. BASALT will process and generate results under current working directory. We suggest users prepare a copy or generate soft links of raw reads and assembly files under the working directory to avoid error occurring.


   
### Cite our article
Z Qiu, L Yuan, C Lian, B Lin, J Chen, R Mu, X Qiao, L Zhang, Z Xu, L Fan, Y Zhang, S Wang, J Li, H Cao, B Li, B Chen, C Song, Y Liu, L Shi, Y Tian, J Ni, T Zhang, J Zhou, W Zhuang, K Yu. BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. Nat. Commun. 2024, 15, 2179. https://doi.org/10.1038/s41467-024-46539-7



### References
1. Uritskiy, G.V., DiRuggiero, J. & Taylor, J. MetaWRAP—a flexible pipeline for genome-resolved metagenomic data analysis. Microbiome 6, 1-13 (2018).
2.	Sieber, C.M. et al. Recovery of genomes from metagenomes via a dereplication, aggregation and scoring strategy. Nature microbiology 3, 836-843 (2018).
3.	Olm, M.R., Brown, C.T., Brooks, B. & Banfield, J.F. dRep: a tool for fast and accurate genomic comparisons that enables improved genome recovery from metagenomes through de-replication. The ISME journal 11, 2864-2868 (2017).
4.	Qiu, Z. et al. BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis. Nature Communications 15, 2179 (2024).

