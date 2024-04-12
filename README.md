## BASALT - Binning Across a Series of AssembLies Toolkit

## Type of input data for BASALT
BASALT is a versatile tool with high efficiency for binning and post-binning refinement. BASALT can generate high quality metagenome-assembled genomes (MAGs) from various input data types including: 1) assembly from short-read sequences (SRS); 2) assembly from long-read sequences (LRS); [Note: only PacBio-HiFi data is supported in v1.0.1 for long-read only assemblies, other types of LRS data will be available in v1.0.2.] 3) hybrid assembly from SRS + LRS. Specific features of BASALT are listed below:

1.	Multiple assemblies as input with dereplication function
BASALT developed a comprehensive method incorporating multiple assembly files, including single assemblies (SAs) and co-assemblies (CAs) in one run. Additionally, a dereplication step is applied after initial binning to efficiently remove redundant bins. Comparatively, prominent binning tools such as metaWRAP 1 and DASTool 2 only support single assembly file as input, where multiple binning processes are required if there are multiple assembly files in a dataset. Moreover, redundant bins generated under SA + CA mode need to be removed using dereplication tools such as dRep 3. Although BASALT takes longer time than metaWRAP and DASTool in one single run, considerable amount of time will be saved when processing multiple assembly files, and a significantly more and better quality MAGs can be generated by BASALT than other tools, based on our assessment 4.
2.	Standalone Refinement module
BASALT can effectively identify and remove potential contamination sequences in Refinement module based on neural networks. Specifically, users can import their bins along with raw sequences to BASALT to run the Refinement module independently without initial binning using BASALT.
3.	High read use efficiency of LRS
BASALT maximized the utilization of LRS in the post-binning refinement steps. Firstly, LRS will be used at Sequence retrieval step by recruiting unused sequences to target bins via pair-end tracking. After processing the Sequence retrieval function, an extra polishing step will be performed in the existence of LRS. Polishing at this step will save ~90% of the computation time than conducting at assembly step with same iterations, as the data size is largely reduced. Furthermore, LRS will be exploited again at reassembly step using the SPAdes Hybrid function (default). [Note: reassembly function is not applicable on LRS-alone dataset in v1.0.1, but will be available in later version.] Although reassembly may take a considerable amount of time, large augmentation of genome quality can be observed after reassembly.
For any issue compiling and running BASALT, as well as bug report, please do not hesitate to contact us (yuke.sz@pku.edu.cn). Thanks for using BASALT!


### SYSTEM REQUIREMENTS

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
Download ‘BASALT_setup.py’ and run:
```
python BASALT_setup.py
```
Please remain patient, as the installation process may take an extended period.
	For users in China mainland who may experience a network issue, please download the alternative script ‘BASALT_setup_China_mainland.py’ and run:
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




#### Option 1: Installation in one step with the conda command:
```
git clone https://github.com/EMBL-PKU/BASALT.git
cd BASALT 
conda env create -n BASALT --file basalt_env.yml
```
Please remain patient, as the installation process may take an extended period.

#### Option 2: Installation with BASALT_setup.py (one step):

1 Download this BASALT_setup.py.

2 You can easily install the software using the command provided below.
```
python BASALT_setup.py
```
Please remain patient, as the installation process may take an extended period.

#### Option 3: Manual installation (this is best, if you are comfortable):

1 Download this basalt_env.yml.

2 (Optional) You can also use mirrors to increase the download speed of BASALT dependent software. For example,if you are in China, you can do this:


```
site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
conda config --add channels ${site}/pkgs/free/ 
conda config --add channels ${site}/pkgs/main/
conda config --add channels ${site}/cloud/conda-forge/
conda config --add channels ${site}/cloud/bioconda/ 
```
 



3 Make a new conda environment to install and manage all dependancies:
```
conda env create -n BASALT --file basalt_env.yml
```
The step 3 will take a while.

4 Download this BASALT_script.zip. Then upload this file to your Linux directory.
```
unzip BASALT_script.zip
chmod -R 777 BASALT_script
```
Move BASALT_script files to your conda BASALT environment. Please note that the conda path may be different on each computer. In general, your conda BASALT environment is located in the subdirectory envs of the conda installation (e.g. /home/anaconda2/envs/).
```
mv BASALT_script/* your_conda/envs/BASALT/bin
```
##### 5 Download the Trained Models for Neural Networks (You can execute this step beforehand, since the neural network model is a foundational requirement for the BASALT operation) 
Download prerequisite script ‘BASALT_models_download.py’ and run the python script:
```
python BASALT_models_download.py
```
If you encounter an error when downloading the models, one possible reason could be the issue of network connections. Please use a VPN to retry or use the manual downloads instead.
For manual downloads, you can download directly from Figshare:
```
https://figshare.com/ndownloader/files/41093033
```
Then rename 41093033 to BASALT.zip.
or Tencent Weiyun from a web browser:
```
https://share.weiyun.com/r33c2gqa
```
When finish downloading, copy downloaded files to .cache folder:
```
mv BASALT.zip ~/.cache
cd ~/.cache
unzip BASALT.zip
```

### TEST FILES
We have also prepared sample demo files (see BASALT demo files) for testing whether the BASALT script can be successfully performed and the bins can be generated. The demo files contain Data.tar.gz,Final_bestbinset.tar.gz and basalt.sh. 
```
Data.tar.gz -> short read and long read raw sequence files and an OPERA-MS assembled contig file.
Final_bestbinset.tar.gz -> expected output of final bins.
basalt.sh -> script running this demo
```
A workstation with a configuration of Intel(R) Xeon(R) Gold 5218 CPU @ 2.30GHz with 32 cores is expected to complete processing of this demo dataset within 6 hours.

### USAGE
1 If you install with Conda or BASALT_setup.py, use the following command to run BASALT:
```
BASALT [-h] [-a ASSEMBLIES] [-s SR_DATASETS] [-l LONG_DATASETS] [-hf HIFI_DATASET] [-c HI_C_DATASET] [-t THREADS] [-m RAM]
       [-e EXTRA_BINNER] [-qc QC_SOFTWARE] [--min-cpn MIN_COMPLETENESS] [--max-ctn MAX_CONTAMINATION] [--mode RUNNING_MODE]
       [--module FUNCTIONAL_MODULE] [--autopara AUTOBINING_PARAMETERS] [--refinepara REFINEMENT_PARAMTER]
```
2 If you use stanalone version of BASALT:
```
python BASALT.py [-h] [-a ASSEMBLIES] [-s SR_DATASETS] [-l LONG_DATASETS] [-hf HIFI_DATASET] [-c HI_C_DATASET] [-t THREADS]
                 [-m RAM] [-qc QC_SOFTWARE] [-e EXTRA_BINNER] [--min-cpn MIN_COMPLETENESS] [--max-ctn MAX_CONTAMINATION]
                 [--mode RUNNING_MODE] [--module FUNCTIONAL_MODULE] [--autopara AUTOBINING_PARAMETERS] [--refinepara REFINEMENT_PARAMTER]
```
BASALT
```
optional arguments:
  -h, --help            show this help message and exit
  -a ASSEMBLIES, --assemblies ASSEMBLIES
                        List of assemblies, e.g.: as1.fa,as2.fa
  -s SR_DATASETS, --shortreads SR_DATASETS
                        List of paired-end reads, e.g.: r1_1.fq,r1_2.fq/r2_1.fq,r2_2.fq (paried_ends reads need '/' to
                        seperate)
  -l LONG_DATASETS, --longreads LONG_DATASETS
                        List of long reads, e.g.: lr1.fq,lr2.fq
  -hf HIFI_DATASETS, --hifi HIFI_DATASETS
                        List of hifi reads, e.g.: hifi1.fq,hifi2.fq
  -c HI_C_DATASET, --HIC HI_C_DATASET
                        List of Hi-C dataset(s), e.g.: hc1.fq,hc2.fq
  -t THREADS, --threads THREADS
                        Number of threads, e.g.: 64
  -m RAM, --ram RAM     Number of ram, minimum ram suggested: 32G
  -qc QC_SOFTWARE, --quality-check QC_SOFTWARE
                        Quality check software, default CHECKM, option CHECKM2. e.g.: -qc checkm2
  -e EXTRA_BINNER, --extra_binner EXTRA_BINNER
                        Extra binner for binning: m: metabinner, v: vamb; for instance: -e m, means BASALT will use
                        metabinner for binning besides metabat2, maxbin2, and concoct
  --min-cpn MIN_COMPLETENESS
                        Min completeness of kept bins (default: 35)
  --max-ctn MAX_CONTAMINATION
                        Max contamination of kept bins (default: 20)
  --mode RUNNING_MODE   Start a new project (new) or continue to run (continue). e.g. --mode continue / --mode new
  --module FUNCTIONAL_MODULE
                        Three modules: 1. autobinning; 2. refinement; 3. reassembly. Default will run all modules. But
                        you could set the only perform modle. e.g. --module reassembly
  --autopara AUTOBINING_PARAMETERS
                        Three parameters to chose: 1. more-sensitive; 2. sensitive; 3. quick. Default: more-sensitive.
                        e.g. --autopara sensitive
  --refinepara REFINEMENT_PARAMTER
                        Two refinement parameters to chose: 1. deep; 2. quick. Default: deep. e.g. --refpara quick
```

Example:

Short-reads is essential for current version of BASALT. We will update BASALT once we found accurate way to use only long-read datasets for binning.

1 If you have both short-read datasets and long-read datasets:
   ```
   BASALT -a assembly1.fa,assembly2.fa,assembly3.fa -s SR1_r1.fq,SR1_r2.fq/SR2_r1.fq,SR2_r2.fq -l lr1.fq,lr2.fq -t 60 -m 250
   ```
   You may put as many assemblies as you have, and as many SR or LR datasets as you have
   
2 If you have both short-read datasets only:
   ```
   BASALT -a assembly1.fa,assembly2.fa,assembly3.fa -s SR1_r1.fq,SR1_r2.fq/SR2_r1.fq,SR2_r2.fq -t 60 -m 250
   ```
3 If you have short-read datasets and Hi-C datasets:
   ```
   BASALT -a assembly1.fa,assembly2.fa,assembly3.fa -s SR1_r1.fq,SR1_r2.fq/SR2_r1.fq,SR2_r2.fq -c hc1.fq,hc2.fq -t 60 -m 250
   ```
4 If you have short-read datasets, long-read datasets and Hi-C datasets:
   ```
   BASALT -a assembly1.fa,assembly2.fa,assembly3.fa -s SR1_r1.fq,SR1_r2.fq/SR2_r1.fq,SR2_r2.fq -l lr1.fq,lr2.fq -c hc1.fq,hc2.fq -t 60 -m 250
   ```
5 If you have hifi dataset only:
   ```
   BASALT -a assembly1.fa -hf hifi1.fq -t 60 -m 250
   ```
6 If you have both hifi dataset and short-reads:
   ```
   BASALT -a assembly1.fa -hf hifi1.fq -s SR1_r1.fq,SR1_r2.fq -t 60 -m 250
   ```
7 If you want to use checkm2 instead of checkm for quality check:
   ```
   BASALT -a assembly1.fa -hf hifi1.fq -s SR1_r1.fq,SR1_r2.fq -t 60 -m 250 -qc checkm2
   ```


### PUBLICATIONS
Qiu, Z.*, Lian, C. A., Yuan, L., ... & Yu, K.*# (2024). BASALT refines binning from metagenomic data and increases resolution of genome-resolved metagenomic analysis, Nature Communications, in press
