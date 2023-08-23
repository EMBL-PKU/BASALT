## BASALT - Binning Across a Series of AssembLies Toolkit

BASALT is a versatile toolkit that recovers, compares and optimizes MAGs across a series of assemblies assembled from short-read, long-read or hybrid strategies. Firstly, BASALT uses high-throughput assembly methods to automatically assemble/co-assemble multiple files in parallel to reduce the manual input; Next, BASALT incorporates self-designed algorithms which automates the separation of redundant bins to elongate and refine best bins and improve contiguity; Further, BASALT facilitates state-of-art refinement tools using third-generation sequencing data to calibrate assembled bins and complement genome gaps that unable to be recalled from bins; Lastly, BASALT is an open frame toolkit that allows multiple integration of bioinformatic tools, which can optimize a wide range of datasets from various of assembly and binning software.

### SYSTEM REQUIREMENTS
Required dependencies:

•	Linux x64 systems

•	Python (>=3) modules: biopython, pandas, numpy, scikit-learn, copy, multiprocessing, collections，pytorch, tensorboardx

•	Perl

•	Java (>= v1.7)

•	Binning tools: MetaBAT2, MaxBin2, CONCOCT

•	Sequences processing tools: Bowtie2, BWA, SAMtools, CheckM, BLAST, Prodigal, HMMER

•	Sequences assembly and polishing tools: SPAdes, IDBA-UD, Pilon, Unicycler (optional)

Please see basalt_env.yml for details. The resource requirements for this pipeline will be based on the amount of data being processed. However, due to large memory requirements of many softwares used (e.g. SPAdes), 8+ cores and 125GB+ RAM are recommended. BASALT officially supports only Linux x64 systems.

### INSTALLATION
As there are many required dependencies needed to install, we highly recommend that you install and manage BASALT with conda. Once you have Miniconda or Anaconda installed, create a conda environment and install the BASALT within the specific environment. There is three ways to install the package. Please refer to the following description and contact us (yuke.sz@pku.edu.cn) if you have an issue when installing BASALT.

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
BASALT [-h] [-a ASSEMBLIES] [-s SR_DATASETS] [-l LONG_DATASETS] [-c HI_C_DATASET] [-t THREADS] [-m RAM]
       [-e EXTRA_BINNER] [--min-cpn MIN_COMPLETENESS] [--max-ctn MAX_CONTAMINATION] [--mode RUNNING_MODE]
       [--module FUNCTIONAL_MODULE] [--autopara AUTOBINING_PARAMETERS] [--refinepara REFINEMENT_PARAMTER]
```
2 If you use stanalone version of BASALT:
```
python BASALT.py [-h] [-a ASSEMBLIES] [-s SR_DATASETS] [-l LONG_DATASETS] [-c HI_C_DATASET] [-t THREADS] [-m RAM]
                 [-e EXTRA_BINNER] [--min-cpn MIN_COMPLETENESS] [--max-ctn MAX_CONTAMINATION] [--mode RUNNING_MODE]
                 [--module FUNCTIONAL_MODULE] [--autopara AUTOBINING_PARAMETERS] [--refinepara REFINEMENT_PARAMTER]
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
  -c HI_C_DATASET, --HIC HI_C_DATASET
                        List of Hi-C dataset(s), e.g.: hc1.fq,hc2.fq
  -t THREADS, --threads THREADS
                        Number of threads, e.g.: 64
  -m RAM, --ram RAM     Number of ram, minimum ram suggested: 32G
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

### PUBLICATIONS
Yu, K., Qiu, Z., Mu, R., Qiao, X., Zhang, L., Lian, C. A., ... & Zhuang, W. (2021). Recovery of high-qualitied Genomes from a deep-inland Salt Lake Using BASALT. bioRxiv.
