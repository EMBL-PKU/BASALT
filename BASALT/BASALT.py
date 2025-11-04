#!/usr/bin/python
# -*- coding: UTF-8 -*-
#coding=utf-8

import time, sys, os
import argparse
from glob import glob

parser = argparse.ArgumentParser(description='BASALT')
parser.add_argument('-a', '--assemblies', type=str, dest='assemblies',
                    help='List of assemblies, e.g.: as1.fa,as2.fa')
parser.add_argument('-s', '--shortreads', type=str, dest='sr_datasets',
                    help='List of paired-end reads, e.g.: r1_1.fq,r1_2.fq/r2_1.fq,r2_2.fq (paried_ends reads need \'/\' to seperate)')
parser.add_argument('-l', '--longreads', type=str, dest='long_datasets',
                    help='Including ont and pb dataset, excluding hifi dataset. List of long reads, e.g.: lr1.fq,lr2.fq')
parser.add_argument('-hf', '--hifi', type=str, dest='hifi_datasets',
                    help='Hifi dataset. List of hifi, e.g.: hf1.fq,hf2.fq')
# parser.add_argument('-c','--HIC', type=str, dest='Hi_C_dataset',
                    # help='List of Hi-C dataset(s), e.g.: hc1.fq,hc2.fq')
parser.add_argument('-t','--threads', type=int, dest='threads', default=4,
                    help='Number of threads, e.g.: 64')
parser.add_argument('-m','--ram', type=int, dest='ram', default=32,
                    help='Number of ram, minimum ram suggested: 32G')
parser.add_argument('-e','--extra_binner', type=str, dest='extra_binner',
                    help='Extra binner for binning: m: metabinner, v: vamb; for instance: -e m, means BASALT will use metabinner for binning besides metabat2, maxbin2, and concoct')
parser.add_argument('-o','--out', type=str, dest='output_folder_name', default='Final_binset',
                    help='Name of the output folder. For binning, E.g. -o Anammox. BASALT would put those bins into folder Anammox_final_binset; for data feeding, e.g. -o Anammox; output files will under the folder of Anammox_data_feeded')
parser.add_argument('-q','--quality-check', type=str, dest='quality_check', default='checkm', 
                    help='Chance checkm version, default: checkm; you may use: \'-q checkm2\' to specify checkm for quality check when running BASALT')
parser.add_argument('--min-cpn', type=int, dest='Min_completeness', default=35,
                    help='Min completeness of kept bins (default: 35)')
parser.add_argument('--max-ctn', type=int, dest='Max_contamination', default=20,
                    help='Max contamination of kept bins (default: 20)')
parser.add_argument('--mode', type=str, dest='running_mode', default='continue',
                    help='Start a new project (new) or continue to run (continue). e.g. --mode continue / --mode new')
parser.add_argument('--module', type=str, dest='functional_module', default='all',
                    help='Modules for binning. Four modules: 1. autobinning; 2. refinement; 3. reassembly; 4. all. Default will run all modules. But you could set the only perform modle. e.g. --module reassembly. In the module, ')
# parser.add_argument('--autopara', type=str, dest='autobining_parameters', default='more-sensitive',
#                     help='Three parameters to chose: 1. more-sensitive; 2. sensitive; 3. quick. Default: more-sensitive. e.g. --autopara sensitive')
parser.add_argument('--refinepara', type=str, dest='refinement_paramter', default='quick',
                    help='Two refinement parameters to chose: 1. deep; 2. quick. Default: quick. e.g. --refpara deep')
parser.add_argument('--sensitive', type=str, dest='binning_sensitive', default='sensitive',
                    help='Three parameters to chose: 1. quick; 2. sensitive; 3. more-sensitive. Default: sensitive. If you want to change the sensitive, use: e.g. --sensitive sensitive')
# parser.add_argument('--hybrid', type=str, dest='hybrid_reassembly', default='no',
#                     help='Use hybrid re-assembly. e.g. --hybrid y / --hybrid n; defalt no')
# parser.add_argument('--compression', type=str, dest='compression', default=0,
#                     help='Two refinement parameters to chose: 1. deep; 2. quick. Default: deep. e.g. --refpara quick')
# parser.add_argument('--hifi-only', action='store_true', default=False, dest='hifi',
#                     help='Only Hifi data')
# parser.add_argument('--ont-only', action='store_true', default=False, dest='ont',
#                     help='Only ont data')
# parser.add_argument('--hybrid', action='store_true', default=False, dest='ont',
#                     help='Hybrid data, including both HTS and ont/pb datasets')
parser.add_argument('-d','--data-feeding-folder', type=str, dest='data_feeding_folder',
                    help='List of folder name of extra binset(s) for data feeding, e.g.: -d binset1_folder_name,binset2_folder_name')
parser.add_argument('--binset-index', type=str, dest='extra_binset_start_index', default=500,
                    help='Optional parameter for data feeding. The start index of the extra binset, e.g.: -bi 5. BASALT already set a default index, but if you already had 4 assemblies, the binset start index could be 5')
# parser.add_argument('--only-refinement', action='store_true', dest='only_refinement',
#                     help='Only carry out refinement, e.g.: --only-refinement')
parser.add_argument('-r', '--refinement-binset', type=str, dest='refinement_binset', default='',
                    help='Specify binset folder name for refinement e.g.: -r Human_gut_microbime_MAGs')
parser.add_argument('-c', '--coverage-list', type=str, dest='coverage_list',
                    help='List of depth file for refinement. Coverage file(s) could be generated from data feeding modole. e.g.: -c Coverage_matrix_for_binning_1_assembly.fa.txt,Coverage_matrix_for_binning_2_assembly.fa.txt')
parser.add_argument('-b', '--binsets-list', type=str, dest='binsets_list',
                    help='List of binsets for de-replication. Binset depth file(s) could be generated from data feeding modole. e.g.: -b 1_assembly_BestBinsSet,2_assembly_BestBinsSet')

args = parser.parse_args()
assemblies=args.assemblies
sr_datasets=args.sr_datasets
long_reads_list=args.long_datasets
hifi_list=args.hifi_datasets
# Hi_C_reads_list=args.Hi_C_dataset
QC_software=args.quality_check
num_threads=args.threads
ram=args.ram
extrabinner=args.extra_binner
min_cpn=args.Min_completeness
max_ctn=args.Max_contamination
continue_mode=new=args.running_mode
functional_module=args.functional_module
# autobining_parameters=args.autobining_parameters
refinement_paramter=args.refinement_paramter
output_folder=args.output_folder_name
sensitivity=args.binning_sensitive
data_feeding_folder=args.data_feeding_folder
binsetindex=args.extra_binset_start_index
# only_refinement=args.only_refinement
refinement_binset=args.refinement_binset
coverage_list=args.coverage_list
binsets_list=args.binsets_list

try:
    datasets_list=sr_datasets.split('/')
    datasets, n = {}, 0
    for item in datasets_list:
        n+=1
        datasets[str(n)]=[]
        pr=str(item).strip().split(',')
        datasets[str(n)].append(pr[0].strip())
        datasets[str(n)].append(pr[1].strip())
except:
    datasets={}

try:
    assembly_list=assemblies.split(',')
except:
    assembly_list=[]

try:
    lr_list2=long_reads_list.split(',')
    lr_list=[]
    for item in lr_list2:
        lr_list.append(item.strip())
except:
    lr_list=[]

try:
    hifi_list2=hifi_list.split(',')
    hifi_list=[]
    for item in hifi_list2:
        hifi_list.append(item.strip())
except:
    hifi_list=[]

try:
    hic_list2=Hi_C_reads_list.split(',')
    hic_list=[]
    for item in hic_list2:
        hic_list.append(item.strip())
except:
    hic_list=[]

try:
    eb_list2=extrabinner.split(',')
    eb_list=[]
    for item in eb_list2:
        eb_list.append(item.strip())
except:
    eb_list=[]

try:
    data_feeding_folder2=data_feeding_folder.split(',')
    data_feeding_folder=[]
    for item in data_feeding_folder2:
        data_feeding_folder.append(item.strip())
except:
    data_feeding_folder=[]

try:
    coverage_list2=coverage_list.split(',')
    coverage_list=[]
    for item in coverage_list2:
        coverage_list.append(item.strip())
except:
    coverage_list=[]

try:
    binsets_list2=binsets_list.split(',')
    binsets_list=[]
    for item in binsets_list2:
        binsets_list.append(item.strip())
except:
    binsets_list=[]

pwd=os.getcwd()
print('Processing assemblies:', str(assembly_list))
print('Processing short-reads:', str(datasets))
print('Processing long-reads:', str(lr_list))
print('Processing hifi-reads:', str(hifi_list))
print('Processing Hi-C reads:', str(hic_list))
print('Output folder name will be:', str(output_folder))
print('Process with extra binner:', str(eb_list))
print('Quality check software:', str(QC_software))
print('Binning sensitivity:', str(sensitivity))
print('Processing with:', str(num_threads), 'threads')
print('Processing with:', str(ram), 'G')
print('Running status:', str(continue_mode))
print('Binning module:', str(functional_module))
print('Min completeness:', str(min_cpn))
print('Max contamination:', str(max_ctn))
# print('Autobinning parameter:', str(autobining_parameters))
print('Refinement parameter:', str(refinement_paramter))
# print('Hybrid reassembly:', str(hybrid_reassembly))
print('Extra binset(s) for data feeding:', str(data_feeding_folder))
# print('Only refinement:', str(only_refinement))
print('Refinement binset:', str(refinement_binset))
print('List of coverage file(s):', str(coverage_list))
print('Binset(s) list:', str(binsets_list))

if continue_mode == 'continue':
    continue_mode = 'last'

if QC_software == 'checkm2':
    if len(data_feeding_folder) != 0:
        from Data_feeding import *
        if output_folder != 'Final_binset':
            output_folder=output_folder+'_data_feeded'
        else:
            output_folder='Data_feeded'
        data_feeding(data_feeding_folder, datasets, binsetindex, num_threads, output_folder, QC_software)
        
    elif len(binsets_list) != 0:
        from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_10242023 import *
        step='initial_drep' 
        # coverage_list: Coverage maxtrix
        # dataset: OK with both PE dataset or original datasets
        multiple_assembly_comparitor_main(assembly_list, binsets_list, coverage_list, datasets, step, num_threads)

    elif refinement_binset != '':
        pwd=os.getcwd()
        from S5_Outlier_remover_DL_11012023 import *
        # coverage_list: Coverage maxtrix
        # dataset: OK with both PE dataset or original datasets
        # assembly_list: data feeded assembly
        outlier_remover_main(refinement_binset, coverage_list, datasets, lr_list, hifi_list, assembly_list, pwd, num_threads)

    else:
        from BASALT_main_d import * 
        BASALT_main_d(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitivity, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)

else:
    if len(binsets_list) != 0:
        from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm import *
        step='initial_drep' 
        # coverage_list: Coverage maxtrix
        # dataset: OK with both PE dataset or original datasets
        multiple_assembly_comparitor_main(assembly_list, binsets_list, coverage_list, datasets, step, num_threads)

    elif refinement_binset != '':
        pwd=os.getcwd()
        from S5_Outlier_remover_DL_checkm import *
        # coverage_list: Coverage maxtrix
        # dataset: OK with both PE dataset or original datasets
        # assembly_list: data feeded assembly
        outlier_remover_main(refinement_binset, coverage_list, datasets, assembly_list, pwd, num_threads)

    else:
        if len(assembly_list) != 0:
            if functional_module == 'autobinning' or functional_module == 'all':
                from BASALT_main_c_autobinning import * 
                BASALT_main_c_autobinning(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitivity, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)

            if functional_module == 'refinement' or functional_module == 'all':  
                from BASALT_main_c_refinement import * 
                BASALT_main_c_refinement(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitivity, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)

            if functional_module == 'reassembly' or functional_module == 'all':
                from BASALT_main_c_re_assembly import * 
                BASALT_main_c_re_assembly(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitivity, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)

            # from BASALT_main_c import * 
            # BASALT_main_c(assembly_list, datasets, num_threads, lr_list, hifi_list, hic_list, eb_list, ram, continue_mode, functional_module, sensitivity, refinement_paramter, max_ctn, min_cpn, pwd, QC_software, output_folder)
        
            if len(data_feeding_folder) != 0:
                pwd=os.getcwd()
                from BASALT_main_c_datafeeding import *
                data_feeding_main(assembly_list, datasets, num_threads, data_feeding_folder, pwd, QC_software, output_folder, binsetindex, continue_mode)
            
            from Cleanup import cleanup
            cleanup(assembly_list)
            print('All accomplish!')

        else:
            if len(data_feeding_folder) != 0:
                from Data_feeding import *
                if output_folder != 'Final_binset':
                    output_folder=output_folder+'_data_feeded'
                else:
                    output_folder='Data_feeded'
                pe='y'
                data_feeding(data_feeding_folder, datasets, binsetindex, num_threads, output_folder, QC_software, pe)