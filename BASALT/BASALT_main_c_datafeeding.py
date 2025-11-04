#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

import time, sys, os, copy
from Bio import SeqIO
from Data_feeding import *
from S4_Multiple_Assembly_Comparitor_multiple_processes_bwt_checkm import *
from S5_Outlier_remover_DL_checkm import *
from glob import glob
from Cleanup import *

def data_feeding_main(assembly_list, datasets, num_threads, data_feeding_folder, pwd, QC_software, output_folder, binsetindex, continue_mode):
    #### Check existence of models
    user_dir = os.path.expanduser('~')
    local_dir = f"{user_dir}/.cache/BASALT"
    os.chdir(local_dir)
    model_list=glob(r'*_ensemble.csv')
    os.chdir(pwd)
    # print(model_list)
    if len(model_list) == 5:
        x=0
    else:
        print('BASALT models lacking. Start download the model')
        # os.system('python BASALT_models_download.py')
        os.system('BASALT_models_download.py')

    #### Program start
    last_step=0
    if continue_mode == 'last':
        try:
            n=0
            for line in open('Basalt_checkpoint.txt', 'r'):
                n+=1

            n1=0
            for line in open('Basalt_checkpoint.txt', 'r'):
                n1+=1
                if n1 == n:
                    ls=str(line)[0]
                    try:
                        ls2=int(str(line)[1])
                        last_step=int(str(ls)+str(ls2))
                    except:
                        last_step=int(ls)
                    # last_step=int(str(line).replace('th','').replace('st','').replace('nd','').replace('rd','').split(' ')[0])
        except:
            f_cp_m=open('Basalt_checkpoint.txt', 'w')
            f_cp_m.close()
    else:
        print('Start a new project')
        cleanup(assembly_list)
        f_cp_m=open('Basalt_checkpoint.txt', 'w')
        f_cp_m.close()

    print('BASALT started from step: '+str(last_step))
    try:
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT started from step: '+str(last_step)+'\n')
        fx.close()
    except:
        fx=open('Basalt_log.txt','w')
        fx.write('BASALT started from step: '+str(last_step)+'\n')
        fx.close()
    fx.close()
    
    try:
        f=open('BASALT_command.txt','a')
        f.write('BASALT restarted from step: '+str(last_step)+'\n')
        f.close()
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT restarted from step: '+str(last_step)+'\n')
        fx.close()
    except:
        f=open('BASALT_command.txt','w')
        f.write('BASALT started from 1st step'+'\n')
        f.write('Assemblies: '+str(str(assembly_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        f.write('Datasets: ')
        fx=open('Basalt_log.txt','a')
        fx.write('BASALT started from 1st step'+'\n')
        fx.write('Assemblies: '+str(str(assembly_list).replace('[','').replace(']','').replace(' ','').replace('\'',''))+'\n')
        fx.write('Datasets: ')
        datasets_list=datasets.split(',')
        x=0
        for ds in datasets_list:
            x+=1
            ds1=str(ds).split('[')[1].replace('\'','').replace(']','').replace(' ','')
            if x == 1:
                f.write(str(ds1))
                fx.write(str(ds1))
            else:
                f.write('/'+str(ds1))
                fx.write('/'+str(ds1))
        f.write('\n')
        fx.write('\n')

    print('Performing data feeding of '+str(data_feeding_folder))
    if output_folder != 'Final_binset':
        output_folder2=output_folder+'_data_feeded'
    else:
        output_folder2='Data_feeded'

    if last_step < 14:
        data_feeding_folder.append(output_folder)
        
        PE_datasets={}
        for item in datasets.keys():
            PE_datasets[item]=[datasets[item][0],datasets[item][1]]
    
        pe='n'
        data_feeding(data_feeding_folder, datasets, binsetindex, num_threads, output_folder2, QC_software, pe)
        f_cp_m=open('Basalt_checkpoint.txt', 'a')
        f_cp_m.write('14th Data feeding finish!'+'\n')
        f_cp_m.close()

    try:
        os.chdir(pwd+'/'+output_folder2)
        extra_condense_connection_list=glob(r'condense_connections_*')
        assembly_list2, binsets_list2, coverage_list2 = [], [], []
        for connection in extra_condense_connection_list:
            folder_name=str(connection).split('condense_connections_')[1].split('.txt')[0]
            index=int(folder_name.split('_')[0])
            binsets_list2.append(str(folder_name))
            assembly_list2.append(str(index)+'_assembly.fa')
            coverage_list2.append('Coverage_matrix_for_binning_'+str(index)+'_assembly.fa.txt')
            os.system('mv '+str(folder_name)+' '+str(index)+'_assembly.fa Coverage_matrix_for_binning_'+str(index)+'_assembly.fa.txt '+pwd)
    except:
        print(output_folder2+' absent')
    os.chdir(pwd)

    if last_step < 15:
        print('Processing de-replication after data feeding')
        step='initial_drep'
        try:
            os.system('mv BestBinset BestBinset_before_datafeeding')
        except:
            xyzzz=1
        multiple_assembly_comparitor_main(assembly_list2, binsets_list2, coverage_list2, PE_datasets, step, num_threads)
        os.system('rm -rf CC_*')
        f_cp_m=open('Basalt_checkpoint.txt', 'a')
        f_cp_m.write('15th Da-replication after data-feeding done!'+'\n')
        f_cp_m.close()

    if last_step < 16:
        print('Processing outlier removal after data feeding')
        try:
            os.system('mv BestBinset_outlier_refined BestBinset_outlier_refined_before_datafeeding')
        except:
            xyzzz=1

        outlier_remover_main('BestBinset', coverage_list2, PE_datasets, assembly_list2, pwd, num_threads)
        f_cp_m=open('Basalt_checkpoint.txt', 'a')
        f_cp_m.write('16th 2nd outlier removal done!'+'\n')
        f_cp_m.close()

    if last_step < 17:
        # for item in binsets_list2:
        #     os.system('rm -rf '+item)
        paired_bins, bins_list = {}, []
        for line in open(pwd+'/BestBinset/Bins_wth_sameQua_difSize.txt','r'):
            bins=str(line).strip().split('\t')[0].split('---')
            bin1,bin2 = bins[0], bins[1]
            bin1_list=bin1.split('.')
            bin1_list.pop()
            bin1_c='.'.join(bin1_list)
            bin2_list=bin2.split('.')
            bin2_list.pop()
            bin2_c='.'.join(bin2_list)
            paired_bins[bin1_c]=bin2_c
            bins_list.append(bins[0])
            bins_list.append(bins[1])
        
        os.system('mv '+output_folder+' '+output_folder+'_before_feeding')
        os.system('mv BestBinset_outlier_refined '+output_folder)

        selected_bin, remove_bin, red = [], [], []
        if len(paired_bins) != 0:
            os.chdir(output_folder)
            checkm={}
            for line in open('bin_stats_ext.tsv','r'):
                bin_id=str(line).strip().split('\t')[0]
                line=line.strip().replace('{','').replace('}','')
                comp=line.split('Completenes')[1].split(':')[1].split(',')[0]
                cont=line.split('Contamination')[1].split(':')[1]
                ml=line.split('marker lineage')[1].split(':')[1].split(',')[0]
                gz=line.split('Genome size')[1].split(':')[1].split(',')[0]
                checkm[bin_id]={}
                checkm[bin_id]['Completenes']=float(comp)
                checkm[bin_id]['Contamination']=float(cont)
                checkm[bin_id]['marker lineage']=str(ml)
                checkm[bin_id]['Genome size']=str(gz)

            for bin1 in paired_bins.keys():
                q1=checkm[bin1]['Completenes']-5*checkm[bin1]['Contamination']
                bin2=paired_bins[bin1]
                q2=checkm[bin2]['Completenes']-5*checkm[bin2]['Contamination']
                if q1 > q2:
                    selected_bin.append(bin1)
                    remove_bin.append(bin2)
                    os.system('rm '+bin2+'.fa')
                elif q1 < q2:
                    selected_bin.append(bin2)
                    remove_bin.append(bin1)
                    os.system('rm '+bin1+'.fa')
                else:
                    red.append(bin2)
                    remove_bin.append(bin2)

            if len(red) != 0:
                os.system('mkdir Potential_redundent_bins')
                for item in red:
                    os.system('mv '+item+'.fa Potential_redundent_bins')
            
            os.system('mv bin_stats_ext.tsv bin_stats_ext_before_de-replication.txt')
            f=open('bin_stats_ext.tsv','w')
            for item in checkm.keys():
                if item not in remove_bin:
                    f.write(item+'\t'+str(checkm[item])+'\n')
            f.close()
            os.chdir('Potential_redundent_bins')
            f=open('Redundant_bin_stats_ext.tsv','w')
            for item in red:
                 f.write(item+'\t'+str(checkm[item])+'\n')
            f.close()
            os.chdir(pwd)

        os.system('rm -rf '+output_folder2)
        f_cp_m=open('Basalt_checkpoint.txt', 'a')
        f_cp_m.write('17th 2nd de-replication done!'+'\n')
        f_cp_m.close()
    print('Data feeding and refinement done!')

if __name__ == '__main__': 
    assembly_list=['8_medium_S001_SPAdes_scaffolds.fasta','10_medium_cat_SPAdes_scaffolds.fasta']
    datasets={'1':['RM2_S001_insert_270_mate1.fq','RM2_S001_insert_270_mate2.fq'], '2':['RM2_S002_insert_270_mate1.fq','RM2_S002_insert_270_mate2.fq']}
    num_threads=30
    pwd=os.getcwd()
    QC_software='checkm' # checkm or checkm2
    output_folder='Final_binset'
    data_feeding_folder=[] # your list of extra binset(s)
    binsetindex=500
    continue_mode='last'
    data_feeding_main(assembly_list, datasets, num_threads, data_feeding_folder, pwd, QC_software, output_folder, binsetindex, continue_mode)
    
