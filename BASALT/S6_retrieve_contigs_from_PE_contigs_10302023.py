#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio import SeqIO
import os, threading, glob, copy, gc
from sklearn.decomposition import PCA
import numpy as np
from multiprocessing import Pool
# from Outlier_remover import *

def lr_connecting_contigs(assembly, lr_connnecting_file, binset, bin_contig, num_threads):
    connecting_contigs = {}
    for bins in bin_contig.keys():
        connecting_contigs[bins]={}

    for line in open(lr_connnecting_file,'r'):
        try:
            contigs_list_l=str(line).strip().split('\t')[1].split(',')
            contigs_list=[]
            for item in contigs_list_l:
                contigs_list.append(item.split('\'')[1].strip())
                
            for contigs in contigs_list:
                for bins in bin_contig.keys():
                    if contigs in bin_contig[bins].keys():
                        for contig2 in contigs_list:
                            if contigs != contig2 and contig2 not in bin_contig[bins].keys():
                                try:
                                    connecting_contigs[bins][contig2]+=1
                                except:
                                    connecting_contigs[bins][contig2]=1
        except:
            xyzza=0
    
    fbc=open('LR_bin_contigs_connection.txt','w')
    connecting_contigs2, connecting_contigs2_level, xyz ={}, {}, 0
    for bins in connecting_contigs.keys():
        xyz+=1
        print('Processing LR connecting '+str(bins))
        bin_contig_sort, n = [], 0
        for contigs in connecting_contigs[bins].keys():
            if connecting_contigs[bins][contigs] >= 3:
                n+=1
                try:
                    connecting_contigs2[bins][contigs]=connecting_contigs[bins][contigs]
                    bin_contig_sort.append(connecting_contigs[bins][contigs])
                except:
                    connecting_contigs2[bins]={}
                    connecting_contigs2_level[bins]={}
                    connecting_contigs2[bins][contigs]=connecting_contigs[bins][contigs]
                    bin_contig_sort.append(connecting_contigs[bins][contigs])

                fbc.write(str(bins)+'\t'+str(contigs)+'\t'+str(connecting_contigs2[bins][contigs])+'\n')

        if n >=8:
            bin_contig_sort.sort(reverse=True)
            # m1, m2, m3, x = int(n/4), 2*int(n/4), 3*int(n/4), 0
            x = 0
            for item in bin_contig_sort:
                x+=1
                if x == 3:
                    Q2=item

            connecting_contigs2_level[bins][1]=[]
            # connecting_contigs2_level[bins][2]=[]
            # connecting_contigs2_level[bins][3]=[]
            # connecting_contigs2_level[bins][4]=[]

            for contigs in connecting_contigs2[bins].keys():
                if connecting_contigs2[bins][contigs] >= Q2:
                    connecting_contigs2_level[bins][1].append(contigs)
                # if connecting_contigs2[bins][contigs] >= Q3:
                #     connecting_contigs2_level[bins][1].append(contigs)
                # elif connecting_contigs2[bins][contigs] >= Q2:
                #     connecting_contigs2_level[bins][2].append(contigs)
                # elif connecting_contigs2[bins][contigs] >= Q1:
                #     connecting_contigs2_level[bins][3].append(contigs)
                # else:
                #     connecting_contigs2_level[bins][4].append(contigs)
        else:
            if bins in connecting_contigs2.keys():
                connecting_contigs2_level[bins]={}
                connecting_contigs2_level[bins][1]=[]
                for contigs in connecting_contigs2[bins].keys():
                    connecting_contigs2_level[bins][1].append(contigs)
    fbc.close()

    connecting_contigs3=copy.deepcopy(connecting_contigs2)
    # f=open(assembly+'_LR_connecting_contigs.txt','w')
    # f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
    for bins in connecting_contigs2.keys():
        if len(connecting_contigs2[bins]) == 0:
        #     f.write(str(bins)+'\t'+str(connecting_contigs2[bins])+'\n')
        # else:
            del connecting_contigs3[bins]
            del connecting_contigs2_level[bins]
    # f.close()

    connecting_contigs3_level=copy.deepcopy(connecting_contigs2_level)
    for bins in connecting_contigs2_level.keys():
        if len(connecting_contigs2_level[bins][1]) == 0:
            del connecting_contigs3_level[bins]

    # if len(connecting_contigs3) == 0:
    #     os.system('rm '+assembly+'_LR_connecting_contigs.txt')

    return connecting_contigs3, connecting_contigs3_level

def coverage_matrix_mpt(coverage_matrix_file, num):
    contig_cov={}
    n=0
    for line in open(str(coverage_matrix_file),'r'):
        n+=1
        if n >= 2:
            ids=str(line).strip().split('\t')[0]
            contig_cov[ids]={}
            for i in range(1,num+1):
                contig_cov[ids][i]=float(str(line).strip().split('\t')[3*i+1])
    return contig_cov, num

def bin_contig_recruite(bin_contig_cov, bin_contigs, bin_contigs_mock, contig_id, contig_cov, num):
    for bin_name in bin_contigs.keys():
        bin_contigs_mock[bin_name][contig_id]=1
        if len(bin_contigs_mock[bin_name]) == len(bin_contigs[bin_name]):
            bin_contig_cov[bin_name][contig_id]={}
            for i in range(1,num+1):
                bin_contig_cov[bin_name][contig_id][i]=contig_cov[contig_id][i]
        else:
            del bin_contigs_mock[bin_name][contig_id]
    return bin_contig_cov, num

def record_bin_coverage(binset, num_threads, coverage_matrix_list):
    pwd=os.getcwd()
    bin_contigs, contig_bin_list, bin_contigs_mock, m = {}, {}, {}, 0
    print('Parsing bins')
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        os.chdir(pwd+'/'+str(binset))
        for file in files:
            hz=file.split('.')[-1]
            if 'fa' in hz or 'fna' in hz:
                m+=1
                bin_contigs[file]={}
                bin_contigs_mock[file]={}
                for record in SeqIO.parse(file, 'fasta'):
                    contig_bin_list[record.id]=[]
                    bin_contigs[file][record.id]=record.seq
                    bin_contigs_mock[file][record.id]=0
                
                for record in SeqIO.parse(file, 'fasta'):
                    contig_bin_list[record.id].append(file)
    print('Parsed', m, 'bins')
    os.chdir(pwd)

    print('Recording the coverage of contigs from bins')
    x = 0 
    for item in coverage_matrix_list:
        # print('Parsing coverage_matrix_for_binning_'+str(item)+'.txt')
        x+=1
        if x == 1:
            n=0
            for line in open(item,'r'):
                n+=1
                if n == 2:
                    ls=str(line).strip().split('\t')
                    num=int((len(ls)-4)/3)

    pool=Pool(processes=num_threads)
    result, x = {}, 0
    for item in coverage_matrix_list:
        x+=1
        print('Processing '+str(item))
        result[item]=pool.apply_async(coverage_matrix_mpt, args=(item, num,))
    pool.close()
    pool.join()

    result2={}
    for item in result:
        result2[item]=result[item].get()

    contig_cov={}
    for item in result2.keys():
        contig_cov.update(result2[item][0])

    print('Writing bin coverage matrix file')
    os.system('mkdir Bin_coverage_after_contamination_removal')
    os.chdir('Bin_coverage_after_contamination_removal')
    f_coverage=open('Coverage_matrix_total.txt','w')
    for item in contig_cov.keys():
        f_coverage.write(str(item)+'\t'+str(contig_cov[item])+'\n')
    f_coverage.close()

    bin_contig_cov={}    
    for bin_name in bin_contigs.keys():
        bin_contig_cov[bin_name]={} 

    for contig_id in contig_bin_list.keys():
        for i in range(0, len(contig_bin_list[contig_id])):
            bin_name=contig_bin_list[contig_id][i]
            bin_contig_cov[bin_name][contig_id]={}
            for i2 in range(1,num+1):
                bin_contig_cov[bin_name][contig_id][i2]=contig_cov[contig_id][i2]

    for bins in bin_contig_cov.keys():
        f=open(bins+'_coverage_matrix.txt', 'w')
        # f.write('Contig'+'\t'+'Coverage'+'\n')
        for contigs in bin_contig_cov[bins].keys():
            f.write(str(contigs)+'\t'+str(bin_contig_cov[bins][contigs])+'\n')
        f.close()
    os.chdir(pwd)
    return bin_contig_cov, bin_contigs, contig_cov

def cycle_mt(bin_connecting_contigs, bin_connecting_contigs2, bin_connecting_contigs3, connections, bins):
    m, m_before, m_after=1, 0, 1
    print('Parsing', bins)
    bin_connecting_contigs_x, num_list = {}, []
    while m <= 1: ### Consider connection level less than 2
        # print(bins, 'cycle', m)
        if m_before != m_after:
            m_before = len(bin_connecting_contigs)
            for contigs in connections.keys():
                for item in connections[contigs].keys():
                    connections_num=int(eval(connections[contigs][item]))
                    if contigs in bin_connecting_contigs.keys() and item not in bin_connecting_contigs.keys():
                        bin_connecting_contigs[item]=m
                        bin_connecting_contigs2[item]=m
                        bin_connecting_contigs_x[item]=connections_num
                        num_list.append(connections_num)
            m_after = len(bin_connecting_contigs)
            m+=1
        else:
            m=6

    if len(bin_connecting_contigs_x) <= 20:
        bin_connecting_contigs3=bin_connecting_contigs_x
    else:
        num_list.sort(reverse=True)
        num_list_top=num_list[0:19]
        ### ranking
        for item in bin_connecting_contigs_x.keys():
            if bin_connecting_contigs_x[item] in num_list_top:
                bin_connecting_contigs3[item]=bin_connecting_contigs_x[item]

    return bin_connecting_contigs, bin_connecting_contigs2, bin_connecting_contigs3

def Parsing_kmer_file(assemblies_list, binset, bin_extract_contig, num_threads):
    print('Processing kmer file of '+str(assemblies_list))
    pwd=os.getcwd()
    bin_contigs, bin_contigs_mock = {}, {}
    os.chdir(pwd+'/'+binset)
    for root, dirs, files in os.walk(pwd+'/'+binset):
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                # qz=file.split('_genomes.')[0]
                if '.fa' in hz or '.fasta' in hz:
                    bin_contigs[file]={}
                    bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_contigs[file][str(record.id)]=0
                        bin_contigs_mock[file][str(record.id)]=0
    os.chdir(pwd)

    contig_bin_kmer, contig_kmer= {}, {}
    for assembly in assemblies_list:
        print('Parsing', assembly, 'kmer file')
        n=0
        for bin in bin_contigs.keys():
            contig_bin_kmer[bin]={}

    contig_kmer= {}
    for assembly in assemblies_list:
        n=0
        kmer_file=str(assembly)+'.kmer.txt'
        print('READing',str(kmer_file))
        for line in open(kmer_file,'r'):
            n+=1
            if n >= 2:
                ids=str(line).strip().split('\t')[0].strip()
                kmer_list=str(line).strip().split('\t')
                kmer_list.remove(kmer_list[0])
                kmer='\t'.join(kmer_list)
                contig_kmer[ids]=kmer
                for bin in bin_contigs.keys():
                    bin_contigs_mock[bin][ids]=0
                    if len(bin_contigs_mock[bin]) == len(bin_contigs[bin]):
                        contig_bin_kmer[bin][ids]=kmer
                    else:
                        del bin_contigs_mock[bin][ids]

    for bins in contig_bin_kmer.keys():
        f=open(bins+'_kmer.txt', 'w')
        for contigs in contig_bin_kmer[bins].keys():
            f.write(str(contigs)+'['+str(contig_bin_kmer[bins][contigs])+'\n')
        f.close()

        if bins in bin_extract_contig.keys():
            f=open(bins+'_connecting_contigs_kmer.txt', 'w')
            for connecting_contig in bin_extract_contig[bins].keys():
                f.write(str(connecting_contig)+'['+str(contig_kmer[connecting_contig])+'\n')
            f.close()

    return contig_bin_kmer, contig_kmer

def PE_connecting_contigs(assembly, PE_connections_file, binset, num_threads):
    pwd=os.getcwd()
    print('Finding common contigs of---', assembly, '---from---', binset, 'using ---', PE_connections_file)
    bin_contigs, bin_contigs_mock, bin_seq, bin_select_contigs={}, {}, {}, {}
    os.chdir(pwd+'/'+binset)
    for root, dirs, files in os.walk(pwd+'/'+binset):
        # os.chdir(pwd+'/'+str(binset_1_assembly))
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                qz=file.split('_genomes.')[0]
                if '.fa' in hz:
                    # print('Parsing', file
                    bin_contigs[file]={}
                    bin_select_contigs[file]={}
                    bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_contigs[file][str(record.id)]=0
                        bin_select_contigs[file][str(record.id)]=0
                        bin_contigs_mock[file][str(record.id)]=0
                        bin_seq[str(record.id)]=str(record.seq)
    os.chdir(pwd)

    ### Finding possible contigs within group of multiple autobinners
    # mirror_bins, n={}, 0
    # try:
    #     os.chdir(pwd+'/'+'comparison_files')
    #     for root, dirs, files in os.walk(pwd+'/'+'comparison_files'):
    #         for file in files:
    #             if 'Best_bin_set_iteration_' in file:
    #                 for line in open(file, 'r'):
    #                     bin1=str(line).strip().split('\t')[1].split('---')[0].strip()
    #                     bin2=str(line).strip().split('\t')[1].split('---')[1].strip()
    #                     if len(mirror_bins) == 0:
    #                         n=1
    #                         mirror_bins[n]={}
    #                         mirror_bins[n][bin1]=''
    #                         mirror_bins[n][bin2]=''
    #                     else:
    #                         for item in mirror_bins.keys():
    #                             if bin1 in mirror_bins[item].keys() or bin2 in mirror_bins[item].keys():
    #                                 mirror_bins[item][bin1]=''
    #                                 mirror_bins[item][bin2]=''
    # except:
    #     print('')

    print('Parsing', assembly, 'PE connections file')
    n, connections=0, {}
    for line in open(PE_connections_file,'r'):
        n+=1
        if n >= 2:
            id1=str(line).strip().split('\t')[0]
            id2=str(line).strip().split('\t')[2]
            connection_num=str(line).strip().split('\t')[3]
            if id1 not in connections.keys():
                connections[id1]={}
                connections[id1][str(id2)]=str(connection_num) ### Number suggests the connecting level
            else:
                connections[id1][str(id2)]=str(connection_num)

            if id2 not in connections.keys():
                connections[id2]={}
                connections[id2][str(id1)]=str(connection_num)
            else:
                connections[id2][str(id1)]=str(connection_num)

    f=open('Connections_'+str(assembly)+'.txt', 'w')
    for item in connections.keys():
        f.write(str(item)+'\t'+str(connections[item])+'\n')
    f.close()

    pool=Pool(processes=num_threads)
    bin_connecting_contigs, bin_connecting_contigs2, bin_connecting_contigs3, result={}, {}, {}, {}
    for bins in bin_contigs.keys():
        bin_connecting_contigs[bins]={}
        bin_connecting_contigs2[bins]={}
        bin_connecting_contigs3[bins]={}
        bin_connecting_contigs[bins]=bin_contigs[bins]
        result[bins]=pool.apply_async(cycle_mt,args=(bin_connecting_contigs[bins], bin_connecting_contigs2[bins], bin_connecting_contigs3[bins], connections, bins,))
    pool.close()
    pool.join()

    result2={}
    for bins in result:
        result2[bins]=result[bins].get()

    for bins in result:
        bin_connecting_contigs[bins].update(result2[bins][0])
        bin_connecting_contigs2[bins].update(result2[bins][1])
        bin_connecting_contigs3[bins].update(result2[bins][2])

    del_bin={}
    for item in bin_connecting_contigs2.keys():
        if len(bin_connecting_contigs2[item]) == 0:
            del_bin[item]=''

    for item in del_bin.keys():
        del bin_connecting_contigs2[item]

    f=open('Bin_connecting_contigs_'+str(assembly)+'.txt', 'w')
    f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
    for item in bin_connecting_contigs2.keys():
        f.write(str(item)+'\t'+str(bin_connecting_contigs2[item])+'\n')
    f.close()

    f=open('Bin_connecting_top20_contigs_'+str(assembly)+'.txt', 'w')
    f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
    for item in bin_connecting_contigs3.keys():
        f.write(str(item)+'\t'+str(bin_connecting_contigs3[item])+'\n')
    f.close()

    bin_connecting_contigs_level={}
    for item in bin_connecting_contigs3.keys():
        bin_connecting_contigs_level[item]={}
        for contigs in bin_connecting_contigs3[item].keys(): ### top contigs
            level=int(bin_connecting_contigs2[item][contigs])
            try:
                bin_connecting_contigs_level[item][level].append(contigs)
            except:
                bin_connecting_contigs_level[item][level]=[contigs]

    os.chdir(pwd)
    return bin_connecting_contigs3, connections, bin_connecting_contigs_level

def test_outlier(connecting_contig, item_data, test_index):
    # print('Judging', connecting_contig
    four = pd.Series(item_data).describe()
    # print(four)
    # print('Q1= {0}, Q2= {1}, Q3={2}'.format(four['25%'],four['50%'],four['75%']))
    Q1 = four['25%']
    Q3 = four['75%']
    IQR = Q3 - Q1
    upper1 = Q3 + 1.5 * IQR
    lower1 = Q1 - 1.5 * IQR

    stat=0
    if item_data[test_index] > float(upper1) or item_data[0] < float(lower1):
        stat=0
    else:
        stat=1
    return stat

def PCA_slector(data_array, num_contig):
    pca = PCA(n_components=1)
    pca.fit(data_array)
    explained_variance_ratio=pca.explained_variance_ratio_
    # print(explained_variance_ratio)
    # num=len(explained_variance_ratio)
    newData=pca.fit_transform(data_array)
    newData2=newData.reshape((1,num_contig))
    # print('Shape', num, num_contig
    newData_list=newData2.tolist()
    n=0
    for item in newData_list:
        n+=1
        if n == 1:
            newData_list_item=item
    return newData_list_item, explained_variance_ratio

def bin_comparison(original_bin_folder, new_bins_checkm, new_bin_folder, refinement_mode, num_threads):
    pwd=os.getcwd()
    print('Comparing bins before and after refining process')
    os.chdir(pwd+'/'+str(original_bin_folder))

    bin_checkm, bin_seq_rec, selected_bin = {}, {}, {}
    for root, dirs, files in os.walk(pwd+'/'+str(original_bin_folder)):
        for file in files:
            if 'quality_report.tsv' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2:
                        binID=str(line).strip().split('\t')[0].strip()
                        selected_bin[binID]=1
                        genome_size=str(line).strip().split('\t')[1].strip()
                        completeness=str(line).strip().split('\t')[2].strip()
                        contamination=str(line).strip().split('\t')[3].strip()
                        N50=str(line).strip().split('\t')[4].strip()

                        bin_checkm[str(binID)]={}
                        bin_checkm[str(binID)]['N50']=int(N50)
                        bin_checkm[str(binID)]['Completeness']=float(completeness)
                        bin_checkm[str(binID)]['Genome size']=float(genome_size)
                        bin_checkm[str(binID)]['Contamination']=float(contamination)
            else:
                bin_id_list=str(file).split('.')
                bin_id_list.remove(bin_id_list[-1])
                bin_id='.'.join(bin_id_list)
                if 'fa' in str(file).split('.')[-1]:
                    bin_seq_rec[bin_id]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_seq_rec[bin_id][record.id]=record.seq
    # os.chdir(pwd)

    os.chdir(pwd+'/'+str(new_bin_folder))
    bin_comparison, bin_comparison_list={}, {}
    for refined_bin in new_bins_checkm.keys():
        bin_id=refined_bin.split('.retrieved')[0]
        if bin_id in bin_checkm.keys():
            if bin_id not in bin_comparison.keys():
                bin_comparison[bin_id]=bin_id+'---'+refined_bin
                bin_comparison_list[bin_id]=[]
                bin_comparison_list[bin_id].append(bin_id)
                bin_comparison_list[bin_id].append(refined_bin)
            else:
                bin_comparison[bin_id]+='---'+refined_bin
                bin_comparison_list[bin_id].append(refined_bin)

    print(str(bin_comparison_list))
    f=open('Refined_bin_comparison.txt','w')
    f2=open('Deep_refine_bins.txt', 'w')
    f.write('Bin'+'\t'+'Related_bin'+'\t'+'checkm'+'\n')
    bestbin, further_refined_bin={}, {}
    for item in bin_comparison_list.keys():
        bestbin[item]={}
        write_out=str(item)+'\t'+(bin_comparison[item])
        for i in range(0, len(bin_comparison_list[item])):
            if i == 0:
                bin_id=bin_comparison_list[item][i]
                write_out+='\t'+str(bin_checkm[bin_id])
                bestbin[item][bin_id]={}
                bestbin[item][bin_id]=bin_checkm[bin_id]
            else:
                refined_bin_id=bin_comparison_list[item][i]
                write_out+='---'+str(new_bins_checkm[refined_bin_id])
                # re_connections=int(new_bins_checkm[refined_bin_id]['Connections'])
                re_cpn=float(new_bins_checkm[refined_bin_id]['Completeness'])
                re_ctn=float(new_bins_checkm[refined_bin_id]['Contamination'])
                # re_taxon=str(new_bins_checkm[refined_bin_id]['marker lineage'])
                # re_genome_size=float(new_bins_checkm[refined_bin_id]['Genome size'])
                re_delta=re_cpn-re_ctn
                re_5delta=re_cpn-5*re_ctn
                for bin_id2 in bestbin[item].keys():
                    ori_bin=bin_id2 
                    # ori_connections=int(bestbin[item][bin_id2]['Connections'])
                    ori_cpn=float(bestbin[item][bin_id2]['Completeness'])
                    ori_ctn=float(bestbin[item][bin_id2]['Contamination'])
                    # ori_taxon=str(bestbin[item][bin_id2]['marker lineage'])
                    # ori_genome_size=float(bestbin[item][bin_id2]['Genome size'])
                    ori_delta=ori_cpn-ori_ctn
                    ori_5delta=ori_cpn-5*ori_ctn

                if re_delta > ori_delta:
                    if re_5delta >= ori_5delta:
                        del bestbin[item][ori_bin]
                        bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]
                    # elif re_5delta < ori_5delta and ori_cpn >= 40 and ori_ctn <= 15 and re_ctn <= 50:
                    else:
                        further_refined_bin[refined_bin_id]=item
                        f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                    # elif re_cpn < ori_cpn and ori_cpn >= 40 and ori_ctn <= 15 and re_ctn <= 50:
                    #     further_refined_bin[refined_bin_id]=item
                    #     f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                elif re_cpn > ori_cpn:
                    if ori_cpn >= 40 and ori_ctn <= 10 and re_ctn <= 50:
                        further_refined_bin[refined_bin_id]=item
                        f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                else:
                    continue
        f.write(str(write_out)+'\n')
    f.close()
    f2.close()

    # print(str(further_refined_bin))

    f=open('Selected_bin.txt', 'w')
    f.write('Original bin'+'\t'+'Selected bin'+'\t'+'Factors'+'\n')
    for item in bestbin.keys():
        for bins in bestbin[item].keys():
            f.write(str(item)+'\t'+str(bins)+'\t'+str(bestbin[item][bins])+'\n')
            del selected_bin[item]
            selected_bin[bins]=1
    f.close()

    n, further_refined_bin_contig = 0, {}
    for root, dirs, files in os.walk(pwd+'/'+str(new_bin_folder)):
        for file in files:
            if '_genomes.' in file and '.fa' in file:
                # if 'retrieve' in file:
                #     ids_list=file.split('.')
                #     ids_list.remove(ids_list[-1])
                #     ids_list.remove(ids_list[-1])
                #     bin_id='.'.join(ids_list)
                # else:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                if bin_id in further_refined_bin.keys():
                    n+=1
                    org_bin_id = further_refined_bin[bin_id]
                    further_refined_bin_contig[org_bin_id]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        if record.id not in bin_seq_rec[org_bin_id].keys():
                            further_refined_bin_contig[org_bin_id][record.id]=record.seq
    print('Record', str(n), 'bins for deep refinement')
    # total_drf_num=n

    fx=open('tst.txt','w')
    fx.write(str(further_refined_bin_contig)+'\n'+str(further_refined_bin))
    fx.close()

    n = 0
    for root, dirs, files in os.walk(pwd+'/'+str(new_bin_folder)):
        for file in files:
            if '_genomes.' in file and '.fa' in file:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                if str(bin_id) not in selected_bin.keys():
                    n+=1
                    os.system('rm '+file)
    print('Removed', n, 'refined bins')

    n=0
    os.chdir(pwd+'/'+str(original_bin_folder))
    for root, dirs, files in os.walk(pwd+'/'+str(original_bin_folder)):
        for file in files:     
            if '_genomes.' in file:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                # if bin_id in bestbin.keys():
                if bin_id in selected_bin.keys():
                    n+=1
                    os.system('cp '+file+' '+pwd+'/'+str(new_bin_folder))
    print('Moved', n, 'to refined bins folder')

    os.chdir(pwd)

    replace_bin = {}
    if refinement_mode == 'deep':
        if len(further_refined_bin_contig) != 0:
            print('Processing deep refinement')
            try:
                f_dr=open('DeepRefinement_checkpoint.txt','a')
            except:
                f_dr=open('DeepRefinement_checkpoint.txt','w')
        #     f_dr.write('6th bins quality check done'+'\n')
            f_dr.close()

            n = 0
            for line in open('DeepRefinement_checkpoint.txt','r'):
                n+=1
            
            lp, n1 = 0, 0
            for line in open('DeepRefinement_checkpoint.txt','r'):
                n1+=1
                if n1 == n:
                    lp=int(eval(str(line)[0]))
            
            contig_neutral, contig_positive, contig_negative, contig_positive_suspect, contig_negative_suspect, contig_num = {}, {}, {}, {}, {}, {}
            for org_bin_id in further_refined_bin_contig.keys():
                contig_neutral[org_bin_id], contig_positive[org_bin_id], contig_negative[org_bin_id] = {}, {}, {}
                contig_positive_suspect[org_bin_id], contig_negative_suspect[org_bin_id], contig_num[org_bin_id] = {}, {}, {}
                n=0
                for contigs in further_refined_bin_contig[org_bin_id].keys():
                    n+=1
                    contig_num[org_bin_id][n]=str(contigs)
            
            try:
                os.system('mkdir Deep_retrieved_bins')
            except:
                print('Deep_retrieved_bins folder existed')
            
            try:
                os.system('mkdir Deep_retrieved_bins_split')
            except:
                print('Deep_retrieved_bins_split folder existed')

            os.chdir(pwd+'/'+str(new_bin_folder))
            if os.path.exists('Bin_contigs_white_black_list.txt'):
                fxxx=open('Bin_contigs_white_black_list.txt','a')
            else:
                fxxx=open('Bin_contigs_white_black_list.txt','w')
            os.chdir(pwd)

            if lp < 1:
                # m=0
                for org_bin_id in further_refined_bin_contig.keys():
                    os.chdir('Deep_retrieved_bins_split')
                    n=0
                    for contigs in further_refined_bin_contig[org_bin_id].keys():
                        n+=1
                        f=open(org_bin_id+'_deep_'+str(n)+'.fa','w')
                        f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                        # contig_num[org_bin_id][n]=str(contigs)
                        for contigs_id in bin_seq_rec[org_bin_id].keys():
                            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                        f.close()

                    os.chdir(pwd)
                
                try:
                    f_dr=open('DeepRefinement_checkpoint.txt','a')
                except:
                    f_dr=open('DeepRefinement_checkpoint.txt','w')
                f_dr.write('1st Writen retrieve bins done!'+'\n')
                f_dr.close()

            if lp < 2:
                # print(str(contig_positive_suspect))
                # os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(org_bin_id)+'_deep_retrieval -x fa -o '+str(org_bin_id)+'_deep_retrieval_checkm')
                # os.chdir(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm/')
                os.system('checkm2 predict -t '+str(num_threads)+' -i Deep_retrieved_bins_split -x fa -o Deep_retrieved_bins_split_checkm')
                # n = 0
                # for line in open('quality_report.tsv','r'):
                #     n+=1
                #     if n >= 2:
                #         bin_id=str(line).strip().split('\t')[0].strip()
                #         if bin_id in further_refined_bin_contig.keys():
                #             org_bin_cpn=float(str(line).strip().split('\t')[1].strip())
                #             org_bin_ctn=float(str(line).strip().split('\t')[2].strip())
                #             ori_delta=org_bin_cpn-5*org_bin_ctn
                #             ori_1delta=org_bin_cpn-org_bin_ctn
                try:
                    f_dr=open('DeepRefinement_checkpoint.txt','a')
                except:
                    f_dr=open('DeepRefinement_checkpoint.txt','w')
                f_dr.write('2nd retrieve bins quality check done!'+'\n')
                f_dr.close()
            
            if lp < 3:
                os.chdir(pwd+'/Deep_retrieved_bins_split_checkm/')
                n=0
                for line in open('quality_report.tsv','r'):
                    n+=1
                    if n >= 2:
                        refined_bin_id=str(line).strip().split('\t')[0].strip()
                        refined_bin_cpn=float(str(line).strip().split('\t')[1].strip())
                        refined_bin_ctn=float(str(line).strip().split('\t')[2].strip())
                        refined_delta=refined_bin_cpn-5*refined_bin_ctn
                        refined_1delta=refined_bin_cpn-refined_bin_ctn

                        org_bin_id=refined_bin_id.split('_deep_')[0]
                        for sel_bin in bestbin[org_bin_id].keys():
                            org_bin_cpn=float(bestbin[org_bin_id][sel_bin]['Completeness'])
                            org_bin_ctn=float(bestbin[org_bin_id][sel_bin]['Contamination'])
                        ori_delta=org_bin_cpn-5*org_bin_ctn
                        ori_1delta=org_bin_cpn-org_bin_ctn
                        # refined_delta=refined_bin_cpn-5*refined_bin_ctn

                        if refined_delta > ori_delta:
                            del bestbin[org_bin_id][sel_bin]
                            replace_bin[sel_bin]=''
                            del selected_bin[sel_bin]
                            selected_bin[refined_bin_id]=''
                            new_bins_checkm[refined_bin_id]={}
                            new_bins_checkm[refined_bin_id]['Completeness']=refined_bin_cpn
                            new_bins_checkm[refined_bin_id]['Contamination']=refined_bin_ctn
                            new_bins_checkm[refined_bin_id]['Genome size']=str(line).strip().split('\t')[8].strip()
                            new_bins_checkm[refined_bin_id]['N50']=int(str(line).strip().split('\t')[6].strip())
                            # new_bins_checkm[refined_bin_id]['marker lineage']=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                            bestbin[org_bin_id][refined_bin_id]=new_bins_checkm[refined_bin_id]
                            os.system('cp '+pwd+'/Deep_retrieved_bins_split/'+refined_bin_id+'.fa '+pwd+'/'+str(new_bin_folder))
                            try:
                                contigs=contig_num[org_bin_id][int(refined_bin_id.split('_deep_')[1].split('.fa')[0])]
                                contig_positive[org_bin_id][str(contigs)]=0
                                fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Positive'+'\n')
                            except:
                                os.chdir(pwd+'/Deep_retrieved_bins_split')
                                for record in SeqIO.parse(refined_bin_id+'.fa','fasta'):
                                    if record.id not in bin_seq_rec[org_bin_id].keys():
                                        try:
                                            contig_positive[org_bin_id][str(record.id)]=0
                                        except:
                                            contig_positive[org_bin_id]={}
                                            contig_positive[org_bin_id][str(record.id)]=str(record.seq)
                                            
                                        try:
                                            further_refined_bin_contig[org_bin_id][str(record.id)]=str(record.seq)
                                        except:
                                            further_refined_bin_contig[org_bin_id]={}
                                            further_refined_bin_contig[org_bin_id][str(record.id)]=str(record.seq)
                                os.chdir(pwd)
                                fxxx.write(str(org_bin_id)+'\t'+str(record.id)+'\t'+'Positive'+'\n')
                            
                        elif refined_delta == ori_delta:
                            try:
                                contigs=contig_num[org_bin_id][int(refined_bin_id.split('_deep_')[1].split('.fa')[0])]
                                contig_neutral[org_bin_id][str(contigs)]=0
                                fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Neutral'+'\n')
                            except:
                                os.chdir(pwd+'/Deep_retrieved_bins_split')
                                for record in SeqIO.parse(refined_bin_id+'.fa','fasta'):
                                    if record.id not in bin_seq_rec[org_bin_id].keys():
                                        try:
                                            contig_neutral[org_bin_id][str(record.id)]=0
                                        except:
                                            contig_neutral[org_bin_id]={}
                                            contig_neutral[org_bin_id][str(record.id)]=str(record.seq)

                                        try:
                                            further_refined_bin_contig[org_bin_id][str(record.id)]=str(record.seq)
                                        except:
                                            further_refined_bin_contig[org_bin_id]={}
                                            further_refined_bin_contig[org_bin_id][str(record.id)]=str(record.seq)
                                os.chdir(pwd)
                                fxxx.write(str(org_bin_id)+'\t'+str(record.id)+'\t'+'Neutral'+'\n')
                        # elif refined_delta < ori_delta:
                        #     if refined_1delta >= ori_1delta:
                        #         contig_positive_suspect[org_bin_id][str(contigs)]=0
                        #         fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Positive suspect'+'\n')
                        #     elif refined_bin_cpn > org_bin_cpn:
                        #         contig_negative_suspect[org_bin_id][str(contigs)]=0
                        #         fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Negative suspect'+'\n')
                        # else:
                        #     contig_negative[org_bin_id][str(contigs)]=0
                        #     fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Negative'+'\n')
                    # os.chdir(pwd)

                os.chdir(pwd+'/Deep_retrieved_bins')
                for org_bin_id in contig_neutral.keys():
                    print('Writing neutral contigs of bins of '+str(org_bin_id))
                    xxxx=0
                    if len(contig_neutral[org_bin_id]) != 0:
                        f=open(org_bin_id+'_RT-D1.fa','w')
                        for contigs in contig_neutral[org_bin_id].keys():
                            try:
                                f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                            except:
                                f.write('>'+str(contigs)+'\n'+str(contig_neutral[org_bin_id][contigs])+'\n')
                                
                        try:
                            for contigs in contig_positive[org_bin_id].keys():
                                try:
                                    f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                                except:
                                    f.write('>'+str(contigs)+'\n'+str(contig_positive[org_bin_id][contigs])+'\n')
                        except:
                            x=0

                        for contigs_id in bin_seq_rec[org_bin_id].keys():
                            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                        f.close()
                        xxxx=1

                for org_bin_id in contig_positive.keys():
                    print('Writing positive contigs of bins of '+str(org_bin_id))
                    if len(contig_positive[org_bin_id]) != 0:
                        # os.chdir(pwd+'/Deep_retrieved_bins')
                        f=open(org_bin_id+'_RT-D2.fa','w')
                        for contigs in contig_positive[org_bin_id].keys():
                            try:
                                f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                            except:
                                f.write('>'+str(contigs)+'\n'+str(contig_positive[org_bin_id][contigs])+'\n')
                        for contigs_id in bin_seq_rec[org_bin_id].keys():
                            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                        f.close()
                        xxxx=1
                
                os.chdir(pwd)

                try:
                    f_dr=open('DeepRefinement_checkpoint.txt','a')
                except:
                    f_dr=open('DeepRefinement_checkpoint.txt','w')
                f_dr.write('3rd collect positive contigs and wrote retrieve bins'+'\n')
                f_dr.close()

                # for org_bin_id in contig_positive_suspect.keys():
                #     if len(contig_positive_suspect[org_bin_id]) != 0:
                #         print('Writing positive_suspect contigs of bins of '+str(org_bin_id))
                #         # os.chdir(pwd+'/Deep_retrieved_bins')
                #         f=open(org_bin_id+'_RT-D3.fa','w')
                #         for contigs in contig_neutral[org_bin_id].keys():
                #             f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                #         for contigs in contig_positive[org_bin_id].keys():
                #             f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                #         for contigs in contig_positive_suspect[org_bin_id].keys():
                #             f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                #         for contigs_id in bin_seq_rec[org_bin_id].keys():
                #             f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                #         f.close()
                #         xxxx=1

                # for org_bin_id in contig_negative_suspect.keys():
                #     if len(contig_negative_suspect[org_bin_id]) != 0:
                #         os.chdir(pwd+'/Deep_retrieved_bins')
                #         f=open(org_bin_id+'_RT-D4.fa','w')
                #         for contigs in contig_negative_suspect[org_bin_id].keys():
                #             f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                #         for contigs in contig_neutral[org_bin_id].keys():
                #             f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                #         for contigs in contig_positive[org_bin_id].keys():
                #             f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                #         for contigs in contig_positive_suspect[org_bin_id].keys():
                #             f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                #         for contigs_id in bin_seq_rec[org_bin_id].keys():
                #             f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                #         f.close()
                #         xxxx=1
            os.chdir(pwd)

            if lp < 4:        
                xxxx=0
                for root, dirs, files in os.walk(pwd+'/Deep_retrieved_bins'):
                    for file in files:
                        xxxx+=1

                if xxxx>=1:
                    os.system('checkm2 predict -t '+str(num_threads)+' -i Deep_retrieved_bins -x fa -o Deep_retrieved_bins_checkm')

                    os.chdir(pwd+'/Deep_retrieved_bins_checkm/')
                    n=0
                    for line in open('quality_report.tsv','r'):
                        n+=1
                        if n >= 2:
                            refined_bin_id=str(line).strip().split('\t')[0].strip()
                            refined_bin_cpn=float(str(line).strip().split('\t')[1].strip())
                            refined_bin_ctn=float(str(line).strip().split('\t')[2].strip())
                            org_bin_id=str(line).strip().split('\t')[0].strip().split('_RT-D')[0]
                            for sel_bin in bestbin[org_bin_id].keys():
                                ori_cpn=float(bestbin[org_bin_id][sel_bin]['Completeness'])
                                ori_ctn=float(bestbin[org_bin_id][sel_bin]['Contamination'])
                            org_delta=ori_cpn-5*ori_ctn
                            refined_delta=refined_bin_cpn-5*refined_bin_ctn

                            if refined_delta > org_delta:
                                del bestbin[org_bin_id][sel_bin]
                                replace_bin[sel_bin]=''
                                del selected_bin[sel_bin]
                                selected_bin[refined_bin_id]=''
                                new_bins_checkm[refined_bin_id]={}
                                new_bins_checkm[refined_bin_id]['Completeness']=refined_bin_cpn
                                new_bins_checkm[refined_bin_id]['Contamination']=refined_bin_ctn
                                new_bins_checkm[refined_bin_id]['Genome size']=str(line).strip().split('\t')[8].strip()
                                new_bins_checkm[refined_bin_id]['N50']=int(str(line).strip().split('\t')[6].strip())
                                # new_bins_checkm[refined_bin_id]['marker lineage']=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                                bestbin[org_bin_id][refined_bin_id]=new_bins_checkm[refined_bin_id]
                                try:
                                    os.system('cp '+pwd+'/Deep_retrieved_bins/'+refined_bin_id+'.fa '+pwd+'/'+str(new_bin_folder))
                                except:
                                    os.system('cp '+pwd+'/Deep_retrieved_bins/'+refined_bin_id+'.fasta '+pwd+'/'+str(new_bin_folder))
                            elif refined_delta == org_delta:
                                if refined_bin_cpn < 85:
                                    if '_RT-D' in sel_bin and '_RT-D' in refined_bin_id:
                                        num1=int(eval(sel_bin.split('_RT-D')[1]))
                                        num2=int(eval(refined_bin_id.split('_RT-D')[1]))
                                        if num1 > num2:
                                            del bestbin[org_bin_id][sel_bin]
                                            replace_bin[sel_bin]=''
                                            del selected_bin[sel_bin]
                                            selected_bin[refined_bin_id]=''
                                            new_bins_checkm[refined_bin_id]={}
                                            new_bins_checkm[refined_bin_id]['Completeness']=refined_bin_cpn
                                            new_bins_checkm[refined_bin_id]['Contamination']=refined_bin_ctn
                                            new_bins_checkm[refined_bin_id]['Genome size']=str(line).strip().split('\t')[8].strip()
                                            new_bins_checkm[refined_bin_id]['N50']=int(str(line).strip().split('\t')[6].strip())
                                            # new_bins_checkm[refined_bin_id]['marker lineage']=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                                            bestbin[org_bin_id][refined_bin_id]=new_bins_checkm[refined_bin_id]
                                            try:
                                                os.system('cp '+pwd+'/Deep_retrieved_bins/'+refined_bin_id+'.fa '+pwd+'/'+str(new_bin_folder))
                                            except:
                                                os.system('cp '+pwd+'/Deep_retrieved_bins/'+refined_bin_id+'.fasta '+pwd+'/'+str(new_bin_folder))
                    os.chdir(pwd)
                    try:
                        f_dr=open('DeepRefinement_checkpoint.txt','a')
                    except:
                        f_dr=open('DeepRefinement_checkpoint.txt','w')
                    f_dr.write('4th collect positive contigs and wrote retrieve bins, quality check done!'+'\n')
                    f_dr.close()

            os.chdir(pwd+'/'+str(new_bin_folder))
        else:
            print('There is no bin in deep refinement list')
        
    os.chdir(pwd+'/'+str(new_bin_folder))
    for root, dirs, files in os.walk(pwd+'/'+str(new_bin_folder)):
        for file in files:
            if '_genomes.' in file and '.fa' in file:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                if str(bin_id) not in selected_bin.keys():
                    os.system('rm '+file)

    total_bin_checkm={}
    total_bin_checkm.update(bin_checkm)
    total_bin_checkm.update(new_bins_checkm)
    print(len(bin_checkm),' total number of orginal bins')
    print(len(new_bins_checkm),'total number of new bins')
    print(len(total_bin_checkm),' total number of bins')
    print(len(selected_bin), 'selected number of bins')

    f_bin_checkm=open(new_bin_folder+'_quality_report.tsv', 'w')
    f_bin_checkm.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completenes'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    for item in selected_bin.keys():
        if '_RT-D' in item:
            item_name=str(item).split('_RT-D')[0]
            f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item]['Genome size'])+'\t'+str(total_bin_checkm[item]['Completeness'])+'\t'+str(total_bin_checkm[item]['Contamination'])+'\t'+str(total_bin_checkm[item]['N50'])+'\n')
        elif '.retrieved_level' in item:
            item_name=str(item).split('.retrieved_level')[0]
            f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item]['Genome size'])+'\t'+str(total_bin_checkm[item]['Completeness'])+'\t'+str(total_bin_checkm[item]['Contamination'])+'\t'+str(total_bin_checkm[item]['N50'])+'\n')
        elif '_deep_' in item:
            item_name=str(item).split('_deep_')[0]
            f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item]['Genome size'])+'\t'+str(total_bin_checkm[item]['Completeness'])+'\t'+str(total_bin_checkm[item]['Contamination'])+'\t'+str(total_bin_checkm[item]['N50'])+'\n')
        else:
            f_bin_checkm.write(str(item)+'\t'+str(total_bin_checkm[item]['Genome size'])+'\t'+str(total_bin_checkm[item]['Completeness'])+'\t'+str(total_bin_checkm[item]['Contamination'])+'\t'+str(total_bin_checkm[item]['N50'])+'\n')
    f_bin_checkm.close()

    if len(replace_bin) != 0:
        f_replace=open(new_bin_folder+'_replace_bin_test.txt', 'w')
        for bins in replace_bin.keys():
            f_replace.write(str(bins)+'\n')
        f_replace.close()

    ### Rename
    f=open('Rename_files.txt','w')
    for root, dirs, files in os.walk(pwd+'/'+str(new_bin_folder)):
        for file in files:
            if '_genomes.' in file and '.fa' in file:
                if '_RT-D' in file or '.retrieved_level' in file or '_deep_' in file:
                    if '_RT-D' in file:
                        file_name=str(file).split('_RT-D')[0]
                    if '.retrieved_level' in file:
                        file_name=str(file).split('.retrieved_level')[0]
                    if '_deep_' in file:
                        file_name=str(file).split('_deep_')[0]
                    os.system('mv '+file+' '+str(file_name)+'.fa')
                    f.write(str(file)+'\t'+str(file_name)+'.fa'+'\n')
    f.close()
    os.chdir(pwd)   
    print('Contig retrieve done!')

def coverage_filtration(bin, m):
    os.system('S6p_coverage_filtration_mpt_06102022.py -b '+str(bin)+' -n '+str(m)+' -p coverage_filtration -c 1.5 -l 1')

def TNF_filtration(bin, m):
    os.system('S6p_coverage_filtration_mpt_06102022.py -b '+str(bin)+' -n '+str(m)+' -p TNF_filtration -c 1.5 -l 1')

def coverage_filtration_contigs(bin_connecting_contigs_total, bin_contig_cov, contig_cov, selected_1st, bin_extract_contig, elemimated_contig, elemimated_contig_total, connecting_contig, bin):
    test_contig_bin_cov={}
    test_contig_bin_cov.update(bin_contig_cov[bin])
    contig_num=len(bin_contig_cov[bin])

    if connecting_contig in contig_cov.keys():
        test_contig_bin_cov[connecting_contig]=contig_cov[connecting_contig]
        num_coverage=len(contig_cov[connecting_contig])

    ### Filtration of contigs
    cov_index_total_num, cov_index={}, {}
    for item in bin_contig_cov[bin].keys():  
        for i in range(0, num_coverage):
            if i not in cov_index_total_num.keys():
                cov_index_total_num[i]=float(bin_contig_cov[bin][item][i+1])
                cov_index[i]=[]
                cov_index[i].append(float(bin_contig_cov[bin][item][i+1]))
            else:
                cov_index_total_num[i]+=float(bin_contig_cov[bin][item][i+1])
                cov_index[i].append(float(bin_contig_cov[bin][item][i+1]))

    judgement, upper, lower = 0, {}, {}
    for item in cov_index.keys():
        cov_index[item].sort() ### sort low to high
        list_key_num=len(cov_index[item])
        p25=int(0.25*list_key_num)
        p75=int(0.75*list_key_num)

        n=0
        for i in cov_index[item]:
            n+=1
            if n == p25:
                Q1=i
            elif n == p75:
                Q3=i
            else:
                continue

        IQR = Q3 - Q1
        upper[item+1] = Q3 + 1.5 * IQR
        lower[item+1] = Q1 - 1.5 * IQR

    for i in range(1, num_coverage+1):
        if contig_cov[connecting_contig][i] <= upper[i] and contig_cov[connecting_contig][i] >= lower[i]:
            # print(str(i), str(contig_cov[connecting_contig][i]), str(upper[i]), str(lower[i])
            judgement+=1

    if judgement == num_coverage:
        if bin not in selected_1st.keys():
            selected_1st[bin]={}
            selected_1st[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
        else:
            selected_1st[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

        coverage_data, contigs_ids, coverage_list, n, test_index={}, [], [], 0, 0
        for item in test_contig_bin_cov.keys():
            n+=1
            contigs_ids.append(item)
            for i in range(1, num_coverage+1):
                if i not in coverage_data.keys():
                    coverage_data[i]=[]
                # print(i
                coverage_data[i].append(test_contig_bin_cov[item][i])
                coverage_list.append(test_contig_bin_cov[item][i])

            if item == connecting_contig:
                test_index=n-1

        coverage_array=np.array(coverage_list).reshape((n,num_coverage))
        A=PCA_slector(coverage_array, n)
        newData=A[0]
        explained_variance_ratio=A[1]
        bin_outlier=test_outlier(connecting_contig, newData, test_index)
        if bin_outlier == 1:
            if bin not in bin_extract_contig.keys():
                bin_extract_contig[bin]={}
                bin_extract_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            else:
                bin_extract_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
        else:
            if bin not in elemimated_contig.keys():
                elemimated_contig[bin]={}
                elemimated_contig_total[bin]={}
                elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            else:
                elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
    else:
        if bin not in elemimated_contig.keys():
            elemimated_contig[bin]={}
            elemimated_contig_total[bin]={}
            elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
        else:
            elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
    return selected_1st, bin_extract_contig, elemimated_contig, elemimated_contig_total

def recal_PE(bin_extract_contig, connection_total, bin_contig, bin):
    new_bin_connections={}
    for contigs in bin_extract_contig[bin].keys():
        if contigs in connection_total.keys():
            for connecting_contigs in connection_total[contigs].keys():
                if connecting_contigs not in bin_extract_contig[bin].keys() and connecting_contigs not in bin_contig[bin].keys():
                    if bin not in new_bin_connections.keys():
                        new_bin_connections[bin]=int(connection_total[contigs][connecting_contigs])
                    else:
                        new_bin_connections[bin]+=int(connection_total[contigs][connecting_contigs])

    for contigs in bin_contig[bin].keys():
        if contigs in connection_total.keys():
            for connecting_contigs in connection_total[contigs].keys():
                if connecting_contigs not in bin_extract_contig[bin].keys() and connecting_contigs not in bin_contig[bin].keys():
                    if bin not in new_bin_connections.keys():
                        new_bin_connections[bin]=int(connection_total[contigs][connecting_contigs])
                    else:
                        new_bin_connections[bin]+=int(connection_total[contigs][connecting_contigs])      
    return new_bin_connections

def kmer_cal(input_file, output_file):
    os.system('calc.kmerfreq.pl -i '+str(input_file)+' -o '+str(output_file))

def parse_dict(related_file):
    contig_level={}
    try:
        for line in open(str(related_file),'r'):
            bin_id=str(line).strip().split('\t')[0]
            contig_level[bin_id]={}
            contig_matrix_list=str(line).strip().replace('{','').replace('}','').split('\t')[1].split(',')
            for item in contig_matrix_list:
                contig_id=item.split('\'')[1].strip()
                level=int(item.split(':')[1].strip())
                contig_level[bin_id][contig_id]=level
    except:
        print('---')
    return contig_level

def parse_bin_in_bestbinset(assemblies_list, binset, outlier_remover_folder, PE_connections_list, lr_connection_list, num_threads, last_step, coverage_matrix_list, refinement_mode):
    pwd=os.getcwd()
    assemblies={}
    for item in assemblies_list:
        assemblies[item]=[]

    ### Record bins from source
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        os.chdir(pwd+'/'+str(binset))
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz or 'fna' in hz:
                for item in assemblies.keys():
                    if item in file:
                        assemblies[item].append(file)
    os.chdir(pwd)

    elemimated_contig_total  = {}
    A=record_bin_coverage(binset, num_threads, coverage_matrix_list)
    bin_contig_cov=A[0]
    bin_contig=A[1]
    contig_cov=A[2]

    if last_step == 0:
        bin_connecting_contigs_total_sr, connection_total, bin_connecting_contigs_total_level_sr = {}, {}, {}
        if len(PE_connections_list) != 0:
            for i in range(0, len(assemblies_list)):
                A=PE_connecting_contigs(assemblies_list[i], PE_connections_list[i], binset, num_threads)
                bin_connecting_contigs=A[0]
                connections_single=A[1]
                bin_connecting_contigs_level=A[2]
                bin_connecting_contigs_total_sr.update(bin_connecting_contigs)
                connection_total.update(connections_single)
                bin_connecting_contigs_total_level_sr.update(bin_connecting_contigs_level)
                print('Parsed', str(assemblies_list[i]))
                print('--------------------')

        # bin_connecting_contigs_total, connection_total, bin_connecting_contigs_total_level = {}, {}, {}
        # for i in range(0, len(assemblies_list)):
        #     A=PE_connecting_contigs(assemblies_list[i], PE_connections_list[i], binset, num_threads)
        #     bin_connecting_contigs=A[0]
        #     connections_single=A[1]
        #     bin_connecting_contigs_level=A[2]
        #     bin_connecting_contigs_total.update(bin_connecting_contigs)
        #     connection_total.update(connections_single)
        #     bin_connecting_contigs_total_level.update(bin_connecting_contigs_level)
        #     print('Parsed', str(assemblies_list[i]))
        #     print('--------------------')

        bin_connecting_contigs_total_lr, bin_connecting_contigs_total_lr_level = {}, {}
        if len(lr_connection_list) != 0:
            for i in range(0, len(assemblies_list)):
                connecting_contigs, connecting_contigs2_level=lr_connecting_contigs(assemblies_list[i], lr_connection_list, binset, bin_contig, num_threads)
                bin_connecting_contigs_total_lr.update(connecting_contigs)
                bin_connecting_contigs_total_lr_level.update(connecting_contigs2_level)

        f=open('LR_connecting_contigs.txt','w')
        f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
        for bins in bin_connecting_contigs_total_lr.keys():
            f.write(str(bins)+'\t'+str(bin_connecting_contigs_total_lr[bins])+'\n')
        f.close()

        f2=open('LR_connecting_contigs_level.txt','w')
        f2.write('Bin'+'\t'+'Level'+'\t'+'Connecting Contigs'+'\n')
        for bins in bin_connecting_contigs_total_lr_level.keys():
            for level in bin_connecting_contigs_total_lr_level[bins].keys():
                f2.write(str(bins)+'\t'+str(level)+'\t'+str(bin_connecting_contigs_total_lr_level[bins][level])+'\n')
        f2.close()

        f=open('SR_connecting_contigs.txt','w')
        f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
        for bins in bin_connecting_contigs_total_sr.keys():
            f.write(str(bins)+'\t'+str(bin_connecting_contigs_total_sr[bins])+'\n')
        f.close()

        f2=open('SR_connecting_contigs_level.txt','w')
        f2.write('Bin'+'\t'+'Level'+'\t'+'Connecting Contigs'+'\n')
        for bins in bin_connecting_contigs_total_level_sr.keys():
            for level in bin_connecting_contigs_total_level_sr[bins].keys():
                f2.write(str(bins)+'\t'+str(level)+'\t'+str(bin_connecting_contigs_total_level_sr[bins][level])+'\n')
        f2.close()

        bin_connecting_contigs_total, bin_connecting_contigs_total_level = {}, {}
        if len(bin_connecting_contigs_total_level_sr) !=0 and len(bin_connecting_contigs_total_lr) != 0:
            print('Merging both LR and SR connecting contigs')
            bin_connecting_contigs_total.update(bin_connecting_contigs_total_sr)
            bin_connecting_contigs_total.update(bin_connecting_contigs_total_lr)
            for bins in bin_connecting_contigs_total_lr_level.keys():
                bin_connecting_contigs_total_level[bins]={}
                if bins in bin_connecting_contigs_total_level_sr.keys():
                    for contigs in bin_connecting_contigs_total_lr_level[bins][1]:
                        try:
                            if contigs in bin_connecting_contigs_total_level_sr[bins][1]:
                                try:
                                    bin_connecting_contigs_total_level[bins][1].append(contigs)
                                except:
                                    bin_connecting_contigs_total_level[bins][1]=[]
                                    bin_connecting_contigs_total_level[bins][1].append(contigs)
                            else:
                                try:
                                    bin_connecting_contigs_total_level[bins][3].append(contigs)
                                except:
                                    bin_connecting_contigs_total_level[bins][3]=[]
                                    bin_connecting_contigs_total_level[bins][3].append(contigs)
                        except:
                            try:
                                bin_connecting_contigs_total_level[bins][3].append(contigs)
                            except:
                                bin_connecting_contigs_total_level[bins][3]=[]
                                bin_connecting_contigs_total_level[bins][3].append(contigs)
                    
                    try:
                        for contigs in bin_connecting_contigs_total_level_sr[bins][1]:
                            try:
                                if contigs not in bin_connecting_contigs_total_lr_level[bins][1]:
                                    try:
                                        bin_connecting_contigs_total_level[bins][2].append(contigs)
                                    except:
                                        bin_connecting_contigs_total_level[bins][2]=[]
                                        bin_connecting_contigs_total_level[bins][2].append(contigs)
                            except:
                                try:
                                    bin_connecting_contigs_total_level[bins][2].append(contigs)
                                except:
                                    bin_connecting_contigs_total_level[bins][2]=[]
                                    bin_connecting_contigs_total_level[bins][2].append(contigs)
                    except:
                        print('There is not connecting sequences from '+bins+' at level 1')

                else:
                    bin_connecting_contigs_total_level[bins].update(bin_connecting_contigs_total_lr_level[bins])
            
            for bins in bin_connecting_contigs_total_level_sr.keys():

                if bins not in bin_connecting_contigs_total_lr_level.keys():
                    bin_connecting_contigs_total_level[bins]={}
                    bin_connecting_contigs_total_level[bins].update(bin_connecting_contigs_total_level_sr[bins])

        elif len(bin_connecting_contigs_total_lr) != 0:
            print('Sorting both LR connecting contigs')
            bin_connecting_contigs_total.update(bin_connecting_contigs_total_lr)
            bin_connecting_contigs_total_level.update(bin_connecting_contigs_total_lr_level)
        elif len(bin_connecting_contigs_total_level_sr) != 0:
            print('Sorting both SR connecting contigs')
            bin_connecting_contigs_total.update(bin_connecting_contigs_total_sr)
            bin_connecting_contigs_total_level.update(bin_connecting_contigs_total_level_sr)

        f=open('Bin_connecting_contigs.txt','w')
        f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
        for bins in bin_connecting_contigs_total.keys():
            f.write(str(bins)+'\t'+str(bin_connecting_contigs_total[bins])+'\n')
        f.close()

        f2=open('Bin_connecting_contigs_level.txt','w')
        f2.write('Bin'+'\t'+'Level'+'\t'+'Connecting Contigs'+'\n')
        for bins in bin_connecting_contigs_total_level.keys():
            for level in bin_connecting_contigs_total_level[bins].keys():
                f2.write(str(bins)+'\t'+str(level)+'\t'+str(bin_connecting_contigs_total_level[bins][level])+'\n')
        f2.close()

        # f=open('Bin_connecting_contigs.txt', 'w')
        # f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
        # for item in bin_connecting_contigs_total.keys():
        #     f.write(str(item)+'\t'+str(bin_connecting_contigs_total[item])+'\n')
        # f.close()

        # f=open('Bin_connecting_contigs_level.txt', 'w')
        # f.write('Bin'+'\t'+'Level'+'\t'+'Connecting contigs'+'\n')
        # for item in bin_connecting_contigs_total_level.keys():
        #     for level in bin_connecting_contigs_total_level[item].keys():
        #         f.write(str(item)+'\t'+str(level)+'\t'+str(bin_connecting_contigs_total_level[item][level])+'\n')
        # f.close()

        f_cp=open('S6_checkpoint.txt','a')
        f_cp.write('1st Bin_connecting done'+'\n')
        f_cp.close()

    if int(last_step) >= 1:
        black_contigs=finding_black_contigs(binset, outlier_remover_folder, pwd)
        bin_connecting_contigs_total, bin_connecting_contigs_total_level, n = {}, {}, 0
        for line in open('Bin_connecting_contigs.txt', 'r'):
            n+=1
            if n >= 2:
                bins=str(line).strip().split('\t')[0]
                bin_name_list=bins.split('.')
                bin_name_list.pop()
                bin_name='.'.join(bin_name_list)
                bin_connecting_contigs_total[bins]={}
                dict_key_list=str(line).strip().replace('{','').replace('}','').split('\t')[1].split(',')
                if bin_name in black_contigs.keys():
                    for item in dict_key_list:
                        contig=item.split('\'')[1]
                        if contig not in black_contigs[bin_name]:
                            bin_connecting_contigs_total[bins][contig]=int(item.split(':')[1].strip())
                else:
                    for item in dict_key_list:
                        try:
                            contig=item.split('\'')[1]
                            bin_connecting_contigs_total[bins][contig]=int(item.split(':')[1].strip())
                        except:
                            xyz=0
                            
        n=0
        for line in open('Bin_connecting_contigs_level.txt', 'r'):
            n+=1
            if n >= 2:
                bins=str(line).strip().split('\t')[0]
                bin_connecting_contigs_total_level[bins]={}
        
        n=0
        for line in open('Bin_connecting_contigs_level.txt', 'r'):
            n+=1
            if n >= 2:
                bins=str(line).strip().split('\t')[0]
                level=int(str(line).strip().split('\t')[1])
                bin_connecting_contigs_total_level[bins][level]=[]

        n=0
        for line in open('Bin_connecting_contigs_level.txt', 'r'):
            n+=1
            if n >= 2:
                bins=str(line).strip().split('\t')[0]
                level=int(str(line).strip().split('\t')[1])
                contig_list=str(line).strip().split('\t')[2].replace('[','').replace(']','').replace('\'','').replace(' ','').split(',')
                for item in contig_list:
                    bin_connecting_contigs_total_level[bins][level].append(item)

    if int(last_step) < 2:
        pool=Pool(processes=num_threads)
        bin_extract_contig, selected_1st, elemimated_contig, m={}, {}, {}, 0
        selected_files, extract_files, elemimated_files = [], [], [] 
        for bin in bin_connecting_contigs_total.keys():
            m+=1
            pool.apply_async(coverage_filtration,args=(bin, m,))
            selected_files.append(str(bin)+'_selected_contigs.txt')
            extract_files.append(str(bin)+'_bin_extract_contig.txt')
            elemimated_files.append(str(bin)+'_elemimated_contig.txt')
        pool.close()
        pool.join()

        try:
            os.system('mkdir S6_coverage_filtration_matrix')
        except:
            print('Folder S6_coverage_filtration_matrix existed. Deleteced and recreated.')
            os.system('rm -rf S6_coverage_filtration_matrix')
            os.system('mkdir S6_coverage_filtration_matrix')

        for item in selected_files:
            selected_1st.update(parse_dict(item))
            os.system('mv '+str(item)+' S6_coverage_filtration_matrix')

        for item in extract_files:
            bin_extract_contig.update(parse_dict(item))
            os.system('mv '+str(item)+' S6_coverage_filtration_matrix')

        for item in elemimated_files:
            elemimated_contig.update(parse_dict(item))
            os.system('mv '+str(item)+' S6_coverage_filtration_matrix')
            
        elemimated_contig_total.update(elemimated_contig)

        print('Coverage filtration done!')
        f=open('Coverage_1st_filtrated_bin_connecting_contigs.txt','w')
        for bin in selected_1st.keys():
            f.write(str(bin)+'\t'+str(selected_1st[bin])+'\n')
        f.close()

        f=open('Bin_extract_contigs_after_coverage_filtration.txt','w')
        for bin in bin_extract_contig.keys():
            f.write(str(bin)+'\t'+str(bin_extract_contig[bin])+'\n')
        f.close()

        f=open('Coverage_eliminated_bin_connecting_contigs.txt','w')
        for bin in elemimated_contig.keys():
            f.write(str(bin)+'\t'+str(elemimated_contig[bin])+'\n')
        f.close()

        f_cp=open('S6_checkpoint.txt','a')
        f_cp.write('2nd Bin_connecting done'+'\n')
        f_cp.close()

    if int(last_step) < 3:
        kmer_list=glob.glob(r'*.kmer.txt')
        if len(kmer_list) == 0:
            print('There is no kmer file in this folder. Re-cal kmer')
            kmer_list = []
            pool=Pool(processes=num_threads)
            for i in range(0, len(assemblies_list)):
                print('Parsing kmer', str(assemblies_list[i]))
                # kmer_list.append(str(assemblies_list[i])+'.kmer.txt')
                pool.apply_async(kmer_cal, args=(str(assemblies_list[i]), str(assemblies_list[i])+'.kmer.txt',))
            pool.close()
            pool.join()

        f_cp=open('S6_checkpoint.txt','a')
        f_cp.write('3rd Kmer calculation done'+'\n')
        f_cp.close()

    if int(last_step) >= 2:
        selected_1st, bin_extract_contig, elemimated_contig_total  = {}, {}, {}
        for line in open('Coverage_1st_filtrated_bin_connecting_contigs.txt','r'):
            bins=str(line).strip().split('\t')[0]
            selected_1st[bins]={}
            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
            for item in dict_key_list:
                selected_1st[bins][item.split('\'')[1]]=int(item.split(':')[1].strip()) ###

        for line in open('Bin_extract_contigs_after_coverage_filtration.txt','r'):
            bins=str(line).strip().split('\t')[0]
            bin_extract_contig[bins]={}
            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
            for item in dict_key_list:
                bin_extract_contig[bins][item.split('\'')[1]]=item.split(':')[1].strip() ###

        for line in open('Coverage_eliminated_bin_connecting_contigs.txt','r'):
            bins=str(line).strip().split('\t')[0]
            elemimated_contig_total[bins]={}
            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
            for item in dict_key_list:
                elemimated_contig_total[bins][item.split('\'')[1]]=item.split(':')[1].strip() ###

    if int(last_step) < 4:
        # for i in range(0, len(assemblies_list)):
        #     Parsing_kmer_file(assemblies_list[i], binset, str(assemblies_list[i])+'.kmer.txt', bin_extract_contig, num_threads)
        Parsing_kmer_file(assemblies_list, binset, bin_extract_contig, num_threads)

        try:
            os.mkdir('Bin_kmer')
        except:
            print('Bin_kmer exists')

        pool=Pool(processes=num_threads)
        bin_extract_contig_TNF, elemimated_contig_TNF, TNFs_exceptional_contigs, bin_extract_contig_TNF_list, elemimated_contig_TNF_list, TNFs_exceptional_contigs_list, m ={}, {}, {}, [], [], [], 0
        for bin in bin_connecting_contigs_total.keys():
            if bin in bin_extract_contig.keys():
                m+=1
                pool.apply_async(TNF_filtration,args=(bin, m,))
                bin_extract_contig_TNF_list.append(str(bin)+'_bin_extract_contig_TNF.txt')
                elemimated_contig_TNF_list.append(str(bin)+'_elemimated_contig_TNF.txt')
                TNFs_exceptional_contigs_list.append(str(bin)+'_TNFs_exceptional_contigs.txt')
        pool.close()
        pool.join()

        try:
            os.system('mkdir S6_TNF_filtration_matrix')
        except:
            # print('Folder S6_TNF_filtration_matrix existed')
            print('Folder S6_TNF_filtration_matrix existed. Deleteced and recreated.')
            os.system('rm -rf S6_TNF_filtration_matrix')
            os.system('mkdir S6_TNF_filtration_matrix')

        for item in bin_extract_contig_TNF_list:
            bin_extract_contig_TNF.update(parse_dict(item))
            os.system('mv '+str(item)+' S6_TNF_filtration_matrix')

        for item in elemimated_contig_TNF_list:
            elemimated_contig_TNF.update(parse_dict(item))
            os.system('mv '+str(item)+' S6_TNF_filtration_matrix')

        elemimated_contig_total.update(elemimated_contig_TNF)

        for item in TNFs_exceptional_contigs_list:
            TNFs_exceptional_contigs.update(parse_dict(item))
            os.system('mv '+str(item)+' S6_TNF_filtration_matrix')    

        f=open('TNF_filtrated_bin_connecting_contigs.txt','w')
        for bin in bin_extract_contig_TNF.keys():
            f.write(str(bin)+'\t'+str(bin_extract_contig_TNF[bin])+'\n')
        f.close()

        f=open('TNF_eliminated_bin_connecting_contigs.txt','w')
        for bin in elemimated_contig_TNF.keys():
            f.write(str(bin)+'\t'+str(elemimated_contig_TNF[bin])+'\n')
        f.close()

        f=open('Total_eliminated_bin_connecting_contigs.txt','w')
        for bin in elemimated_contig_total.keys():
            f.write(str(bin)+'\t'+str(elemimated_contig_total[bin])+'\n')
        f.close()

        f=open('TNFs_exceptional_contigs.txt','w')
        for bins in TNFs_exceptional_contigs.keys():
            for contigs in TNFs_exceptional_contigs[bins].keys():
                f.write(str(bins)+'\t'+str(contigs)+'\t'+str(TNFs_exceptional_contigs[bins][contigs])+'\n')
        f.close()

        f_cp=open('S6_checkpoint.txt','a')
        f_cp.write('4th TNF filtration done'+'\n')
        f_cp.close()

    if int(last_step) == 4:
        bin_extract_contig_TNF={}
        for line in open('TNF_filtrated_bin_connecting_contigs.txt','r'):
            bins=str(line).strip().split('\t')[0]
            bin_extract_contig_TNF[bins]={}
            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
            for item in dict_key_list:
                 bin_extract_contig_TNF[bins][item.split('\'')[1]]=int(item.split(':')[1].strip())

    if int(last_step) < 5:
        print('Writing retrieved bins')
        try:
            os.system('mkdir '+binset+'_retrieved')
        except:
            print(str(binset)+'_retrieved existed. Deleted and recreated')
            os.system('rm -rf '+binset+'_retrieved')
            os.system('mkdir '+binset+'_retrieved')
        
        total_bin_contigs={}
        for assemblies in assemblies_list:
            for record in SeqIO.parse(assemblies, 'fasta'):
                total_bin_contigs[record.id]=record.seq

        os.chdir(binset+'_retrieved')
        bin_level_contigs={}
        for item in bin_extract_contig_TNF.keys():
        # for item in bin_extract_contig.keys():
            bin_name_list=str(item).split('.')
            bin_name_list.remove(bin_name_list[-1])
            bin_name='.'.join(bin_name_list)
            for level in bin_connecting_contigs_total_level[item].keys():
                f=open(bin_name+'.retrieved_level'+str(level)+'.fa','w')
                bin_level_contigs[bin_name]={}
                bin_level_contigs[bin_name]['.retrieved_level'+str(level)+'.fa']={}
                if item in bin_contig.keys():
                    for ids in bin_contig[item].keys():
                        f.write('>'+str(ids)+'\n'+str(bin_contig[item][ids])+'\n')
                        bin_level_contigs[bin_name]['.retrieved_level'+str(level)+'.fa'][str(ids)]=''

                for contigs in bin_extract_contig_TNF[item].keys():
                # for contigs in bin_extract_contig[item].keys():
                    f.write('>'+str(contigs)+'\n'+str(total_bin_contigs[contigs])+'\n')
                    bin_level_contigs[bin_name]['.retrieved_level'+str(level)+'.fa'][str(contigs)]=''
                f.close()

        del total_bin_contigs ### Release ram
        gc.collect() ### Cleanup

        for bins in bin_level_contigs.keys():
            m=len(bin_level_contigs[bins])
            for i in range(0,m-1):
                m1=m-i
                m2=m-i-1
                if bin_level_contigs[bins]['.retrieved_level'+str(m1)+'.fa'] == bin_level_contigs[bins]['.retrieved_level'+str(m2)+'.fa']:
                    os.system('rm '+str(bins)+'.retrieved_level'+str(m1)+'.fa')
        os.chdir(pwd)

        # connection_total={}
        # for assembly in assemblies_list:
        #     for line in open('Connections_'+str(assembly)+'.txt', 'r'):
        #         contig=str(line).strip().split('\t')[0]
        #         connection_total[contig]={}
        #         dict_key_list=str(line).strip().split('\t')[1].replace('{','').replace('}','').split(',')
        #         for item in dict_key_list:
        #             connection_total[contig][item.split('\'')[1]]=item.split('\'')[1].strip()

        # print('Recalculation of connections file of retrieved bins')
        # pool=Pool(processes=num_threads)
        # result, new_bin_connections={}, {}
        # for bin in bin_extract_contig_TNF.keys():
        #     result[bins]=pool.apply_async(recal_PE,args=(bin_extract_contig_TNF, connection_total, bin_contig, bin))
        # pool.close()
        # pool.join()

        # for bins in result:
        #     new_bin_connections[bins]=result[bins].get()
        
        # os.chdir(binset+'_retrieved')
        # new_bin_connections2, new_bin_connections3={}, {}
        # f=open('Retrieved_bins_connections.txt','w')
        # for item in new_bin_connections.keys():
        #     bin_name_list=str(item).split('.')
        #     bin_name_list.remove(bin_name_list[-1])
        #     bin_name='.'.join(bin_name_list)
        #     bin_name_new=bin_name+'.retrieved.fa'
        #     bin_name_new2=bin_name+'.retrieved'
        #     f.write(str(bin_name_new)+'\t'+str(new_bin_connections[item])+'\n')
        #     new_bin_connections2[str(bin_name_new)]=str(new_bin_connections[item])
        #     new_bin_connections3[str(bin_name_new2)]=str(new_bin_connections[item])
        # f.close()

        f_cp=open('S6_checkpoint.txt','a')
        f_cp.write('5th bins retrieved done'+'\n')
        f_cp.close()

    if int(last_step) < 6:
        os.chdir(pwd)
        print('Checking quality of retrieved bins')
        # checkm(str(binset)+'_retrieved', num_threads)
        # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset)+'_retrieved '+str(binset)+'_retrieved_checkm')
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(binset)+'_retrieved -x fa -o '+str(binset)+'_retrieved_checkm')
        f_cp=open('S6_checkpoint.txt','a')
        f_cp.write('6th bins quality check done'+'\n')
        f_cp.close()

    if int(last_step) < 7:
        print('Comparing bins with retrieved bins')
        os.chdir(str(binset)+'_retrieved_checkm/')
        print('Parsing '+str(binset)+' checkm output')
    
        # try:
        bins_checkm, n = {}, 0
        try:
            for line in open('quality_report.tsv','r'):
                n+=1
                if n >= 2:
                    binID=str(line).strip().split('\t')[0].strip()
                    genome_size=str(line).strip().split('\t')[8].strip()
                    completeness=str(line).strip().split('\t')[1].strip()
                    contamination=str(line).strip().split('\t')[2].strip()
                    N50=str(line).strip().split('\t')[6].strip()

                    bins_checkm[str(binID)]={}
                    bins_checkm[str(binID)]['N50']=int(N50)
                    bins_checkm[str(binID)]['Completeness']=float(completeness)
                    bins_checkm[str(binID)]['Genome size']=int(genome_size)
                    bins_checkm[str(binID)]['Contamination']=float(contamination)
        except:
            print('There is no quality_report.tsv file')
        os.chdir(pwd)
        bin_comparison(str(binset), bins_checkm, str(binset)+'_retrieved', refinement_mode, num_threads)
        # except:
        #     print('Contig retrieve done! There is not retrieved bin')
        #     os.chdir(pwd)
        #     os.system('rm -rf '+str(binset)+'_retrieved_checkm')

    os.chdir(pwd)
    os.system('mv *_filtrated_bin_connecting_contigs.txt *_eliminated_bin_connecting_contigs.txt Bin_connecting_contigs* '+binset+'_retrieved')
    os.system('rm -rf *_deep_retrieval *_deep_retrieval_checkm')
    os.system('mkdir '+binset+'_retrieved_kmer')
    os.system('mv *a_kmer.txt '+binset+'_retrieved_kmer')

def binset_filtration(binset, pwd, cpn_cutoff, ctn_cutoff):
    filtrated_checkm={}
    os.chdir(pwd+'/'+str(binset))
    
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        for file in files:
            if 'quality_report.tsv' in file:
                n=0
                for line in open('quality_report.tsv','r'):
                    n+=1
                    if n >= 2:
                        binID=str(line).strip().split('\t')[0].strip()
                        completeness=float(str(line).strip().split('\t')[2].strip())
                        contamination=float(str(line).strip().split('\t')[3].strip())
                        # N50=float(str(line).strip().split('\t')[4].strip())
                        # genome_size=int(str(line).strip().split('\t')[1].strip())

                        if float(completeness) >= int(cpn_cutoff) and float(contamination) < int(ctn_cutoff):
                            filtrated_checkm[str(binID)]=str(line)
    
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        for file in files:
            if '_remove_com_outlier.fa' in str(file):
                file1=str(file).split('_remove_com_outlier.fa')[0]
            else:
                file1=str(file)
            file_name_list=str(file1).split('.')
            hz='.'+str(file_name_list[-1])
            file_name_list.remove(file_name_list[-1])
            file_name='.'.join(file_name_list)
            if file_name in filtrated_checkm.keys():
                os.system('cp '+file+' '+pwd+'/'+str(binset)+'_filtrated/'+str(file_name)+str(hz))

    os.chdir(pwd+'/'+str(binset)+'_filtrated')
    f=open('quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completenes'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    for item in filtrated_checkm.keys():
        f.write(str(filtrated_checkm[item]))
    f.close()
    os.chdir(pwd)

def finding_black_contigs(binset, outlier_remover_folder, pwd):
    black_contigs={}
    try:
        n=0
        for line in open(pwd+'/'+outlier_remover_folder+'/Contig_list.txt','r'):
            n+=1
            if n >= 2:
                bins_name=str(line).strip().split('\t')[0].strip()
                contig_id=str(line).strip().split('\t')[1].strip()
                stats=str(line).strip().split('\t')[2].strip()
                if 'Negative' in stats:
                    try:
                        black_contigs[bins_name][contig_id]=''
                    except:
                        black_contigs[bins_name]={}
                        black_contigs[bins_name][contig_id]=''
    except:
        xyz=0

    black_contigs2={}
    for item in black_contigs.keys():
        black_contigs2[item]=[]
        for item2 in black_contigs[item].keys():
            black_contigs2[item].append(item2)
    return black_contigs2

def Contig_recruiter_main(binset, outlier_remover_folder, num_threads, parameter, cpn_cutoff, ctn_cutoff, assemblies_list, PE_connections_list, lr_connection_list, coverage_matrix_list, refinement_mode, pwd):
    
    print('--------------------------------------')
    print('Processing contigs retrieving process')
    print('Binset: '+str(binset))
    print('Binset after outlier removal: '+str(outlier_remover_folder))
    print('The minimun completeness to keep: '+str(cpn_cutoff))
    print('The maximal contaminaition to keep: '+str(ctn_cutoff))
    print('Assemblies: '+str(assemblies_list))
    print('Connections: '+str(PE_connections_list))
    print('Long-reads connections: '+str(lr_connection_list))
    print('Coverage matrix: '+str(coverage_matrix_list))
    print('.....................................')

    try:
        os.mkdir(binset+'_filtrated')
    except:
        print(binset+'_filtrated existed. Deleted and recreated')
        os.system('rm -rf '+binset+'_filtrated')
        os.mkdir(binset+'_filtrated')
    binset_filtration(binset, pwd, cpn_cutoff, ctn_cutoff)

    last_step=0
    if parameter == 'last':
        try:
            n=0
            for line in open('S6_checkpoint.txt', 'r'):
                n+=1

            n1=0
            for line in open('S6_checkpoint.txt', 'r'):
                n1+=1
                if n1 == n:
                    last_step=int(str(line)[0])
        except:
            last_step=0
            f_cp=open('S6_checkpoint.txt', 'w')
            f_cp.close()
    else:
        last_step=0
        f_cp=open('S6_checkpoint.txt', 'w')
        f_cp.close()
    start_step=last_step+1
    print('Start from step '+str(start_step))

    parse_bin_in_bestbinset(assemblies_list, binset+'_filtrated', outlier_remover_folder, PE_connections_list, lr_connection_list, num_threads, last_step, coverage_matrix_list, refinement_mode)
    os.system('rm -rf Deep_retrieved_bins_split')
    os.chdir(pwd)

if __name__ == '__main__': 
    best_binset_from_multi_assemblies='BestBinset_outlier_refined'
    outlier_remover_folder='BestBinset_outlier_refined'
    # assemblies_list=['1_CAMI_cat_opera_contigs_polished.fasta','2_S001_opera_contigs_polished.fasta','3_S002_opera_contigs.fasta',
    # '4_S003_opera_contigs.fasta','5_S004_opera_contigs.fasta','6_S005_opera_contigs.fasta']
    assemblies_list=['1_HumanGut_cat.fasta']
    coverage_matrix_list=['Coverage_matrix_for_binning_1_HumanGut_cat.fasta.txt']
    PE_connections_list=[]
    lr_connection_list='Long_reads_connecting_contigs.txt'  ### e.g. lr_connection_list='Long_reads_connecting_contigs.txt'; if there is not, let it blank, e.g. lr_connection_list=''
    num_threads=60
    refinement_mode='deep' ### 'quick'; 'deep'
    cpn_cutoff=35 ### The minimun completeness to keep
    ctn_cutoff=20 ### The maximal contaminaition to keep
    pwd=os.getcwd()
    parameter='last' ### 'last' or 'none'. checkpoint: using 'last', the script will continue to run the previous incompleted job. 'none' suggests that the job will start from the beginning. 
    Contig_recruiter_main(best_binset_from_multi_assemblies, outlier_remover_folder, num_threads, parameter, cpn_cutoff, ctn_cutoff, assemblies_list, PE_connections_list, lr_connection_list, coverage_matrix_list, refinement_mode, pwd)