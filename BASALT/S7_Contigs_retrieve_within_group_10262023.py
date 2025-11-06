#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio import SeqIO
import os, threading, copy
from sklearn.decomposition import PCA
import numpy as np
from multiprocessing import Pool
# from Outlier_remover import *
from glob import glob

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

def record_bin_coverage(binset, num_threads, assembly_list, coverage_matrix_list, level_num):
    pwd=os.getcwd()
    bin_contigs, contig_bin_list, bin_contigs_mock, total_bin_contigs, m = {}, {}, {}, {}, 0
    print('Parsing bins in retrieval_'+str(level_num))
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        os.chdir(pwd+'/'+str(binset))
        for file in files:
            if '_genomes.' in file:
                hz=file.split('_genomes.')[1]
                # qz=file.split('_genomes.')[0]
                # assembly_name_list=qz.split('_')
                # if len(assembly_name_list) >= 3:
                #     assembly_name_list.remove(assembly_name_list[-1])
                #     assembly_name_list.remove(assembly_name_list[-1])
                #     assembly_name='_'.join(assembly_name_list)
                # else:
                #     assembly_name=assembly_name_list[0]
                # assembly_list[assembly_name]=1
                if '.fasta' in hz or '.fa' in hz:
                    m+=1
                    bin_contigs[file]={}
                    bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        contig_bin_list[record.id]=[]
                        bin_contigs[file][record.id]=record.seq
                        bin_contigs_mock[file][record.id]=0
                    
                    for record in SeqIO.parse(file, 'fasta'):
                        contig_bin_list[record.id].append(file)
    print('Parsed', m, 'bins in', str(binset))

    os.chdir(pwd)
    for assemblies in assembly_list:
        for record in SeqIO.parse(assemblies, 'fasta'):
            total_bin_contigs[record.id]=record.seq

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
        result[item]=pool.apply_async(coverage_matrix_mpt, args=(item, num))
    pool.close()
    pool.join()

    result2={}
    for item in result:
        result2[item]=result[item].get()

    contig_cov={}
    for item in result2.keys():
        contig_cov.update(result2[item][0])

    print('Writing bin coverage matrix file')
    try:
        os.system('mkdir Bin_coverage_after_contamination_removal')
    except:
        print('Bin_coverage_after_contamination_removal presented. Re-created.')
        os.system('rm -rf Bin_coverage_after_contamination_removal')
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
    return bin_contig_cov, bin_contigs, contig_cov, total_bin_contigs

def cycle_mt(bin_connecting_contigs, bin_connecting_contigs2, connections, bins):
    m, m_before, m_after=1, 0, 1
    print('Parsing', bins)
    while m <= 3: ### Consider connection level less than 2
        # print(bins, 'cycle', m)
        if m_before != m_after:
            m_before = len(bin_connecting_contigs)
            for contigs in connections.keys():
                for item in connections[contigs].keys():
                    if contigs in bin_connecting_contigs.keys() and item not in bin_connecting_contigs.keys():
                        bin_connecting_contigs[item]=m
                        bin_connecting_contigs2[item]=m
            m_after = len(bin_connecting_contigs)
            m+=1
        else:
            m=6
    return bin_connecting_contigs, bin_connecting_contigs2

def Parsing_kmer_file(assemblies_list, binset, bin_extract_contig, num_threads):
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
    bin_connecting_contigs, bin_connecting_contigs2, result={}, {}, {}
    for bins in bin_contigs.keys():
        bin_connecting_contigs[bins]={}
        bin_connecting_contigs2[bins]={}
        bin_connecting_contigs[bins]=bin_contigs[bins]
        result[bins]=pool.apply_async(cycle_mt,args=(bin_connecting_contigs[bins], bin_connecting_contigs2[bins], connections, bins))
    pool.close()
    pool.join()

    result2={}
    for bins in result:
        result2[bins]=result[bins].get()

    for bins in result:
        bin_connecting_contigs[bins].update(result2[bins][0])
        bin_connecting_contigs2[bins].update(result2[bins][1])

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

    bin_connecting_contigs_level={}
    for item in bin_connecting_contigs2.keys():
        bin_connecting_contigs_level[item]={}
        for contigs in bin_connecting_contigs2[item].keys():
            level=int(bin_connecting_contigs2[item][contigs])
            try:
                bin_connecting_contigs_level[item][level].append(contigs)
            except:
                bin_connecting_contigs_level[item][level]=[contigs]

    os.chdir(pwd)
    return bin_connecting_contigs2, connections, bin_connecting_contigs_level

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

def checkm(bin_folder, num_threads):
    pwd=os.getcwd()
    # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(bin_folder)+' '+str(bin_folder)+'_checkm')
    os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(bin_folder)+' -x fa -o '+str(bin_folder)+'_checkm')

    print('Parsing '+bin_folder+' checkm output')
    refined_checkm={}
    try:
        os.chdir(str(bin_folder)+'_checkm/')
        n=0
        for line in open('quality_report.tsv','r'):
            n+=1
            if n >= 2:
                binID=str(line).strip().split('\t')[0].strip()
                refined_checkm[str(binID)]={}

                try:
                    genome_size=str(line).strip().split('\t')[8].strip()
                    completeness=str(line).strip().split('\t')[1].strip()
                    contamination=str(line).strip().split('\t')[2].strip()
                    N50=str(line).strip().split('\t')[6].strip()
                except:
                    genome_size=str(line).strip().split('\t')[1].strip()
                    completeness=str(line).strip().split('\t')[2].strip()
                    contamination=str(line).strip().split('\t')[3].strip()
                    N50=str(line).strip().split('\t')[4].strip()

                refined_checkm[str(binID)]['N50']=int(eval(N50))
                refined_checkm[str(binID)]['Completeness']=float(eval(completeness))
                refined_checkm[str(binID)]['Genome size']=int(eval(genome_size))
                refined_checkm[str(binID)]['Contamination']=float(eval(contamination))
    except:
        print('There is not bins found in the retrieved folder')
    os.chdir(pwd)
    return refined_checkm

def bin_comparison(original_bin_folder, new_bins_checkm, new_bin_folder, num_threads, total_contigs, level_num):
    pwd=os.getcwd()
    print('Comparing bins before and after refining process')
    os.chdir(pwd+'/'+str(original_bin_folder))
    bin_checkm, bin_seq_rec, selected_bin = {}, {}, {}
    for root, dirs, files in os.walk(pwd+'/'+str(original_bin_folder)):
        for file in files:
            if 'quality_report.tsv'  in file:
                n=0
                for line in open(file, 'r'):
                    n+=1
                    if n >= 2:
                        binID=str(line).strip().split('\t')[0].strip()
                        bin_checkm[binID]={}
                        selected_bin[binID]=1

                        try:
                            genome_size=str(line).strip().split('\t')[8].strip()
                            completeness=str(line).strip().split('\t')[1].strip()
                            contamination=str(line).strip().split('\t')[2].strip()
                            N50=str(line).strip().split('\t')[6].strip()
                        except:
                            genome_size=str(line).strip().split('\t')[1].strip()
                            completeness=str(line).strip().split('\t')[2].strip()
                            contamination=str(line).strip().split('\t')[3].strip()
                            N50=str(line).strip().split('\t')[4].strip()

                        bin_checkm[binID]['N50']=int(eval(N50))
                        bin_checkm[binID]['Completeness']=float(eval(completeness))
                        bin_checkm[binID]['Genome size']=int(eval(genome_size))
                        bin_checkm[binID]['Contamination']=float(eval(contamination))
            else:
                bin_id_list=str(file).split('.')
                bin_id_list.remove(bin_id_list[-1])
                bin_id='.'.join(bin_id_list)
                if 'fa' in str(file).split('.')[-1]:
                    bin_seq_rec[bin_id]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_seq_rec[bin_id][record.id]=record.seq

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

    # print(str(bin_comparison_list))
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
                re_taxon=str(new_bins_checkm[refined_bin_id]['N50'])
                re_genome_size=float(new_bins_checkm[refined_bin_id]['Genome size'])
                re_delta=re_cpn-re_ctn
                re_5delta=re_cpn-5*re_ctn
                for bin_id2 in bestbin[item].keys():
                    ori_bin=bin_id2 
                    # ori_connections=int(bestbin[item][bin_id2]['Connections'])
                    ori_cpn=float(bestbin[item][bin_id2]['Completeness'])
                    ori_ctn=float(bestbin[item][bin_id2]['Contamination'])
                    ori_taxon=str(bestbin[item][bin_id2]['N50'])
                    ori_genome_size=float(bestbin[item][bin_id2]['Genome size'])
                    ori_delta=ori_cpn-ori_ctn
                    ori_5delta=ori_cpn-5*ori_ctn
                    
                if re_5delta > ori_5delta:
                # if re_delta > ori_delta:
                    del bestbin[item][ori_bin]
                    bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id] 
                elif re_5delta == ori_5delta:
                # elif re_delta == ori_delta:
                    if '.retrieved' in ori_bin:
                        cutoff1=float(ori_bin.split('.retrieved_level')[1].split('_')[1].split('.fa')[0])
                        if '.retrieved' in refined_bin_id:
                            cutoff2=float(refined_bin_id.split('.retrieved_level')[1].split('_')[1].split('.fa')[0])
                            if cutoff2 == 1:
                                del bestbin[item][ori_bin]
                                bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]
                        else:
                            if cutoff1 != 1:
                                del bestbin[item][ori_bin]
                                bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]
                    else:
                        if '.retrieved' in refined_bin_id:
                            cutoff2=float(refined_bin_id.split('.retrieved_level')[1].split('_')[1].split('.fa')[0])
                            if cutoff2 == 1:
                                del bestbin[item][ori_bin]
                                bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]

                    # if '.retrieved' in refined_bin_id:
                    #     cutoff2=float(refined_bin_id.split('.retrieved_level')[1].split('_')[1].split('.fa')[0])
                    #     if '.retrieved' not in ori_bin:
                    #         if cutoff2 == 1:
                    #             del bestbin[item][ori_bin]
                    #             bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]
                    #     else:
                    #         cutoff1=float(ori_bin.split('.retrieved_level')[1].split('_')[1].split('.fa')[0])
                    #         if cutoff1 > cutoff2:
                    #             del bestbin[item][ori_bin]
                    #             bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]
                    # else: ### refined bin is not the original bin
                    #     if '.retrieved' in ori_bin:
                    #         cutoff1=float(ori_bin.split('.retrieved_level')[1].split('_')[1].split('.fa')[0])
                    #         if cutoff1 != 1:
                    #             del bestbin[item][ori_bin]
                    #             bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id]


                # if re_delta > ori_delta:
                #     del bestbin[item][ori_bin]
                #     bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id] 
                #     if re_5delta < ori_5delta and ori_cpn >= 40 and ori_ctn <= 15 and re_ctn <= 50:
                #         further_refined_bin[refined_bin_id]=item
                #         f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                #     elif re_cpn < ori_cpn and ori_cpn >= 40 and ori_ctn <= 15 and re_ctn <= 50:
                #         further_refined_bin[refined_bin_id]=item
                #         f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                # elif re_cpn > ori_cpn:
                #     if ori_cpn >= 40 and ori_ctn <= 10 and re_ctn <= 50:
                #         further_refined_bin[refined_bin_id]=item
                #         f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                # else:
                #     continue
        f.write(str(write_out)+'\n')
    f.close()
    f2.close()

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
    total_drf_num=n

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
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                ids_list=file.split('.')
                ids_list.remove(ids_list[-1])
                bin_id='.'.join(ids_list)
                # if bin_id in bestbin.keys():
                if bin_id in selected_bin.keys():
                    n+=1
                    os.system('cp '+file+' '+pwd+'/'+str(new_bin_folder))
    print('Moved', n, 'bins to refined bins folder')

    os.chdir(pwd)
    print('Processing deep refinement')
    try:
        os.system('mkdir Deep_retrieved_bins_'+str(level_num))
    except:
        print('Deep_retrieved_bins folder existed')
        os.system('rm -rf Deep_retrieved_bins_'+str(level_num))
        os.system('mkdir Deep_retrieved_bins_'+str(level_num))

    os.chdir(pwd+'/'+str(new_bin_folder))
    if os.path.exists('Bin_contigs_white_black_list.txt'):
        fxxx=open('Bin_contigs_white_black_list.txt','a')
    else:
        fxxx=open('Bin_contigs_white_black_list.txt','w')
    os.chdir(pwd)

    further_refined_bin_contig={}
    m=0
    for org_bin_id in further_refined_bin_contig.keys():
        m+=1
        pct=str(round(100*m/total_drf_num,2))+'%'
        print('Processing bin-'+str(m)+' '+str(org_bin_id)+', accounting of '+str(pct))
        try:
            os.system('mkdir '+org_bin_id+'_deep_retrieval')
        except:
            os.system('rm -rf '+org_bin_id+'_deep_retrieval')
            os.system('mkdir '+org_bin_id+'_deep_retrieval')
        os.chdir(org_bin_id+'_deep_retrieval')
        f=open(org_bin_id+'.fa','w')
        for contigs_id in bin_seq_rec[org_bin_id].keys():
            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
        f.close()

        contig_neutral, contig_positive, contig_negative, contig_positive_suspect, contig_negative_suspect, contig_num, n = {}, {}, {}, {}, {}, {}, 0
        for contigs in further_refined_bin_contig[org_bin_id].keys():
            if org_bin_id in total_contigs.keys():
                if contigs not in total_contigs[org_bin_id]:
                    n+=1
                    f=open(org_bin_id+'_deep_'+str(n)+'.fa','w')
                    f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                    contig_num[n]=str(contigs)
                    for contigs_id in bin_seq_rec[org_bin_id].keys():
                        f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                    f.close()
            else:
                n+=1
                f=open(org_bin_id+'_deep_'+str(n)+'.fa','w')
                f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                contig_num[n]=str(contigs)
                for contigs_id in bin_seq_rec[org_bin_id].keys():
                    f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                f.close()

        os.chdir(pwd)
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(org_bin_id)+'_deep_retrieval -x fa -o '+str(org_bin_id)+'_deep_retrieval_checkm')
        os.chdir(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm')
        n=0
        for line in open('quality_report.tsv', 'r'):
            n+=1
            if n >= 2:
                bin_id=str(line).strip().split('\t')[0].strip()
                if bin_id == org_bin_id:
                    org_bin_cpn=float(str(line).strip().split('\t')[1].strip())
                    org_bin_ctn=float(str(line).strip().split('\t')[2].strip())
                    ori_delta=org_bin_cpn-5*org_bin_ctn
                    ori_1delta=org_bin_cpn-org_bin_ctn
            n=0
            for line in open('quality_report.tsv', 'r'):
                n+=1
                if n >= 2:
                    refined_bin_id=str(line).strip().split('\t')[0].strip()
                    refined_bin_cpn=float(str(line).strip().split('\t')[1].strip())
                    refined_bin_ctn=float(str(line).strip().split('\t')[2].strip())
                    refined_delta=refined_bin_cpn-5*refined_bin_ctn
                    refined_1delta=refined_bin_cpn-refined_bin_ctn

                    if refined_bin_id != org_bin_id:
                        if refined_delta > ori_delta:
                            contigs=contig_num[int(refined_bin_id.split('_deep_')[1].split('.fa')[0])]
                            contig_positive[str(contigs)]=0
                            fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Positive'+'\n')
                        elif refined_delta == ori_delta:
                            contigs=contig_num[int(refined_bin_id.split('_deep_')[1].split('.fa')[0])]
                            contig_neutral[str(contigs)]=0
                            fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Neutral'+'\n')
                        elif refined_delta < ori_delta:
                            if refined_1delta >= ori_1delta:
                                contig_positive_suspect[str(contigs)]=0
                                fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Positive suspect'+'\n')
                            elif refined_bin_cpn > org_bin_cpn:
                                contig_negative_suspect[str(contigs)]=0
                                fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Negative suspect'+'\n')
                        else:
                            contig_negative[str(contigs)]=0
                            fxxx.write(str(org_bin_id)+'\t'+str(contigs)+'\t'+'Negative'+'\n')
            os.system(pwd)

        xxxx=0
        if len(contig_positive) != 0:
            os.chdir(pwd+'/Deep_retrieved_bins_'+str(level_num))
            f=open(org_bin_id+'_RT-D3.fa','w')
            for contigs in contig_positive.keys():
                f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
            for contigs_id in bin_seq_rec[org_bin_id].keys():
                f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
            f.close()
            xxxx=1

        if len(contig_positive_suspect) != 0:
            os.chdir(pwd+'/Deep_retrieved_bins_'+str(level_num))
            f=open(org_bin_id+'_RT-D4.fa','w')
            for contigs in contig_positive.keys():
                f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
            for contigs in contig_positive_suspect.keys():
                f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

            for contigs_id in bin_seq_rec[org_bin_id].keys():
                f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
            f.close()
            xxxx=1
        os.chdir(pwd)
    os.chdir(pwd)

    xxxx=0
    for root, dirs, files in os.walk(pwd+'/Deep_retrieved_bins_'+str(level_num)):
        for file in files:
            xxxx+=1

    replace_bin = {}
    if xxxx>=1:
        os.system('checkm2 predict -t '+str(num_threads)+' -i Deep_retrieved_bins_'+str(level_num)+' -x fa -o Deep_retrieved_bins_'+str(level_num)+'_checkm')
        os.chdir(pwd+'/Deep_retrieved_bins_checkm/')
        n=0
        for line in open('quality_report.tsv', 'r'):
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
                    new_bins_checkm[refined_bin_id]['Genome size']=float(str(line).strip().split('\t')[8].strip())
                    new_bins_checkm[refined_bin_id]['N50']=float(str(line).strip().split('\t')[6].strip())
                    bestbin[org_bin_id][refined_bin_id]=new_bins_checkm[refined_bin_id]
                    os.system('cp '+pwd+'/Deep_retrieved_bins_'+str(level_num)+'/'+refined_bin_id+'.fa '+pwd+'/'+str(new_bin_folder))
                elif refined_delta == org_delta:
                # else:
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
                                new_bins_checkm[refined_bin_id]['Genome size']=float(str(line).strip().split('\t')[8].strip())
                                new_bins_checkm[refined_bin_id]['N50']=float(str(line).strip().split('\t')[6].strip())
                                bestbin[org_bin_id][refined_bin_id]=new_bins_checkm[refined_bin_id]
                                os.system('cp '+pwd+'/Deep_retrieved_bins_'+str(level_num)+'/'+refined_bin_id+'.fa '+pwd+'/'+str(new_bin_folder))
                                
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

    f_bin_checkm=open(new_bin_folder+'_retrieved_quality_report.tsv', 'w')
    f_bin_checkm.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    for item in selected_bin.keys():
        if '_RT-D' in item:
            item_name=str(item).split('_RT-D')[0]
            # f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item])+'\n')
            f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item]['Genome size'])+'\t'+str(total_bin_checkm[item]['Completeness'])+'\t'+str(total_bin_checkm[item]['Contamination'])+'\t'+str(total_bin_checkm[item]['N50'])+'\n')
        elif '.retrieved_level' in item:
            item_name=str(item).split('.retrieved_level')[0]
            # f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item])+'\n')
            f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item]['Genome size'])+'\t'+str(total_bin_checkm[item]['Completeness'])+'\t'+str(total_bin_checkm[item]['Contamination'])+'\t'+str(total_bin_checkm[item]['N50'])+'\n')
        else:
            # f_bin_checkm.write(str(item)+'\t'+str(total_bin_checkm[item])+'\n')
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
                if '_RT-D' in file or '.retrieved_level' in file:
                    if '_RT-D' in file:
                        file_name=str(file).split('_RT-D')[0]
                    if '.retrieved_level' in file:
                        file_name=str(file).split('.retrieved_level')[0]
                    # if 'metabat' in file_name or 'metabinner' in file_name or 'vamb' in file_name:
                    #     os.system('mv '+file+' '+str(file_name)+'.fa')
                    #     f.write(str(file)+'\t'+str(file_name)+'.fa'+'\n')
                    else:
                        os.system('mv '+file+' '+str(file_name)+'.fa')
                        f.write(str(file)+'\t'+str(file_name)+'.fa'+'\n')
    f.close()
    os.chdir(pwd)   
    # os.system('rm -rf *_deep_retrieval_checkm *_deep_retrieval')
    print('Contig retrieve '+str(level_num)+' done!')

def coverage_filtration(bin, m, level_num, cutoff):
    os.system('S6p_coverage_filtration_mpt_06102022.py -b '+str(bin)+' -n '+str(m)+' -p coverage_filtration -c '+str(cutoff)+' -l '+str(level_num))

def TNF_filtration(bin, m, level_num, cutoff):
    os.system('S6p_coverage_filtration_mpt_06102022.py -b '+str(bin)+' -n '+str(m)+' -p TNF_filtration -c '+str(cutoff)+' -l '+str(level_num))

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

def combat_within_group(bestbinset, assembly_list, pwd, cpn_cutoff, ctn_cutoff):
    # bestbinset='BestBinset'
    # assembly_list=['16_high_cat_spades.fasta','17_high_S001_spades.fasta','18_high_S002_spades.fasta','19_high_S003_spades.fasta','20_high_S004_spades.fasta','21_high_S005_spades.fasta']
    # # assembly_list=['16_high_cat_spades.fasta']
    # pwd=os.getcwd()
    bin_contigs, bin_contigs_mock, bin_seq, bin_select_contigs={}, {}, {}, {}
    os.chdir(pwd+'/'+bestbinset)
    for root, dirs, files in os.walk(pwd+'/'+bestbinset):
            # os.chdir(pwd+'/'+str(binset_1_assembly))
        for file in files:
            if '_genomes.' in file:
                hz=file.split('.')[-1]
                qz=file.split('_genomes.')[0]
                if 'fa' in hz or 'fna' in hz:
                        # print('Parsing', file
                    bin_contigs[file]={}
                    # bin_select_contigs[file]={}
                    # bin_contigs_mock[file]={}
                    for record in SeqIO.parse(file, 'fasta'):
                        bin_contigs[file][str(record.id)]=0
                        # bin_select_contigs[file][str(record.id)]=0
                        # bin_contigs_mock[file][str(record.id)]=0
                        # bin_seq[str(record.id)]=str(record.seq)
    os.chdir(pwd)

    print('Trying to read combat bins from the same assembly')
    paired, mirror_bins, n={}, {}, 0
    try:
        for assembly in assembly_list:
            paired[assembly]={}
            # print(str(assembly))
            os.chdir(pwd+'/'+str(assembly)+'_comparison_files')
            for root, dirs, files in os.walk(pwd+'/'+str(assembly)+'_comparison_files'):
                for file in files:
                    if 'Best_bin_set_iteration_' in file:
                        # ids=int(str(file).split('.')[0].split('Best_bin_set_iteration_')[1])
                        # print(str(ids))
                        # paired[assembly][ids]={}
                        # print(str(file))
                        n1=0
                        for line in open(file, 'r'):
                            n1+=1
                            if '---' in line:
                                bin1=str(line).strip().split('\t')[1].split('---')[0].strip()
                                bin2=str(line).strip().split('\t')[1].split('---')[1].strip()
                                bin1_checkm=str(line).strip().split('\t')[4]
                                bin1_checkm_cpn=float(bin1_checkm.split('Completeness\':')[1].split(',')[0].strip())
                                bin1_checkm_ctn=float(bin1_checkm.split('Contamination\':')[1].split('}')[0].split(',')[0].strip())
                                bin2_checkm=str(line).strip().split('\t')[5]
                                bin2_checkm_cpn=float(bin2_checkm.split('Completeness\':')[1].split(',')[0].strip())
                                bin2_checkm_ctn=float(bin2_checkm.split('Contamination\':')[1].split('}')[0].split(',')[0].strip())
                                if bin1_checkm_cpn >= float(cpn_cutoff) and bin1_checkm_ctn <= float(ctn_cutoff) and bin2_checkm_cpn >= float(cpn_cutoff) and bin2_checkm_ctn <= float(ctn_cutoff):                          
                                    if bin1 not in paired[assembly]:
                                        paired[assembly][bin1]=[]
                                    paired[assembly][bin1].append(bin2)

                                    if bin2 not in paired[assembly]:
                                        paired[assembly][bin2]=[]
                                    paired[assembly][bin2].append(bin1)
            os.chdir(pwd)
    except:
        print('Within group comparison not appliable')

    bin_withingroup={}
    for bins in bin_contigs.keys():
        bin_withingroup[bins]={}

    for assembly in assembly_list:
        if len(paired[assembly]) != 0:
            for bins in bin_contigs.keys():
                if bins in paired[assembly].keys():
                    for bins2 in paired[assembly][bins]:
                        n+=1
                        bin_withingroup[bins][bins2]=''
                        if bins2 in paired[assembly]:
                            for bins3 in paired[assembly][bins2]:
                                n+=1
                                bin_withingroup[bins][bins3]=''
                                if bins3 in paired[assembly]:
                                    for bins4 in paired[assembly][bins3]:
                                        n+=1
                                        bin_withingroup[bins][bins4]=''

    print(len(bin_withingroup))
    for item in bin_withingroup.keys():
        n=0
        for item2 in bin_withingroup[item].keys():
            n+=1
            bin_withingroup[item][item2]=n
    # os.chdir(pwd)

    bin_withingroup2=copy.deepcopy(bin_withingroup)
    f=open('Combat_within_group.txt','w')
    for bins in bin_withingroup2.keys():
        if len(bin_withingroup2[bins]) != 0:
            f.write(str(bins)+'\t'+str(bin_withingroup2[bins])+'\n')
        else:
            del bin_withingroup[bins]
    f.close()                        
    print(len(bin_withingroup))

    max_iteration=0
    for bins in bin_withingroup.keys():
        for bins2 in bin_withingroup[bins].keys():
            if int(bin_withingroup[bins][bins2]) > max_iteration:
                max_iteration = int(bin_withingroup[bins][bins2])

    bin_withingroup_level={}
    for i in range(1, max_iteration+1):
        bin_withingroup_level[i]={}

    # extra_contigs={}
    for t_bin in bin_withingroup.keys():
        # extra_contigs[t_bin]={}
        m=0
        for e_bin in bin_withingroup[t_bin].keys():
            if e_bin != t_bin:
                gennome_folder_list=str(e_bin).split('.')
                gennome_folder_list.pop()
                gennome_folder_list.pop()
                gennome_folder='.'.join(gennome_folder_list)
                m+=1
                e_bin_contigs={}
                os.chdir(pwd+'/'+str(gennome_folder))
                level_num=bin_withingroup[t_bin][e_bin]
                # for record in SeqIO.parse(e_bin,'fasta'): 
                #     e_bin_contigs[record.id]=level_num
                if t_bin not in bin_withingroup_level[level_num].keys():
                    bin_withingroup_level[level_num][t_bin]={}
                bin_withingroup_level[level_num][t_bin][e_bin]=''
                    # bin_withingroup_level[level_num][t_bin]={}
                    # bin_withingroup_level[level_num][t_bin][record.id]=''
                os.chdir(pwd)

                # for item in e_bin_contigs.keys():
                #     if item not in bin_contigs[t_bin].keys():
                #         extra_contigs[t_bin][item]=level_num

    # f=open('Combat_extract_contigs.txt','w')
    # for item in extra_contigs.keys():
    #     f.write(str(item)+'\t'+str(extra_contigs[item])+'\n')
    # f.close()

    f=open('Combat_extract_contigs_level.txt','w')
    max_iteration2 = 0
    for item in bin_withingroup_level.keys():
        for t_bin in bin_withingroup_level[item].keys():
            if len(bin_withingroup_level[item][t_bin]) == 0:
                del bin_withingroup_level[item][t_bin]

    bin_withingroup_level2=copy.deepcopy(bin_withingroup_level)
    for item in bin_withingroup_level2.keys():
        if len(bin_withingroup_level2[item]) == 0:
            del bin_withingroup_level[item]

    n, comfired_level = 0, {}
    for item in bin_withingroup_level.keys():
        if item not in comfired_level.keys():
            n+=1
            comfired_level[item]=n

        for t_bin in bin_withingroup_level[item].keys():
            f.write(str(comfired_level[item])+'\t'+str(t_bin)+'\t'+str(bin_withingroup_level[item][t_bin])+'\n')
            if int(comfired_level[item]) > max_iteration2:
                max_iteration2=int(comfired_level[item])
    f.close()
    return max_iteration2

def Contig_retrieve_within_group(assemblies_list, binset, outlier_remover_folder, PE_connections_list, num_threads, last_step, coverage_matrix_list, pwd, cpn_cutoff, ctn_cutoff):
    # pwd=os.getcwd()
    assemblies={}
    for item in assemblies_list:
        assemblies[item]=[]

    ### Record bins from source
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        os.chdir(pwd+'/'+str(binset))
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if '.fa' in hz:
                    for item in assemblies.keys():
                        if item in file:
                            assemblies[item].append(file)
    os.chdir(pwd)

    elemimated_contig_total  = {}
    fxy=open('Combat_extract_contigs_level2.txt','w')
    fxy.close()

    if int(last_step) < 1:
        ### Finding combats
        max_iteration=combat_within_group(binset, assemblies_list, pwd, cpn_cutoff, ctn_cutoff)
        # f_cp=open('Retrieve_from_contigs_within_group.txt','a')
        # f_cp.write('1st Bin_connecting done'+'\n')
        # f_cp.close()

        os.chdir(pwd)
        bin_connecting_contigs_total_level = {}
        for line in open('Combat_extract_contigs_level.txt','r'):
            level_n=int(str(line).strip().split('\t')[0].strip())
            bin_connecting_contigs_total_level[level_n]={}
        
        print('Total '+str(max_iteration)+' iterations')
        print('Iteration starts')
        try:
            f_checkpoit=open('S7_contigs_retriev_within_group_checkpoint.txt','a')
        except:
            f_checkpoit=open('S7_contigs_retriev_within_group_checkpoint.txt','w')
        f_checkpoit.close()

        latest_iteration, xyzzz = 0, 0
        for line in open('S7_contigs_retriev_within_group_checkpoint.txt','r'):
            xyzzz+=1

        xyzzz2=0
        for line in open('S7_contigs_retriev_within_group_checkpoint.txt','r'):
            xyzzz2+=1
            if xyzzz2 == xyzzz:
                try:
                    latest_iteration=int(str(line).strip().split(' ')[0])
                except:
                    latest_iteration=0
                    
        # comparied_level, latest_iteration, comparied_level2={},0 ,{}
        for level_num in range(1, int(max_iteration)+1):
            if level_num > latest_iteration:
                if level_num == 1:
                    print('Start form 1st iteration')
                else:
                    print('Continue with iteration '+str(level_num))

                if level_num in bin_connecting_contigs_total_level.keys():
                    # actual_level_num+=1
                    # comparied_level[actual_level_num]=level_num
                    # comparied_level2[level_num]=actual_level_num

                    if level_num == 1:
                        A=record_bin_coverage(binset, num_threads, assemblies_list, coverage_matrix_list, level_num)
                        binset2=binset
                    else:
                        # level_num2=level_num-1
                        A=record_bin_coverage(binset+'_retrieved_'+str(latest_iteration), num_threads, assemblies_list, coverage_matrix_list, level_num)
                        binset2=binset+'_retrieved_'+str(latest_iteration)
                        os.system('cp '+pwd+'/'+binset+'/Bin_contigs_white_black_list.txt '+pwd+'/'+binset2+'/Bin_contigs_white_black_list.txt')
                    bin_contig_cov=A[0]
                    bin_contig=A[1]
                    contig_cov=A[2]
                    total_bin_contigs=A[3]

                    black_contigs, grey_contigs, total_contigs=finding_black_contigs(binset, binset2, outlier_remover_folder, pwd, level_num)
                    
                    os.chdir(pwd)
                    
                    fxy=open('Combat_extract_contigs_level2.txt','a')
                    fxyzz=open('Record_error_bin.txt','w')
                    for line in open('Combat_extract_contigs_level.txt','r'):
                        level_n=int(str(line).strip().split('\t')[0].strip())
                        if level_n == level_num:
                            target_bin=str(line).strip().split('\t')[1].strip()
                            bin_connecting_contigs_total_level[level_n][target_bin]={}
                            bin_name_list=target_bin.split('.')
                            bin_name_list.pop()
                            bin_name='.'.join(bin_name_list) 
                            ebin_list=str(line).strip().split('\t')[2].strip().replace('{','').replace('}','').replace('\'','').replace(':','').strip().split(',')
                            for ebin in ebin_list:
                                gennome_folder_list=str(ebin).strip().split('.')
                                gennome_folder_list.pop()
                                gennome_folder_list.pop()
                                gennome_folder='.'.join(gennome_folder_list)
                                os.chdir(pwd+'/'+str(gennome_folder))
                                try:
                                    for record in SeqIO.parse(ebin,'fasta'):
                                        if record.id not in bin_contig_cov[target_bin].keys():
                                            bin_connecting_contigs_total_level[level_n][target_bin][record.id]=''
                                except:
                                    fxyzz.write(str(ebin)+' dose not exist in folder '+pwd+'/'+str(gennome_folder)+'\n')
                                os.chdir(pwd)
                                if len(bin_connecting_contigs_total_level[level_n][target_bin]) != 0:
                                    fxy.write(str(level_n)+'\t'+str(target_bin)+'\t'+str(bin_connecting_contigs_total_level[level_n][target_bin])+'\n')
                                else:
                                    del bin_connecting_contigs_total_level[level_n][target_bin]
                    fxy.close()
                    fxyzz.close()

                    # print('Iteration-'+str(level_num))
                    thresholds=[1, 1.2, 1.5]
                    ##########
                    for cutoff in thresholds:
                        print('Processing cutoff '+str(cutoff)+' at level '+str(level_num))
                        bin_connecting_contigs_total, bin_connecting_contigs_total_x = {}, {}
                        bin_connecting_contigs_total=copy.deepcopy(bin_connecting_contigs_total_level[level_num])
                        bin_connecting_contigs_total_x=copy.deepcopy(bin_connecting_contigs_total_level[level_num])
                        # print(bin_connecting_contigs_total)
                        for bins in bin_connecting_contigs_total_x.keys():
                            for contigs in bin_connecting_contigs_total_x[bins].keys():
                                if contigs in bin_contig_cov[bins].keys():
                                    del bin_connecting_contigs_total[bins][contigs]
                        
                        bin_connecting_contigs_total_x={}
                        bin_connecting_contigs_total_x=copy.deepcopy(bin_connecting_contigs_total)

                        for bins in bin_connecting_contigs_total_x.keys():
                            if len(bin_connecting_contigs_total_x[bins]) == 0:
                                del bin_connecting_contigs_total[bins]

                        fct=open('bin_connecting_contigs_total_'+str(level_num)+'.txt','w')
                        for item in bin_connecting_contigs_total.keys():
                            fct.write(str(item)+'\t'+str(bin_connecting_contigs_total[item])+'\n')
                        fct.close()

                        pool=Pool(processes=num_threads)
                        bin_extract_contig, selected_1st, elemimated_contig, m= {}, {}, {}, 0
                        selected_files, extract_files, elemimated_files = [], [], []
                        print(str(len(bin_connecting_contigs_total)), 'bins are waiting for coverage filtration')
                        print('-----------------------------------------------')
                        for bin in bin_connecting_contigs_total.keys():
                            m+=1
                            pool.apply_async(coverage_filtration,args=(bin, m, level_num, cutoff))
                            selected_files.append(str(bin)+'_'+str(level_num)+'_selected_contigs.txt')
                            extract_files.append(str(bin)+'_'+str(level_num)+'_bin_extract_contig.txt')
                            elemimated_files.append(str(bin)+'_'+str(level_num)+'_elemimated_contig.txt')
                        pool.close()
                        pool.join()

                        for bin in bin_connecting_contigs_total.keys():
                            os.system('mv '+str(bin)+'_selected_contigs.txt '+str(bin)+'_'+str(level_num)+'_selected_contigs.txt')
                            os.system('mv '+str(bin)+'_bin_extract_contig.txt '+str(bin)+'_'+str(level_num)+'_bin_extract_contig.txt')
                            os.system('mv '+str(bin)+'_elemimated_contig.txt '+str(bin)+'_'+str(level_num)+'_elemimated_contig.txt')

                        try:
                            os.system('mkdir coverage_filtration_matrix_'+str(level_num))
                        except:
                            print('Folder coverage_filtration_matrix_'+str(level_num)+' existed. Deleteced and recreated.')
                            os.system('rm -rf coverage_filtration_matrix_'+str(level_num))
                            os.system('mkdir coverage_filtration_matrix_'+str(level_num))

                        for item in selected_files:
                            selected_1st.update(parse_dict(item))
                            os.system('mv '+str(item)+' coverage_filtration_matrix_'+str(level_num))

                        for item in extract_files:
                            bin_extract_contig.update(parse_dict(item))
                            os.system('mv '+str(item)+' coverage_filtration_matrix_'+str(level_num))

                        for item in elemimated_files:
                            elemimated_contig.update(parse_dict(item))
                            os.system('mv '+str(item)+' coverage_filtration_matrix_'+str(level_num))
                            
                        elemimated_contig_total.update(elemimated_contig)

                        print(str(level_num)+' Coverage filtration done!')
                        f=open('Coverage_'+str(level_num)+'_filtrated_bin_connecting_contigs.txt','w')
                        for bin in selected_1st.keys():
                            f.write(str(bin)+'\t'+str(selected_1st[bin])+'\n')
                        f.close()

                        f=open('Bin_extract_contigs_after_coverage_filtration_'+str(level_num)+'.txt','w')
                        for bin in bin_extract_contig.keys():
                            f.write(str(bin)+'\t'+str(bin_extract_contig[bin])+'\n')
                        f.close()

                        f=open('Coverage_eliminated_bin_connecting_contigs_'+str(level_num)+'.txt','w')
                        for bin in elemimated_contig.keys():
                            f.write(str(bin)+'\t'+str(elemimated_contig[bin])+'\n')
                        f.close()

                        selected_1st, bin_extract_contig, elemimated_contig_total  = {}, {}, {}
                        for line in open('Coverage_'+str(level_num)+'_filtrated_bin_connecting_contigs.txt','r'):
                            bins=str(line).strip().split('\t')[0]
                            selected_1st[bins]={}
                            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
                            for item in dict_key_list:
                                selected_1st[bins][item.split('\'')[1]]=int(item.split(':')[1].strip()) ###

                        for line in open('Bin_extract_contigs_after_coverage_filtration_'+str(level_num)+'.txt','r'):
                            bins=str(line).strip().split('\t')[0]
                            bin_extract_contig[bins]={}
                            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
                            for item in dict_key_list:
                                bin_extract_contig[bins][item.split('\'')[1]]=item.split(':')[1].strip() ###

                        for line in open('Coverage_eliminated_bin_connecting_contigs_'+str(level_num)+'.txt','r'):
                            bins=str(line).strip().split('\t')[0]
                            elemimated_contig_total[bins]={}
                            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
                            for item in dict_key_list:
                                elemimated_contig_total[bins][item.split('\'')[1]]=item.split(':')[1].strip() ###

                        Parsing_kmer_file(assemblies_list, binset, bin_extract_contig, num_threads)

                        try:
                            os.mkdir('Bin_kmer')
                        except:
                            os.system('rm -rf Bin_kmer')
                            os.mkdir('Bin_kmer')
                            print('Bin_kmer exists')

                        print(str(len(bin_extract_contig)), 'bins are waiting for TNFs filtration')

                        pool=Pool(processes=num_threads)
                        bin_extract_contig_TNF, elemimated_contig_TNF, TNFs_exceptional_contigs, bin_extract_contig_TNF_list, elemimated_contig_TNF_list, TNFs_exceptional_contigs_list, m ={}, {}, {}, [], [], [], 0
                        for bin in bin_connecting_contigs_total.keys():
                            if bin in bin_extract_contig.keys():
                                m+=1
                                pool.apply_async(TNF_filtration,args=(bin, m, level_num, cutoff))
                                bin_extract_contig_TNF_list.append(str(bin)+'_'+str(level_num)+'_bin_extract_contig_TNF.txt')
                                elemimated_contig_TNF_list.append(str(bin)+'_'+str(level_num)+'_elemimated_contig_TNF.txt')
                                TNFs_exceptional_contigs_list.append(str(bin)+'_'+str(level_num)+'_TNFs_exceptional_contigs.txt')
                        pool.close()
                        pool.join()

                        for bin in bin_connecting_contigs_total.keys():
                            if bin in bin_extract_contig.keys():
                                os.system('mv '+str(bin)+'_bin_extract_contig_TNF.txt '+str(bin)+'_'+str(level_num)+'_bin_extract_contig_TNF.txt')
                                os.system('mv '+str(bin)+'_elemimated_contig_TNF.txt '+str(bin)+'_'+str(level_num)+'_elemimated_contig_TNF.txt')
                                os.system('mv '+str(bin)+'_TNFs_exceptional_contigs.txt '+str(bin)+'_'+str(level_num)+'_TNFs_exceptional_contigs.txt')

                        try:
                            os.system('mkdir TNF_filtration_matrix_'+str(level_num))
                        except:
                            # print('Folder TNF_filtration_matrix_2 existed')
                            print('Folder TNF_filtration_matrix_'+str(level_num)+' existed. Deleteced and recreated.')
                            os.system('rm -rf TNF_filtration_matrix_'+str(level_num))
                            os.system('mkdir TNF_filtration_matrix_'+str(level_num))

                        for item in bin_extract_contig_TNF_list:
                            bin_extract_contig_TNF.update(parse_dict(item))
                            os.system('mv '+str(item)+' TNF_filtration_matrix_'+str(level_num))

                        for item in elemimated_contig_TNF_list:
                            elemimated_contig_TNF.update(parse_dict(item))
                            os.system('mv '+str(item)+' TNF_filtration_matrix_'+str(level_num))

                        elemimated_contig_total.update(elemimated_contig_TNF)

                        for item in TNFs_exceptional_contigs_list:
                            TNFs_exceptional_contigs.update(parse_dict(item))
                            os.system('mv '+str(item)+' TNF_filtration_matrix_'+str(level_num))    

                        f=open('TNF_filtrated_bin_connecting_contigs_'+str(level_num)+'.txt','w')
                        for bin in bin_extract_contig_TNF.keys():
                            f.write(str(bin)+'\t'+str(bin_extract_contig_TNF[bin])+'\n')
                        f.close()

                        f=open('TNF_eliminated_bin_connecting_contigs_'+str(level_num)+'.txt','w')
                        for bin in elemimated_contig_TNF.keys():
                            f.write(str(bin)+'\t'+str(elemimated_contig_TNF[bin])+'\n')
                        f.close()

                        f=open('Total_eliminated_bin_connecting_contigs_'+str(level_num)+'.txt','w')
                        for bin in elemimated_contig_total.keys():
                            f.write(str(bin)+'\t'+str(elemimated_contig_total[bin])+'\n')
                        f.close()

                        f=open('TNFs_exceptional_contigs_'+str(level_num)+'.txt','w')
                        for bins in TNFs_exceptional_contigs.keys():
                            for contigs in TNFs_exceptional_contigs[bins].keys():
                                f.write(str(bins)+'\t'+str(contigs)+'\t'+str(TNFs_exceptional_contigs[bins][contigs])+'\n')
                        f.close()

                        bin_extract_contig_TNF={}
                        for line in open('TNF_filtrated_bin_connecting_contigs_'+str(level_num)+'.txt','r'):
                            bins=str(line).strip().split('\t')[0]
                            bin_extract_contig_TNF[bins]={}
                            dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
                            for item in dict_key_list:
                                bin_extract_contig_TNF[bins][item.split('\'')[1]]=int(item.split(':')[1].strip())

                        print('Writing retrieved bins')
                        try:
                            os.system('mkdir '+binset+'_retrieved_'+str(level_num))
                        except:
                            print(str(binset)+'_retrieved_'+str(level_num)+' existed. Deleted and recreated')
                            # os.system('rm -rf '+binset+'_retrieved_'+str(level_num))
                            # os.system('mkdir '+binset+'_retrieved_'+str(level_num))
                    
                        os.chdir(binset+'_retrieved_'+str(level_num))
                        bin_level_contigs={}
                        for item in bin_extract_contig_TNF.keys():
                        # for item in bin_extract_contig.keys():
                            bin_name_list=str(item).split('.')
                            bin_name_list.remove(bin_name_list[-1])
                            bin_name='.'.join(bin_name_list)
                            # for level in bin_connecting_contigs_total_level[item].keys():
                            f=open(bin_name+'.retrieved_level'+str(level_num)+'_'+str(cutoff)+'.fa','w')
                            bin_level_contigs[bin_name]={}
                            bin_level_contigs[bin_name]['.retrieved_level'+str(level_num)+'_'+str(cutoff)+'.fa']={}
                            bin_contigs_wrote, xxx, xxx2 = {}, 0, 0
                            if item in bin_contig.keys():
                                for ids in bin_contig[item].keys():
                                    xxx += 1
                                    xxx2 += 1
                                    f.write('>'+str(ids)+'\n'+str(bin_contig[item][ids])+'\n')
                                    bin_contigs_wrote[str(ids)]=''
                                    bin_level_contigs[bin_name]['.retrieved_level'+str(level_num)+'_'+str(cutoff)+'.fa'][str(ids)]=''

                            for contigs in bin_extract_contig_TNF[item].keys():
                            # for contigs in bin_extract_contig[item].keys():
                                if contigs not in bin_contigs_wrote.keys():
                                    xxx2 += 1
                                    f.write('>'+str(contigs)+'\n'+str(total_bin_contigs[contigs])+'\n')
                                    bin_level_contigs[bin_name]['.retrieved_level'+str(level_num)+'_'+str(cutoff)+'.fa'][str(contigs)]=''
                            f.close()

                            if xxx == xxx2:
                                os.system('rm -rf '+bin_name+'.retrieved_level'+str(level_num)+'_'+str(cutoff)+'.fa')
                                del bin_level_contigs[bin_name]
                        os.chdir(pwd)
                        #########

                    os.chdir(pwd)
                    grey_contigs, xxyy={}, 0
                    for line in open('Total_contig_grey_list.txt','r'):
                        xxyy+=1
                        try:
                            bin_name=str(line).strip().split('\t')[0].strip()
                            grey_contigs[bin_name]={}
                            item_dict=str(line).strip().split('\t')[1].strip().replace('\'','').replace('{','').replace('}','').split(',')
                            for item in item_dict:
                                contig=item.split(':')[0].strip()
                                level_id=item.split(':')[1].strip()
                                grey_contigs[bin_name][contig]=level_id
                        except:
                            print(str(xxyy)+' '+str(line)+' error')

                    os.chdir(pwd+'/'+str(binset)+'_retrieved_'+str(level_num))
                    bin_re, bin_re_list = {}, {}
                    for root, dirs, files in os.walk(pwd+'/'+str(binset)+'_retrieved_'+str(level_num)):
                        for file in files:
                            hz=str(file).split('.')[-1]
                            if 'fa' in file:
                                if '.retrieved_level' in file:
                                    bin_name=str(file).split('.retrieved_level')[0]
                                    cutoff=float(str(file).split('.retrieved_level')[1].split('_')[1].split('.')[0])
                                    if bin_name not in bin_re.keys():
                                        bin_re[bin_name]={}
                                        bin_re_list[bin_name]=[]
                                    bin_re[bin_name][cutoff]=file
                                    if cutoff not in bin_re_list[bin_name]:
                                        bin_re_list[bin_name].append(cutoff)

                    bin_remove={}
                    for bin_name in bin_re.keys():
                        if len(bin_re[bin_name]) >= 2:
                            for i in range(1, len(bin_re[bin_name])):
                                i0=i-1
                                if cutoff0 < cutoff1:
                                    cutoff0=bin_re_list[bin_name][i0]
                                    cutoff1=bin_re_list[bin_name][i]
                                else:
                                    cutoff0=bin_re_list[bin_name][i]
                                    cutoff1=bin_re_list[bin_name][i0]

                                bin0=bin_re[bin_name][cutoff0]
                                bin1=bin_re[bin_name][cutoff1]

                            bin0_id={}
                            for record in SeqIO.parse(bin0, 'fasta'):
                                bin0_id[record.id]=''

                            xxx = 0
                            for record in SeqIO.parse(bin1, 'fasta'):
                                if record.id not in bin0_id.keys():
                                    xxx=1

                        if xxx == 0:
                            bin_remove[bin1]=''

                    for item in bin_remove.keys():
                        os.system('rm '+item)

                print('Checking quality of retrieved bins')
                xyyzzz=0
                for root, dirs, files in os.walk(pwd+'/'+str(binset)+'_retrieved_'+str(level_num)):
                    os.chdir(pwd+'/'+str(binset)+'_retrieved_'+str(level_num))
                    for file in files:
                        hz=str(file).split('.')[-1]
                        if 'fa' in hz:
                            xyyzzz+=1
                os.chdir(pwd)

                if xyyzzz >= 1:
                    bins_checkm=checkm(str(binset)+'_retrieved_'+str(level_num), num_threads)

                    if int(level_num) == 1:   
                        bin_comparison(str(binset), bins_checkm, str(binset)+'_retrieved_'+str(level_num), num_threads, total_contigs, str(level_num))
                        latest_iteration=1
                        # fa_list=glob(r'./'+str(binset)+'_retrieved_1/*.fa')
                        # fa_list2=glob(r'./'+str(binset)+'_retrieved_1/*.fasta')
                        # t=len(fa_list)+len(fa_list2)
                        # if t == 0:
                        #     o_binset=binset
                        # else:
                        #     o_binset=+str(binset)+'_retrieved_1'
                    else:
                        
                        if level_num in bin_connecting_contigs_total_level.keys():
                            # actual_num=comparied_level2[level_num]
                            # level_num2=comparied_level[actual_num-1]
                            # level_num2=int(level_num)-1
                            
                            xyz=0
                            while latest_iteration != 0:
                                xyz+=1
                                latest_iteration=level_num-xyz
                                if latest_iteration >= 1:
                                    fa_list=glob(r'./'+str(binset)+'_retrieved_'+str(latest_iteration)+'/*.fa')
                                    fa_list2=glob(r'./'+str(binset)+'_retrieved_'+str(latest_iteration)+'/*.fasta')
                                    t=len(fa_list)+len(fa_list2)
                                    if t == 0:
                                        latest_iteration=1
                                    else:
                                        o_binset=str(binset)+'_retrieved_'+str(latest_iteration)
                                        latest_iteration=0
                                else:
                                    o_binset=binset
                            bin_comparison(o_binset, bins_checkm, str(binset)+'_retrieved_'+str(level_num), num_threads, total_contigs, str(level_num))
                left_level=int(max_iteration)-int(level_num)
                print('Iteration-'+str(level_num)+' accomplished. Left '+str(left_level)+' iteration(s).')
                f_checkpoit=open('S7_contigs_retriev_within_group_checkpoint.txt','a')
                f_checkpoit.write(str(level_num)+' iteration accomplished'+'\n')
                f_checkpoit.close()
            
        os.chdir(pwd)
        os.system('mv *_filtrated_bin_connecting_contigs_'+str(level_num)+'.txt *_eliminated_bin_connecting_contigs_'+str(level_num)+'.txt '+binset+'_retrieved_'+str(level_num))
        os.system('rm -rf coverage_filtration_matrix_* TNF_filtration_matrix_*')
    ei = 0
    for i in range(1, max_iteration+1):
        try:
            fa_list=glob(r'./'+str(binset)+'_retrieved_'+str(i)+'/*.fa')
            fa_list2=glob(r'./'+str(binset)+'_retrieved_'+str(i)+'/*.fasta')
            t=len(fa_list)+len(fa_list2)
            if t != 0:
                ei=i
        except:
            fx=open('BASALT_log.txt','a')
            fx.write('There is not folder:'+str(binset)+'_retrieved_'+str(i)+'\n')
            fx.close()
            print('There is not folder:'+str(binset)+'_retrieved_'+str(i))
    if ei != 0:
        os.system('mv '+binset+'_retrieved_'+str(ei)+' '+binset+'_retrieved')
        fx=open('BASALT_log.txt','a')
        fx.write('Rename '+binset+'_retrieved_'+str(ei)+' to '+binset+'_retrieved'+'\n')
        fx.close()
        print('Rename '+binset+'_retrieved_'+str(ei)+' to '+binset+'_retrieved')
    else:
        fx=open('BASALT_log.txt','a')
        fx.write('Within group retrieve did not find any extract contig'+'\n')
        fx.close()  
        print('Within group retrieve did not find any extract contig')
    os.system('mv Bin_connecting_contigs* Total_contig_black_list.txt Total_contig_grey_list.txt '+binset+'_retrieved')
    for i in range(1, max_iteration+1):
        os.system('rm -rf Deep_retrieved_bins_'+str(i))
        os.system('rm -rf '+str(binset)+'_retrieved_'+str(i))
        os.system('rm -rf '+str(binset)+'_retrieved_'+str(i)+'_checkm')
        os.system('rm -rf coverage_filtration_matrix_'+str(i))
        os.system('rm -rf TNF_filtration_matrix_'+str(i))
    try:
        os.mkdir('bin_extract-eleminated-selected_contig')
    except:
        print('bin_extract-eleminated-selected_contig already existed')
    os.system('mv *_selected_contigs.txt *_elemimated_contig.txt *_bin_extract_contig.txt TNF_filtrated_bin_connecting_contigs_* TNF_eliminated_bin_connecting_contigs_* TNFs_exceptional_contigs_* Total_eliminated_bin_connecting_contigs_* bin_extract-eleminated-selected_contig')
    os.system('mv bin_connecting_contigs_total_* Bin_extract_contigs_after_coverage_filtration_* *_filtrated_bin_connecting_contigs.txt Coverage_eliminated_bin_connecting_contigs_* bin_extract-eleminated-selected_contig')

def binset_filtration(binset, pwd, cpn_cutoff, ctn_cutoff):
    filtrated_checkm={}
    os.chdir(pwd+'/'+str(binset))
    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        for file in files:
            if 'quality_report.tsv' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2:
                        binID=str(line).strip().split('\t')[0].strip()
                        try:
                            genome_size=str(line).strip().split('\t')[8].strip()
                            completeness=str(line).strip().split('\t')[1].strip()
                            contamination=str(line).strip().split('\t')[2].strip()
                            N50=str(line).strip().split('\t')[6].strip()
                        except:
                            genome_size=str(line).strip().split('\t')[1].strip()
                            completeness=str(line).strip().split('\t')[2].strip()
                            contamination=str(line).strip().split('\t')[3].strip()
                            N50=str(line).strip().split('\t')[4].strip()
                        
                        if float(completeness) >= int(cpn_cutoff) and float(contamination) < int(ctn_cutoff):
                            filtrated_checkm[str(binID)]=str(line)

                        # bins_checkm[genome_ids]['N50']=int(N50)
                        # bins_checkm[genome_ids]['Completeness']=float(completeness)
                        # bins_checkm[genome_ids]['Genome size']=int(genome_size)
                        # bins_checkm[genome_ids]['Contamination']=float(contamination)


    for root, dirs, files in os.walk(pwd+'/'+str(binset)):
        for file in files:
            file_name_list=str(file).split('.')
            file_name_list.remove(file_name_list[-1])
            file_name='.'.join(file_name_list)
            if file_name in filtrated_checkm.keys():
                os.system('cp '+file+' '+pwd+'/'+str(binset)+'_filtrated')

    os.chdir(pwd+'/'+str(binset)+'_filtrated')
    f=open('quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    for item in filtrated_checkm.keys():
        f.write(str(filtrated_checkm[item]))
    f.close()
    os.chdir(pwd)

def finding_black_contigs(binset, binset2, outlier_remover_folder, pwd, level_num):
    black_contigs, grey_contigs, total_contigs = {}, {}, {}
    try:
        for line in open(pwd+'/'+binset+'/Coverage_eliminated_bin_connecting_contigs.txt','r'):
            bins=str(line).strip().split('\t')[0]
            bins_list=bins.split('.')
            bins_list.pop()
            bins_name='.'.join(bins_list)
            black_contigs[bins_name]={}
            contigs_list=str(line).strip().split('\t')[1].replace('\'','').replace('{','').replace('}','').split(',')
            for item in contigs_list:
                contig_id=item.split(':')[0].strip()
                black_contigs[bins_name][contig_id]=''
    except:
        print('Coverage_eliminated_bin_connecting_contigs.txt not found')

    try:
        for line in open(pwd+'/'+binset+'/TNF_eliminated_bin_connecting_contigs.txt','r'):
            bins=str(line).strip().split('\t')[0]
            bins_list=bins.split('.')
            bins_list.pop()
            bins_name='.'.join(bins_list)
            if bins_name not in black_contigs.keys():
                black_contigs[bins_name]={}
            contigs_list=str(line).strip().split('\t')[1].replace('\'','').replace('{','').replace('}','').split(',')
            for item in contigs_list:
                contig_id=item.split(':')[0].strip()
                black_contigs[bins_name][contig_id]=''
    except:
        print('TNF_eliminated_bin_connecting_contigs.txt not found')

    try:
        n=0
        for line in open(pwd+'/'+binset2+'/Bin_contigs_white_black_list.txt','r'):
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
                if bins_name not in total_contigs.keys():
                    total_contigs[bins_name]={}
                total_contigs[bins_name][contig_id]=''
    except:
        print('Bin_contigs_white_black_list.txt not found')

    try:
        for line in open(pwd+'/'+binset+'/Total_eliminated_bin_connecting_contigs.txt','r'):
            bins=str(line).strip().split('\t')[0]
            bins_list=bins.split('.')
            bins_list.pop()
            bins_name='.'.join(bins_list)
            if bins_name not in black_contigs.keys():
                black_contigs[bins_name]={}
            contigs_list=str(line).strip().split('\t')[1].replace('\'','').replace('{','').replace('}','').split(',')
            for item in contigs_list:
                contig_id=item.split(':')[0].strip()
                black_contigs[bins_name][contig_id]=''
    except:
        print('Total_eliminated_bin_connecting_contigs.txt not found')

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
                
                if bins_name not in total_contigs.keys():
                    total_contigs[bins_name]={}
                total_contigs[bins_name][contig_id]=''
    except:
        print('Contig_list.txt not found')

    try:
        for line in open(pwd+'/'+binset+'/TNF_filtrated_bin_connecting_contigs.txt','r'):
            bins=str(line).strip().split('\t')[0]
            bins_list=bins.split('.')
            bins_list.pop()
            bins_name='.'.join(bins_list)
            if bins_name not in grey_contigs.keys():
                grey_contigs[bins_name]={}

            # if bins_name not in total_contigs.keys():
            #     total_contigs[bins_name]={}

            contigs_list=str(line).strip().split('\t')[1].replace('\'','').replace('{','').replace('}','').split(',')
            if bins_name not in black_contigs.keys():
                for item in contigs_list:
                    contig_id=item.split(':')[0].strip()
                    level_id=int(item.split(':')[1].strip())
                    grey_contigs[bins_name][contig_id]=level_id
                    # total_contigs[bins_name][contig_id]=level_id
            else:
                for item in contigs_list:
                    contig_id=item.split(':')[0].strip()
                    if contig_id not in black_contigs[bins_name].keys():
                        level_id=int(item.split(':')[1].strip())
                        grey_contigs[bins_name][contig_id]=level_id
                        # total_contigs[bins_name][contig_id]=level_id
    except:
        print('TNF_filtrated_bin_connecting_contigs.txt not found')


    black_contigs2={}
    for item in black_contigs.keys():
        black_contigs2[item]=[]
        for item2 in black_contigs[item].keys():
            black_contigs2[item].append(item2)
    
    total_contigs2={}
    for item in total_contigs.keys():
        total_contigs2[item]=[]
        for item2 in total_contigs[item].keys():
            total_contigs2[item].append(item2)    

    f=open('Total_contig_black_list.txt','w')
    for item in black_contigs2.keys():
        f.write(str(item)+'\t'+str(black_contigs2[item])+'\n')
    f.close()

    f=open('Total_contig_grey_list.txt','w')
    for item in grey_contigs.keys():
        f.write(str(item)+'\t'+str(grey_contigs[item])+'\n')
    f.close()

    f=open('Total_contig_eliminated_from_deep-refinment.txt','w')
    for item in total_contigs2.keys():
        f.write(str(item)+'\t'+str(total_contigs2[item])+'\n')
    f.close()

    return black_contigs2, grey_contigs, total_contigs2

def Contig_retrieve_within_group_main(binset, outlier_remover_folder, num_threads, parameter, cpn_cutoff, ctn_cutoff, assemblies_list, PE_connections_list, coverage_matrix_list):
    pwd=os.getcwd()

    last_step=0
    print('--------------------------------------')
    print('Processing contigs retrieving within group process')
    print('Binset: '+str(binset))
    print('Binset after outlier removal: '+str(outlier_remover_folder))
    print('The minimun completeness to keep: '+str(cpn_cutoff))
    print('The maximal contaminaition to keep: '+str(ctn_cutoff))
    print('Assemblies: '+str(assemblies_list))
    print('Connections: '+str(PE_connections_list))
    print('Coverage matrix: '+str(coverage_matrix_list))
    print('.....................................')
    try:
        flog=open('BASALT_log.txt','w')
    except:
        flog=open('BASALT_log.txt','a')

    flog.write('Processing contigs retrieving within group process'+'\n'+'Binset: '+str(binset)+'\n'+'Binset after outlier removal: '+str(outlier_remover_folder)+'\n')
    flog.write('The minimun completeness to keep: '+str(cpn_cutoff)+'\n'+'The maximal contaminaition to keep: '+str(ctn_cutoff)+'\n'+'Assemblies: '+str(assemblies_list)+'\n')
    flog.write('Connections: '+str(PE_connections_list)+'\n'+'Coverage matrix: '+str(coverage_matrix_list)+'\n')
    flog.close()

    # if parameter == 'last':
    #     try:
    #         n=0
    #         for line in open('Retrieve_from_contigs_within_group.txt', 'r'):
    #             n+=1

    #         n1=0
    #         for line in open('Retrieve_from_contigs_within_group.txt', 'r'):
    #             n1+=1
    #             if n1 == n:
    #                 last_step=int(str(line)[0])
    #     except:
    #         last_step=0
    #         f_cp=open('Retrieve_from_contigs_within_group.txt', 'w')
    #         f_cp.close()
    # else:
    #     last_step=0
    #     f_cp=open('Retrieve_from_contigs_within_group.txt', 'w')
    #     f_cp.close()
    # start_step=last_step+1
    # print('Start from step '+str(start_step))
    os.chdir(pwd)

    Contig_retrieve_within_group(assemblies_list, binset, outlier_remover_folder, PE_connections_list, num_threads, last_step, coverage_matrix_list, pwd, cpn_cutoff, ctn_cutoff)
    os.chdir(pwd)
    try:
        os.system('mkdir '+binset+'_retrieved_kmer')
    except:
        print(binset+'_retrieved_kmer already existed')    
    os.system('mv *a_kmer.txt '+binset+'_retrieved_kmer')
    os.system('rm -rf *a_kmer.txt '+binset+'_retrieved_kmer')
    # os.system('rm -rf *deep_refinement *deep_refinement_checkm')
    print('Contig retrieve within group done')

if __name__ == '__main__':
    best_binset_from_multi_assemblies='BestBinset_outlier_refined_filtrated_retrieved'
    outlier_remover_folder='BestBinset_outlier_refined'
    assemblies_list=['1_RH_insert_cat_270_final.contigs.fa','2_RH_S001_insert_270_final.contigs.fa','3_RH_S002_insert_270_final.contigs.fa','4_RH_S003_insert_270_final.contigs.fa','5_RH_S004_insert_270_final.contigs.fa','6_RH_S005_insert_270_final.contigs.fa']
    coverage_matrix_list=['Coverage_matrix_for_binning_1_RH_insert_cat_270_final.contigs.fa.txt','Coverage_matrix_for_binning_2_RH_S001_insert_270_final.contigs.fa.txt','Coverage_matrix_for_binning_3_RH_S002_insert_270_final.contigs.fa.txt','Coverage_matrix_for_binning_4_RH_S003_insert_270_final.contigs.fa.txt','Coverage_matrix_for_binning_5_RH_S004_insert_270_final.contigs.fa.txt','Coverage_matrix_for_binning_6_RH_S005_insert_270_final.contigs.fa.txt']
    PE_connections_list=['condense_connections_1_RH_insert_cat_270_final.contigs.fa.txt','condense_connections_2_RH_S001_insert_270_final.contigs.fa.txt','condense_connections_3_RH_S002_insert_270_final.contigs.fa.txt','condense_connections_4_RH_S003_insert_270_final.contigs.fa.txt','condense_connections_5_RH_S004_insert_270_final.contigs.fa.txt','condense_connections_6_RH_S005_insert_270_final.contigs.fa.txt']
    num_threads=60
    cpn_cutoff=35 ### The minimun completeness to keep
    ctn_cutoff=5 ### The maximal contaminaition to keep
    parameter='last' ### 'last' or 'none'. checkpoint: using 'last', the script will continue to run the previous incompleted job. 'none' suggests that the job will start from the beginning. 
    Contig_retrieve_within_group_main(best_binset_from_multi_assemblies, outlier_remover_folder, num_threads, parameter, cpn_cutoff, ctn_cutoff, assemblies_list, PE_connections_list, coverage_matrix_list)
