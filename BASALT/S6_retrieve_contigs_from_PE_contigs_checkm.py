#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio import SeqIO
import os, threading, glob, gc
from sklearn.decomposition import PCA
import numpy as np
from multiprocessing import Pool
# from Outlier_remover import *

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

def cycle_mt(bin_connecting_contigs, bin_connecting_contigs2, connections, bins):
    m, m_before, m_after=1, 0, 1
    print('Parsing', bins)
    while m <= 1: ### Consider connection level less than 2
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
    bin_connecting_contigs, bin_connecting_contigs2, result={}, {}, {}
    for bins in bin_contigs.keys():
        bin_connecting_contigs[bins]={}
        bin_connecting_contigs2[bins]={}
        bin_connecting_contigs[bins]=bin_contigs[bins]
        result[bins]=pool.apply_async(cycle_mt,args=(bin_connecting_contigs[bins], bin_connecting_contigs2[bins], connections, bins,))
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

# def checkm(bin_folder, num_threads):
#     pwd=os.getcwd()
#     os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(bin_folder)+' '+str(bin_folder)+'_checkm')
#     # os.chdir(str(bin_folder)+'_checkm/storage/')
    # print('Parsing '+bin_folder+' checkm output')
    # refined_checkm={}
    # for line in open('bin_stats_ext.tsv','r'):
    #     binID=str(line).strip().split('{\'')[0].strip()
    #     genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
    #     taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
    #     completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
    #     contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].split('}')[0].strip()
    #     # GC=round(float(str(line).strip().split('\'GC\':')[1].split(', \'GCN4\'')[0].strip())*100, 1)
    #     refined_checkm[str(binID)]={}
    #     refined_checkm[str(binID)]['marker lineage']=str(taxon)
    #     refined_checkm[str(binID)]['Completeness']=float(completeness)
    #     refined_checkm[str(binID)]['Genome size']=float(eval(genome_size))
    #     refined_checkm[str(binID)]['Contamination']=float(contamination)
    os.chdir(pwd)
    return refined_checkm

def bin_comparison(original_bin_folder, new_bins_checkm, new_bin_folder, refinement_mode, num_threads):
    pwd=os.getcwd()
    print('Comparing bins before and after refining process')
    os.chdir(pwd+'/'+str(original_bin_folder))
    bin_checkm, bin_seq_rec, selected_bin = {}, {}, {}
    for root, dirs, files in os.walk(pwd+'/'+str(original_bin_folder)):
        for file in files:
            if 'bin_stats_ext.tsv' in file:
                for line in open(file, 'r'):
                    binID=str(line).strip().split('{\'')[0].strip()
                    selected_bin[binID]=1
                    taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                    completeness=float(str(line).strip().split('Completeness\':')[1].split(',')[0].strip())
                    contamination=float(str(line).strip().split('Contamination\':')[1].split(',')[0].split('}')[0].strip())
                    genome_size=float(eval(str(line).strip().split('Genome size\':')[1].split(',')[0].strip()))
                    bin_checkm[str(binID)]={}
                    # bin_checkm[str(binID)]['Connections']=int(connections)
                    bin_checkm[str(binID)]['marker lineage']=str(taxon)
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
                re_taxon=str(new_bins_checkm[refined_bin_id]['marker lineage'])
                re_genome_size=float(new_bins_checkm[refined_bin_id]['Genome size'])
                re_delta=re_cpn-re_ctn
                re_5delta=re_cpn-5*re_ctn
                for bin_id2 in bestbin[item].keys():
                    ori_bin=bin_id2 
                    # ori_connections=int(bestbin[item][bin_id2]['Connections'])
                    ori_cpn=float(bestbin[item][bin_id2]['Completeness'])
                    ori_ctn=float(bestbin[item][bin_id2]['Contamination'])
                    ori_taxon=str(bestbin[item][bin_id2]['marker lineage'])
                    ori_genome_size=float(bestbin[item][bin_id2]['Genome size'])
                    ori_delta=ori_cpn-ori_ctn
                    ori_5delta=ori_cpn-5*ori_ctn

                if re_delta > ori_delta:
                    del bestbin[item][ori_bin]
                    bestbin[item][refined_bin_id]=new_bins_checkm[refined_bin_id] 
                    if re_5delta < ori_5delta and ori_cpn >= 40 and ori_ctn <= 15 and re_ctn <= 50:
                        further_refined_bin[refined_bin_id]=item
                        f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                    elif re_cpn < ori_cpn and ori_cpn >= 40 and ori_ctn <= 15 and re_ctn <= 50:
                        further_refined_bin[refined_bin_id]=item
                        f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                elif re_cpn > ori_cpn:
                    if ori_cpn >= 40 and ori_ctn <= 10 and re_ctn <= 50:
                        further_refined_bin[refined_bin_id]=item
                        f2.write(str(item)+'\t'+str(refined_bin_id)+'\n')
                else:
                    continue
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
                os.system('mkdir Deep_retrieved_bins')
            except:
                print('Deep_retrieved_bins folder existed')

            os.chdir(pwd+'/'+str(new_bin_folder))
            if os.path.exists('Bin_contigs_white_black_list.txt'):
                fxxx=open('Bin_contigs_white_black_list.txt','a')
            else:
                fxxx=open('Bin_contigs_white_black_list.txt','w')
            os.chdir(pwd)

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
                    n+=1
                    f=open(org_bin_id+'_deep_'+str(n)+'.fa','w')
                    f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                    contig_num[n]=str(contigs)
                    for contigs_id in bin_seq_rec[org_bin_id].keys():
                        f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                    f.close()
                
                # f_contig_num=open(org_bin_id+'_contig_id.txt','w')
                # for item in contig_num.keys():
                #     f_contig_num.write(str(item)+'\t'+str(contig_num[item])+'\n')
                # f_contig_num.close()

                marker_lineage=bin_checkm[str(org_bin_id)]['marker lineage']
                if '__' in marker_lineage:
                    marker_lineage_level=marker_lineage.split('__')[0].strip()
                    marker_lineage_taxon=marker_lineage.split('__')[1].strip()
                    if marker_lineage_level == 'k':
                        marker_lineage_level_name='domain'
                    elif marker_lineage_level == 'p':
                        marker_lineage_level_name='phylum'
                    elif marker_lineage_level == 'f':
                        marker_lineage_level_name='family'
                    elif marker_lineage_level == 'o':
                        marker_lineage_level_name='order'
                    elif marker_lineage_level == 'c':
                        marker_lineage_level_name='class'
                    elif marker_lineage_level == 's':
                        marker_lineage_level_name='species'
                    elif marker_lineage_level == 'g':
                        marker_lineage_level_name='genus'
                elif 'root' in marker_lineage:
                    marker_lineage_level_name='life'

                split_folder=1
                if n%500 == 0:
                    split_folder=n/500
                    num_bin_split=n/split_folder
                else:
                    split_folder=int(n/500)+1
                    num_bin_split=int(n/split_folder)+1

                os.chdir(pwd)
                if split_folder != 1:
                    m, i, num_moved, floder_name=0, 1, 0, {}
                    for ix in range(1, split_folder+1):
                        os.system('mkdir '+org_bin_id+'_deep_retrieval_'+str(ix))
                        floder_name[org_bin_id+'_deep_retrieval_'+str(ix)]=''

                    os.chdir(pwd+'/'+org_bin_id+'_deep_retrieval')
                    for root, dirs, files in os.walk(pwd+'/'+org_bin_id+'_deep_retrieval'):
                        for file in files:
                            if '.fa' in file:
                                m+=1
                                if m < num_bin_split:
                                    os.system('mv '+str(file)+' '+pwd+'/'+org_bin_id+'_deep_retrieval_'+str(i))
                                elif m == num_bin_split:
                                    os.system('mv '+str(file)+' '+pwd+'/'+org_bin_id+'_deep_retrieval_'+str(i))
                                    m=0
                                    i+=1
                os.chdir(pwd)

                # try:
                xxxx=0
                if os.path.exists(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm/storage/bin_stats_ext.tsv'):
                    xxxx=1
                else:
                    xxxx=0
                    os.system('rm -rf '+str(org_bin_id)+'_deep_retrieval_checkm')
                print(str(marker_lineage_level_name)+': '+str(bin_checkm[str(org_bin_id)]['marker lineage']))

                if xxxx == 0:
                    if split_folder != 1:
                        bin_stats_record=[]
                        for f_name in floder_name.keys():
                            os.system('S6p_Run_checkm_taxonomy_wf.py -l '+str(marker_lineage_level_name)+' -t '+str(num_threads)+' -c '+str(marker_lineage_taxon)+' -f '+str(f_name))
                            if os.path.exists(pwd+'/'+str(f_name)+'_checkm/storage/bin_stats_ext.tsv'):
                                print(str(f_name)+' taxonomy wf done')
                            else:
                                os.system('rm -rf '+str(f_name)+'_checkm')
                                os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(f_name)+' '+str(f_name)+'_checkm')
                            os.chdir(pwd+'/'+str(f_name)+'_checkm/storage')
                            for line in open('bin_stats_ext.tsv', 'r'):
                                bin_stats_record.append(str(line).strip())
                            os.chdir(pwd)
                        os.system('mkdir '+str(org_bin_id)+'_deep_retrieval_checkm')
                        os.system('mkdir '+str(org_bin_id)+'_deep_retrieval_checkm/storage')
                        os.chdir(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm/storage')
                        f_checkm_total=open('bin_stats_ext.tsv','w')
                        for item in bin_stats_record:
                            f_checkm_total.write(item+'\n')
                        f_checkm_total.close()
                        xxxx=1
                        os.chdir(pwd)
                    else:
                        if marker_lineage_taxon != 'algicola':
                            os.system('S6p_Run_checkm_taxonomy_wf.py -l '+str(marker_lineage_level_name)+' -t '+str(num_threads)+' -c '+str(marker_lineage_taxon)+' -f '+str(org_bin_id)+'_deep_retrieval')
                            if os.path.exists(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm/storage/bin_stats_ext.tsv'):
                                print(str(org_bin_id)+' taxonomy wf done')
                                xxxx=1
                            else:
                                os.system('rm -rf '+str(org_bin_id)+'_deep_retrieval_checkm')
                                os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(org_bin_id)+'_deep_retrieval '+str(org_bin_id)+'_deep_retrieval_checkm')
                                if os.path.exists(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm/storage/bin_stats_ext.tsv'):
                                    print(str(org_bin_id)+' lineage wf done')
                                    xxxx=1
                        else:
                            os.system('rm -rf '+str(org_bin_id)+'_deep_retrieval_checkm')
                            os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(org_bin_id)+'_deep_retrieval '+str(org_bin_id)+'_deep_retrieval_checkm')
                            if os.path.exists(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm/storage/bin_stats_ext.tsv'):
                                print(str(org_bin_id)+' lineage wf done')
                                xxxx=1
                try:
                    if xxxx == 1:
                        os.chdir(pwd+'/'+str(org_bin_id)+'_deep_retrieval_checkm/storage/')
                        for line in open('bin_stats_ext.tsv', 'r'):
                            bin_id=str(line).strip().split('\t')[0].strip()
                            if bin_id == org_bin_id:
                                org_bin_cpn=float(str(line).strip().split('Completeness\':')[1].split(',')[0].split('}')[0].strip())
                                org_bin_ctn=float(str(line).strip().split('Contamination\':')[1].split(',')[0].split('}')[0].strip())
                                ori_delta=org_bin_cpn-5*org_bin_ctn
                                ori_1delta=org_bin_cpn-org_bin_ctn

                        for line in open('bin_stats_ext.tsv', 'r'):
                            refined_bin_id=str(line).strip().split('\t')[0].strip()
                            refined_bin_cpn=float(str(line).strip().split('Completeness\':')[1].split(',')[0].strip())
                            refined_bin_ctn=float(str(line).strip().split('Contamination\':')[1].split(',')[0].split('}')[0].strip())
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
                        os.chdir(pwd)

                    xxxx=0
                    if len(contig_neutral) != 0:
                        os.chdir(pwd+'/Deep_retrieved_bins')
                        f=open(org_bin_id+'_RT-D1.fa','w')
                        for contigs in contig_neutral.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                        for contigs in contig_positive.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                        for contigs_id in bin_seq_rec[org_bin_id].keys():
                            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                        f.close()
                        xxxx=1

                    if len(contig_positive) != 0:
                        os.chdir(pwd+'/Deep_retrieved_bins')
                        f=open(org_bin_id+'_RT-D2.fa','w')
                        for contigs in contig_positive.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                        for contigs_id in bin_seq_rec[org_bin_id].keys():
                            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                        f.close()
                        xxxx=1

                    if len(contig_positive_suspect) != 0:
                        os.chdir(pwd+'/Deep_retrieved_bins')
                        f=open(org_bin_id+'_RT-D3.fa','w')
                        for contigs in contig_neutral.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                        for contigs in contig_positive.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')
                        for contigs in contig_positive_suspect.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                        for contigs_id in bin_seq_rec[org_bin_id].keys():
                            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                        f.close()
                        xxxx=1

                    if len(contig_negative_suspect) != 0:
                        os.chdir(pwd+'/Deep_retrieved_bins')
                        f=open(org_bin_id+'_RT-D4.fa','w')
                        for contigs in contig_negative_suspect.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                        for contigs in contig_neutral.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                        for contigs in contig_positive.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                        for contigs in contig_positive_suspect.keys():
                            f.write('>'+str(contigs)+'\n'+str(further_refined_bin_contig[org_bin_id][contigs])+'\n')

                        for contigs_id in bin_seq_rec[org_bin_id].keys():
                            f.write('>'+str(contigs_id)+'\n'+str(bin_seq_rec[org_bin_id][contigs_id])+'\n')
                        f.close()
                        xxxx=1
                    os.chdir(pwd)
                except:
                    print(str(org_bin_id)+'_deep_retrieval did not performed correctly')
            os.chdir(pwd)

            xxxx=0
            for root, dirs, files in os.walk(pwd+'/Deep_retrieved_bins'):
                for file in files:
                    xxxx+=1

            if xxxx>=1:
                os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa Deep_retrieved_bins Deep_retrieved_bins_checkm')
                os.chdir(pwd+'/Deep_retrieved_bins_checkm/storage/')
                for line in open('bin_stats_ext.tsv', 'r'):
                    refined_bin_id=str(line).strip().split('\t')[0].strip()
                    refined_bin_cpn=float(str(line).strip().split('Completeness\':')[1].split(',')[0].split('}')[0].strip())
                    refined_bin_ctn=float(str(line).strip().split('Contamination\':')[1].split(',')[0].split('}')[0].strip())
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
                        new_bins_checkm[refined_bin_id]['Genome size']=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
                        new_bins_checkm[refined_bin_id]['marker lineage']=str(line).strip().split('lineage')[1].split('\'')[2].strip()
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
                                    new_bins_checkm[refined_bin_id]['Genome size']=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
                                    new_bins_checkm[refined_bin_id]['marker lineage']=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                                    bestbin[org_bin_id][refined_bin_id]=new_bins_checkm[refined_bin_id]
                                    try:
                                        os.system('cp '+pwd+'/Deep_retrieved_bins/'+refined_bin_id+'.fa '+pwd+'/'+str(new_bin_folder))
                                    except:
                                        os.system('cp '+pwd+'/Deep_retrieved_bins/'+refined_bin_id+'.fasta '+pwd+'/'+str(new_bin_folder))

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

    f_bin_checkm=open(new_bin_folder+'_bin_stats_ext.tsv', 'w')
    for item in selected_bin.keys():
        if '_RT-D' in item:
            item_name=str(item).split('_RT-D')[0]
            f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item])+'\n')
        elif '.retrieved_level' in item:
            item_name=str(item).split('.retrieved_level')[0]
            f_bin_checkm.write(str(item_name)+'\t'+str(total_bin_checkm[item])+'\n')
        else:
            f_bin_checkm.write(str(item)+'\t'+str(total_bin_checkm[item])+'\n')
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
                    if 'metabat' in file_name or 'metabinner' in file_name or 'vamb' in file_name:
                        os.system('mv '+file+' '+str(file_name)+'.fa')
                        f.write(str(file)+'\t'+str(file_name)+'.fa'+'\n')
                    else:
                        os.system('mv '+file+' '+str(file_name)+'.fasta')
                        f.write(str(file)+'\t'+str(file_name)+'.fasta'+'\n')
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

def parse_bin_in_bestbinset(assemblies_list, binset, outlier_remover_folder, PE_connections_list, num_threads, last_step, coverage_matrix_list, refinement_mode):
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
        bin_connecting_contigs_total, connection_total, bin_connecting_contigs_total_level = {}, {}, {}
        for i in range(0, len(assemblies_list)):
            A=PE_connecting_contigs(assemblies_list[i], PE_connections_list[i], binset, num_threads)
            bin_connecting_contigs=A[0]
            connections_single=A[1]
            bin_connecting_contigs_level=A[2]
            bin_connecting_contigs_total.update(bin_connecting_contigs)
            connection_total.update(connections_single)
            bin_connecting_contigs_total_level.update(bin_connecting_contigs_level)
            print('Parsed', str(assemblies_list[i]))
            print('--------------------')

        f=open('Bin_connecting_contigs.txt', 'w')
        f.write('Bin'+'\t'+'Connecting Contigs'+'\n')
        for item in bin_connecting_contigs_total.keys():
            f.write(str(item)+'\t'+str(bin_connecting_contigs_total[item])+'\n')
        f.close()

        f=open('Bin_connecting_contigs_level.txt', 'w')
        f.write('Bin'+'\t'+'Level'+'\t'+'Connecting contigs'+'\n')
        for item in bin_connecting_contigs_total_level.keys():
            for level in bin_connecting_contigs_total_level[item].keys():
                f.write(str(item)+'\t'+str(level)+'\t'+str(bin_connecting_contigs_total_level[item][level])+'\n')
        f.close()

        # f=open('Bin_kmer_total.txt', 'w')
        # f.write('Contig'+'\t'+'Kmer'+'\n')
        # for bins in bin_kmer_total.keys():
        #     for contigs in bin_kmer_total[bins].keys():
        #         f.write(str(contigs)+'\t'+str(contigs)+'\t'+'['+str(bin_kmer_total[bins][contigs])+'\n')
        # f.close()

        # f=open('Kmer_total.txt', 'w')
        # f.write('Contig'+'\t'+'Kmer'+'\n')
        # for item in kmer_total.keys():
        #     f.write(str(item)+'\t'+'['+str(kmer_total[item])+'\n')
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
                        contig=item.split('\'')[1]
                        bin_connecting_contigs_total[bins][contig]=int(item.split(':')[1].strip())
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
        gc.collect()

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
        f_cp.write('5th bins retrived done'+'\n')
        f_cp.close()

    if int(last_step) < 6:
        os.chdir(pwd)
        print('Checking quality of retrieved bins')
        # checkm(str(binset)+'_retrieved', num_threads)
        os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset)+'_retrieved '+str(binset)+'_retrieved_checkm')
        f_cp=open('S6_checkpoint.txt','a')
        f_cp.write('6th bins quality check done'+'\n')
        f_cp.close()

    if int(last_step) < 7:
        print('Comparing bins with retrieved bins')
        os.chdir(str(binset)+'_retrieved_checkm/storage/')
        print('Parsing '+str(binset)+'_retrieved_checkm output')
        bins_checkm={}
        try:
            for line in open('bin_stats_ext.tsv','r'):
                binID=str(line).strip().split('{\'')[0].strip()
                genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
                taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
                contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].split('}')[0].strip()
                # GC=round(float(str(line).strip().split('\'GC\':')[1].split(', \'GCN4\'')[0].strip())*100, 1)
                bins_checkm[str(binID)]={}
                bins_checkm[str(binID)]['marker lineage']=str(taxon)
                bins_checkm[str(binID)]['Completeness']=float(completeness)
                bins_checkm[str(binID)]['Genome size']=float(eval(genome_size))
                bins_checkm[str(binID)]['Contamination']=float(contamination)
        except:
            print('There is no retrieved bin')
        os.chdir(pwd)
        bin_comparison(str(binset), bins_checkm, str(binset)+'_retrieved', refinement_mode, num_threads)

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
            if 'bin_stats_ext.tsv' in file:
                for line in open(file,'r'):
                    #print(str(line).strip())
                    binID=str(line).strip().split('{\'')[0].strip()
                    genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].strip()
                    taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                    completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
                    contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].split('}')[0].strip()
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
    f=open('bin_stats_ext.tsv','w')
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

def Contig_recruiter_main(binset, outlier_remover_folder, num_threads, parameter, cpn_cutoff, ctn_cutoff, assemblies_list, PE_connections_list, coverage_matrix_list, refinement_mode, pwd):
    
    print('--------------------------------------')
    print('Processing contigs retrieving process')
    print('Binset: '+str(binset))
    print('Binset after outlier removal: '+str(outlier_remover_folder))
    print('The minimun completeness to keep: '+str(cpn_cutoff))
    print('The maximal contaminaition to keep: '+str(ctn_cutoff))
    print('Assemblies: '+str(assemblies_list))
    print('Connections: '+str(PE_connections_list))
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
#### Data feeding problematic
    # os.chdir(pwd+'/'+str(binset)+'_filtrated')
    # assembly_dict, PE_connections_list={}, []
    # for root, dirs, files in os.walk(pwd+'/'+str(binset)+'_filtrated'):
    #     for file in files:
    #         if '_genomes.' in file:
    #             qz=file.split('_genomes.')[0]
    #             assembly_name_list=qz.split('_')
    #             assembly_name_list.remove(assembly_name_list[-1])
    #             assembly_name_list.remove(assembly_name_list[-1])
    #             assembly_name='_'.join(assembly_name_list)
    #             assembly_dict[assembly_name]=1

    # os.chdir(pwd)
    # for item in assembly_dict.keys():
    #     ids_list=item.split('_')
    #     ids_list.remove(ids_list[0])
    #     connection_name='_'.join(ids_list)
    #     PE_connections_list.append('condense_connections_'+connection_name+'.txt')
#### Data feeding problematic

    parse_bin_in_bestbinset(assemblies_list, binset+'_filtrated', outlier_remover_folder, PE_connections_list, num_threads, last_step, coverage_matrix_list, refinement_mode)
    os.chdir(pwd)

if __name__ == '__main__': 
    best_binset_from_multi_assemblies='BestBinset_outlier_refined'
    outlier_remover_folder='BestBinset_outlier_refined'
    assemblies_list=['1_CAMI_cat_opera_contigs_polished.fasta','2_S001_opera_contigs_polished.fasta','3_S002_opera_contigs.fasta',
    '4_S003_opera_contigs.fasta','5_S004_opera_contigs.fasta','6_S005_opera_contigs.fasta']
    coverage_matrix_list=['Coverage_matrix_for_binning_1_CAMI_cat_opera_contigs_polished.fasta.txt','Coverage_matrix_for_binning_2_S001_opera_contigs_polished.fasta.txt',
                   'Coverage_matrix_for_binning_3_S002_opera_contigs.fasta.txt','Coverage_matrix_for_binning_4_S003_opera_contigs.fasta.txt',
                   'Coverage_matrix_for_binning_5_S004_opera_contigs.fasta.txt','Coverage_matrix_for_binning_6_S005_opera_contigs.fasta.txt']
    PE_connections_list=['condense_connections_CAMI_cat_opera_contigs_polished.fasta.txt','condense_connections_S001_opera_contigs_polished.fasta.txt',
                      'condense_connections_S002_opera_contigs.fasta.txt','condense_connections_S003_opera_contigs.fasta.txt',
                      'condense_connections_S004_opera_contigs.fasta.txt','condense_connections_S005_opera_contigs.fasta.txt']
    num_threads=60
    refinement_mode='deep' ### 'quick'; 'deep'
    cpn_cutoff=35 ### The minimun completeness to keep
    ctn_cutoff=20 ### The maximal contaminaition to keep
    pwd=os.getcwd()
    parameter='last' ### 'last' or 'none'. checkpoint: using 'last', the script will continue to run the previous incompleted job. 'none' suggests that the job will start from the beginning. 
    Contig_recruiter_main(best_binset_from_multi_assemblies, outlier_remover_folder, num_threads, parameter, cpn_cutoff, ctn_cutoff, assemblies_list, PE_connections_list, coverage_matrix_list, refinement_mode, pwd)
