#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio import SeqIO
import os, sys, threading, copy
from sklearn.decomposition import PCA
import numpy as np
from multiprocessing import Pool
# from Outlier_remover import *

def test_outlier(connecting_contig, item_data, test_index, coff):
    # print('Judging', connecting_contig
    four = pd.Series(item_data).describe()
    # print(four)
    # print('Q1= {0}, Q2= {1}, Q3={2}'.format(four['25%'],four['50%'],four['75%']))
    Q1 = four['25%']
    Q3 = four['75%']
    IQR = Q3 - Q1
    upper1 = Q3 + coff*IQR
    lower1 = Q1 - coff*IQR

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

def coverage_filtration_bin_mpt(bin_connecting_contigs_total, bin_contig_cov, contig_cov, selected_1st, bin_extract_contig, elemimated_contig, bin, m, coff):
    print('Bin-'+str(m)+' '+bin, 'coverage filtration of single copy contigs: bins with mt')
    for connecting_contig in bin_connecting_contigs_total[bin].keys():
        test_contig_bin_cov={}
        test_contig_bin_cov.update(bin_contig_cov[bin])
        contig_num=len(bin_contig_cov[bin])
        # print(str(contig_num)+' number')

        if connecting_contig in contig_cov.keys():
            test_contig_bin_cov[connecting_contig]=contig_cov[connecting_contig]
            num_coverage=len(contig_cov[connecting_contig])

        ### Find the possible multiple copies gene with ranking, transforming the data matrix
        test_contig_bin_cov2={}
        if num_coverage >= 2:
            ### Find and set the coordinate; average coverage >= 1
            total_depth, average_depth, max_depth, base_line_d = {}, 0, 0, 0
            for i in range(1, num_coverage+1):
                total_depth[i]=0

            xyzxyz=0
            for contig in test_contig_bin_cov.keys():
                xyzxyz+=1
                for i in range(1, num_coverage+1):
                    total_depth[i]+=test_contig_bin_cov[contig][i]
            
            ## Triming coverage
            passed_quality_coverage_index=[]
            for i in range(1, num_coverage+1):
                average_depth=total_depth[i]/xyzxyz
                if average_depth > max_depth and average_depth >= 1: ### Tuning parameter. This one is very import
                    max_depth=average_depth
                    base_line_d=i
                    passed_quality_coverage_index.append(i)
            
            if len(passed_quality_coverage_index) >= 2:
                for contig in bin_contig_cov[bin_id].keys():
                    test_contig_bin_cov2[contig]={}
                    base_line=test_contig_bin_cov[contig][base_line_d]
                    if base_line == 0:
                        base_line = 0.1

                    i2=0
                    for i in passed_quality_coverage_index:
                        if i != base_line_d:
                            i2+=1
                            normalized_data=test_contig_bin_cov[contig][i]/base_line
                            test_contig_bin_cov2[contig][i2]=normalized_data
            else:
                test_contig_bin_cov2=copy.deepcopy(test_contig_bin_cov)
        else:
            test_contig_bin_cov2=copy.deepcopy(test_contig_bin_cov)

        ### Filtration of contigs
        n1=0
        for contigs in test_contig_bin_cov2.keys():
            n1+=1
            if n1 == 1:
                num_coverage=len(test_contig_bin_cov2[contigs])
            else:
                break

        cov_index_total_num, cov_index={}, {}
        for item in test_contig_bin_cov2.keys():  
            for i in range(0, num_coverage):
                if i not in cov_index_total_num.keys():
                    cov_index_total_num[i]=float(test_contig_bin_cov2[item][i+1])
                    cov_index[i]=[]
                    cov_index[i].append(float(test_contig_bin_cov2[item][i+1]))
                else:
                    cov_index_total_num[i]+=float(test_contig_bin_cov2[item][i+1])
                    cov_index[i].append(float(test_contig_bin_cov2[item][i+1]))

        judgement, upper, lower = 0, {}, {}
        for item in cov_index.keys():
            cov_index[item].sort() ### sort low to high
            ### if contig_num >= 5:
            list_key_num=len(cov_index[item])
            p25=int(0.25*list_key_num)
            p75=int(0.75*list_key_num)

            n, Q3, Q1 =0, 0, 0
            for i in cov_index[item]:
                n+=1
                if n == p25:
                    Q1=i
                elif n == p75:
                    Q3=i
                else:
                    continue

            IQR = Q3 - Q1
            upper[item+1] = Q3 + coff*IQR
            lower[item+1] = Q1 - coff*IQR

        # print(str(len(upper))+' upper number')
        # print(str(len(contig_cov[connecting_contig]))+' contigs number')

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
            bin_outlier=test_outlier(connecting_contig, newData, test_index, coff)
            if bin_outlier == 1:
                if bin not in bin_extract_contig.keys():
                    bin_extract_contig[bin]={}
                    bin_extract_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    bin_extract_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            else:
                if bin not in elemimated_contig.keys():
                    elemimated_contig[bin]={}
                    # elemimated_contig_total[bin]={}
                    elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    # elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                    # elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
        else:
            if bin not in elemimated_contig.keys():
                elemimated_contig[bin]={}
                # elemimated_contig_total[bin]={}
                elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                # elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            else:
                elemimated_contig[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                # elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
    return selected_1st, bin_extract_contig, elemimated_contig

def TNF_filtration(bin_connecting_contigs_total, bin_extract_contig, bin_extract_contig_TNF, elemimated_contig_TNF, TNFs_exceptional_contigs, bin, m, coff):
    print('Bin-'+str(m)+' '+bin+' TNFs filtration of contigs')
    bin_kmer_total, kmer_total = {}, {}
    bin_kmer_total[bin]={}
    for line in open(bin+'_kmer.txt', 'r'):
        contig=str(line).strip().split('[')[0].strip()
        kmer=str(line).strip().split('[')[1].strip()
        bin_kmer_total[bin][contig]=kmer

    for line in open(bin+'_connecting_contigs_kmer.txt', 'r'): 
        kmer_total[str(line).strip().split('[')[0].strip()]=str(line).strip().split('[')[1].strip()

    for connecting_contig in bin_extract_contig[bin].keys():
        Bins_TNFs_test=[]
        test_contig_bin_kmer={}
        test_contig_bin_kmer.update(bin_kmer_total[bin])
        test_contig_bin_kmer[connecting_contig]=kmer_total[connecting_contig]
        num_contig=len(test_contig_bin_kmer)
        n1, n_connecting_contig=0, 0
        for item in test_contig_bin_kmer.keys():
            n1+=1
            if item == connecting_contig:
                n_connecting_contig=n1-1
            lis=test_contig_bin_kmer[item].split('\t')
            for i in range(0, len(lis)):
                Bins_TNFs_test.append(lis[i])

        TNF_array=np.array(Bins_TNFs_test).reshape((num_contig, 256))
        try:
            A=PCA_slector(TNF_array, num_contig)
            newData=A[0]
            explained_variance_ratio=A[1]
            bin_outlier=test_outlier(connecting_contig, newData, n_connecting_contig, coff)

            if bin_outlier == 1:
                if bin not in bin_extract_contig_TNF.keys():
                    bin_extract_contig_TNF[bin]={}
                    bin_extract_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    bin_extract_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            else:
                if bin not in elemimated_contig_TNF.keys():
                    elemimated_contig_TNF[bin]={}
                    elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                else:
                    elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

                # if bin not in elemimated_contig_total.keys():
                #     elemimated_contig_total[bin]={}
                #     elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
                # else:
                #     elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
        except:
            # print(bin, 'contig-', connecting_contig, 'TNFs clustering error')
            # print('Adding', connecting_contig, 'to eliminated contigs list')
            if bin not in TNFs_exceptional_contigs.keys():
                TNFs_exceptional_contigs[bin]={}
                TNFs_exceptional_contigs[bin][connecting_contig]=kmer_total[connecting_contig]
            else:
                TNFs_exceptional_contigs[bin][connecting_contig]=kmer_total[connecting_contig]

            if bin not in elemimated_contig_TNF.keys():
                elemimated_contig_TNF[bin]={}
                elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            else:
                elemimated_contig_TNF[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]

            # if bin not in elemimated_contig_total.keys():
            #     elemimated_contig_total[bin]={}
            #     elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
            # else:
            #     elemimated_contig_total[bin][connecting_contig]=bin_connecting_contigs_total[bin][connecting_contig]
    os.system('mv '+str(bin)+'_kmer.txt '+str(bin)+'_connecting_contigs_kmer.txt Bin_kmer')
    return bin_extract_contig_TNF, elemimated_contig_TNF, TNFs_exceptional_contigs
        # return bin_extract_contig_TNF, elemimated_contig_TNF, elemimated_contig_total, TNFs_exceptional_contigs

def parse_coverage_matrix(pwd, coverage_file):
    contig_cov={}
    for line in open(pwd+'/Bin_coverage_after_contamination_removal/'+coverage_file,'r'):
        contig_id=str(line).strip().split('\t')[0]
        contig_cov[contig_id]={}
        coverage_matrix_list=str(line).strip().replace('{','').replace('}','').split('\t')[1].split(',')
        for item in coverage_matrix_list:
            coverage_id=int(item.split(':')[0].strip())
            coverage=float(item.split(':')[1].strip())
            contig_cov[contig_id][coverage_id]=coverage
    return contig_cov

def contig_filtration(bin_id, parameter, m, coff, level_n, pwd):
    # print(bin_contig_cov)
    bin_connecting_contigs_total, connecting_contigs, n = {}, {}, 0
    try:
        for line in open('Bin_connecting_contigs.txt', 'r'):
            n+=1
            if n >= 2:
                bins=str(line).strip().split('\t')[0]
                if bins == bin_id:
                    bin_connecting_contigs_total[bins]={}
                    dict_key_list=str(line).strip().replace('{','').replace('}','').split('\t')[1].split(',')
                    for item in dict_key_list:
                        bin_connecting_contigs_total[bins][item.split('\'')[1]]=int(item.split(':')[1].strip())
                        connecting_contigs[item.split('\'')[1]]=0
    except:
        for line in open('Combat_extract_contigs_level2.txt','r'):
            level_num=int(str(line).strip().split('\t')[0].strip())
            if level_num == int(level_n):
                target_bin=str(line).strip().split('\t')[1].strip()
                if target_bin == bin_id:
                    bin_connecting_contigs_total[target_bin]={}
                    bin_name_list=target_bin.split('.')
                    bin_name_list.pop()
                    bin_name='.'.join(bin_name_list) 
                    contigs_list=str(line).strip().split('\t')[2].strip().replace('{','').replace('}','').replace('\'','').replace(':','').strip().split(',')
                    for contig in contigs_list:
                        bin_connecting_contigs_total[target_bin][contig.strip()]=0
                        connecting_contigs[contig.strip()]=0
    # else:
    #     print('No file found. Error!')

    bin_extract_contig, selected_1st, elemimated_contig ={}, {}, {}
    if parameter == 'coverage_filtration':
        bin_contig_cov, contig_cov = {}, {}
        for line in open(pwd+'/Bin_coverage_after_contamination_removal/Coverage_matrix_total.txt','r'):
            contig_id=str(line).strip().split('\t')[0]
            if contig_id in connecting_contigs.keys():
                contig_cov[contig_id]={}
                coverage_matrix_list=str(line).strip().replace('{','').replace('}','').split('\t')[1].split(',')
                for item in coverage_matrix_list:
                    coverage_id=int(item.split(':')[0].strip())
                    coverage=float(item.split(':')[1].strip())
                    contig_cov[contig_id][coverage_id]=coverage

        # contig_cov.update(parse_coverage_matrix(pwd, 'Coverage_matrix_total.txt'))
        bin_contig_cov[bin_id]={}
        bin_contig_cov[bin_id].update(parse_coverage_matrix(pwd, bin_id+'_coverage_matrix.txt'))
        A=coverage_filtration_bin_mpt(bin_connecting_contigs_total, bin_contig_cov, contig_cov, selected_1st, bin_extract_contig, elemimated_contig, bin_id, m, coff)
        selected_1st.update(A[0])
        bin_extract_contig.update(A[1])
        elemimated_contig.update(A[2])
        # elemimated_contig_total.update(A[3])

        f=open(bin_id+'_selected_contigs.txt', 'w')
        f1=open(bin_id+'_bin_extract_contig.txt', 'w')
        f2=open(bin_id+'_elemimated_contig.txt', 'w')
        for item in selected_1st.keys():
            f.write(str(item)+'\t'+str(selected_1st[item])+'\n')
        f.close()

        for item in bin_extract_contig.keys():
            f1.write(str(item)+'\t'+str(bin_extract_contig[item])+'\n')
        f1.close()

        for item in elemimated_contig.keys():
            f2.write(str(item)+'\t'+str(elemimated_contig[item])+'\n')
        f2.close()
        bin_contig_cov, contig_cov = {}, {}

    elif parameter == 'TNF_filtration':
        bin_extract_contig={}
        try:
            for line in open('Bin_extract_contigs_after_coverage_filtration_'+str(level_n)+'.txt','r'):
                bins=str(line).strip().split('\t')[0]
                bin_extract_contig[bins]={}
                dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
                for item in dict_key_list:
                    bin_extract_contig[bins][item.split('\'')[1]]=item.split(':')[1].strip() ###
        except:
            for line in open('Bin_extract_contigs_after_coverage_filtration.txt','r'):
                bins=str(line).strip().split('\t')[0]
                bin_extract_contig[bins]={}
                dict_key_list=str(line).strip().split('\t')[1].replace('\"','').replace('{','').replace('}','').split(',')
                for item in dict_key_list:
                    bin_extract_contig[bins][item.split('\'')[1]]=item.split(':')[1].strip() ###

        bin_extract_contig_TNF, elemimated_contig_TNF, TNFs_exceptional_contigs={}, {}, {}
        A=TNF_filtration(bin_connecting_contigs_total, bin_extract_contig, bin_extract_contig_TNF, elemimated_contig_TNF, TNFs_exceptional_contigs, bin_id, m, coff)

        bin_extract_contig_TNF.update(A[0])
        elemimated_contig_TNF.update(A[1])
        TNFs_exceptional_contigs.update(A[2])

        f=open(bin_id+'_bin_extract_contig_TNF.txt', 'w')
        f1=open(bin_id+'_elemimated_contig_TNF.txt', 'w')
        f2=open(bin_id+'_TNFs_exceptional_contigs.txt', 'w')
        for item in bin_extract_contig_TNF.keys():
            f.write(str(item)+'\t'+str(bin_extract_contig_TNF[item])+'\n')
        f.close()

        for item in elemimated_contig_TNF.keys():
            f1.write(str(item)+'\t'+str(elemimated_contig_TNF[item])+'\n')
        f1.close()

        for item in TNFs_exceptional_contigs.keys():
            f2.write(str(item)+'\t'+str(TNFs_exceptional_contigs[item])+'\n')
        f2.close()
        bin_extract_contig={}
    bin_connecting_contigs_total, connecting_contigs, n = {}, {}, 0
    bin_extract_contig, selected_1st, elemimated_contig ={}, {}, {}

if __name__ == '__main__': 
    binID_num, bin_num, para_num=0, 0, 0
    for i in range(1, len(sys.argv)):
        if '-b' in str(sys.argv[i]):
            binID_num=int(i)+1
        elif '-n' in str(sys.argv[i]):
            bin_num=int(i)+1
        elif '-p' in str(sys.argv[i]):
            para_num=int(i)+1
        elif '-c' in str(sys.argv[i]):
            cutoff_num=int(i)+1
        elif '-l' in str(sys.argv[i]):
            level_num=int(i)+1
        else:
            continue

    bin_id=str(sys.argv[binID_num])
    num=int(sys.argv[bin_num])
    ### 'coverage_filtration' or 'TNF_filtration'
    parameter=str(sys.argv[para_num])
    coff=float(sys.argv[cutoff_num])
    level_n=int(sys.argv[level_num])
    pwd=os.getcwd()
        
    # binID=bin_name
    # num=m
    # parameter='coverage_filtration_bin_mpt' 
    contig_filtration(bin_id, parameter, num, coff, level_n, pwd)
