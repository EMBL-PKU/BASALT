#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

from ast import excepthandler
from stat import S_ISBLK
from Bio import SeqIO
from Bio.Seq import Seq
try:
    from Bio.Alphabet import generic_dna
    generic_dna_t=1
except ImportError:
    generic_dna_t=0
import os, copy
from sklearn.decomposition import PCA
import numpy as np
import pandas as pd
from multiprocessing import Pool

def self_connecting(bins, bin_folder, bin_contigs_length, pwd, num_threads):
    os.system('makeblastdb -dbtype nucl -in '+str(bins)+' -logfile makeblastdb_log.txt')
    os.system('blastn -query '+str(bins)+' -db '+str(bins)+' -outfmt 6 -evalue 1e-5 -max_target_seqs 50 -out '+str(bins)+'_selfblast.txt -num_threads 1')
    os.system('rm '+str(bins)+'.nsq '+str(bins)+'.nin '+str(bins)+'.nhr')
    total_seq={}
    for record in SeqIO.parse(str(bins),'fasta'):
        total_seq[record.id]=record.seq
    yy=len(total_seq)

    fil = {}
    for line in open(str(bins)+'_selfblast.txt','r'):
        que=str(line).strip().split('\t')[0]
        sub=str(line).strip().split('\t')[1]
        sim=float(str(line).strip().split('\t')[2])
        leng=int(str(line).strip().split('\t')[3])
        if sim >= 99 and leng >= 30:
            if que != sub:
                que_s=int(str(line).strip().split('\t')[6])
                que_e=int(str(line).strip().split('\t')[7])
                sub_s=int(str(line).strip().split('\t')[8])
                sub_e=int(str(line).strip().split('\t')[9])
                if que_s == 1 or que_e == bin_contigs_length[bins][que]:
                    if que_s == 1:
                        if sub_e > sub_s:
                            if sub_e == bin_contigs_length[bins][sub]:
                                if que > sub:
                                    fil[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                else:
                                    fil[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line
                        else:
                            if sub_e == 1:
                                if que > sub:
                                    fil[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                else:
                                    fil[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line
                    elif que_e == bin_contigs_length[bins][que]:
                        if sub_e > sub_s:
                            if sub_s == 1:
                                if que > sub:
                                    fil[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                else:
                                    fil[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line
                        else:
                            if sub_s == bin_contigs_length[bins][sub]:
                                if que > sub:
                                    fil[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                else:
                                    fil[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line
                                  
    branch, branch2 = {}, {}               
    if len(fil) != 0:
        f=open('Filtrated_'+str(bins)+'_selfblast.txt','w')
        for item in fil.keys():
            f.write(str(fil[item]))
            que=str(fil[item]).strip().split('\t')[0]
            sub=str(fil[item]).strip().split('\t')[1]
            branch[que]={}
            branch[sub]={}
            branch2[que]={}
            branch2[sub]={}

        for item in fil.keys():
            que=str(fil[item]).strip().split('\t')[0]
            sub=str(fil[item]).strip().split('\t')[1]
            que_s=int(str(fil[item]).strip().split('\t')[6])
            que_e=int(str(fil[item]).strip().split('\t')[7])
            sub_s=int(str(fil[item]).strip().split('\t')[8])
            sub_e=int(str(fil[item]).strip().split('\t')[9])    
            if que_s == 1:
                try:
                    branch[que][1]+=1
                except:
                    branch[que][1]=1

                try:
                    branch2[que]+=1
                except:
                    branch2[que]=1
            elif que_e == bin_contigs_length[bins][que]:
                try:
                    branch[que][2]+=1
                except:
                    branch[que][2]=1
                
                try:
                    branch2[que]+=1
                except:
                    branch2[que]=1
            
            if sub_s == 1:
                try:
                    branch[sub][1]+=1
                except:
                    branch[sub][1]=1
                
                try:
                    branch2[sub]+=1
                except:
                    branch2[sub]=1
            elif sub_e == bin_contigs_length[bins][sub]:
                try:
                    branch[sub][2]+=1
                except:
                    branch[sub][2]=1
                
                try:
                    branch2[sub]+=1
                except:
                    branch2[sub]=1
            elif sub_s == bin_contigs_length[bins][sub]:
                try:
                    branch[sub][1]+=1
                except:
                    branch[sub][1]=1

                try:
                    branch2[sub]+=1
                except:
                    branch2[sub]=1
            elif sub_e == 1:
                try:
                    branch[sub][2]+=1
                except:
                    branch[sub][2]=1
                
                try:
                    branch2[sub]+=1
                except:
                    branch2[sub]=1
        f.close()
        
        f1=open('Filtrated_'+str(bins)+'_groups.txt','w')
        f2=open('Filtrated_'+str(bins)+'_straight.txt','w')
        uniq_con, mul_con, merged_group = [], [], {}
        for item in fil.keys():
            uniq_con.append(item)
            que=str(fil[item]).strip().split('\t')[0]
            sub=str(fil[item]).strip().split('\t')[1]
            que_s=int(str(fil[item]).strip().split('\t')[6])
            que_e=int(str(fil[item]).strip().split('\t')[7])
            sub_s=int(str(fil[item]).strip().split('\t')[8])
            sub_e=int(str(fil[item]).strip().split('\t')[9])
            if branch2[que] == 1 and branch2[sub] == 1:
                f2.write(str(fil[item]))
            else:
                mul_con.append(str(fil[item]))
                f1.write(str(fil[item]))
        f1.close()
        f2.close()

        total_contigs={}
        for item in mul_con:
            que=str(item).split('\t')[0]
            sub=str(item).split('\t')[1]
            total_contigs[que]=0
            total_contigs[sub]=0

        mul_con2=copy.deepcopy(mul_con)
        x, processed = 0, {}
        for item in mul_con:
            x+=1
            que=str(item).split('\t')[0]
            sub=str(item).split('\t')[1]
            y, y2=len(processed), 1
            if x == 1 and len(processed) != len(total_contigs):
                if que not in processed.keys():
                    merged_group[que]={}
                    merged_group[que][que]=0
                    merged_group[que][sub]=0
                    processed[que]=0
                    processed[sub]=0
                    while y != y2:
                        y=len(processed)
                        for item2 in mul_con2:
                            if item2 != item:
                                que2=str(item2).split('\t')[0]
                                sub2=str(item2).split('\t')[1]
                                if que2 in merged_group[que].keys() or sub2 in merged_group[que].keys():
                                    merged_group[que][que2]=0
                                    merged_group[que][sub2]=0
                                    processed[que2]=0
                                    processed[sub2]=0
                        y2=len(processed)
                    x=0
                # elif sub not in processed.keys():
                #     merged_group[sub]={}
                #     merged_group[sub][que]=0
                #     merged_group[sub][sub]=0
                #     processed[que]=0
                #     processed[sub]=0
                #     while y != y2:
                #         y=len(processed)
                #         for item2 in mul_con2:
                #             if item2 != item:
                #                 que2=str(item2).split('\t')[0]
                #                 sub2=str(item2).split('\t')[1]
                #                 if que2 in merged_group[sub].keys() or sub2 in merged_group[sub].keys():
                #                     merged_group[sub][que2]=0
                #                     merged_group[sub][sub2]=0
                #                     processed[que2]=0
                #                     processed[sub2]=0
                #         y2=len(processed)
                #     x=0
        # print(str(merged_group))
        group_blast={}
        for line in open('Filtrated_'+str(bins)+'_groups.txt','r'):
            query=str(line).strip().split('\t')[0]
            subject=str(line).strip().split('\t')[1]
            if query in merged_group.keys():
                try:
                    group_blast[query][str(line).strip()]=0
                except:
                    group_blast[query]={}
                    group_blast[query][str(line).strip()]=0
            else:
                for que in merged_group.keys():
                    if query in merged_group[que].keys():
                        try:
                            group_blast[que][str(line).strip()]=0
                        except:
                            group_blast[que]={}
                            group_blast[que][str(line).strip()]=0
            # print(str(que)+'\t'+str(group_blast[que]))

        num_seq, processed_seq_num, processed_contigs, merged_seq, merged_file_name = 0, 0, {}, [], []
        for line in open ('Filtrated_'+str(bins)+'_straight.txt', 'r'):
            num_seq+=1
            query=str(line).strip().split('\t')[0]
            subject=str(line).strip().split('\t')[1]
            sim=float(str(line).strip().split('\t')[2])
            leng=int(str(line).strip().split('\t')[3])
            query_start=int(str(line).strip().split('\t')[6])
            query_end=int(str(line).strip().split('\t')[7])
            subject_start=int(str(line).strip().split('\t')[8])
            subject_end=int(str(line).strip().split('\t')[9])

            A=seq_merge(total_seq, query, subject, query_start, query_end, subject_start, subject_end, num_seq, bins)
            total_seq.update(A[0])
            # processed_contigs.update(A[1])
            processed_contigs[query]=1
            processed_contigs[subject]=1
            merged_seq.append(A[2])
            merged_file_name.append(A[3])
            processed_seq_num+=1
        
        # f=open('Filtrated_'+str(bins)+'_with_branches.txt','w')
        branches, branches_blast_output = {}, {}
        for group_id in group_blast.keys():
            for item in group_blast[group_id].keys():
                query=str(item).strip().split('\t')[0]
                subject=str(item).strip().split('\t')[1]
                branches[query], branches[subject] = {}, {}
                branches_blast_output[query], branches_blast_output[subject] = {}, {}
                branches_blast_output[query]['end']={}
                branches_blast_output[query]['start']={}
                branches_blast_output[subject]['end']={}
                branches_blast_output[subject]['start']={}
        
        for group_id in group_blast.keys():
            for item in group_blast[group_id].keys():
                query=str(item).strip().split('\t')[0]
                subject=str(item).strip().split('\t')[1]
                query_start=int(str(item).strip().split('\t')[6])
                query_end=int(str(item).strip().split('\t')[7])
                subject_start=int(str(item).strip().split('\t')[8])
                subject_end=int(str(item).strip().split('\t')[9])

                if query_end == bin_contigs_length[bins][query]:
                    try:
                        branches[query]['end']+=1
                    except:
                        branches[query]['end']=1
                    branches_blast_output[query]['end'][item]=0
                elif query_start == 1:
                    try:
                        branches[query]['start']+=1
                    except:
                        branches[query]['start']=1
                    branches_blast_output[query]['start'][item]=0
                if subject_end > subject_start:
                    if subject_start == 1:
                        try:
                            branches[subject]['start']+=1
                        except:    
                            branches[subject]['start']=1
                        branches_blast_output[subject]['start'][item]=0
                    else:
                        try:
                            branches[subject]['end']+=1
                        except:    
                            branches[subject]['end']=1
                        branches_blast_output[subject]['end'][item]=0
                else:
                    if subject_end == 1:
                        try:
                            branches[subject]['end']+=1
                        except:    
                            branches[subject]['end']=1
                        branches_blast_output[subject]['end'][item]=0
                    else:
                        try:
                            branches[subject]['start']+=1
                        except:    
                            branches[subject]['start']=1
                        branches_blast_output[subject]['start'][item]=0

        branches_blastoutput2, branches_blastoutput3 = {}, []
        f=open('Filtrated_'+str(bins)+'_with_branches.txt', 'w')
        for item in branches.keys():
            if 'start' in branches[item].keys():
                if branches[item]['start'] == 2:
                    branches_blastoutput2[item]={}
                    branches_blastoutput2[item]['start']=[]
                    f.write(str(item)+'\t'+'start')
                    for item2 in branches_blast_output[item]['start'].keys():
                        f.write('\t'+str(item2))
                        branches_blastoutput2[item]['start'].append(item2)
                        branches_blastoutput3.append(item2)
                    f.write('\n')

            if 'end' in branches[item].keys():
                if branches[item]['end'] == 2:
                    branches_blastoutput2[item]={}
                    branches_blastoutput2[item]['end']=[]
                    f.write(str(item))
                    for item2 in branches_blast_output[item]['end'].keys():
                        f.write('\t'+str(item2))
                        branches_blastoutput2[item]['end'].append(item2)
                        branches_blastoutput3.append(item2)
                    f.write('\n')
        f.close()

        group_blast2=copy.deepcopy(group_blast)
        for group_id in group_blast2.keys():
            for item in group_blast2[group_id].keys():
                if item in branches_blastoutput3:
                    try:
                        del group_blast[group_id][item]
                    except:
                        xyz=0

        for group_id in group_blast.keys():
            f=open('Filtrated_group_'+str(group_id)+'_blastoutput_without_branches.txt','w')
            for item in group_blast[group_id].keys():
                f.write(str(item)+'\n')
            f.close()
            # os.system('mv Filtrated_group_'+str(group_id)+'_blastoutput_without_branches.txt Selfblast_output')
        
        for group_id in group_blast.keys():
            print('Processing group '+str(group_id))
            contigs={}
            for item in group_blast[group_id].keys():
                query=item.strip().split('\t')[0]
                subject=item.strip().split('\t')[1]
                contigs[query]=0
                contigs[subject]=0
            print(str(contigs))

            if len(contigs) != 0:
                num, processed_item, xyok, processed_contigs2, new_contigs = 0, {}, 0, {}, {}
                # while len(processed_item) != len(group_blast[group_id]) and xyok <= len(group_blast[group_id]):
                for item in group_blast[group_id].keys():
                    num+=1
                    if num == 1:
                        query=str(item).strip().split('\t')[0]
                        subject=str(item).strip().split('\t')[1]
                        query_start=int(str(item).strip().split('\t')[6])
                        query_end=int(str(item).strip().split('\t')[7])
                        subject_start=int(str(item).strip().split('\t')[8])
                        subject_end=int(str(item).strip().split('\t')[9])
                        sim=float(str(item).strip().split('\t')[2])
                        leng=int(str(item).strip().split('\t')[3])

                        A=seq_merge(total_seq, query, subject, query_start, query_end, subject_start, subject_end, num_seq, bins) 
                        processed_contigs[query]=1
                        processed_contigs[subject]=1
                        processed_contigs2[query]=1
                        processed_contigs2[subject]=1
                        merged_seq_temp=A[2]
                        # merged_seq.append(A[2])
                        merged_file_name_temp=A[3]
                        merged_seq_id=A[4]
                        contigs[merged_seq_id]=0
                        merged_file_name.append(merged_file_name_temp)
                        new_contigs[merged_seq_id]=0
                        total_seq.update(A[0])
                        # processed_item[item]=0
                        processed_seq_num+=1

                num1=len(processed_contigs2)
                num2=0
                while num1 != num2 and xyok <= len(contigs):
                    xyok+=1
                    # if len(processed_contigs2) != len(contigs):
                    for contig in contigs.keys():
                        if contig not in processed_contigs2.keys():
                            temp_file_name=str(bins)+'_'+str(contig)+'.fa'
                            seq_id=contig

                            contigs_len = {}
                            for record in SeqIO.parse(str(merged_file_name_temp),'fasta'):
                                contigs_len[record.id]=len(record.seq)
                                # total_seq2[record.id]=record.seq

                            f=open(temp_file_name,'w')
                            if seq_id in total_seq.keys():
                                f.write('>'+str(seq_id)+'\n'+str(total_seq[seq_id])+'\n')
                                contigs_len[seq_id]=len(total_seq[seq_id])
                                # total_seq2[seq_id]=total_seq[seq_id]
                            f.close()
                            os.system('makeblastdb -dbtype nucl -in '+str(merged_file_name_temp)+' -logfile makeblastdb_log.txt')
                            os.system('blastn -db '+str(merged_file_name_temp)+' -query '+str(temp_file_name)+' -outfmt 6 -evalue 1e-3 -max_target_seqs 50 -out '+str(merged_file_name_temp)+'_blast.txt  -num_threads 1')
                            os.system('rm '+str(merged_file_name_temp)+'.nsq '+str(merged_file_name_temp)+'.nin '+str(merged_file_name_temp)+'.nhr')

                            fil2 = {}
                            for line in open(str(merged_file_name_temp)+'_blast.txt','r'):
                                sim=float(str(line).strip().split('\t')[2])
                                leng=int(str(line).strip().split('\t')[3])
                                if sim >= 99:
                                    que=str(line).strip().split('\t')[0]
                                    sub=str(line).strip().split('\t')[1]
                                    if que != sub:
                                        que_s=int(str(line).strip().split('\t')[6])
                                        que_e=int(str(line).strip().split('\t')[7])
                                        sub_s=int(str(line).strip().split('\t')[8])
                                        sub_e=int(str(line).strip().split('\t')[9])
                                        if que_s == 1 or que_e == contigs_len[que]:
                                            if sub_e > sub_s:
                                                if sub_s == 1 or sub_e == contigs_len[sub]:
                                                    if que > sub:
                                                        fil2[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                                    else:
                                                        fil2[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line

                                            else:
                                                if sub_e == 1 or sub_s == contigs_len[sub]:
                                                    if que > sub:
                                                        fil2[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                                    else:
                                                        fil2[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line
                            # os.system('rm '+str(temp_file_name)+' '+str(temp_file_name)+'_'+str(merged_seq_temp)+'.txt')
                            for item in fil2.keys():
                                # query=item.split('\t')[1]
                                # subject=item.split('\t')[0]
                                line=fil2[item].strip()
                                query=str(line).strip().split('\t')[0]
                                subject=str(line).strip().split('\t')[1]
                                query_start=int(str(line).strip().split('\t')[6])
                                query_end=int(str(line).strip().split('\t')[7])
                                subject_start=int(str(line).strip().split('\t')[8])
                                subject_end=int(str(line).strip().split('\t')[9])

                            A=seq_merge(total_seq, query, subject, query_start, query_end, subject_start, subject_end, num_seq, bins) 
                            processed_contigs[query]=1
                            processed_contigs[subject]=1
                            processed_contigs2[query]=1
                            processed_contigs2[subject]=1
                            merged_seq_temp=A[2]
                            # merged_seq.append(A[2])
                            merged_file_name_temp=A[3]
                            merged_seq_id=A[4]
                            # contigs[merged_seq_id]=0
                            merged_file_name.append(merged_file_name_temp)
                            total_seq.update(A[0])
                            os.system('rm '+str(temp_file_name)+' '+str(temp_file_name)+'_'+str(merged_file_name_temp)+'_blast.txt')
                            num2=len(processed_contigs2)
                            processed_seq_num+=1

        os.system('mv Filtrated_'+str(bins)+'_selfblast.txt '+str(target_bin_folder)+'_selfblast_output')
        os.system('mv Filtrated_'+str(bins)+'_groups.txt '+str(target_bin_folder)+'_selfblast_output')
        os.system('mv Filtrated_'+str(bins)+'_straight.txt '+str(target_bin_folder)+'_selfblast_output')
        os.system('mv Filtrated_'+str(bins)+'_with_branches.txt '+str(target_bin_folder)+'_selfblast_output')
        # os.system('mv '+str(merged_file_name_temp)+'_blast.txt '+str(target_bin_folder)+'_selfblast_output')

    os.system('mkdir '+str(bins)+'_merged_files')
    # print(str(merged_file_name))
    # print(processed_contigs)

    try:  
        for file_id in merged_file_name:
            os.system('mv '+str(file_id)+' '+str(bins)+'_merged_files')
        os.system('mv '+str(bins)+'_selfblast.txt '+str(target_bin_folder)+'_selfblast_output')
    except:
        xyzzz=0

    # if processed_seq_num != 0:
    try:
        if len(processed_contigs) != 0:
            fx=open(str(target_bin_folder)+'_sc_bins_summary.txt','a')
            bin_id=str(bins).split('.fa')[0]
            bin_id2=bin_id+'_sc.fa'
            f=open(bin_id2,'w')
            x=0
            for seq_id in total_seq.keys():
                if seq_id not in processed_contigs.keys():
                    x+=1
                    f.write('>'+str(seq_id)+'\n'+str(total_seq[seq_id])+'\n')
            f.close()
            fx.write(str(bin_id2)+'\t'+str(yy)+'\t'+str(x)+'\n')
            fx.close()
            os.system('mv '+bin_id2+' '+str(target_bin_folder)+'_self_connected_bins')
            os.system('rm '+str(bins))
        else:
            os.system('mv '+str(bins)+' '+str(target_bin_folder)+'_self_connected_bins')
    except:
        os.system('mv '+str(bins)+' '+str(target_bin_folder)+'_self_connected_bins')

def seq_merge(total_seq, query, subject, query_start, query_end, subject_start, subject_end, num_seq, target_bin):
    print('Merging '+str(target_bin)+' sequences: '+str(query)+' and '+str(subject))
    processed_contigs={}
    delta_target_seq_alignment=abs(query_end-query_start)+1
    delta_vs_seq_alignment=abs(subject_end-subject_start)+1
    if delta_vs_seq_alignment == len(total_seq[subject]): ### In this case, subject sequence will totally covered by query seq
        merged_seq=total_seq[query]
    else:
        # print(str(subject))
        # print(str(subject_start))
        # print(str(total_seq[subject]))
        target_seq_start_nucl=total_seq[query][query_start-1]
        elong_seq_start_nucl=total_seq[subject][subject_start-1]
        if delta_target_seq_alignment == len(total_seq[query]):
            if subject_end > subject_start: ### Positive direction
                if subject_start == 1:
                    if target_seq_start_nucl == elong_seq_start_nucl:
                        merged_seq=total_seq[query]+total_seq[subject][subject_end:]
                    else:
                        # c_seq=str(total_seq[subject]).complement()
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        merged_seq=total_seq[query]+c_seq[subject_end:]
                else: ### Alignment start from Suject End
                    if target_seq_start_nucl == elong_seq_start_nucl:
                        merged_seq=total_seq[subject][:subject_start-1]+total_seq[query] ### Be careful about the subject start
                    else:
                        # c_seq=str(total_seq[subject]).complement()
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        merged_seq=c_seq[:subject_start-1]+total_seq[query]
            else: ### Reverse direction
                # r_seq=str(total_seq[subject])[::-1]
                # subject_reverse_start_position=len(subject)-subject_start
                if subject_end == 1:
                    # if total_seq[query][-1] == r_seq[-1]:
                    if total_seq[query][-1] == total_seq[subject][0]:
                        # merged_seq=r_seq[:subject_reverse_start_position]+total_seq[query]
                        merged_seq=total_seq[subject][subject_start:][::-1]+total_seq[query]
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1])).complement()
                        merged_seq=c_seq+total_seq[query]
                else: ### Alignment start from Suject End
                    if total_seq[query][0] == total_seq[subject][-1]:
                        # merged_seq=total_seq[query]+r_seq[:subject_start-1] ### ?
                        # merged_seq=total_seq[query]+r_seq[delta_vs_seq_alignment:]
                        merged_seq=total_seq[query]+total_seq[subject][:subject_end][::-1]
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][:subject_end][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][:subject_end][::-1])).complement()
                        # merged_seq=total_seq[query]+c_seq[:subject_start-1]
                        merged_seq=total_seq[query]+c_seq
        else: ### part of query and subject seqs merge together
            if subject_end > subject_start: 
                if subject_start == 1: ### Positive direction
                    if total_seq[query][query_start-1] == total_seq[subject][0]:
                        # merged_seq=total_seq[query][:query_end-1]+total_seq[subject][subject_end:]
                        merged_seq=total_seq[query]+total_seq[subject][subject_end:] #check
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        # merged_seq=total_seq[query][:query_end-1]+c_seq[subject_end:]
                        merged_seq=total_seq[query]+c_seq[subject_end:]
                else: ### Subject End
                    if total_seq[query][0] == total_seq[subject][subject_start-1]:
                        merged_seq=total_seq[subject][:subject_start-1]+total_seq[query] #check
                    else:
                        # c_seq=str(total_seq[subject]).complement()
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject])).complement()
                        merged_seq=c_seq[:subject_start-1]+total_seq[query]
            else: ### Reverse seq
                if subject_end == 1:
                    if total_seq[query][query_end-1] == total_seq[subject][0]:
                        merged_seq=total_seq[subject][subject_start:][::-1]+total_seq[query] ### Checked
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][subject_start:][::-1])).complement()
                        merged_seq=c_seq+total_seq[query]
                else:
                    if total_seq[query][query_end-1] == total_seq[subject][subject_end-1]:
                        merged_seq=total_seq[query]+total_seq[subject][:subject_end-1][::-1] ### Checked
                    else:
                        if generic_dna_t == 1:
                            c_seq=Seq(str(total_seq[subject][:subject_end-1][::-1]), generic_dna).complement()
                        else:
                            c_seq=Seq(str(total_seq[subject][:subject_end-1][::-1])).complement()
                        merged_seq=total_seq[query]+c_seq

    if '--' in subject:
        subject1=subject.split('--')[0]
        subject2=subject.split('--')[-1]
        merged_subject=subject1+'--'+subject2
    else:
        merged_subject=subject

    if '--' in str(query):
        query_ori=str(query).split('--')[0]
        merged_file_name=target_bin+'--'+query_ori+'--'+merged_subject+'_'+str(num_seq)+'_merged_seq.txt'
        merged_seq_file=open(merged_file_name, 'w')
    else:
        merged_file_name=target_bin+'--'+query+'--'+merged_subject+'_'+str(num_seq)+'_merged_seq.txt'
        merged_seq_file=open(merged_file_name, 'w')

    merged_seq_file.write('>'+str(query)+'--'+str(subject)+'\n'+str(merged_seq)+'\n')
    merged_seq_file.close()
    merged_seq_id=str(query)+'--'+str(subject)
    total_seq[query+'--'+subject]=merged_seq
    processed_contigs[query]=1
    processed_contigs[subject]=1
    return total_seq, processed_contigs, merged_seq, merged_file_name, merged_seq_id

# def short_reads_connecting(target_bin_selfconecting_folder, short_reads_folder, num_threads, pwd):
#     os.chdir(pwd+'/'+str(short_reads_folder))
#     bin_sr={}
#     for root, dirs, files in os.walk(pwd+'/'+str(short_reads_folder)):
#         for file in files:
#             if '.fq' in file:
#                 bin_id=str(file).split('_')[0]
#                 sr=bin_id+'_sr.fa'
#                 bin_sr[bin_id]=sr
#                 try:
#                     f=open(sr,'a')
#                 except:
#                     f=open(sr,'w')

#                 for line in open(file,'r'):

def sr_blast(target_bin_folder, bins, sc_bins, bin_sr, self_connected_bins_folder, short_reads_fa_folder, pwd):
    bins_sr_fa=bin_sr[bins]
    sc_bins2=sc_bins[bins]
    print('BLASTing '+str(bins)+' against reads')
    os.system('cp '+pwd+'/'+str(self_connected_bins_folder)+'/'+str(sc_bins2)+' '+str(sc_bins2))
    bin_contig_len={}
    for record in SeqIO.parse(str(sc_bins2),'fasta'):
        bin_contig_len[record.id]=len(record.seq)

    xyz=0
    for line in open(pwd+'/'+str(short_reads_fa_folder)+'/'+str(bins_sr_fa),'r'):
        xyz+=1
        if xyz == 2:
            read_len=len(str(line).strip())
            break
        
    os.system('makeblastdb -in '+str(sc_bins2)+' -dbtype nucl -logfile makeblastdb_log.txt')
    os.system('blastn -query '+pwd+'/'+str(short_reads_fa_folder)+'/'+str(bins_sr_fa)+' -db '+str(sc_bins2)+' -num_threads 1 -evalue 1e-10 -outfmt 6 -out '+str(bins)+'_sr_blast.txt -max_target_seqs 5')
    filtrated, filtrated2, filtrated3 = {}, {}, {}
    for line in open(str(bins)+'_sr_blast.txt','r'):
        query=line.strip().split('\t')[0]
        query2=line.strip().split('\t')[0].split('/')[0]
        sub=line.strip().split('\t')[1]
        sim=float(line.strip().split('\t')[2])
        leng=int(line.strip().split('\t')[3])
        query_s=int(str(line).strip().split('\t')[6])
        query_e=int(str(line).strip().split('\t')[7])
        sub_s=int(str(line).strip().split('\t')[8])
        sub_e=int(str(line).strip().split('\t')[9])
        score=float(str(line).strip().split('\t')[-1])
        if sim >= 99 and leng >= 30:
            max_len=bin_contig_len[sub]
            # max_len=bin_contig_len[sub]-3
            # if sub_s <= 3 or sub_e <= 3 or sub_s >= max_len or sub_e >= max_len:
            if query_s == 1 or query_e == read_len:
                if sub_s == 1 or sub_e == 1 or sub_s == max_len or sub_e == max_len:
                    try:
                        filtrated[query][sub]=str(line).strip()
                    except:
                        filtrated[query]={}
                        filtrated[query][sub]=str(line).strip()
                    
                    try:
                        score1=filtrated3[query2][query]
                        if score > score1:
                            filtrated3[query2][query]=score
                            filtrated2[query2][query]=str(line).strip()
                    except:
                        filtrated3[query2]={}
                        filtrated3[query2][query]=score
                        filtrated2[query2]={}
                        filtrated2[query2][query]=str(line).strip()

    f=open('Filtrated_'+str(bins)+'_sr_blast.txt', 'w')
    f1=open('Filtrated1_'+str(bins)+'_sr_blast.txt', 'w')
    for query in filtrated.keys():
        ali, ali2 = {}, {}
        if len(filtrated[query]) >= 2:
            # n=0
            for sub in filtrated[query].keys():
                f1.write(str(filtrated[query][sub])+'\n')
                # if n == 1:
                line=str(filtrated[query][sub])
                sub1_s=str(line).strip().split('\t')[6]
                sub1_e=str(line).strip().split('\t')[7]
                score=str(line).strip().split('\t')[-1]
                # if sub1_s == 1 and sub1_e == read_len:  ### Negalect of this because of contigs may already connected to each other in previously step
                #     n+=1
                #     alix[n]=line
                # else:
                if sub1_s == 1:
                    try:
                        score1=ali[sub1_s]
                        if score > score1:
                            ali[sub1_s]=score
                            ali2[sub1_s]=line
                    except:
                        ali[sub1_s]=score
                        ali2[sub1_s]=line
                elif sub1_e == read_len:
                    try:
                        score1=ali[sub1_e]
                        if score > score1:
                            ali[sub1_e]=score
                            ali2[sub1_e]=line
                    except:
                        ali[sub1_e]=score
                        ali2[sub1_e]=line
        
        if len(ali) == 2:
            f.write(str(ali2[1])+'\n')
            f.write(str(ali2[read_len])+'\n')
    f.close()
    f1.close()

    f=open('Filtrated_'+str(bins)+'_sr_blast_pe.txt', 'w')
    for query2 in filtrated2.keys():
        if len(filtrated2[query2]) == 2:
            for query in filtrated2[query2].keys():
                f.write(str(filtrated2[query2][query])+'\n')
    f.close()

    n=0
    for line in open('Filtrated_'+str(bins)+'_sr_blast.txt','r'):
        n+=1
        if n == 3:
            break

    if n >= 2:
        os.system('mv Filtrated_'+str(bins)+'_sr_blast.txt '+str(target_bin_folder)+'_filtrated_read_blast_output')
    else:
        os.system('rm Filtrated_'+str(bins)+'_sr_blast.txt')

    n=0
    for line in open('Filtrated_'+str(bins)+'_sr_blast_pe.txt','r'):
        n+=1
        if n == 3:
            break
    
    if n >= 2:
        os.system('mv Filtrated_'+str(bins)+'_sr_blast_pe.txt '+str(target_bin_folder)+'_filtrated_pe_read_blast_output')
    else:
        os.system('rm Filtrated_'+str(bins)+'_sr_blast_pe.txt')
            
    os.system('mv Filtrated1_* '+str(target_bin_folder)+'_all_filtrated_read_blast_output')
    os.system('mv '+str(bins)+'_sr_blast.txt '+str(target_bin_folder)+'_SR_blast_output')
    os.system('rm '+str(sc_bins2)+' '+str(sc_bins2)+'.nsq '+str(sc_bins2)+'.nin '+str(sc_bins2)+'.nhr')

def short_reads_blast(target_bin_folder, self_connected_bins_folder, short_reads_fa_folder, num_threads, pwd):
    os.chdir(pwd+'/'+str(self_connected_bins_folder))
    bin_sr, sc_bins = {}, {}
    for root, dirs, files in os.walk(pwd+'/'+str(self_connected_bins_folder)):
        for file in files:
            bin_id=str(file).split('_')[0]
            bin_sr[bin_id]=bin_id+'_sr.fa'
            sc_bins[bin_id]=file
    os.chdir(pwd)

    os.system('mkdir '+str(target_bin_folder)+'_filtrated_read_blast_output')
    os.system('mkdir '+str(target_bin_folder)+'_filtrated_pe_read_blast_output')
    os.system('mkdir '+str(target_bin_folder)+'_all_filtrated_read_blast_output')
    os.system('mkdir '+str(target_bin_folder)+'_SR_blast_output')
    # for bins in bin_sr.keys():
    #     sr_blast(bins, sc_bins, bin_sr, self_connected_bins_folder, short_reads_fa_folder, pwd)
    num_projects=int(num_threads)
    pool=Pool(processes=num_projects)
    for bins in bin_sr.keys():
        # sr_blast(bins, sc_bins, bin_sr, self_connected_bins_folder, short_reads_fa_folder, pwd)
        pool.apply_async(sr_blast, args=(target_bin_folder, bins, sc_bins, bin_sr, self_connected_bins_folder, short_reads_fa_folder, pwd))
    pool.close()
    pool.join()

def lr_blast(bins, sc_bins, bin_lr, self_connected_bins_folder, long_reads_folder, pwd):
    bins_lr_fa=bin_lr[bins]
    sc_bins2=sc_bins[bins]
    print('BLASTing '+str(bins)+' against long-reads')
    os.system('cp '+pwd+'/'+str(self_connected_bins_folder)+'/'+str(sc_bins2)+' '+str(sc_bins2))
    bin_contig_len={}
    for record in SeqIO.parse(str(sc_bins2),'fasta'):
        bin_contig_len[record.id]=len(record.seq)

    lr_len={}
    for record in SeqIO.parse(pwd+'/'+str(long_reads_folder)+'/'+str(bins_lr_fa),'fasta'):
        lr_len[record.id]=len(record.seq)
        
    os.system('makeblastdb -in '+str(sc_bins2)+' -dbtype nucl -logfile makeblastdb_log.txt')
    os.system('blastn -query '+pwd+'/'+str(long_reads_folder)+'/'+str(bins_lr_fa)+' -db '+str(sc_bins2)+' -num_threads 1 -evalue 1e-10 -outfmt 6 -out '+str(bins)+'_lr_blast.txt -max_target_seqs 5')
    os.system('rm '+str(sc_bins2)+'.nin '+str(sc_bins2)+'.nhr '+str(sc_bins2)+'.nsq')

    ali, ali2={}, {}
    seq_map_road, seq_map_road2 = {}, {}
    for line in open(str(bins)+'_lr_blast.txt','r'):
        lr_id=str(line).strip().split('\t')[0]
        bin_contigs_id=str(line).strip().split('\t')[1]
        simi=float(str(line).strip().split('\t')[2])
        leng=int(str(line).strip().split('\t')[3])
        lr_id_s=int(str(line).strip().split('\t')[6])
        lr_id_e=int(str(line).strip().split('\t')[7])
        bin_contigs_id_s=int(str(line).strip().split('\t')[8])
        bin_contigs_id_e=int(str(line).strip().split('\t')[9])
        evalue=int(str(line).strip().split('\t')[-1])
        if simi >= 85 and leng >= 200:
            lr_len_2=lr_len[lr_id]-50
            if lr_id_s <= 50 or lr_id_e >= lr_len_2:
                bin_contig_len2=bin_contig_len[bin_contigs_id]-50

                if bin_contigs_id_e > bin_contigs_id_s:
                    bin_contigs_id_s2 = bin_contigs_id_s
                    bin_contigs_id_e2 = bin_contigs_id_e
                else:
                    bin_contigs_id_s2 = bin_contigs_id_e
                    bin_contigs_id_e2 = bin_contigs_id_s

                if bin_contigs_id_s2 <= 50 or bin_contigs_id_e2 >= bin_contig_len2:
                    bin_contigs_id_max_len=bin_contig_len[bin_contigs_id]
                    lr_max_len=lr_len[lr_id]

                    # if lr_id in ali2.keys():
                    try:
                        if bin_contigs_id in ali2[lr_id].keys():
                            evalue2=ali2[lr_id][bin_contigs_id]
                            if evalue > evalue2:
                                ali2[lr_id][bin_contigs_id]=evalue
                                ali[lr_id][bin_contigs_id]=line.strip()
                        else:
                            ali2[lr_id][bin_contigs_id]=evalue
                            ali[lr_id][bin_contigs_id]=line.strip()
                    # else:
                    except:
                        ali2[lr_id]={}
                        ali[lr_id]={}
                        ali2[lr_id][bin_contigs_id]=evalue
                        ali[lr_id][bin_contigs_id]=line.strip()
            else:
                aligned_contig_len=abs(bin_contigs_id_e-bin_contigs_id_s)
                bin_contig_len2=bin_contig_len[bin_contigs_id]-100
                if aligned_contig_len >= bin_contig_len2:
                    if lr_id not in seq_map_road.keys():
                        seq_map_road[lr_id]={}
                        seq_map_road2[lr_id]={}
                        seq_map_road[lr_id]['media']=str(lr_id_s)+'|'+str(lr_id_e)+';'+str(bin_contigs_id)+':'+str(bin_contigs_id_s)+'|'+str(bin_contigs_id_e)
                        seq_map_road2[lr_id]['media']=evalue
                    else:
                        if 'media' in seq_map_road2[lr_id].keys():
                            # if evalue > seq_map_road2[lr_id]['media']:
                            seq_map_road[lr_id]['media']+='|$|'+str(lr_id_s)+'|'+str(lr_id_e)+';'+str(bin_contigs_id)+':'+str(bin_contigs_id_s)+'|'+str(bin_contigs_id_e)
                            # seq_map_road2[lr_id]['media']=evalue
                        else:
                            seq_map_road[lr_id]['media']=str(lr_id_s)+'|'+str(lr_id_e)+';'+str(bin_contigs_id)+':'+str(bin_contigs_id_s)+'|'+str(bin_contigs_id_e)
                            # seq_map_road2[lr_id]['media']=evalue
    # os.system('rm '+str(bins)+'_lr_blast.txt')

    f=open('Filtrated_'+str(bins)+'_lr_blast.txt','w')
    f1=open('Roadmap_'+str(bins)+'_lr_blast.txt','w')
    n = 0
    for lr_id in ali.keys():
        if len(ali[lr_id]) >= 2:
            if lr_id not in seq_map_road.keys():
                seq_map_road[lr_id]={}
                seq_map_road2[lr_id]={}
            blast_output, start, end, total, x, record = {}, 0, 0, 0, 0, {}
            lr_len_2=lr_len[lr_id]-50
            for bin_contigs_id in ali[lr_id].keys():
                bop=str(ali[lr_id][bin_contigs_id])
                lr_s=int(bop.split('\t')[6])
                lr_e=int(bop.split('\t')[7])
                contig_s=int(bop.split('\t')[8])
                contig_e=int(bop.split('\t')[9])
                evalue=int(bop.split('\t')[-1])
                # seq_map_road[lr_id]['lr||'+str(lr_id)+'']
                if lr_e >= lr_len_2 and lr_s <= 50:
                    xyz=0
                elif lr_e >= lr_len_2 or lr_s <= 50:
                    x+=1
                    if x == 1:
                        record[bop]=x
                        blast_output[x]={}
                        blast_output[x]['end']=lr_e
                        blast_output[x]['start']=lr_s
                        if lr_e >= lr_len_2:
                            seq_map_road[lr_id]['end']=str(lr_s)+'|'+str(lr_e)+';'+str(bin_contigs_id)+':'+str(contig_s)+'|'+str(contig_e)
                            seq_map_road2[lr_id]['end']=evalue
                        elif lr_s <= 50:
                            seq_map_road[lr_id]['start']=str(lr_s)+'|'+str(lr_e)+';'+str(bin_contigs_id)+':'+str(contig_s)+'|'+str(contig_e)
                            seq_map_road2[lr_id]['start']=evalue
                        # print(str(bop))
                    elif x >= 2:
                        y=0
                        # print(str(bop)) 
                        for x2 in blast_output.keys():
                            s=blast_output[x2]['start']
                            e=blast_output[x2]['end']
                            if s <= 50:
                                delta = e - lr_s
                                # print(str(delta))
                                if delta >= 0:
                                    y = 1
                            elif e >= lr_len_2:
                                delta = lr_e - s
                                # print(str(delta)) 
                                if delta >= 0:
                                    y = 1
                        # print('Y: '+str(y))
                        # print(str(len(record))+' seqs')
                        if y == 0:
                            record[bop]=x
                            if lr_e >= lr_len_2:
                                if 'end' in seq_map_road2[lr_id].keys():
                                    if evalue > seq_map_road2[lr_id]['end']:
                                        seq_map_road[lr_id]['end']=str(lr_s)+'|'+str(lr_e)+';'+str(bin_contigs_id)+':'+str(contig_s)+'|'+str(contig_e)
                                        seq_map_road2[lr_id]['end']=evalue
                                else:
                                    seq_map_road[lr_id]['end']=str(lr_s)+'|'+str(lr_e)+';'+str(bin_contigs_id)+':'+str(contig_s)+'|'+str(contig_e)
                                    seq_map_road2[lr_id]['end']=evalue
                            elif lr_s <= 50:
                                if 'start' in seq_map_road2[lr_id].keys():
                                    if evalue > seq_map_road2[lr_id]['start']:
                                        seq_map_road[lr_id]['start']=str(lr_s)+'|'+str(lr_e)+';'+str(bin_contigs_id)+':'+str(contig_s)+'|'+str(contig_e)
                                        seq_map_road2[lr_id]['start']=evalue
                                else:
                                    seq_map_road[lr_id]['start']=str(lr_s)+'|'+str(lr_e)+';'+str(bin_contigs_id)+':'+str(contig_s)+'|'+str(contig_e)
                                    seq_map_road2[lr_id]['start']=evalue
                        # print(str(len(record))+' seqs')
                        blast_output[x]={}
                        blast_output[x]['end']=lr_e
                        blast_output[x]['start']=lr_s
                # else:
                #     record[bop]=x

            # print('F '+str(len(record))+' seqs')
            # print(str(record))
            if len(record) >= 2:
                # print(str(record))
                for bop in record.keys():
                    f.write(bop.strip()+'\n')
                    n+=1
    f.close()

    n1=0
    for lr_id in seq_map_road.keys():
        for position in seq_map_road[lr_id].keys():
            n1+=1
            f1.write(str(lr_id)+'\t'+str(position)+'\t'+str(seq_map_road[lr_id][position])+'\n')
    f1.close()

    # if n >= 2:
    #     group, x= {}, 0
    #     for line in open('Filtrated_'+str(bins)+'_lr_blast.txt','r'):
    #         n+=1
    #         lr_id=str(line).strip().split('\t')[0]
    #         contig_id=str(line).strip().split('\t')[1]

    if n < 2:
        os.system('rm Filtrated_'+str(bins)+'_lr_blast.txt')
    else:
        os.system('mv Filtrated_'+str(bins)+'_lr_blast.txt '+pwd+'/'+str(target_bin_folder)+'_LR_blast_output')
    
    try:
        os.system('mv '+str(bins)+'_lr_blast.txt '+pwd+'/'+str(target_bin_folder)+'_LR_blast_output')
    except:
        xyz=0
    
    if n1 >= 1:
        os.system('mv Roadmap_'+str(bins)+'_lr_blast.txt '+pwd+'/'+str(target_bin_folder)+'_LR_roadmap')
    else:
        os.system('rm Roadmap_'+str(bins)+'_lr_blast.txt')
    os.system('rm '+sc_bins2)

def seq_merge_lr(contig1, contig2, bin_seq_total, bin_id, lr_id):
    contig1_seq=bin_seq_total[bin_id][contig1]
    contig2_seq=bin_seq_total[bin_id][contig2]
    ### Writing seqs
    contig_len={}
    f_contig1=open(bin_id+'_'+contig1+'.fa','w')
    f_contig1.write('>'+str(contig1)+'\n'+str(contig1_seq)+'\n')
    contig_len[str(contig1)]=len(contig1_seq)
    print(str(contig1)+' seqs len is '+str(len(contig1_seq)))
    f_contig1.close()
    f_contig2=open(bin_id+'_'+contig2+'.fa','w')
    f_contig2.write('>'+str(contig2)+'\n'+str(contig2_seq)+'\n')
    contig_len[str(contig2)]=len(contig2_seq)
    print(str(contig2)+' seqs len is '+str(len(contig2_seq)))
    f_contig2.close()
    total_len=len(contig1_seq)+len(contig2_seq)
    ### BLAST
    os.system('makeblastdb -in '+bin_id+'_'+contig2+'.fa -dbtype nucl -logfile makeblastdb_log.txt')
    os.system('blastn -query '+bin_id+'_'+contig1+'.fa -db '+bin_id+'_'+contig2+'.fa -outfmt 6 -evalue 1e-5 -num_threads 1 -max_target_seqs 5 -out '+bin_id+'_'+str(contig1)+'_VS_'+str(contig2)+'.txt')
    evalue = 0
    for line in open(bin_id+'_'+str(contig1)+'_VS_'+str(contig2)+'.txt','r'):
        evalue2=float(str(line).strip().split('\t')[-1])
        if evalue2 > evalue:
            filtrated_connection=line.strip()
            evalue=evalue2

    os.system('rm '+bin_id+'_'+contig2+'.fa.nin '+bin_id+'_'+contig2+'.fa.nsq '+bin_id+'_'+contig2+'.fa.nhr '+bin_id+'_'+str(contig1)+'_VS_'+str(contig2)+'.txt')
    subject_s=int(filtrated_connection.split('\t')[8])
    subject_e=int(filtrated_connection.split('\t')[9])
    if subject_s > subject_e:
        status_r=1
    else:
        status_r=0
    
    ### Merge seq
    if status_r == 1:
        contig2_seq2=contig2_seq[::-1]

        if generic_dna_t == 1:
            contig2_seq3=Seq(str(contig2_seq2), generic_dna).complement()
        else:
            contig2_seq3=Seq(str(contig2_seq2)).complement()
        
        print('Reverse or complementary merge')
        f=open('CR_'+bin_id+'_'+contig2+'.fa','w')
        f.write('>'+contig2+'\n'+str(contig2_seq3)+'\n')
        contig_len=len(contig2_seq)
        f.close()
        os.system('makeblastdb -in CR_'+bin_id+'_'+contig2+'.fa -dbtype nucl -logfile makeblastdb_log.txt')
        os.system('blastn -query '+bin_id+'_'+contig1+'.fa -db CR_'+bin_id+'_'+contig2+'.fa -outfmt 6 -evalue 1e-5 -num_threads 1 -max_target_seqs 5 -out '+bin_id+'_'+str(contig1)+'_VS_CR_'+str(contig2)+'.txt')
        evalue = 0
        for line in open(bin_id+'_'+str(contig1)+'_VS_CR_'+str(contig2)+'.txt','r'):
            evalue2=float(str(line).strip().split('\t')[-1])
            if evalue2 > evalue:
                filtrated_connection=line.strip()
                evalue=evalue2
        os.system('rm CR_'+bin_id+'_'+contig2+'.fa.nin CR_'+bin_id+'_'+contig2+'.fa.nsq CR_'+bin_id+'_'+contig2+'.fa.nhr '+bin_id+'_'+str(contig1)+'_VS_CR_'+str(contig2)+'.txt CR_'+bin_id+'_'+contig2+'.fa')
    else:
        contig2_seq3=contig2_seq
    
    print('BLAST output: '+str(filtrated_connection))
    query_s=int(filtrated_connection.split('\t')[6])
    query_e=int(filtrated_connection.split('\t')[7])
    subject_s=int(filtrated_connection.split('\t')[8])
    subject_e=int(filtrated_connection.split('\t')[9])
    delta=len(contig1_seq)-query_e
    delta2=len(contig2_seq3)-subject_e
    if query_s <= 50 and delta2 <= 50:
        merged_seq=contig2_seq3+contig1_seq[query_e:]
        print('Merged '+str(contig2)+' with '+str(contig1)+' ['+str(query_e)+':]')
    elif delta <= 50 and subject_s <= 50:
        merged_seq=contig1_seq[:query_s]+contig2_seq3
        print('Merged '+str(contig1)+' [:'+str(query_s)+'] with '+str(contig2))
    else:
        merged_seq=str(contig1_seq).replace(str(contig1_seq[query_s-1:query_e]), str(contig2_seq3))
        print('Replace '+str(contig1)+' ['+str(query_s-1)+':'+str(query_e)+'] with '+str(contig2))

    if len(merged_seq) <= total_len:
        try:
            new_id=contig1.split('_')[0]+'-'+contig2.split('_')[1]
        except:
            new_id=contig1.strip()+'-'+contig2.split('_')[1]
        extra_seq, merged_contigs = {}, {}
        f=open(bin_id+'_'+str(new_id)+'.fa','w')
        f.write('>'+new_id+'\n'+str(merged_seq)+'\n')
        extra_seq[new_id]=merged_seq
        f.close()
        merged_contigs[contig2]=0
        os.system('rm '+bin_id+'_'+contig1+'.fa '+bin_id+'_'+contig2+'.fa')
    else:
        print('Error total length of contig1 and contig2 merged seq ('+str(len(merged_seq))+'bp) is longer than total length of contig1 and contig2 ('+str(total_len)+'bp)')

    return bin_id+'_'+str(new_id)+'.fa', new_id, extra_seq, merged_contigs, merged_seq

def seq_merge_lr_main(select_roadmap_order, bin_id, bin_seq_org, bin_seq_total, target_bin_folder, pwd):
    merged_contigs={}
    bin_new_seq={}
    # n += 1
    # if n <= 5:
    for lr_id in select_roadmap_order[bin_id]:
        num=len(select_roadmap_order[bin_id][lr_id])
        for i in range(0, num-1):
            if i == 0:
                contig1=select_roadmap_order[bin_id][lr_id][i]
                contig2=select_roadmap_order[bin_id][lr_id][i+1]
                print(str(bin_id)+' group '+str(lr_id)+' run'+str(i+1)+' merging. Merging '+str(contig1)+' and '+str(contig2))
                A=seq_merge_lr(contig1, contig2, bin_seq_total, bin_id, lr_id)
                new_merged_seq=A[0]
                new_id=A[1]
                bin_seq_total[bin_id].update(A[2])
                merged_contigs.update(A[3])
            else:
                contig2=select_roadmap_order[bin_id][lr_id][i+1]
                print(str(bin_id)+' group '+str(lr_id)+' run'+str(i+1)+' merging. Merging '+str(contig1)+' and '+str(new_id))
                A=seq_merge_lr(new_id, contig2, bin_seq_total, bin_id, lr_id)
                new_merged_seq=A[0]
                new_id=A[1]
                bin_seq_total[bin_id].update(A[2])
                merged_contigs.update(A[3])
                merged_seq=A[4]
        print(str(bin_id)+' group '+str(lr_id)+' merging: '+str(new_merged_seq))
        bin_new_seq[new_id]=merged_seq
        os.system('rm '+str(new_merged_seq))

    print('Forming new bin of '+str(bin_id))
    f=open(str(bin_id)+'_lr_elongated.fa','w')
    f2=open(str(bin_id)+'_lr_being_merged_contigs.txt','w')
    # for contig in bin_seq_org[bin_id].keys():
    #     if contig not in merged_contigs.keys():
    #         f.write('>'+str(contig)+'\n'+str(bin_seq_org[bin_id][contig])+'\n')
    for contig in merged_contigs.keys():
        f2.write(str(contig)+'\n')
    f2.close()

    for new_id in bin_new_seq.keys():
        f.write('>'+str(new_id)+'\n'+str(bin_new_seq[new_id])+'\n')
    f.close()
    os.system('mv '+str(bin_id)+'_lr_elongated.fa '+pwd+'/'+str(target_bin_folder)+'_LR_elongation')
    os.system('mv '+str(bin_id)+'_lr_being_merged_contigs.txt '+pwd+'/'+str(target_bin_folder)+'_LR_elongation')

def long_reads_merge(target_bin_folder, last_connected_bins_folder, long_reads_folder, num_threads, pwd):
    os.chdir(pwd+'/'+str(last_connected_bins_folder))
    bin_lr, sc_bins = {}, {}
    for root, dirs, files in os.walk(pwd+'/'+str(last_connected_bins_folder)):
        for file in files:
            bin_id=str(file).split('_')[0]
            bin_lr[bin_id]=bin_id+'_lr_polished.fa'
            sc_bins[bin_id]=file
    os.chdir(pwd)

    os.system('mkdir '+str(target_bin_folder)+'_filtrated_long_read_blast_output')
    # os.system('mkdir '+str(target_bin_folder)+'_all_filtrated_read_blast_output')
    os.system('mkdir '+str(target_bin_folder)+'_LR_blast_output')
    os.system('mkdir '+str(target_bin_folder)+'_LR_roadmap')
    num_projects=int(num_threads)
    pool=Pool(processes=num_projects)
    for bins in bin_lr.keys():
        # lr_blast(bins, sc_bins, bin_lr, last_connected_bins_folder, long_reads_folder, pwd)
        pool.apply_async(lr_blast, args=(bins, sc_bins, bin_lr, last_connected_bins_folder, long_reads_folder, pwd))
    pool.close()
    pool.join()

    os.chdir(str(target_bin_folder)+'_LR_roadmap')
    roadmap={}
    for root, dirs, files in os.walk(pwd+'/'+str(target_bin_folder)+'_LR_roadmap'):
        for file in files:
            bin_id=str(file).split('_')[1]
            roadmap[bin_id]={}
            for line in open(file,'r'):
                lr_id=str(line).strip().split('\t')[0]
                position=str(line).strip().split('\t')[1]
                rm=str(line).strip().split('\t')[2]
                if lr_id not in roadmap[bin_id].keys():
                    roadmap[bin_id][lr_id]={}
                    roadmap[bin_id][lr_id][position]=rm
                else:
                    roadmap[bin_id][lr_id][position]=rm
    os.chdir(str(pwd))

    roadmap2={}
    for bin_id in roadmap.keys():
        roadmap2[bin_id]={}
        for lr_id in roadmap[bin_id].keys():
            roadmap2[bin_id][lr_id]=''
            if 'start' in roadmap[bin_id][lr_id].keys():
                rm=roadmap[bin_id][lr_id]['start']
                contig_info=rm.strip().split(';')[1]
                contig_id=contig_info.split(':')[0]
                contig_s=contig_info.split(':')[1].split('|')[0]
                contig_e=contig_info.split(':')[1].split('|')[1]
                if contig_e > contig_s:
                    contig_info2=contig_id+'[0:]'
                else:
                    contig_info2=contig_id+'[r:0]'
                lr_info=rm.strip().split(';')[0]
                # lr_s=lr_info.split('|')[0]
                lr_e=lr_info.split('|')[1]
                # contig_id_seq=rm.strip().split(';')[1].split(':')[1]
                roadmap2[bin_id][lr_id]='$|S|'+contig_info2+'|'+lr_id+'['+lr_e+':'
            if 'media' in roadmap[bin_id][lr_id].keys():
                rm=roadmap[bin_id][lr_id]['media']
                media_list=rm.strip().split('|$|')
                seq_sort_d, seq_sort_l = {}, []
                for item in media_list:
                    contig_info=item.split(';')[1]
                    contig_id=contig_info.split(':')[0]
                    contig_s=contig_info.split(':')[1].split('|')[0]
                    contig_e=contig_info.split(':')[1].split('|')[1]
                    lr_info=item.split(';')[0]
                    lr_s=lr_info.split('|')[0]
                    lr_e=lr_info.split('|')[1]
                    seq_sort_l.append(lr_s)
                    seq_sort_d[lr_s]=str(lr_e)
                seq_sort_l.sort()
                for lr_s in seq_sort_l:
                    if contig_e > contig_s:
                        roadmap2[bin_id][lr_id]+=str(lr_s)+']$|M|'+contig_id+'[0:]'+'|'+lr_id+'['+str(seq_sort_d[lr_s])+':'
                    else:
                        roadmap2[bin_id][lr_id]+=str(lr_s)+']$|M|'+contig_id+'[r:0]'+'|'+lr_id+'['+str(seq_sort_d[lr_s])+':'

            if 'end' in roadmap[bin_id][lr_id].keys():
                rm=roadmap[bin_id][lr_id]['end']
                contig_info=rm.strip().split(';')[1]
                contig_id=contig_info.split(':')[0]
                contig_s=contig_info.split(':')[1].split('|')[0]
                contig_e=contig_info.split(':')[1].split('|')[1]
                lr_info=rm.strip().split(';')[0]
                lr_s=lr_info.split('|')[0]
                # lr_e=lr_info.split('|')[1]
                if contig_e > contig_s:
                    roadmap2[bin_id][lr_id]+=str(lr_s)+']$|E|'+contig_id+'[0:]'
                else:
                    roadmap2[bin_id][lr_id]+=str(lr_s)+']$|E|'+contig_id+'[r:0]'

    f=open('LR_total_roadmap.txt','w')
    f2=open('LR_total_roadmap2.txt','w')
    select_roadmap={}
    for bin_id in roadmap2.keys():
        select_roadmap[bin_id]={}
        for lr_id in roadmap2[bin_id].keys():
            fz=roadmap2[bin_id][lr_id][-1]
            num=roadmap2[bin_id][lr_id].count('$|')
            if num >= 2:
                # select_roadmap[bin_id][lr_id]=0
                if fz != ']':
                    f.write(str(bin_id)+'\t'+str(lr_id)+'\t'+str(roadmap2[bin_id][lr_id])+']'+'\n')
                    f2.write(str(bin_id)+'\t'+str(lr_id)+'\t'+str(roadmap2[bin_id][lr_id])+']'+'\n')
                    select_roadmap[bin_id][lr_id]=str(roadmap2[bin_id][lr_id])+']'
                else:
                    f.write(str(bin_id)+'\t'+str(lr_id)+'\t'+str(roadmap2[bin_id][lr_id])+'\n')
                    f2.write(str(bin_id)+'\t'+str(lr_id)+'\t'+str(roadmap2[bin_id][lr_id])+'\n')
                    select_roadmap[bin_id][lr_id]=str(roadmap2[bin_id][lr_id])
            else:
                if fz != ']':
                    f.write(str(bin_id)+'\t'+str(lr_id)+'\t'+str(roadmap2[bin_id][lr_id])+']'+'\n')
                else:
                    f.write(str(bin_id)+'\t'+str(lr_id)+'\t'+str(roadmap2[bin_id][lr_id])+'\n')
    f.close()
    f2.close()
    
    select_roadmap2=copy.deepcopy(select_roadmap)
    for bin_id in select_roadmap2.keys():
        if len(select_roadmap2[bin_id]) == 0:
            del select_roadmap[bin_id]

    select_roadmap_order, select_roadmap_order2={},{}
    for bin_id in select_roadmap.keys():
        select_roadmap_order[bin_id]={}
        select_roadmap_order2[bin_id]={}
        for lr_id in select_roadmap[bin_id].keys():
            select_roadmap_order[bin_id][lr_id]=[lr_id]
            select_roadmap_order2[bin_id][lr_id]=0
            road=select_roadmap[bin_id][lr_id]
            contig_list=road.split('$|')
            contig_list.remove(contig_list[0])
            for item in contig_list:
                contigs=item[2:]
                contig_list=contigs.split('|')
                for contig in contig_list:
                    contig2=contig.split('[')[0]
                    if contig2 not in select_roadmap_order[bin_id][lr_id]:
                        select_roadmap_order[bin_id][lr_id].append(contig2)
                        select_roadmap_order2[bin_id][contig2]=0

    print('Recording potential merge seqs')
    bin_seq_total, bin_seq_org = {}, {}
    for bin_id in select_roadmap_order.keys():
        bin_seq_total[bin_id]={}
        bin_seq_org[bin_id]={}

    for bin_id in select_roadmap_order.keys():
        os.chdir(pwd+'/'+str(last_connected_bins_folder))
        for root, dirs, files in os.walk(pwd+'/'+str(last_connected_bins_folder)):
            for file in files:
                hz=file.split('.')[-1]
                bin_id2=bin_id+'_'
                bin_id3=bin_id+'.fa'
                if bin_id2 in file and 'fa' in hz:
                    for record in SeqIO.parse(file,'fasta'):
                        bin_seq_org[bin_id][record.id]=record.seq
                        # if record.id in select_roadmap_order2[bin_id]:
                        try:
                            select_roadmap_order2[bin_id][record.id]+=1
                            bin_seq_total[bin_id][record.id]=record.seq
                        except:
                            xyz=0
                elif bin_id3 in file:
                    for record in SeqIO.parse(file,'fasta'):
                        bin_seq_org[bin_id][record.id]=record.seq
                        # if record.id in select_roadmap_order2[bin_id]:
                        try:
                            select_roadmap_order2[bin_id][record.id]+=1
                            bin_seq_total[bin_id][record.id]=record.seq
                        except:
                            xyz=0
                    
        os.chdir(pwd+'/'+str(long_reads_folder))
        for root, dirs, files in os.walk(pwd+'/'+str(long_reads_folder)):
            for file in files:
                hz=file.split('.')[-1]
                bin_id2=bin_id+'_'
                if bin_id2 in file and 'fa' in hz:
                    for record in SeqIO.parse(file,'fasta'):
                        try:
                        # if record.id in select_roadmap_order2[bin_id]:
                            select_roadmap_order2[bin_id][record.id]+=1
                            bin_seq_total[bin_id][record.id]=record.seq
                        except:
                            xyz=0
        os.chdir(pwd)

    os.chdir(pwd)
    os.system('mkdir '+str(target_bin_folder)+'_LR_elongation')
    print('Recorded potential merge seqs. Start to merge seqs')
    # bin_seq, record_seq, merged_contigs, n ={}, [], {}, 0
    num_projects=int(num_threads)
    pool=Pool(processes=num_projects)
    for bin_id in select_roadmap_order.keys():
        # seq_merge_lr_main(select_roadmap_order, bin_id, bin_seq_org, bin_seq_total, target_bin_folder, pwd)
        pool.apply_async(seq_merge_lr_main, args=(select_roadmap_order, bin_id, bin_seq_org, bin_seq_total, target_bin_folder, pwd))
    pool.close()
    pool.join()
    os.system('mv LR_total_roadmap2.txt LR_total_roadmap.txt '+pwd+'/'+str(target_bin_folder)+'_LR_elongation')
    os.system('rm makeblastdb_log.txt')

    os.chdir(pwd+'/'+str(target_bin_folder)+'_LR_elongation')
    f=open('Total_been_merged_contigs.txt','w')
    for root, dirs, files in os.walk(pwd+'/'+str(target_bin_folder)+'_LR_elongation'):
        for file in files:
            if '_lr_being_merged_contigs.txt' in file:
                bin_id=file.split('_lr_being_merged_contigs.txt')[0]
                for line in open(file, 'r'):
                    f.write(bin_id+'\t'+str(line))
    f.close()
    os.chdir(pwd)    

def read_connecting(bin_id, bin_containing_folder, short_read_fa_folder, pwd, num_threads):
    total_seq, read_seq={}, {}
    try:
        for record in SeqIO.parse(pwd+'/'+bin_containing_folder+'/'+str(bin_id)+'_mag_polished_sc.fa','fasta'):
            total_seq[record.id]=record.seq
    except:
        x3=0

    try:
        for record in SeqIO.parse(pwd+'/'+bin_containing_folder+'/'+str(bin_id)+'_mag_polished.fa','fasta'):
            total_seq[record.id]=record.seq
    except:
        x3=0

    # seq_num1=len(total_seq)

    try:
        for record in SeqIO.parse(pwd+'/'+str(short_read_fa_folder)+'/'+str(bin_id)+'_sr.fa','fasta'):
            total_seq[str(record.id).replace('/','_')]=record.seq
            read_seq[str(record.id).replace('/','_')]=0
    except:
        x4=0

    group_blast, group, num1, num2, n = {}, {}, 0, 1, 0
    while num1 != num2:
        num1=len(group_blast)
        for line in open (pwd+'/'+str(target_bin_folder)+'_reads_merged/'+bin_id+'_filtated_SR_blast_de-rep_single_path.txt', 'r'):
            n+=1
            query=str(line).strip().split('\t')[0].replace('/','_')
            subject=str(line).strip().split('\t')[1]
            if n == 1:
                group_blast[subject]={}
                group_blast[subject][line.strip()]=0
                group[subject]={}
                group[subject][subject]=0
                group[subject][query]=0
            else:
                try:
                    group_blast[subject][line.strip()]=0
                    group[subject][query]=0
                except:
                    x=10
                    group2=copy.deepcopy(group)
                    for item in group2.keys():
                        if query in group2[item].keys():
                            group_blast[item][line.strip()]=0
                            group[item][subject]=0
                            group[item][query]=0
                            x=1
                    if x == 10:
                        group_blast[subject]={}
                        group_blast[subject][line.strip()]=0
                        group[subject]={}
                        group[subject][subject]=0
                        group[subject][query]=0
        
        numx1=len(group_blast)
        numx2=0
        while numx1 != numx2:
            numx1=len(group_blast)
            group2=copy.deepcopy(group)
            group3=copy.deepcopy(group)
            for sub2 in group2.keys():
                for sub3 in group3.keys():
                    if sub2 != sub3:
                        if sub2 in group3[sub3].keys():
                            try:
                                group_blast[sub3].update(group_blast[sub2])
                                group[sub3].update(group[sub2])
                                del group[sub2]
                                del group_blast[sub2]
                            except:
                                print(str(sub2)+' has been removed')
                        else:
                            for item in group2[sub2].keys():
                                que=item.split('\t')[0]
                                if que in group3[sub3].keys():
                                    try:
                                        group_blast[sub3].update(group_blast[sub2])
                                        group[sub3].update(group[sub2])
                                        del group[sub2]
                                        del group_blast[sub2]
                                    except:
                                        print(str(sub2)+' has been removed') 
            numx2=len(group_blast)
        num2=len(group_blast)

    num_seq, processed_seq_num, processed_contigs, merged_seq, merged_file_name = 0, 0, {}, [], []
    for group_id in group_blast.keys():
        print('Processing group '+str(group_id))
        contigs = {}
        for item in group_blast[group_id].keys():
            query=item.strip().split('\t')[0].replace('/','_')
            subject=item.strip().split('\t')[1]
            contigs[query]=0
            contigs[subject]=0

        if len(contigs) != 0:
            num, processed_item, xyok, processed_contigs2, new_contigs = 0, {}, 0, {}, {}
            for item in group_blast[group_id].keys():
                num+=1
                if num == 1:
                    num_seq+=1
                    query=str(item).strip().split('\t')[0].replace('/','_')
                    subject=str(item).strip().split('\t')[1]
                    query_start=int(str(item).strip().split('\t')[6])
                    query_end=int(str(item).strip().split('\t')[7])
                    subject_start=int(str(item).strip().split('\t')[8])
                    subject_end=int(str(item).strip().split('\t')[9])
                    sim=float(str(item).strip().split('\t')[2])
                    leng=int(str(item).strip().split('\t')[3])

                    A=seq_merge(total_seq, query, subject, query_start, query_end, subject_start, subject_end, num_seq, bin_id) 
                    processed_contigs[query]=1
                    processed_contigs[subject]=1
                    processed_contigs2[query]=1
                    processed_contigs2[subject]=1
                    merged_seq_temp=A[2]
                    merged_file_name_temp=A[3]
                    merged_seq_id=A[4]
                    contigs[merged_seq_id]=0
                    merged_file_name.append(merged_file_name_temp)
                    new_contigs[merged_seq_id]=0
                    total_seq.update(A[0])
                    processed_seq_num+=1

            num1=len(processed_contigs2)
            num2=0
            while num1 != num2 and xyok <= len(contigs):
                xyok+=1
                for contig in contigs.keys():
                    if contig not in processed_contigs2.keys():
                        temp_file_name=str(bin_id)+'_'+str(contig)+'.fa'
                        seq_id=contig

                        contigs_len = {}
                        for record in SeqIO.parse(str(merged_file_name_temp),'fasta'):
                            contigs_len[record.id]=len(record.seq)
                            # total_seq2[record.id]=record.seq

                        f=open(temp_file_name,'w')
                        if seq_id in total_seq.keys():
                            f.write('>'+str(seq_id)+'\n'+str(total_seq[seq_id])+'\n')
                            contigs_len[seq_id]=len(total_seq[seq_id])
                            # total_seq2[seq_id]=total_seq[seq_id]
                        f.close()
                        os.system('makeblastdb -dbtype nucl -in '+str(merged_file_name_temp)+' -logfile makeblastdb_log.txt')
                        os.system('blastn -db '+str(merged_file_name_temp)+' -query '+str(temp_file_name)+' -outfmt 6 -evalue 1e-3 -max_target_seqs 50 -out '+str(merged_file_name_temp)+'_blast.txt  -num_threads 1')
                        os.system('rm '+str(merged_file_name_temp)+'.nsq '+str(merged_file_name_temp)+'.nin '+str(merged_file_name_temp)+'.nhr')

                        fil2 = {}
                        for line in open(str(merged_file_name_temp)+'_blast.txt','r'):
                            sim=float(str(line).strip().split('\t')[2])
                            leng=int(str(line).strip().split('\t')[3])
                            if sim >= 99:
                                que=str(line).strip().split('\t')[0].replace('/','_')
                                sub=str(line).strip().split('\t')[1]
                                if que != sub:
                                    que_s=int(str(line).strip().split('\t')[6])
                                    que_e=int(str(line).strip().split('\t')[7])
                                    sub_s=int(str(line).strip().split('\t')[8])
                                    sub_e=int(str(line).strip().split('\t')[9])
                                    if que_s == 1 or que_e == contigs_len[que]:
                                        if sub_e > sub_s:
                                            if sub_s == 1 or sub_e == contigs_len[sub]:
                                                if que > sub:
                                                    fil2[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                                else:
                                                    fil2[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line

                                        else:
                                            if sub_e == 1 or sub_s == contigs_len[sub]:
                                                if que > sub:
                                                    fil2[que+'\t'+sub+'\t'+str(sim)+'\t'+str(leng)]=line
                                                else:
                                                    fil2[sub+'\t'+que+'\t'+str(sim)+'\t'+str(leng)]=line
                        # os.system('rm '+str(temp_file_name)+' '+str(temp_file_name)+'_'+str(merged_seq_temp)+'.txt')
                        for item in fil2.keys():
                            # query=item.split('\t')[1]
                            # subject=item.split('\t')[0]
                            line=fil2[item].strip()
                            query=str(line).strip().split('\t')[0].replace('/','_')
                            subject=str(line).strip().split('\t')[1]
                            query_start=int(str(line).strip().split('\t')[6])
                            query_end=int(str(line).strip().split('\t')[7])
                            subject_start=int(str(line).strip().split('\t')[8])
                            subject_end=int(str(line).strip().split('\t')[9])

                        A=seq_merge(total_seq, query, subject, query_start, query_end, subject_start, subject_end, num_seq, bin_id) 
                        processed_contigs[query]=1
                        processed_contigs[subject]=1
                        processed_contigs2[query]=1
                        processed_contigs2[subject]=1
                        merged_seq_temp=A[2]
                        # merged_seq.append(A[2])
                        merged_file_name_temp=A[3]
                        merged_seq_id=A[4]
                        # contigs[merged_seq_id]=0
                        merged_file_name.append(merged_file_name_temp)
                        total_seq.update(A[0])
                        os.system('rm '+str(temp_file_name)+' '+str(merged_file_name_temp)+'_blast.txt')
                        num2=len(processed_contigs2)
                        processed_seq_num+=1

        # os.chdir(pwd)
        os.system('mkdir '+str(bin_id)+'_read_merged_files')

        try:  
            for file_id in merged_file_name:
                os.system('mv '+str(file_id)+' '+str(bin_id)+'_read_merged_files')
        except:
            xyzzz=0

        no_containning={}
        # if processed_seq_num != 0:
        try:
            if len(processed_contigs) != 0:
                # fx=open(str(target_bin_folder)+'_bins_read_merged_summary.txt','a')
                # bin_id=str(bin_id).split('.fa')[0]
                bin_id2=bin_id+'_rm.fa'
                f=open(bin_id2,'w')
                x=0
                for seq_id in total_seq.keys():
                    if seq_id not in processed_contigs.keys() and seq_id not in read_seq.keys():
                        x+=1
                        f.write('>'+str(seq_id)+'\n'+str(total_seq[seq_id])+'\n')
                f.close()
                # fx.write(str(bin_id2)+'\t'+str(yy)+'\t'+str(x)+'\n')
                # fx.close()
                os.system('mv '+bin_id2+' '+str(target_bin_folder)+'_read_connected_bins')
                # os.system('rm '+str(bin_id))
            else:
                no_containning[bin_id]=0
                # os.system('mv '+str(bin_id)+' '+str(target_bin_folder)+'_read_connected_bins')
        except:
            # os.system('mv '+str(bin_id)+' '+str(target_bin_folder)+'_read_connected_bins')
            no_containning[bin_id]=0
    # return no_containning

def short_reads_merging(target_bin_folder, filtrate_blast_output_folder, short_read_fa_folder, bin_containing_folder, num_threads, pwd):
    xyz=0
    os.chdir(pwd+'/'+str(short_read_fa_folder))
    for root, dirs, files in os.walk(pwd+'/'+str(short_read_fa_folder)):
        for file in files:
            xyz+=1
            if xyz == 1:
                xyzzz=0
                for line in open(file,'r'):
                    xyzzz+=1
                    if xyzzz == 2:
                        read_len=len(str(line).strip())
                        break
    os.chdir(pwd)

    filtrated_seq={}
    os.chdir(pwd+'/'+str(filtrate_blast_output_folder))
    for root, dirs, files in os.walk(pwd+'/'+str(filtrate_blast_output_folder)):
        for file in files:
            bin_id=str(file).strip().split('Filtrated1_')[1].split('_sr_blast.txt')[0]
            # print(bin_id)
            filtrated_seq[bin_id]={}
            temp={}
            for line in open(file,'r'):
                query=str(line).strip().split('\t')[0]
                subject=str(line).strip().split('\t')[1]
                q_s=int(str(line).strip().split('\t')[6])
                q_e=int(str(line).strip().split('\t')[7])
                # if query not in filtrated_seq[bin_id].keys():
                #     filtrated_seq[bin_id][query]={}

                if q_s == 1 and q_e != read_len:
                    if query not in temp.keys():
                        temp[query]={}
                        temp[query][1]=1
                    else:
                        temp[query][1]=1
                elif q_e == read_len and q_s != 1:
                    if query not in temp.keys():
                        temp[query]={}
                        temp[query][read_len]=1
                    else:
                        temp[query][read_len]=1

            for item in temp.keys():
                if len(temp[item]) == 2:
                    filtrated_seq[bin_id][item]=[]
            
            for line in open(file,'r'):
                query=str(line).strip().split('\t')[0]
                if query in filtrated_seq[bin_id].keys():
                    filtrated_seq[bin_id][query].append(line)
            
    os.chdir(pwd)
    # fx=open(str(target_bin_folder)+'_bins_read_merged_summary.txt','w')
    os.system('mkdir '+str(target_bin_folder)+'_SR_blast_final_filtrtated')
    os.system('mkdir '+str(target_bin_folder)+'_reads_merged')

    bin_waiting_for_read_merging={}
    os.chdir(str(target_bin_folder)+'_SR_blast_final_filtrtated')
    for bin_id in filtrated_seq.keys():
        if len(filtrated_seq[bin_id]) != 0:
            filtrated_seq2 = {}
            f=open(bin_id+'_filtated_SR_blast.txt','w')
            f2=open(bin_id+'_filtated_SR_blast_de-rep.txt','w')
            de_rep, de_rep2, de_rep_score, de_rep3 = {}, {}, {}, {}
            for query in filtrated_seq[bin_id].keys():
                for item in filtrated_seq[bin_id][query]:
                    f.write(item)
                    subject=str(item).strip().split('\t')[1]
                    score=float(str(item).strip().split('\t')[-1])
                    try:
                        # de_rep[query]+='\t'+str(subject)
                        de_rep[query].append(str(subject))
                        de_rep_score[query]+=score
                    except:
                        de_rep[query]=[str(subject)]
                        de_rep_score[query]=score
            f.close()

            for item in de_rep.keys():
                ss=de_rep[item]
                ss.sort()
                try:
                    score1=de_rep2[str(ss)]
                    if de_rep_score[item] > score1:
                        de_rep2[str(ss)]=de_rep_score[item]
                        de_rep3[str(ss)]=item
                except:
                    de_rep2[str(ss)]=de_rep_score[item]
                    de_rep3[str(ss)]=item

            selected_query={}
            for ss in de_rep2.keys():
                selected_query[de_rep3[ss]]=0
                
            for query in filtrated_seq[bin_id].keys():
                if query in selected_query.keys():
                    for item in filtrated_seq[bin_id][query]:
                        f2.write(item)
                        filtrated_seq2[item]=0
            f2.close()

            subject, subject2, query, query2, total_que, total_sub = {}, {}, {}, {}, {}, {}
            for item in filtrated_seq2.keys():
                que=item.split('\t')[0]
                total_que[que]=0
                query[que]={}
                query2[que]={}
                sub=item.split('\t')[1]
                total_sub[sub]=0
                subject[sub]={}
                subject2[sub]={}

            for item in filtrated_seq2.keys():
                que=item.split('\t')[0]
                sub=item.split('\t')[1]
                q_s=int(item.split('\t')[6])
                q_e=int(item.split('\t')[7])
                s_s=int(item.split('\t')[8])
                s_e=int(item.split('\t')[9])
                if q_s == 1 or q_e == 1:
                    try:
                        query[que]['s']+=1
                        query2[que]['s'].append(sub)
                    except:
                        query[que]['s']=1
                        query2[que]['s']=[sub]
                elif q_s != 1 and q_e != 1:
                    try:
                        query[que]['e']+=1
                        query2[que]['e'].append(sub)
                    except:
                        query[que]['e']=1
                        query2[que]['e']=[sub]
                
                if s_s == 1 or s_e == 1:
                    try:
                        subject[sub]['s']+=1
                        subject2[sub]['s'].append(que)
                    except:
                        subject[sub]['s']=1
                        subject2[sub]['s']=[que]
                elif s_s != 1 and s_e != 1:
                    try:
                        subject[sub]['e']+=1
                        subject2[sub]['e'].append(que)
                    except:
                        subject[sub]['e']=1
                        subject2[sub]['e']=[que]

            del_que, del_sub, del_item, num1, num2, iteration = {}, {}, {}, 0, 1, 0
            while num1 != num2:
                iteration+=1
                num1=len(del_item)
                for item in filtrated_seq2.keys():
                    if item not in del_item.keys():
                        que=item.split('\t')[0]
                        sub=item.split('\t')[1]
                        q_s=int(item.split('\t')[6])
                        q_e=int(item.split('\t')[7])
                        s_s=int(item.split('\t')[8])
                        s_e=int(item.split('\t')[9])
                        
                        if s_s == 1 or s_e == 1:
                            if subject[sub]['s'] >= 2:
                                del_item[item]=0
                        elif s_s != 1 and s_e != 1:
                            if subject[sub]['e'] >= 2:
                                del_item[item]=0

                        if q_s == 1 or q_e == 1:
                            if query[que]['s'] >= 2:
                                del_item[item]=0
                        elif q_s != 1 and q_e != 1:
                            if query[que]['e'] >= 2:
                                del_item[item]=0

                remained_que={}
                for item in filtrated_seq2.keys():
                    if item not in del_item.keys():
                        que=item.split('\t')[0]
                        try:
                            remained_que[que]+=1
                        except:
                            remained_que[que]=1
                
                remained_item={}
                for item in filtrated_seq2.keys():
                    if item not in del_item.keys():
                        que=item.split('\t')[0]
                        if remained_que[que] == 2:
                            remained_item[item]=0
                num2=len(del_item)

            f=open(bin_id+'_filtated_SR_blast_de-rep_multi_paths.txt','w')
            f2=open(bin_id+'_filtated_SR_blast_de-rep_single_path.txt','w')
            for item in filtrated_seq2.keys():
                if item in remained_item.keys():
                    f2.write(item)
                else:
                    f.write(item)
            f.close()
            f2.close()

            xyz=0
            for line in open(bin_id+'_filtated_SR_blast_de-rep_single_path.txt','r'):
                xyz+=1
            
            if xyz >= 1:
                os.system('cp '+bin_id+'_filtated_SR_blast_de-rep_single_path.txt '+pwd+'/'+str(target_bin_folder)+'_reads_merged')
                bin_waiting_for_read_merging[bin_id]=0
    os.chdir(pwd)


    # for bins in bin_contigs_length.keys():
    #     print('Processing '+str(bins))
    #     # self_connecting(bins, bin_folder, short_read_fa_folder, bin_contigs_length, pwd, num_threads)
    #     pool.apply_async(self_connecting, args=(bins, bin_folder, bin_contigs_length, pwd, num_threads))


    # no_containning={}
    pool=Pool(processes=num_threads)
    for bin_id in bin_waiting_for_read_merging.keys():
        print('Processing contig-read-contig merging of: '+str(bin_id))
        # read_connecting(bin_id, bin_containing_folder, short_read_fa_folder, pwd, num_threads)
        pool.apply_async(read_connecting, args=(bin_id, bin_containing_folder, short_read_fa_folder, pwd, num_threads))
    pool.close()
    pool.join()

    os.chdir(pwd)
    bin_read_connected={}
    for root, dirs, files in os.walk(pwd+'/'+str(target_bin_folder)+'_read_connected_bins'):
        for file in files:
            bin_id=str(file).split('_')[0]
            bin_read_connected[bin_id]=0

    os.chdir(pwd+'/'+bin_containing_folder)
    for root, dirs, files in os.walk(pwd+'/'+str(bin_containing_folder)):
        for file in files:
            bin_id2=str(file).split('_')[0]
            if bin_id2 not in bin_read_connected.keys():
                os.system('cp '+str(file)+' '+pwd+'/'+str(target_bin_folder)+'_read_connected_bins')
    os.chdir(pwd)

def parse_checkm(checkm_containing_folder,pwd):
    #pwd=os.getcwd()
    bins_checkm, bin_name = {}, {}
    os.chdir(pwd+'/'+checkm_containing_folder+'/storage')
    for root, dirs, files in os.walk(pwd+'/'+checkm_containing_folder+'/storage'):
        for file in files:        
            if 'bin_stats_ext.tsv' in file: 
                for line in open(file,'r'):
                    bin_id=str(line).strip().split('_')[0]
                    bins_checkm[bin_id]={}

                    try:
                        marker_lineage=str(line).strip().split('\'marker lineage\': \'')[1].strip().split('\'')[0]
                        bins_checkm[bin_id]['marker lineage']=marker_lineage
                    except:
                        print('marker lineage error')
                        bins_checkm[bin_id]['marker lineage']=root

                    try:    
                        completeness=str(line).strip().split('\'Completeness\': ')[1].split(', ')[0]
                        bins_checkm[bin_id]['Completeness']=float(completeness)
                    except:
                        bins_checkm[bin_id]['Completeness']=0

                    try:
                        genome_size=str(line).strip().split('\'Genome size\':')[1].strip().split(', ')[0]
                        bins_checkm[bin_id]['Genome size']=int(genome_size)
                    except:
                        bins_checkm[bin_id]['Genome size']=0

                    bins_checkm[bin_id]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0])

                    try:
                        bins_checkm[bin_id]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split(',')[0].strip())
                    except:
                        bins_checkm[bin_id]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split('}')[0].strip())
    os.chdir(pwd)
    return bins_checkm

def bin_evaluation(target_bin_folder, pre_bin_folder, gap_filling_final_bin_folder, num_threads, pwd):
    # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(pre_bin_folder)+' '+str(pre_bin_folder)+'_checkm')
    bin_checkm1=parse_checkm(str(pre_bin_folder)+'_checkm', pwd)
    os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(gap_filling_final_bin_folder)+' '+str(gap_filling_final_bin_folder)+'_checkm')
    bin_checkm2=parse_checkm(str(gap_filling_final_bin_folder)+'_checkm', pwd)
    # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(target_bin_folder)+' '+str(target_bin_folder)+'_checkm')
    bin_checkm3=parse_checkm(str(target_bin_folder)+'_checkm', pwd)

    f=open('Comparison_test.txt','w')
    move_bin={}
    for bin_id in bin_checkm1.keys():
        # if bin_id not in bin_checkm2.keys():
        #     os.system('mv '+pwd+'/'+pre_bin_folder)
        if bin_id in bin_checkm2.keys() and bin_id in bin_checkm3.keys():
            f.write(str(bin_id)+'\n'+'Target_folder: '+str(bin_checkm3[bin_id])+'\n'+'SR_conncected_folder: '+str(bin_checkm1[bin_id])+'\n'+'LR_conncected_folder: '+str(bin_checkm2[bin_id])+'\n')
        else:
            move_bin[bin_id]=0
    f.close()

def parse_sam_bwa(sam_file, fq, pair, n, mp_run):
    print('Reading reads id')
    try:
        f_not_mapped_reads=open('Not_mapped_reads'+str(mp_run)+'.txt','a')
    except:
        f_not_mapped_reads=open('Not_mapped_reads'+str(mp_run)+'.txt','w')
    m, m1, m2 = 0, 0, 0
    for line in open(sam_file,'r'):
        m1+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id_o=flist[2]
            if bin_id_o != '*':
                bin_id=flist[2].split('_')[0]
                read_id=flist[0]
                read_id_name=str(n)+'_'+read_id
                try:
                    pair[bin_id][read_id_name]+=1
                except:
                    pair[bin_id][read_id_name]=1
        if m1 % 1000000 == 0:
            print('Read', m1,'lines')

    for bin_id in pair.keys():
        for read_id_name in pair[bin_id].keys():
            if pair[bin_id][read_id_name] == 2:
                fq[bin_id][read_id_name]=0

    pair={}
    try:
        f_summary=open('Bin_reads_summary'+str(mp_run)+'.txt','a')
    except:
        f_summary=open('Bin_reads_summary'+str(mp_run)+'.txt','w')
    for item in fq.keys():
        f_summary.write(str(item)+' SEQ number:'+str(len(fq[item]))+'\n')
    f_summary.close()

    print('Parsing', sam_file)
    for line in open(sam_file,'r'):
        m+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id_o=flist[2]
            if bin_id_o != '*':
                bin_id=flist[2].split('_')[0]
                read_id=flist[0]
                read_id_name=str(n)+'_'+read_id #
                fq_seq=flist[9]+'\n'+'+'+'\n'+flist[10]+'\n'
                try:
                    fq[bin_id][read_id_name]+=1
                    i=fq[bin_id][read_id_name]
                    f1=open(str(bin_id)+'_seq_R'+str(i)+'.fq','a')
                    f1.write('@'+str(read_id_name)+' '+str(i)+'\n'+str(fq_seq))
                    # f1.write('@'+str(read_id_name)[:-2]+' '+str(i)+'\n'+str(fq_seq))
                    f1.close()
                    if i == 2:
                        del fq[bin_id][read_id_name]
                except:
                    f_not_mapped_reads.write(str(read_id_name)+'\n')
    
        if m % 1000000 == 0:
            print('Parsed', m,'lines')
    f_not_mapped_reads.close()

def mapping_extraction(lr_elongated_bin_folder, datasets_list, num_threads, pwd, mp_run):
    print('Polishing seqs run '+str(mp_run))
    os.system('mkdir Extract_SR_'+str(mp_run))
    f=open('Total_seq.fa','w')
    record_bin_seq, fq, pair, bin_seq = {}, {}, {}, {}
    os.chdir(pwd+'/'+lr_elongated_bin_folder)
    for root, dirs, files in os.walk(pwd+'/'+lr_elongated_bin_folder):
        for file in files:
            hz=file.split('.')[-1]
            if 'fa' in hz:
                bin_id=str(file).split('_')[0]
                fq[str(bin_id)]={}
                pair[str(bin_id)]={}
                bin_seq[str(bin_id)]=[]
                bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R1.fq')
                bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R2.fq')
                record_bin_seq[bin_id]={}
                for record in SeqIO.parse(file,'fasta'):
                    record_bin_seq[bin_id][record.id]=record.seq
                    f.write('>'+str(bin_id)+'_'+str(record.id)+'\n'+str(record.seq)+'\n')
    os.chdir(pwd)
    f.close()

    for bin_id in bin_seq.keys():
        f1=open(str(bin_id)+'_seq_R1.fq','w')
        f2=open(str(bin_id)+'_seq_R2.fq','w')
        f1.close()
        f2.close()

    ### bwa
    os.system('bwa index Total_seq.fa')
    n, mp_run = 0, 1
    for item in datasets_list.keys():
        print('Mapping '+str(datasets_list[item]))
        n+=1
        os.system('bwa mem -t '+str(num_threads)+' Total_seq.fa '+str(datasets_list[item][0])+' '+str(datasets_list[item][1])+' > '+str(item)+'.sam') ### bwa
        parse_sam_bwa(str(item)+'.sam', fq, pair, n, mp_run)
        os.system('rm '+str(item)+'.sam')
    
    for bin_id in bin_seq.keys():
        for item in bin_seq[bin_id]:
            os.system('mv '+str(item)+' '+pwd+'/Extract_SR_'+str(mp_run))
    os.system('rm Total_seq.fa')

def sr_polishing(bin_id, bin_name, bin_seq, LR_merged_folder, extra_SR_folder, num_threads, pwd):
    r1=bin_seq[bin_id][0]
    r2=bin_seq[bin_id][1]
    
    schedual_polish_iteration, prun, t = 10, 1, 0
    while prun <= schedual_polish_iteration and t == 0:
        print('Processing '+str(bin_id)+' in polishing '+str(prun))
        flog=open('BASALT_log.txt','a')
        flog.write('Processing '+str(bin_id)+' in polishing '+str(prun)+'\n')
        flog.close()
        output_name=bin_name+'_polished'+str(prun)
        if prun == 1:
            bin_name2=bin_name
            os.system('cp '+str(pwd)+'/'+str(LR_merged_folder)+'/'+str(bin_name)+' '+str(bin_name))
            ###
            bin_contigs={}
            try:
                for record in SeqIO.parse(str(bin_name),'fasta'):
                    bin_contigs[record.seq]=0
            except:
                xyzzz=0
        else:
            bin_name2=bin_name+'_polished'+str(prun-1)+'.fasta'
        # num=len(bin_contigs)
        ### bwa
        os.system('bwa index '+str(bin_name2))
        os.system('bwa mem -t '+str(num_threads)+' '+str(bin_name2)+' '+str(pwd)+'/'+str(extra_SR_folder)+'/'+str(r1)+' '+str(pwd)+'/'+str(extra_SR_folder)+'/'+str(r2)+' | samtools sort -@ '+str(num_threads)+' -O bam -o '+str(bin_id)+'.bam')
        
        os.system('samtools view -@ '+str(num_threads)+' -q 30 -b '+str(bin_id)+'.bam > '+str(bin_id)+'_sort.bam')
        os.system('samtools index -@  '+str(num_threads)+' '+str(bin_id)+'_sort.bam')
        # short read consensus call
        os.system('pilon --genome '+str(bin_name2)+' --frags '+str(bin_id)+'_sort.bam --fix all --output '+str(output_name))
        xy, redo_run = 0, 0
        while xy == 0 and redo_run <= 3:
            try:
                xy=0
                for line in open(str(output_name)+'.fasta','r'):
                    xy+=1
                    if xy == 2:
                        break
            except:
                xy=0
            
            if xy == 0:
                print('Redo pilon polishing of processing '+str(bin_id)+' in polishing '+str(prun))
                redo_run += 1
                os.system('pilon --genome '+str(bin_name2)+' --frags '+str(bin_id)+'_sort.bam --fix all --output '+str(output_name))

            try:
                xy=0
                for line in open(str(output_name)+'.fasta','r'):
                    xy+=1
                    if xy == 2:
                        break
            except:
                xy=0

        if xy == 0 and redo_run == 4:
            print('Re-polish the bin '+str(bin_id)+' later')
        else:
            if prun == schedual_polish_iteration:
                print(str(bin_id)+' polishing finished at '+str(prun))
                t=1
                flog=open('BASALT_log.txt','a')
                flog.write(str(bin_id)+' polishing finished at '+str(prun)+'\n')
                flog.close()
                os.system('mv '+str(output_name)+'.fasta '+pwd+'/'+LR_merged_folder+'_polished')
                for i in range(1, prun):
                    os.system('rm '+bin_name+'_polished'+str(i)+'.fasta')
            else:
                for record in SeqIO.parse(str(output_name)+'.fasta','fasta'):
                    try:
                        bin_contigs[record.seq]+=1
                        t=1
                        print(str(bin_id)+' polishing finished at '+str(prun))
                        flog=open('BASALT_log.txt','a')
                        flog.write(str(bin_id)+' polishing finished at '+str(prun)+'\n')
                        flog.close()
                        os.system('mv '+str(output_name)+'.fasta '+pwd+'/'+LR_merged_folder+'_polished')
                        for i in range(1, prun):
                            os.system('rm '+bin_name+'_polished'+str(i)+'.fasta')
                        
                    except:
                        bin_contigs[record.seq]=0
                        t=0

            os.system('rm '+str(bin_name2)+'.sa '+str(bin_name2)+'.pac '+str(bin_name2)+'.ann '+str(bin_name2)+'.amb '+str(bin_name2)+'.bwt '+str(bin_id)+'_sort.bam.bai '+str(bin_id)+'_sort.bam '+str(bin_id)+'.bam')
            prun+=1

def sr_polishing_main(LR_merged_folder, extra_SR_folder, num_threads, ram, pwd):
    bin_lr_merged, bin_seq, x, unpolished_bin = {}, {}, 0, {}
    print('Recording seqs')
    os.chdir(pwd+'/'+str(LR_merged_folder))
    for root, dirs, files in os.walk(pwd+'/'+str(LR_merged_folder)):
        for file in files:
            hz=file.split('.')[-1]
            if 'fa' in hz:
                x+=1
                bin_id=str(file).split('_')[0]
                bin_lr_merged[bin_id]=file
                bin_seq[bin_id]=[bin_id+'_seq_R1.fq',bin_id+'_seq_R2.fq']
                unpolished_bin[bin_id]=0
    os.chdir(pwd)
    os.system('mkdir '+LR_merged_folder+'_polished')

    y, polished_bin = 0, {}
    os.chdir(pwd+'/'+LR_merged_folder+'_polished')
    for root, dirs, files in os.walk(pwd+'/'+LR_merged_folder+'_polished'):
        for file in files:
            hz=file.split('.')[-1]
            if 'fasta' in hz:
                y+=1
                bin_id=str(file).split('_')[0]
                polished_bin[bin_id]=0
    os.chdir(pwd)

    unpolished_bin={}
    for bin_id in bin_lr_merged.keys():
        if bin_id not in polished_bin.keys():
            unpolished_bin[bin_id]=0

    ### Est. meomery usage
    print('Est. meomery usage')
    bin_mem, xyz ={}, 0
    os.chdir(pwd+'/'+extra_SR_folder)
    for root, dirs, files in os.walk(pwd+'/'+extra_SR_folder):
        for file in files:
            if '_seq_R1.fq' in file:
                xyz+=1
                bin_id=file.split('_')[0]
                if xyz == 1:
                    xyzz=0
                    for line in open(file,'r'):
                        xyzz+=1
                        if xyzz == 2:
                            read_len=len(line.strip())
                else:
                    xyzz=0
                    for line in open(file,'r'):
                        xyzz+=1

                bin_mem[bin_id]=int(xyzz * read_len * 150 / 2000000000)
    os.chdir(pwd)

    ### Meomery strategy
    mem_bin, ram_t, project_list, ram_rpp, tpp = {}, {}, {}, [], 6
    max_project=int(num_threads/tpp)
    rpp=int(ram/max_project)
    iteration=int(max_project/2)
    for i in range(1,iteration+1):
        mem_bin[rpp*i]=[]
        projs=int(num_threads/tpp/i)
        tpp=int(num_threads/projs)
        ram_t[rpp*i]=tpp
        ram_rpp.append(rpp*i)
        project_list[rpp*i]=projs
    if ram != rpp*i:
        mem_bin[ram]=[]
        ram_rpp.append(ram)
        ram_t[ram]=num_threads
        project_list[ram]=1
    # biggest_reads=rpp/100

    fm=open('Bin_polish_esb_mem.txt','w')
    fm2=open('Mem_Bin.txt','w')
    processed={}
    for bin_id in bin_mem.keys():
        mem=bin_mem[bin_id]
        fm.write(str(bin_id)+'\t'+str(mem)+'\n')
        for rpp in ram_rpp:
            if mem <= rpp and bin_id not in processed.keys():
                mem_bin[rpp].append(bin_id)
                processed[bin_id]=0
        if mem >= ram:
            mem_bin[ram].append(bin_id)
    fm.close()

    mem_bin2=copy.deepcopy(mem_bin)
    for rpp in mem_bin2.keys():
        fm2.write(str(rpp)+'\t'+str(mem_bin2[rpp])+'\n')
        if len(mem_bin2[rpp]) == 0:
            del mem_bin[rpp]
    fm2.close()

    print('Start to polish')
    # y=0
    # while x != y:
    #     for bin_id in bin_lr_merged.keys():
    #         if bin_id in unpolished_bin.keys():
    #             bin_name=bin_lr_merged[bin_id]
    #             sr_polishing(bin_id, bin_name, bin_seq, LR_merged_folder, extra_SR_folder, num_threads, pwd)
        
    #     y, polished_bin = 0, {}
    #     os.chdir(pwd+'/'+LR_merged_folder+'_polished')
    #     for root, dirs, files in os.walk(pwd+'/'+LR_merged_folder+'_polished'):
    #         for file in files:
    #             hz=file.split('.')[-1]
    #             if 'fasta' in hz:
    #                 y+=1
    #                 bin_id=str(file).split('_')[0]
    #                 polished_bin[bin_id]=0
    #     os.chdir(pwd)

    #     unpolished_bin={}
    #     for bin_id in bin_lr_merged.keys():
    #         if bin_id not in polished_bin.keys():
    #             unpolished_bin[bin_id]=0

    y=0
    while x != y:
        for rpp in mem_bin.keys():
            bin_list=mem_bin[rpp]
            pros=project_list[rpp]
            tpp=ram_t[rpp]
            print(str(tpp)+' threads; '+str(rpp)+'G for each project')
            pool=Pool(processes=pros)
            for bin_id in bin_list:
                if bin_id in unpolished_bin.keys():
                    bin_name=bin_lr_merged[bin_id]
                    pool.apply_async(sr_polishing, args=(bin_id, bin_name, bin_seq, LR_merged_folder, extra_SR_folder, tpp, pwd))
            pool.close()
            pool.join()

        y, polished_bin = 0, {}
        os.chdir(pwd+'/'+LR_merged_folder+'_polished')
        for root, dirs, files in os.walk(pwd+'/'+LR_merged_folder+'_polished'):
            for file in files:
                hz=file.split('.')[-1]
                if 'fasta' in hz:
                    y+=1
                    bin_id=str(file).split('_')[0]
                    polished_bin[bin_id]=0
        os.chdir(pwd)

        unpolished_bin={}
        for bin_id in bin_lr_merged.keys():
            if bin_id not in polished_bin.keys():
                unpolished_bin[bin_id]=0

def Generation_new_bins(LR_elongation_polished_folder, LR_elongation_folder, Read_connected_bins_folder, target_bin_folder, pwd):
    bin_merged_contigs={}
    os.chdir(pwd+'/'+LR_elongation_folder)
    for root, dirs, files in os.walk(pwd+'/'+LR_elongation_folder):
        for file in files:
            if '_lr_being_merged_contigs.txt' in file:
                bin_id=str(file).split('_')[0]
                bin_merged_contigs[bin_id]={}
                for line in open(file,'r'):
                    bin_merged_contigs[bin_id][line.strip()]=0

    bin_polished_contigs={}
    os.chdir(pwd+'/'+LR_elongation_polished_folder)
    for root, dirs, files in os.walk(pwd+'/'+LR_elongation_polished_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                bin_id=str(file).split('_')[0]
                bin_polished_contigs[bin_id]={}
                for record in SeqIO.parse(file,'fasta'):
                    bin_polished_contigs[bin_id][record.id]=record.seq

    os.chdir(pwd)
    os.system('mkdir '+target_bin_folder+'_gap_filling_final')

    os.chdir(pwd+'/'+Read_connected_bins_folder)
    for root, dirs, files in os.walk(pwd+'/'+Read_connected_bins_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                bin_id=str(file).split('_')[0]
                if bin_id not in bin_polished_contigs.keys():
                    os.system('mv '+str(file)+' '+pwd+'/'+target_bin_folder+'_gap_filling_final')
                else:
                    f=open(bin_id+'_gf.fa','w')
                    for record in SeqIO.parse(file,'fasta'):
                        if record.id not in bin_merged_contigs[bin_id].keys():
                            f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    for item in bin_polished_contigs[bin_id].keys():
                        f.write('>'+str(item)+'\n'+str(bin_polished_contigs[bin_id][item])+'\n')
                    f.close()
                    os.system('mv '+str(bin_id)+'_gf.fa '+pwd+'/'+target_bin_folder+'_gap_filling_final')
    os.chdir(pwd)

def gap_filling_main(bin_folder, short_reads_folder, long_reads_folder, datasets_list, num_threads, ram):
    pwd=os.getcwd()
    try:
        f=open('BASALT_log.txt','a')
    except:
        f=open('BASALT_log.txt','w')
    f.close()       

    try:
        f=open('Gap_filling_status.txt','a')
        f.close()
    except:
        f=open('Gap_filling_status.txt','w')
        f.close()

    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Self connecting done!' in line:
            x=1

    if x == 0:
        os.chdir(pwd+'/'+str(bin_folder))
        bin_contigs_length={}
        for root, dirs, files in os.walk(pwd+'/'+str(bin_folder)):
            for file in files:
                # try:
                bin_contigs_length[file]={}
                f=open('sc_'+file, 'w')
                for record in SeqIO.parse(file, 'fasta'):
                    f.write('>'+str(record.id).replace('_pilon','')+'\n'+str(record.seq)+'\n')
                    bin_contigs_length[file][str(record.id).replace('_pilon','')]=len(record.seq)
                f.close()
                os.system('mv sc_'+file+' '+file)
                os.system('cp '+pwd+'/'+str(bin_folder)+'/'+str(file)+' '+pwd)
                # except:
                #     xyzzz=0
        os.chdir(pwd)
        # print(str(len(bin_contigs_length)))

        os.system('mkdir '+str(target_bin_folder)+'_selfblast_output')
        os.system('mkdir '+str(target_bin_folder)+'_self_connected_bins')
        fx=open(str(target_bin_folder)+'_sc_bins_summary.txt','w')
        fx.close()
        flog=open('BASALT_log.txt','a')
        pool=Pool(processes=num_threads)
        for bins in bin_contigs_length.keys():
            print('Processing '+str(bins))
            flog.write('Processing '+str(bins)+'\n')
            # self_connecting(bins, bin_folder, bin_contigs_length, pwd, num_threads)
            pool.apply_async(self_connecting, args=(bins, bin_folder, bin_contigs_length, pwd, num_threads))
        pool.close()
        pool.join()

        os.system('rm *selfblast.txt *_blast.txt *.nhr *.nin *.nsq')
        print('Self connecting done!')
        flog.write('Self connecting done!'+'\n')
        f=open('Gap_filling_status.txt','a')
        f.write('Self connecting done!'+'\n')
        f.close()
        flog.close()

    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Short-reads fq2fa done!' in line:
            x=1

    if x == 0:
        flog=open('BASALT_log.txt','a')
        os.system('mkdir '+str(target_bin_folder)+'_SR_fa')
        os.chdir(pwd+'/'+str(short_reads_folder))
        bin_sr={}
        for root, dirs, files in os.walk(pwd+'/'+str(short_reads_folder)):
            for file in files:
                if '.fq' in file:
                    print('Converting '+str(file))
                    flog.write('Converting '+str(file)+'\n')
                    bin_id=str(file).split('_')[0]
                    sr=bin_id+'_sr.fa'
                    bin_sr[bin_id]=sr
                    try:
                        f=open(sr,'a')
                    except:
                        f=open(sr,'w')
                        
                    n=0
                    for line in open(file,'r'):
                        n+=1
                        if n % 4 == 1:
                            f.write('>'+str(line.strip().replace(' ','/'))[1:]+'\n')
                        elif n % 4 == 2:
                            f.write(str(line))
                    f.close()

        os.system('mv *_sr.fa '+pwd+'/'+str(target_bin_folder)+'_SR_fa')
        os.chdir(pwd)
        f=open('Gap_filling_status.txt','a')
        f.write('Short-reads fq2fa done!'+'\n')
        flog.write('Short-reads fq2fa done!'+'\n')
        f.close()
        flog.close()

    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Short-reads blast done!' in line:
            x=1

    if x == 0:
        flog=open('BASALT_log.txt','a')
        short_reads_blast(str(target_bin_folder), str(target_bin_folder)+'_self_connected_bins', str(target_bin_folder)+'_SR_fa', num_threads, pwd)
        f=open('Gap_filling_status.txt','a')
        f.write('Short-reads blast done!'+'\n')
        flog.write('Short-reads blast done!'+'\n')
        f.close()
        flog.close()

    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Short-reads merging done!' in line:
            x=1

    if x == 0:
        os.system('mkdir '+str(target_bin_folder)+'_read_connected_bins')
        short_reads_merging(str(target_bin_folder), str(target_bin_folder)+'_all_filtrated_read_blast_output', str(target_bin_folder)+'_SR_fa', str(target_bin_folder)+'_self_connected_bins', num_threads, pwd)
        f=open('Gap_filling_status.txt','a')
        f.write('Short-reads merging done!'+'\n')
        f.close()

    ### Shutdown this part ----------------------
    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Long-reads merge done!' in line:
            x=1

    if x == 0:
        long_reads_merge(str(target_bin_folder), str(target_bin_folder)+'_read_connected_bins', long_reads_folder, num_threads, pwd)
        f=open('Gap_filling_status.txt','a')
        f.write('Long-reads merge done!'+'\n')
        f.close()

    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Mapping_extraction done!' in line:
            x=1

    mp_run=1
    if x == 0:
        print('Start mapping')
        mapping_extraction(str(target_bin_folder)+'_LR_elongation', datasets_list, num_threads, pwd, mp_run)
        f=open('Gap_filling_status.txt','a')
        f.write('Mapping_extraction done!'+'\n')    
        f.close()

    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Polishing done!' in line:
            x=1

    if x == 0:
        print('Start polishing')
        sr_polishing_main(str(target_bin_folder)+'_LR_elongation', 'Extract_SR_'+str(mp_run), num_threads, ram, pwd)
        f=open('Gap_filling_status.txt','a')
        f.write('Polishing done!'+'\n')    
        f.close()
    
    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Generation of new bins done!' in line:
            x=1

    if x == 0:
        print('Generating new bins')
        Generation_new_bins(str(target_bin_folder)+'_LR_elongation_polished', str(target_bin_folder)+'_LR_elongation', str(target_bin_folder)+'_read_connected_bins', str(target_bin_folder), pwd)
        f=open('Gap_filling_status.txt','a')
        f.write('Generation of new bins done!'+'\n')    
        f.close()

    x=0
    for line in open('Gap_filling_status.txt','r'):
        if 'Bin 1st evaluation done!' in line:
            x=1

    if x == 0:
        bin_evaluation(target_bin_folder, str(target_bin_folder)+'_read_connected_bins', str(target_bin_folder)+'_gap_filling_final', num_threads, pwd)
        f=open('Gap_filling_status.txt','a')
        f.write('Bin 1st evaluation done!'+'\n')
        f.close()

    ### Shutdown this part ----------------------

if __name__ == '__main__': 
    # target_bin_folder='MAG_polished_t2'
    # target_bin_folder='bin97_test'
    target_bin_folder='1_Opera_unpolished_cat_contigs.fasta_BestBinsSet_outlier_refined_filtrated_retrieved_MAGs_polished'
    short_reads_folder='1_Opera_unpolished_cat_contigs.fasta_BestBinsSet_outlier_refined_filtrated_retrieved_gf_lr_sr_bins_seq'
    long_reads_folder='1_Opera_unpolished_cat_contigs.fasta_BestBinsSet_outlier_refined_filtrated_retrieved_polished_lr_gf_strict'
    datasets_list={'1':['RH_S001_insert_270_1.fastq','RH_S001_insert_270_2.fastq'], '2':['RH_S002_insert_270_1.fastq','RH_S002_insert_270_2.fastq'], '3':['RH_S003_insert_270_1.fastq','RH_S003_insert_270_2.fastq'],'4':['RH_S004_insert_270_1.fastq','RH_S004_insert_270_2.fastq'],'5':['RH_S005_insert_270_1.fastq','RH_S005_insert_270_2.fastq']}
    # '1_Opera_unpolished_cat_contigs.fasta_BestBinsSet_outlier_refined_filtrated_retrieved_gf_lr_polished'
    num_threads=30
    ram=250
    aligned_len_cutoff=200
    similarity_cutoff=99
    coverage_extension=95
    gap_filling_main(target_bin_folder, short_reads_folder, long_reads_folder, datasets_list, num_threads, ram)