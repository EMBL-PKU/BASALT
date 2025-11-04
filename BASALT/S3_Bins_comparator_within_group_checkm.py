#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

from tempfile import TemporaryFile
from Bio import SeqIO
from S2_BinsAbundance_PE_connections_multiple_processes_pool_checkm import *
import os

def contig_id_recorder(genome_folder):
    genomes_sum={}
    pwd=os.getcwd()
    n=0
    for item in genome_folder:
        n+=1
        print('Parsing '+item+' group of bin-set')
        print('--------------------------')
        if n == 1:
            genomes_sum[str(item)]={}
            relation, bin_len, bins1={}, {}, []
            frist=item
            for root, dirs, files in os.walk(item+'_genomes'):
                os.chdir(item+'_genomes')
                for file in files:
                    if '_genomes.' in file:
                        hz=file.split('.')[-1]
                        if 'fasta' in hz and '_maxbin2_genomes' in file:
                            # print('Parsing', file
                            bins1.append(file)
                            genomes_sum[item][file]={}
                            genomes_sum[item][file]['contig']={}
                            bin_len[file]=0
                            for record in SeqIO.parse(file, 'fasta'):
                                bin_len[file]+=len(record.seq)
                                genomes_sum[item][file]['contig'][record.id]=len(record.seq)
                        elif 'fasta' in hz and '_concoct_genomes' in file:
                            # print('Parsing', file
                            bins1.append(file)
                            genomes_sum[item][file]={}
                            genomes_sum[item][file]['contig']={}
                            bin_len[file]=0
                            for record in SeqIO.parse(file, 'fasta'):
                                bin_len[file]+=len(record.seq)
                                genomes_sum[item][file]['contig'][record.id]=len(record.seq)
                        elif 'fa' in hz and '_metabat_genomes' in file and 'unbinned.txt' not in file:
                            # print('Parsing', file
                            bins1.append(file)
                            genomes_sum[item][file]={}
                            genomes_sum[item][file]['contig']={}
                            bin_len[file]=0
                            for record in SeqIO.parse(file, 'fasta'):
                                bin_len[file]+=len(record.seq)
                                genomes_sum[item][file]['contig'][record.id]=len(record.seq)
                        elif 'fa' in hz and '_metabinner_genomes' in file:
                            # print('Parsing', file
                            bins1.append(file)
                            genomes_sum[item][file]={}
                            genomes_sum[item][file]['contig']={}
                            bin_len[file]=0
                            for record in SeqIO.parse(file, 'fasta'):
                                bin_len[file]+=len(record.seq)
                                genomes_sum[item][file]['contig'][record.id]=len(record.seq)
                        elif 'fa' in hz and '_vamb_genomes' in file:
                            # print('Parsing', file
                            bins1.append(file)
                            genomes_sum[item][file]={}
                            genomes_sum[item][file]['contig']={}
                            bin_len[file]=0
                            for record in SeqIO.parse(file, 'fasta'):
                                bin_len[file]+=len(record.seq)
                                genomes_sum[item][file]['contig'][record.id]=len(record.seq)  
                        elif 'fa' in hz and '_semibin_genomes' in file:
                            # print('Parsing', file
                            bins1.append(file)
                            genomes_sum[item][file]={}
                            genomes_sum[item][file]['contig']={}
                            bin_len[file]=0
                            for record in SeqIO.parse(file, 'fasta'):
                                bin_len[file]+=len(record.seq)
                                genomes_sum[item][file]['contig'][record.id]=len(record.seq)
                        elif 'fa' in hz and '_SingleContig_genomes' in file:
                            # print('Parsing', file
                            bins1.append(file)
                            genomes_sum[item][file]={}
                            genomes_sum[item][file]['contig']={}
                            bin_len[file]=0
                            for record in SeqIO.parse(file, 'fasta'):
                                bin_len[file]+=len(record.seq)
                                genomes_sum[item][file]['contig'][record.id]=len(record.seq)
                        else:
                            continue
            os.chdir(pwd)
        else:
            for root, dirs, files in os.walk(item+'_genomes'):
                os.chdir(item+'_genomes')
                bins2=[]
                for file in files:
                    hz=file.split('.')[-1]
                    if 'fasta' in hz and '_maxbin2_genomes' in file:
                        # print('Parsing', file
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue
                    
                    elif 'fasta' in hz and '_concoct_genomes' in file:
                        # print('Parsing', file
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue

                    elif 'fa' in hz and '_metabat_genomes' in file and 'unbinned.txt' not in file:
                        # print('Parsing', file
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue

                    elif 'fa' in hz and '_metabinner_genomes' in file:
                        # print('Parsing', file
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue

                    elif 'fa' in hz and '_vamb_genomes' in file:
                        # print('Parsing', file
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue
                    
                    elif 'fa' in hz and '_semibin_genomes' in file:
                        # print('Parsing', file
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue

                    elif 'fa' in hz and '_SingleContig_genomes' in file:
                        # print('Parsing', file
                        bins2.append(file)
                        bin_len[file]=0
                        for record in SeqIO.parse(file, 'fasta'):
                            bin_len[file]+=len(record.seq)
                            for item2 in genomes_sum[frist].keys():
                                if record.id in genomes_sum[frist][item2]['contig']:
                                    if str(item2+'---'+file) not in relation.keys():
                                        relation[str(item2+'---'+file)]=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                    else:
                                        relation[str(item2+'---'+file)]+=round(100*(float(genomes_sum[frist][item2]['contig'][record.id])/float(bin_len[item2])),2)
                                else:
                                    continue

                    else:
                        continue

                os.chdir(pwd)

    a, b, best_hit_genome, bins1_selected, bins2_selected, bins_extract={}, {}, {}, {}, {}, []
    for item in relation.keys():
        bin_set1=item.split('---')[0]
        bin_set2=item.split('---')[1]
        # print(bin_set1
        if bin_set1 not in a.keys():
            a[bin_set1]=relation[item]
            b[bin_set1]=item+'\t'+str(relation[item])+'\t'+str(round(float(relation[item])*float(bin_len[bin_set1])/float(bin_len[bin_set2]), 2))
        elif float(relation[item]) > float(a[bin_set1]):
            # print(bin_set1, item, float(relation[item]), float(a[bin_set1])
            a[bin_set1]=relation[item]
            b[bin_set1]=item+'\t'+str(relation[item])+'\t'+str(round(float(relation[item])*float(bin_len[bin_set1])/float(bin_len[bin_set2]), 2))
        else:
            continue

    for item in b.keys():
        bin_set1=b[item].split('\t')[0].split('---')[0]
        bin_set2=b[item].split('\t')[0].split('---')[1]
        bin1_score=float(b[item].split('\t')[1])
        bin2_score=float(b[item].split('\t')[2])
        total_score=bin1_score+bin2_score
        # if total_score >= 120 or bin1_score >= 80 or bin2_score >= 80:
        if total_score >= 100 or bin1_score >= 50 or bin2_score >= 50:
            best_hit_genome[item]=str(b[item])
            bins1_selected[bin_set1]=0
            bins2_selected[bin_set2]=0
    
    for item in bins1:
        if item not in bins1_selected.keys():
            bins_extract.append(item)
            # if 'maxbin2' in item:
            #     bins_extract.append(item+'.fasta')
            # elif 'metabat' in item:
            #     bins_extract.append(item+'.fa')

    for item in bins2:
        if item not in bins2_selected.keys():
            bins_extract.append(item)

    return relation, best_hit_genome, bins_extract

def checkm_connections(genome_folder):
    pwd=os.getcwd()
    bins_checkm={}
    print(pwd)
    for item in genome_folder:
        print('Reading checkm results of '+item)
        print('-------------------------')
        try:
            f=open(item+'_genomes/'+item+'_bin_stats_ext.tsv', 'r')
            print('Reading '+item+'_bin_stats_ext.tsv')
            for line in f:
                genome_ids=str(line).strip().split('\t')[0]
                if '_genomes.0' not in genome_ids:
                    bins_checkm[genome_ids]={}
                    #try:
                    marker_lineage=str(line).strip().split('\'marker lineage\': \'')[1].strip().split('\'')[0]
                    bins_checkm[genome_ids]['marker lineage']=marker_lineage
                    #except:
                        #print('marker lineage error'
                        #bins_checkm[genome_ids]['marker lineage']=root
                    #try:    
                    completeness=str(line).strip().split('\'Completeness\': ')[1].split(',')[0]
                    bins_checkm[genome_ids]['Completeness']=float(completeness)
                    #except:
                        #bins_checkm[genome_ids]['Completeness']=0

                    #try:
                    genome_size=str(line).strip().split('\'Genome size\':')[1].strip().split(',')[0]
                    bins_checkm[genome_ids]['Genome size']=int(genome_size)
                    #except:
                        #bins_checkm[genome_ids]['Genome size']=0

                    bins_checkm[genome_ids]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split('}')[0].split(',')[0].strip())
                    bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())
                    #except:
                    #    bins_checkm[genome_ids]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split(',')[0].strip())
                    #    bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0])	
                elif '_semibin_genomes.0' in genome_ids:
                    bins_checkm[genome_ids]={}
                    #try:
                    marker_lineage=str(line).strip().split('\'marker lineage\': \'')[1].strip().split('\'')[0]
                    bins_checkm[genome_ids]['marker lineage']=marker_lineage
                    #except:
                        #print('marker lineage error'
                        #bins_checkm[genome_ids]['marker lineage']=root
                    #try:    
                    completeness=str(line).strip().split('\'Completeness\': ')[1].split(',')[0]
                    bins_checkm[genome_ids]['Completeness']=float(completeness)
                    #except:
                        #bins_checkm[genome_ids]['Completeness']=0

                    #try:
                    genome_size=str(line).strip().split('\'Genome size\':')[1].strip().split(',')[0]
                    bins_checkm[genome_ids]['Genome size']=int(genome_size)
                    #except:
                        #bins_checkm[genome_ids]['Genome size']=0

                    bins_checkm[genome_ids]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split('}')[0].split(',')[0].strip())
                    bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0].strip())
                    #except:
                    #    bins_checkm[genome_ids]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split(',')[0].strip())
                    #    bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0])	
        except:
            print('CheckM output file reading error!')

        # print('Reading dataframe-connections of', item
        # print('-------------------------'
        
        try:
            f=open(item+'_genomes/Bins_total_connections_'+str(item)+'.txt', 'r')
            print('Reading dataframe-connections of Bins_total_connections_'+item+'.txt')
            n=0
            for line in f:
                n+=1
                if n >= 2 and '_genomes.0' not in str(line).strip().split('\t')[0]:
                    bins_checkm[str(line).strip().split('\t')[0]]['Connections']=int(str(line).strip().split('\t')[1])
        except:
            print('Please make sure Bins_total_connections_'+str(item)+'.txt under the folder.')

        for item2 in bins_checkm.keys():
            if 'Connections' not in bins_checkm[item2].keys():
                bins_checkm[item2]['Connections']=0

    return bins_checkm

def genome_selector(best_hit_genome, bin_set_checkm):
    print('Selecting bin-set')
    print('------------------')
    bin_selected={}
    marker_score={'root':0, 'k':1, 'p':1.5, 'c':2.3, 'o':3.4, 'f':5.1, 'g':7.6, 's':11.4}
    for item in best_hit_genome.keys():
        set1=str(best_hit_genome[item]).split('\t')[0].split('---')[0]
        set2=str(best_hit_genome[item]).split('\t')[0].split('---')[1]
        genome_name_list1=set1.split('.')
        genome_name_list2=set2.split('.')
        genome_name_list1.remove(genome_name_list1[-1])
        genome_name_list2.remove(genome_name_list2[-1])
        set1='.'.join(genome_name_list1)
        set2='.'.join(genome_name_list2)
#	print(set1, set2

        if '__' in str(bin_set_checkm[set1]['marker lineage']):
            set1_marker=bin_set_checkm[set1]['marker lineage'].split('__')[0]
            set1_marker_score=marker_score[set1_marker]
        else:
            set1_marker_score=0
        
        if '__' in str(bin_set_checkm[set2]['marker lineage']):
            set2_marker=bin_set_checkm[set2]['marker lineage'].split('__')[0]
            set2_marker_score=marker_score[set2_marker]
        else:
            set2_marker_score=0

        set1_cpn=bin_set_checkm[set1]['Completeness']
        set2_cpn=bin_set_checkm[set2]['Completeness']
        set1_ctn=bin_set_checkm[set1]['Contamination']

        set2_ctn=bin_set_checkm[set2]['Contamination']

        if set1_marker_score == set2_marker_score:
            set1_cpn_ctn=float(set1_cpn)-float(set1_ctn)
            set2_cpn_ctn=float(set2_cpn)-float(set2_ctn)
            if float(set1_cpn_ctn) == float(set2_cpn_ctn):
                if float(bin_set_checkm[set1]['Connections']) != 0 and float(bin_set_checkm[set2]['Connections']) != 0:
                    if float(bin_set_checkm[set1]['Genome size'])/float(bin_set_checkm[set1]['Connections']) >= float(bin_set_checkm[set2]['Genome size'])/float(bin_set_checkm[set2]['Connections']):
                        bin_selected[set1]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
                        # print('Adding', set1, 'to the selected bins set'
                    else:
                        bin_selected[set2]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
                elif float(bin_set_checkm[set1]['Connections']) == 0 and float(bin_set_checkm[set2]['Connections']) != 0:
                    bin_selected[set1]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
                elif float(bin_set_checkm[set1]['Connections']) != 0 and float(bin_set_checkm[set2]['Connections']) == 0:
                    bin_selected[set2]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
                else:
                    if float(bin_set_checkm[set1]['Genome size']) > float(bin_set_checkm[set2]['Genome size']):
                        bin_selected[set1]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
                        # print('Adding', set1, 'to the selected bins set'
                    else:
                        bin_selected[set2]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])

            elif float(set1_cpn_ctn) > float(set2_cpn_ctn):
                bin_selected[set1]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
                # print('Adding', set1, 'to the selected bins set'
            else:
                bin_selected[set2]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
                # print('Adding', set2, 'to the selected bins set'
        elif float(set1_marker_score) > float(set2_marker_score):
            bin_selected[set1]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
            # print('Adding', set1, 'to the selected bins set'
        else:
            bin_selected[set2]=best_hit_genome[item]+'\t'+str(bin_set_checkm[set1])+'\t'+str(bin_set_checkm[set2])
            # print('Adding', set2, 'to the selected bins set'
    print('bin-selecting accomplished!')
    print('---------------------------')
    return bin_selected

def two_groups_comparator(assembly, binset1, binset2, num):
    pwd=os.getcwd()
    try:
        fx=open('Basalt_log.txt','a')
    except:
        ttttx=0
    print('Iteration '+str(num))
    print('Comparing '+binset1+' and '+binset2)
    fx.write('Iteration '+str(num)+'\n')
    fx.write('Comparing '+binset1+' and '+binset2+'\n')
    f=open('All_possible_bin_sets_iteration_'+str(num)+'.txt','w')
    f2=open('Best_hit_bin_sets_iteration_'+str(num)+'.txt','w')

    genome_folder=[binset1, binset2]

    a=contig_id_recorder(genome_folder)
    all_hit_genome=a[0]
    best_hit_genome=a[1]
    bins_extract=a[2]

    for item in all_hit_genome.keys():
        f.write(item+'\t'+str(all_hit_genome[item])+'\n')
    f.close()

    for item in best_hit_genome.keys():
        f2.write(str(best_hit_genome[item])+'\n')
    f2.close()

    f=open('Bins_marker_lineage_completeness_contamination_iteration_'+str(num)+'.txt', 'w')
    
    bin_set_checkm=checkm_connections(genome_folder)

    for item in bin_set_checkm.keys():
        try:
            f.write(item+'\t'+str(bin_set_checkm[item]['marker lineage'])+'\t'+'Completeness:'+str(bin_set_checkm[item]['Completeness'])+'\t'+'Contaomination:'+str(bin_set_checkm[item]['Contamination'])+'\t'+'Genome size:'+str(bin_set_checkm[item]['Genome size'])+'\t'+'Connections:'+str(bin_set_checkm[item]['Connections'])+'\n')
        except:
            print(item+' error')
    f.close()

    bin_selected=genome_selector(best_hit_genome, bin_set_checkm)
    
    f=open('Extract_bins_in_iteration_'+str(num)+'.txt', 'w')
    for item in bins_extract:
        f.write(item+'\n')
    f.close()

    f=open('Best_bin_set_iteration_'+str(num)+'.txt','w')
    for item in bin_selected.keys():
        if 'maxbin2' in item or 'concoct' in item:
            f.write(item+'.fasta'+'\t'+str(bin_selected[item])+'\n')
        elif 'metabat' in item or 'metabinner' in item or 'vamb' in item or 'semibin' in item or 'SingleContig' in item:
            f.write(item+'.fa'+'\t'+str(bin_selected[item])+'\n')
        else:
            continue

    if len(bins_extract) >= 1:
        for item in bins_extract:
            item_list=item .split('.')
            item_list.remove(item_list[-1])
            name='.'.join(item_list)

            f.write(item+'\t'+str(bin_set_checkm[name])+'\n')
            bin_selected[name]='unique genome in', binset2
            print(item+' unique genome in '+binset2)
            print('----------------')
            fx.write(item+' unique genome in '+binset2+'\n')
            fx.write('----------------'+'\n')
    f.close()

    try:
        os.mkdir('Iteration_'+str(num)+'_genomes')
    except:
        os.system('rm -rf Iteration_'+str(num)+'_genomes')
        os.mkdir('Iteration_'+str(num)+'_genomes')
        print('Iteration_'+str(num)+'_genomes exist')
        print('Re-created folder of Iteration_'+str(num)+'_genomes')
        fx.write('Iteration_'+str(num)+'_genomes exist'+'\n')
        fx.write('Re-created folder of Iteration_'+str(num)+'_genomes'+'\n')
    
    f=open('Iteration_'+str(num)+'_genomes/Iteration_'+str(num)+'_bin_stats_ext.tsv','w')
    f2=open('Iteration_'+str(num)+'_genomes/Bins_total_connections_Iteration_'+str(num)+'.txt','w')
    f2.write('Bin'+'\t'+'Total_connections'+'\n')

    for item in bin_selected.keys():
        f.write(item+'\t'+str(bin_set_checkm[item])+'\n')
        f2.write(item+'\t'+str(bin_set_checkm[item]['Connections'])+'\n')
        if 'maxbin2' in item or 'concoct' in item:
            bin=item+'.fasta'
            try:
                folder=item.split('_genomes.')[0]
                os.chdir(pwd+'/'+folder+'_genomes')
                os.system('cp '+bin+' '+pwd+'/Iteration_'+str(num)+'_genomes')
            except:
                print('Copy bin-set error!')
        elif 'metabat' in item or 'metabinner' in item or 'vamb' in item or 'semibin' in item or 'SingleContig' in item:
            bin=item+'.fa'
            try:
                folder=item.split('_genomes.')[0]
                os.chdir(pwd+'/'+folder+'_genomes')
                os.system('cp '+bin+' '+pwd+'/Iteration_'+str(num)+'_genomes')
            except:
                print('Copy bin-set error!')
        else:
            continue
    f.close()  
    f2.close()

    os.chdir(pwd)
    fx.write(binset1+' and '+binset2+' comparison done'+'\n')
    try:
        f_s3=open('S3_checkpoint.txt','a')
    except:
        xyzzzz=0
    f_s3.write(str(assembly)+'\t'+str(num)+'\t'+binset1+' and '+binset2+' comparison done'+'\n')

def bin_within_a_group_comparitor(binset, assembly, num):
    pwd=os.getcwd()
    print('Comparing bins in final iteration')
    print('Parsing bins')
    print('---------------------------------')
    bins_sum={}
    for root, dirs, files in os.walk(binset):
        os.chdir(binset)
        for file in files:
            if '.fasta' in file and '_maxbin2_genomes' in file:
                # print('Parsing', file
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)
            elif '.fasta' in file and '_concoct_genomes' in file:
                # print('Parsing', file
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)
            elif '.fa' in file and '_metabat_genomes' in file and 'unbinned.txt' not in file:
                # print('Parsing', file
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)
                else:
                    continue
            elif '.fa' in file and '_metabinner_genomes' in file:
                # print('Parsing', file
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)
                else:
                    continue
            elif '.fa' in file and '_vamb_genomes' in file:
                # print('Parsing', file
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)
                else:
                    continue
            elif '.fa' in file and '_semibin_genomes' in file:
                # print('Parsing', file
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)
                else:
                    continue
            elif '.fa' in file and '_SingleContig_genomes' in file:
                # print('Parsing', file
                bins_sum[file]={}
                bins_sum[file]['contig']={}
                bins_sum[file]['totallen']=0
                for record in SeqIO.parse(file, 'fasta'):
                    bins_sum[file]['totallen']+=len(record.seq)
                    bins_sum[file]['contig'][record.id]=len(record.seq)
                else:
                    continue

    relation, relation2, processed={}, {}, {}
    for item in bins_sum.keys():
        processed[item]=''
        for contig in bins_sum[item]['contig'].keys():
            for item2 in bins_sum.keys():
                if item2 != item and item2 not in processed.keys():
                     if contig in bins_sum[item2]['contig'].keys():
                        if str(item+'---'+item2) not in relation.keys():
                            relation[str(item+'---'+item2)]=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item]['totallen'])),2)
                            relation2[str(item+'---'+item2)]=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item2]['totallen'])),2)
                        else:
                            relation[str(item+'---'+item2)]+=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item]['totallen'])),2)
                            relation2[str(item+'---'+item2)]+=round(100*(float(bins_sum[item2]['contig'][contig])/float(bins_sum[item2]['totallen'])),2)
                else:
                    continue
        #del bins_sum[item]

    final_iteration_checkm={}
    f=open('Iteration_'+str(num)+'_bin_stats_ext.tsv', 'r')
    for line in f:
        ids=str(line).strip().split('\t')[0]
        connections=int(str(line).strip().split('Connections\': ')[1].split('}')[0].split(',')[0])
        marker_lineage=str(line).strip().split('marker lineage\': \'')[1].split('}')[0].split('\'')[0]
        completeness=float(str(line).strip().split('\'Completeness\': ')[1].split('}')[0].split(',')[0])
        # try:
        contamination=float(str(line).strip().split('\'Contamination\': ')[1].split('}')[0].split(',')[0])
        # except:
        #     contamination=float(str(line).strip().split('\'Contamination\': ')[1].split(',')[0])
        # try:
        Mean_scaffold_length=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].split(',')[0].strip())
        # except:
            # Mean_scaffold_length=float(str(line).strip().split('Mean scaffold length\': ')[1].split('}')[0].strip())
        genome_size=int(str(line).strip().split('\'Genome size\': ')[1].split('}')[0].split(',')[0])
        if '_maxbin2_genomes.' in ids or '_concoct_genomes.' in ids:
            bin_id=ids+'.fasta'
        elif '_metabat_genomes.' in ids or '_semibin_genomes.' in ids or '_SingleContig_genomes.' in ids:
            bin_id=ids+'.fa'
        elif '_metabinner_genomes.' in ids:
            bin_id=ids+'.fa'
        elif '_vamb_genomes.' in ids:
            bin_id=ids+'.fa'
        else:
            continue 
        final_iteration_checkm[bin_id]={'Connections': connections, 'marker lineage': marker_lineage, 'Completeness': completeness, 'Genome size': genome_size, 'Mean scaffold length': Mean_scaffold_length, 'Contamination': contamination}

    os.chdir(pwd)
    try:
        os.system('mkdir '+str(assembly)+'_BestBinsSet')
    except:
        os.system('rm -rf '+str(assembly)+'_BestBinsSet')
        os.system('mkdir '+str(assembly)+'_BestBinsSet')
        print(str(assembly)+'_BestBinsSet exist')
        print('Re-created folder of '+str(assembly)+'_BestBinsSet')

    os.chdir(pwd+'/'+str(assembly)+'_BestBinsSet')
    pos_bins={}
    f=open('Possible_same_bin.txt','w')
    f2=open('Highly_possible_same_bin.txt','w')
    for item in relation.keys():
        # if float(relation[item]) >= 50:
        f.write(item+'\t'+str(relation[item])+'\t'+str(relation2[item])+'\n')
        sim=float(relation[item]) + float(relation2[item])
        if sim >= 100 or float(relation[item]) >= 50 or float(relation2[item]) >= 50:
        # if sim >= 120 or float(relation[item]) >= 80 or float(relation2[item]) >= 80:
            f2.write(item+'\t'+str(relation[item])+'\t'+str(relation2[item])+'\t'+str(final_iteration_checkm[str(item).split('---')[0]])+'\t'+str(final_iteration_checkm[str(item).split('---')[1]])+'\n')
            pos_bins[item]=str(relation[item])+'\t'+str(relation2[item])

    f.close()
    f2.close()

    marker_score={'root':0, 'k':1, 'p':1.5, 'c':2.3, 'o':3.4, 'f':5.1, 'g':7.6, 's':11.4}
    remain_bin, del_bin={}, {}
    for item in pos_bins.keys():
        set1=str(item).split('---')[0]
        set2=str(item).split('---')[1]

        if '__' in str(final_iteration_checkm[set1]['marker lineage']):
            set1_marker=final_iteration_checkm[set1]['marker lineage'].split('__')[0]
            set1_marker_score=marker_score[set1_marker]
        else:
            set1_marker_score=0
        
        if '__' in str(final_iteration_checkm[set2]['marker lineage']):
            set2_marker=final_iteration_checkm[set2]['marker lineage'].split('__')[0]
            set2_marker_score=marker_score[set2_marker]
        else:
            set2_marker_score=0

        set1_cpn=final_iteration_checkm[set1]['Completeness']
        set2_cpn=final_iteration_checkm[set2]['Completeness']
        set1_ctn=final_iteration_checkm[set1]['Contamination']
        set2_ctn=final_iteration_checkm[set2]['Contamination']

        set1_cpn_ctn=float(set1_cpn)-float(set1_ctn)
        set2_cpn_ctn=float(set2_cpn)-float(set2_ctn)

        if set1_marker_score == set2_marker_score:
            if float(set1_cpn_ctn) == float(set2_cpn_ctn):
                if float(final_iteration_checkm[set1]['Connections']) != 0 and float(final_iteration_checkm[set2]['Connections']) != 0:
                    if float(final_iteration_checkm[set1]['Genome size'])/float(final_iteration_checkm[set1]['Connections']) >= float(final_iteration_checkm[set2]['Genome size'])/float(final_iteration_checkm[set2]['Connections']):
                        remain_bin[set1]=0
                        del_bin[set2]=0
                    else:
                        remain_bin[set2]=0
                        del_bin[set1]=0
                elif float(final_iteration_checkm[set1]['Connections']) == 0 and float(final_iteration_checkm[set2]['Connections']) != 0:
                    remain_bin[set1]=0
                    del_bin[set2]=0
                elif float(final_iteration_checkm[set1]['Connections']) != 0 and float(final_iteration_checkm[set2]['Connections']) == 0:
                    remain_bin[set2]=0
                    del_bin[set1]=0
                else:
                    if float(final_iteration_checkm[set1]['Genome size']) > float(final_iteration_checkm[set2]['Genome size']):
                        remain_bin[set1]=0
                        del_bin[set2]=0
                    else:
                        remain_bin[set2]=0
                        del_bin[set1]=0
            elif float(set1_cpn_ctn) > float(set2_cpn_ctn):
                remain_bin[set1]=0
                del_bin[set2]=0
            else:
                remain_bin[set2]=0
                del_bin[set1]=0
        else:
            if float(set1_marker_score)*float(set1_cpn_ctn) >= float(set2_marker_score)*set2_cpn_ctn:
                remain_bin[set1]=0
                del_bin[set2]=0
            else:
                remain_bin[set2]=0
                del_bin[set1]=0

    #for item in del_bin.keys():
    #    if item in remain_bin.keys():
    #        del del_bin[item]
    #    else:
    #        continue
    
    #for item in final_iteration_checkm.keys():
    #    if item in del_bin.keys():
    for item in del_bin.keys():
        if item in final_iteration_checkm.keys():
            del final_iteration_checkm[item]
            print('Deleted '+item+' from the selected bins set')
        else:
            continue

    bin_selected={}
    f=open(assembly+'_BestBinSet_bin_stats_ext.tsv','w')
    for item in final_iteration_checkm.keys():
        checkm_id_list=item.split('.')
        checkm_id_list.remove(checkm_id_list[-1])
        checkm_id='.'.join(checkm_id_list)
        f.write(checkm_id+'\t'+str(final_iteration_checkm[item])+'\n')
        bin_selected[checkm_id]=0
    f.close()

    os.chdir(pwd+'/'+'Iteration_'+str(num)+'_genomes')
    for root, dirs, files in os.walk(pwd+'/'+'Iteration_'+str(num)+'_genomes'):
        for file in files:
            if file in final_iteration_checkm.keys():
                os.system('cp '+file+' '+pwd+'/'+str(assembly)+'_BestBinsSet')

    os.chdir(pwd+'/'+str(assembly)+'_BestBinsSet')
    f3=open(str(assembly)+'_BestBinsSet.depth.txt','w')
    f4=open('prebinned_genomes_output_for_dataframe_'+str(assembly)+'_BestBinsSet.txt','w')
    f5=open('Genome_group_all_list_'+str(assembly)+'_BestBinsSet.txt','w')

    binset_folders={}
    for item in final_iteration_checkm.keys():
        binset_folders[item.split('_genomes.')[0]]=0

    print('----------------------------')
    print('Parsing bins in best bin-set')
    contig_bin={}
    for root, dirs, files in os.walk(pwd+'/'+str(assembly)+'_BestBinsSet'):
        for file in files:
            if '_genomes.' in file:
                hz=str(file).split('_genomes.')[1]
                if '.fasta' in hz or '.fa' in hz:
                    for record in SeqIO.parse(file, 'fasta'):
                        contig_bin[record.id]=str(file)

    title1, title2, title3=[], [], []
    for item in binset_folders.keys():
        print('Parsing files in folder '+item)
        for root, dirs, files in os.walk(pwd+'/'+str(item)+'_genomes'):
            os.chdir(pwd+'/'+str(item)+'_genomes')
            for file in files:
                if '.depth.txt' in file:
                    n=0
                    for line in open(file, 'r'):
                        n+=1
                        if n == 1 and len(title1) == 0:
                            title1.append(str(line))
                            f3.write(str(line))
                        else:
                            if str(line).strip().split('\t')[0] in contig_bin.keys():
                                f3.write(str(line))

                if 'prebinned_genomes_output_for_dataframe_' in file:
                    n=0
                    for line in open(file, 'r'):
                        n+=1
                        if n == 1 and len(title2) == 0:
                            title2.append(str(line))
                            f4.write(str(line))
                        else:
                            if str(line).strip().split('\t')[1] in bin_selected.keys():
                                f4.write(str(line))
                
                if 'Genome_group_all_list_' in file:
                    n=0
                    for line in open(file, 'r'): 
                        n+=1
                        if n == 1 and len(title3) == 0:
                            f5.write(str(line))
                            title3.append(str(line))
                        else:
                            if str(line).strip().split('\t')[0] in bin_selected.keys():
                                f5.write(str(line))

    f3.close()
    f4.close()
    f5.close()
    os.chdir(pwd)
    try:
        f_s3=open('S3_checkpoint.txt','a')
    except:
        xyzzzz=0
    num2=num+1
    f_s3.write(str(assembly)+'\t'+str(num2)+'\t'+str(assembly)+'_BestBinsSet done'+'\n')
    return str(assembly)+'_BestBinsSet'

def binset_filtration(binset):
    print('Parsing '+binset)
    pwd=os.getcwd()
    os.chdir(pwd+'/'+binset)
    del_bin, bin_checkm=[], []
    for root, dirs, files in os.walk(pwd+'/'+binset):
        for file in files:
            if '_bin_stats_ext.tsv' in file:
                #xxx=0
                for line in open(file, 'r'):
                    #xxx+=1
                    bin_id=str(line).strip().split('\t')[0]
                    if '_maxbin2_genomes.' in bin_id or '_concoct_genomes.' in bin_id:
                        bin_id_f=bin_id+'.fasta'
                        bin_checkm.append(bin_id_f)
                    elif '_metabat_genomes.' in bin_id or '_metabinner_genomes.' in bin_id or '_vamb_genomes.' in bin_id  or 'semibin' in bin_id or 'SingleContig' in bin_id:
                        bin_id_f=bin_id+'.fa'
                        bin_checkm.append(bin_id_f)
                    else:
                        continue
                    #print(str(xxx))
                    #print(str(line).strip().split('marker lineage\':')[1])
                    try:
                        marker_lineage=str(line).strip().split('marker lineage\': \'')[1].split('\',')[0]
                        Completeness=float(str(line).strip().split('Completeness\': ')[1].split(',')[0])
                        genome_size=int(str(line).strip().split('Genome size\':')[1].split(',')[0].strip())
                    except:
                        marker_lineage='root'
                        Completeness=0
                        genome_size=0
 
                    if marker_lineage == 'root':
                        del_bin.append(bin_id_f)

                    if Completeness <= 5:
                        del_bin.append(bin_id_f)
                
                    if genome_size <= 200000:
                        del_bin.append(bin_id_f)
            
        for file in files:
            hz=file.split('.')[-1]
            if 'fasta' in hz or 'fa' in hz:
                if file not in bin_checkm:
                    del_bin.append(file)
            else:
                continue

    # fx=open('test.txt','w')
    # for item in bin_checkm:
    #     fx.write(str(item)+'\n')
    # fx.close()

    os.mkdir('Remove_bins')
    for root, dirs, files in os.walk(pwd+'/'+binset):
        for file in files:
            if file in del_bin:
                print(file)
                os.system('mv '+file+' Remove_bins')
    
    os.system('tar zcvf Remove_bins.tar.gz Remove_bins')
    os.system('rm -rf Remove_bins')

    os.chdir(pwd)

def bins_comparator_multiple_groups(genome_folder, assembly):
    try:
        f_s3=open('S3_checkpoint.txt','a')
    except:
        f_s3=open('S3_checkpoint.txt','w')
    
    finished_step=[]
    for line in open('S3_checkpoint.txt','r'):
        ass=str(line).strip().split('\t')[0].strip()
        if ass == assembly:
            finished_step.append(str(line).strip().split('\t')[1].strip())

    for item in genome_folder:
        binset_filtration(item)

    num_binset=len(genome_folder)
    num=0
    for num in range(0,num_binset-1):
        num+=1
        if num == 1:
            if str(num) not in finished_step:
                two_groups_comparator(assembly, str(genome_folder[0]).split('_genomes')[0], str(genome_folder[-1]).split('_genomes')[0], num)
        else:
            genome_folder.remove(genome_folder[0])
            genome_folder.remove(genome_folder[-1])
            genome_folder.append('Iteration_'+str(num-1))
            if str(num) not in finished_step:
                two_groups_comparator(assembly, str(genome_folder[-1]).split('_genomes')[0], str(genome_folder[0]).split('_genomes')[0], num)
    i=0
    while i < num_binset-2:
        i+=1
        os.system('rm -rf Iteration_'+str(i)+'_genomes')

    num2=num+1
    if str(num2) not in finished_step:
        bestbinset=bin_within_a_group_comparitor('Iteration_'+str(num)+'_genomes', str(assembly), num)    
    
    os.system('rm -rf Iteration_'+str(num)+'_genomes')
    os.system('mkdir '+str(assembly)+'_comparison_files')
    os.system('mv Bins_marker_lineage_completeness_contamination_iteration_* All_possible_bin_sets_iteration_* Best_hit_bin_sets_iteration_* Extract_bins_in_iteration_* Best_bin_set_iteration_* '+str(assembly)+'_comparison_files')

    # try:
    #     f_s3=open('S3_checkpoint.txt','a')
    # except:
    #     xyztt=0
    # f_s3.write('3rd bin selection within multiple groups done!')
    print('Done!')
    return bestbinset

if __name__ == '__main__': 
    bins_folders_name_list=['1_adh_dn1_contigs.fasta_0.3_maxbin2_genomes', '1_adh_dn1_contigs.fasta_0.5_maxbin2_genomes', '1_adh_dn1_contigs.fasta_0.7_maxbin2_genomes', '1_adh_dn1_contigs.fasta_0.9_maxbin2_genomes', '1_adh_dn1_contigs.fasta_200_metabat_genomes', '1_adh_dn1_contigs.fasta_300_metabat_genomes', '1_adh_dn1_contigs.fasta_400_metabat_genomes', '1_adh_dn1_contigs.fasta_500_metabat_genomes', '1_adh_dn1_contigs.fasta_300_concoct_genomes', '1_adh_dn1_contigs.fasta_400_concoct_genomes']
    # depth_file='assembly.depth.txt'
    assembly='1_adh_dn1_contigs.fasta'
    # PE_connections_file='condensed.cytoscape.connections.tab'

    # binsabundance_pe_connections(bins_folders_name_list, depth_file, PE_connections_file)
    bins_comparator_multiple_groups(bins_folders_name_list, assembly)
