#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

from lib2to3.fixes import fix_buffer
from Bio import SeqIO
import sys, os, time
from collections import Counter
from multiprocessing import Pool

def metabinner(assembly_file, depth_file, num_threads, ram, pwd, QC_software):
    condaenv=os.popen('conda info --envs').read()
    a=condaenv.split('\n')
    # print(str(a))
    for item in a:
        a=str(item).split('/')
        # name=a[0].strip()
        if 'BASALT' in item:
            path='/'.join(a[1:len(a)])
            path='/'+path+'/bin/'

    os.system('cat '+str(depth_file)+' | cut -f -1,4- > '+str(assembly_file)+'_coverage_profile.tsv')
    # os.system('python gen_kmer.py '+pwd+'/'+str(assembly_file)+' 1000 4')
    os.system(str(path)+'scripts/gen_kmer.py '+pwd+'/'+str(assembly_file)+' 1000 4')
    assembly_name_list=assembly_file.split('.')
    assembly_name_list.remove(assembly_name_list[-1])
    assembly_name='.'.join(assembly_name_list)

    kmer_file=str(assembly_name)+'_kmer_4_f1000.csv'
    n, contig_id = 0, {}
    for line in open(kmer_file,'r'):
       n+=1
       if n >= 2:
           contig_id[str(line).split(',')[0]]=0
           
    fcovout=open(str(assembly_file)+'_coverage_profile2.tsv','w')
    n=0
    for line in open(str(assembly_file)+'_coverage_profile.tsv','r'):
        n+=1
        if n == 1:
            fcovout.write(line)
        else:
            c_id=str(line).split('\t')[0]
            if c_id in contig_id.keys():
                fcovout.write(line)
    fcovout.close()
    os.system('mv '+str(assembly_file)+'_coverage_profile2.tsv '+str(assembly_file)+'_coverage_profile.tsv')
    os.system(path+'/run_metabinner.sh -a '+pwd+'/'+str(assembly_file)+' -o '+pwd+'/'+str(assembly_file)+'_metabinner -d '+pwd+'/'+str(assembly_file)+'_coverage_profile.tsv -k '+pwd+'/'+str(assembly_name)+'_kmer_4_f1000.csv -t '+str(num_threads)+' -p '+path)
    metabinner_bin_contig, mbn={}, {}
    for line in open(pwd+'/'+str(assembly_file)+'_metabinner/metabinner_res/metabinner_result.tsv','r'):
        bin_id=str(line).strip().split('\t')[1]
        contig=str(line).strip().split('\t')[0]
        metabinner_bin_contig[contig]=bin_id
        mbn[str(assembly_file)+'_100_metabinner_genomes.'+str(bin_id)+'.fa']=1

    for record in SeqIO.parse(assembly_file,'fasta'):
        try:
            fmbn=open(str(assembly_file)+'_100_metabinner_genomes.'+str(metabinner_bin_contig[record.id])+'.fa','a')
            fmbn.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            fmbn.close()
        except:
            if str(record.id) in metabinner_bin_contig.keys():
                fmbn=open(str(assembly_file)+'_100_metabinner_genomes.'+str(metabinner_bin_contig[record.id])+'.fa','w')
                fmbn.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                fmbn.close()          
    
    os.system('mkdir '+str(assembly_file)+'_100_metabinner_genomes')
    for item in mbn.keys():
        # print(item)
        os.system('mv '+str(item)+' '+pwd+'/'+str(assembly_file)+'_100_metabinner_genomes')
    if QC_software == 'checkm':
        os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(assembly_file)+'_100_metabinner_genomes '+str(assembly_file)+'_100_metabinner_checkm')
    elif QC_software == 'checkm2':
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(assembly_file)+'_100_metabinner_genomes  -x fa -o '+str(assembly_file)+'_100_metabinner_checkm')
    os.system('rm -rf '+str(assembly_file)+'_metabinner '+str(assembly_file)+'_coverage_profile.tsv '+str(assembly_name)+'_kmer_4_f1000.csv')

def vamb(assembly_file, datasets, num_threads, pwd, QC_software):
    assembly_num=str(assembly_file).split('_')[0]
    for i in range(1, len(datasets)+1):
        if i == 1:
            bam_list=str(assembly_num)+'_DNA-'+str(i)+'.bam'
        else:
            bam_list+=' '+str(assembly_num)+'_DNA-'+str(i)+'.bam'
    os.system('vamb --outdir '+str(assembly_file)+'_vamb --fasta '+str(assembly_file)+' --bamfiles '+str(bam_list)+' --minfasta 500000')
    # os.system('vamb --outdir '+str(assembly_file)+'_100_vamb_genomes --fasta '+str(assembly_file)+' --bamfiles '+str(bam_list)+' -o C')

    vamb_bin_contig, vbn={}, {}
    for line in open(pwd+'/'+str(assembly_file)+'_vamb/clusters.tsv','r'):
        bin_id=str(line).strip().split('\t')[0]
        contig=str(line).strip().split('\t')[1]
        vamb_bin_contig[contig]=bin_id
        vbn[str(assembly_file)+'_100_vamb_genomes.'+str(bin_id)+'.fa']=1

    for record in SeqIO.parse(assembly_file,'fasta'):
        try:
            fmbn=open(str(assembly_file)+'_100_vamb_genomes.'+str(vamb_bin_contig[record.id])+'.fa','a')
            fmbn.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            fmbn.close()
        except:
            if str(record.id) in vamb_bin_contig.keys():
                fmbn=open(str(assembly_file)+'_100_metabinner_genomes.'+str(vamb_bin_contig[record.id])+'.fa','w')
                fmbn.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                fmbn.close()       

    os.system('mkdir '+str(assembly_file)+'_100_vamb_genomes')
    for item in vbn.keys():
        # print(item)
        os.system('mv '+str(item)+' '+pwd+'/'+str(assembly_file)+'_100_vamb_genomes')
    
    os.chdir(pwd+'/'+str(assembly_file)+'_100_vamb_genomes')
    for root, dirs, files in os.walk(pwd+'/'+str(assembly_file)+'_100_vamb_genomes'):
        for file in files:
            bin_len=0
            for record in SeqIO.parse(file,'fasta'):
                bin_len+=len(record.seq)
            if bin_len < 500000:
                os.system('rm '+str(file))
    os.chdir(pwd)

    if QC_software == 'checkm':
        os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(assembly_file)+'_100_vamb_genomes '+str(assembly_file)+'_100_vamb_checkm')
    elif QC_software == 'checkm2':
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(assembly_file)+'_100_vamb_genomes  -x fa -o '+str(assembly_file)+'_100_vamb_checkm')
    # os.system('rm *.seed *.out *.err *.nto *.gff *.ffn *.faa *.ndb *.njs *.not *.ntf')
    os.system('rm -rf '+str(assembly_file)+'_vamb')

def extra_binner(binner, datasets, assembly_file, depth_file, num_threads, ram, pwd, QC_software):
    extra_bin_folder=[]
    if binner == 'm':
        metabinner(assembly_file, depth_file, num_threads, ram, pwd, QC_software)
        extra_bin_folder.append(str(assembly_file)+'_100_metabinner_genomes')
    elif binner == 'v':
        vamb(assembly_file, datasets, num_threads, pwd, QC_software)
        extra_bin_folder.append(str(assembly_file)+'_100_vamb_genomes')
    os.system('rm *.seed *.out *.err *.nto *.gff *.ffn *.faa *.ndb *.njs *.not *.ntf')
    return extra_bin_folder

if __name__ == '__main__': 
    num_threads=20
    ram=250
    pwd=os.getcwd()
    assembly_file='1_assembly_sample1.fa'
    depth_file='1_assembly.depth.txt'
    datasets={'1':['sample1.R1.fq','sample1.R2.fq'], '2':['sample2.R1.fq','sample2.R2.fq']}
    binner='v' ### 'm': metabinner; 'v': vamb
    QC_software='checkm2' ### checkm or checkm2
    extra_binner(binner, datasets, assembly_file, depth_file, num_threads, ram, pwd, QC_software)