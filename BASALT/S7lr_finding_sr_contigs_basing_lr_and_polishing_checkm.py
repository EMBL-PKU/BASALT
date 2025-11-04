#!/usr/bin/env python
from ast import excepthandler
from imp import NullImporter
from selectors import EpollSelector
from telnetlib import SE
from xmlrpc.server import SimpleXMLRPCRequestHandler
from Bio import SeqIO
import sys, os, threading, copy, math
from multiprocessing import Pool


def record_bin(binset_folder):
    pwd=os.getcwd()
    os.chdir(pwd+'/'+binset_folder)
    record_seq, record_bin_seq = {}, {}
    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            if 'fa' in hz:
                record_bin_seq[file]={}
                for record in SeqIO.parse(file, 'fasta'):
                    record_bin_seq[file]
                    try:
                        record_seq[record.seq]+='|'+str(record.id)
                    except:
                        record_seq[record.seq]=str(record.id)

    os.chdir(pwd)
    f=open('Total_bins.fa','w')
    for seqs in record_seq.keys():
        f.write('>'+str(record_seq[seqs])+'\n'+str(seqs)+'\n')
    f.close()

def mod_bin(binset_folder):
    pwd=os.getcwd()
    bins_checkm={}
    try:
        os.mkdir(str(binset_folder)+'_mod')
    except:
        print(str(binset_folder)+'_mod is exist. Re-create the folder')
        os.system('rm -rf '+str(binset_folder)+'_mod')
        os.mkdir(str(binset_folder)+'_mod')
    
    os.chdir(pwd+'/'+binset_folder)
    f=open('Mod_contig_id.txt','w')
    f1=open('Bin_name_mod.txt','w')
    f2=open('Bin_name_mod_bin_stats_ext.tsv','w')
    n, record_seq, mod_bin_list, mod_bin_dict, bin_contig_num, bin_name_change = 0, {}, [], {}, {}, {}
    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            hz=str(file).split('.')[-1]
            # print(hz
            if 'fa' in hz:
                n+=1
                m=0
                record_seq['bin'+str(n)]={}
                mod_bin_list.append('bin'+str(n))
                checkm_name_list=str(file).split('.')
                checkm_name_list.remove(checkm_name_list[-1])
                checkm_name='.'.join(checkm_name_list)
                mod_bin_dict[checkm_name]='bin'+str(n)
                bin_name_change[str(file).split('_')[0]]='bin'+str(n)
                f1.write(str(file)+'\t'+'bin'+str(n)+'.fa'+'\n')
                for record in SeqIO.parse(file, 'fasta'):
                    m+=1
                    f.write('bin'+str(n)+'_'+str(m)+'\t'+str(record.id)+'\n')
                    record_seq['bin'+str(n)]['bin'+str(n)+'_'+str(m)]=str(record.seq)
                bin_contig_num['bin'+str(n)]=m

    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            if 'bin_stats_ext.tsv' in str(file) and str(file) != 'Bin_name_mod_bin_stats_ext.tsv':
                for line in open(str(file),'r'):
                    checkm_name=str(line).strip().split('\t')[0]
                    mod_bin_checkm_name=str(mod_bin_dict[checkm_name])
                    bins_checkm[mod_bin_checkm_name]={}
                    bins_checkm[mod_bin_checkm_name]['marker lineage']=str(line).strip().split('\'marker lineage\': \'')[1].strip().split('\'')[0]
                    bins_checkm[mod_bin_checkm_name]['Completeness']=float(str(line).strip().split('\'Completeness\': ')[1].split(', ')[0])
                    bins_checkm[mod_bin_checkm_name]['Genome size']=float(str(line).strip().split('\'Genome size\':')[1].strip().split(', ')[0].replace('\'','').replace('\"','').replace('}',''))
                    try:
                        #bins_checkm[mod_bin_checkm_name]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split(',')[0].strip())
                        bins_checkm[mod_bin_checkm_name]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0])
                    except:
                        #bins_checkm[mod_bin_checkm_name]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split('}')[0].strip())
                        bins_checkm[mod_bin_checkm_name]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split(',')[0].strip())
                    bins_checkm[mod_bin_checkm_name]['Mean scaffold length']=float(bins_checkm[mod_bin_checkm_name]['Genome size'])/bin_contig_num[mod_bin_checkm_name]
                    f2.write(mod_bin_checkm_name+'\t'+str(line).strip().split('\t')[1]+'\n')
    f.close()
    f1.close()
    f2.close()

    os.system('mv Mod_contig_id.txt Bin_name_mod.txt Bin_name_mod_bin_stats_ext.tsv '+pwd+'/'+str(binset_folder)+'_mod')
    
    os.chdir(pwd+'/'+str(binset_folder)+'_mod')
    f1=open('Total_bins.fa','w')
    for bin_id in record_seq.keys():
        f=open(str(bin_id)+'.fa','w')
        for contigs in record_seq[bin_id].keys():
            f.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
            f1.write('>'+str(contigs)+'\n'+str(record_seq[bin_id][contigs])+'\n')
        f.close()
    f1.close()
    os.system('mv Total_bins.fa '+pwd)
    os.chdir(pwd)
    return str(binset_folder)+'_mod', 'Total_bins.fa', mod_bin_list, bins_checkm, bin_name_change

def parse_sam_bwa(sam_file, fq, pair, n, batch, mp_run):
    print('Reading reads id')
    f_not_mapped_reads=open('Not_mapped_reads_'+str(batch)+'_'+str(mp_run)+'.txt','a')
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
    f_summary=open('Bin_reads_summary.txt','a')
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
    # return fq

def parse_sam(sam_file, fq, pair, n, batch, mp_run):
    print('Reading reads id')
    f_not_mapped_reads=open('Not_mapped_reads_'+str(batch)+'_'+str(mp_run)+'.txt','a')
    m, m1, m2 = 0, 0, 0
    for line in open(sam_file,'r'):
        m1+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id=flist[2].split('_')[0]
            read_id=flist[0]
            read_id_name=str(n)+'_'+read_id.split('_')[0]
            pair[bin_id][read_id_name]=0
        if m1 % 1000000 == 0:
            print('Read', m1,'lines')

    print('Collecting paired reads')
    for line in open(sam_file,'r'):
        m2+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id=flist[2].split('_')[0] ###
            read_id=flist[0] ###
            read_id_name=str(n)+'_'+read_id.split('_')[0]
            pair[bin_id][read_id_name]+=1

            if pair[bin_id][read_id_name] == 2:
                fq[bin_id][read_id_name]=0
                del pair[bin_id][read_id_name]

        if m2 % 1000000 == 0:
            print('Read', m2,'lines')

    f_summary=open('Bin_reads_summary.txt','w')
    for item in fq.keys():
        f_summary.write(str(item)+' SEQ number:'+str(len(fq[item]))+'\n')
    f_summary.close()

    print('Parsing', sam_file)
    for line in open(sam_file,'r'):
        m+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id=flist[2].split('_')[0]
            read_id=flist[0]
            reads=str(n)+'_'+read_id #
            # read_num=reads.split('_')[-1] #
            read_id_name=str(n)+'_'+read_id.split('_')[0] #
            fq_seq=flist[9]+'\n'+'+'+'\n'+flist[10]+'\n'
            try:
                fq[bin_id][read_id_name]+=1
                i=fq[bin_id][read_id_name]
                f1=open(str(bin_id)+'_seq_R'+str(i)+'.fq','a')
                f1.write('@'+str(reads)[:-2]+' '+str(i)+'\n'+str(fq_seq))
                f1.close()
                if i == 2:
                    del fq[bin_id][read_id_name]
            except:
                f_not_mapped_reads.write(str(read_id_name)+'\n')
    
        if m % 1000000 == 0:
            print('Parsed', m,'lines')
    f_not_mapped_reads.close()
    # return fq

def parse_lr_sam(sam_file, long_read, sn):
    print('Reading long reads id '+str(long_read))
    # f_not_mapped_reads=open('Not_mapped_reads.txt','w')
    bin_lr, bin_lr2, lr_bin, lr_bin2 = {}, {}, {}, {}
    m, m1, m2 = 0, 0, 0
    for line in open(sam_file,'r'):
        m1+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            read_id=flist[0]
            bin_id=flist[2].split('_')[0]
            bin_lr[bin_id]={}
            lr_bin[read_id]={}

        if m1 % 1000000 == 0:
            print('Read', m1,'lines')

    print('Collecting long reads '+str(long_read))
    for line in open(sam_file,'r'):
        m2+=1
        flist=str(line).split('\t')
        if len(flist) >= 12:
            bin_id=flist[2].split('_')[0]
            seqs=flist[9]
            if 'bin' in bin_id and len(seqs) >= 100:
                read_id=flist[0]
                bin_lr[bin_id][read_id]=''
                lr_bin[read_id][bin_id]=''

        if m2 % 1000000 == 0:
            print('Read', m2,'lines')

    bin_lr2=copy.deepcopy(bin_lr)
    lr_bin2=copy.deepcopy(lr_bin)

    for item in bin_lr2.keys():
        if len(bin_lr2[item]) == 0:
            del bin_lr[item]

    for item in lr_bin2.keys():
        if len(lr_bin2[item]) == 0:
            del lr_bin[item]

    f=open('Bin_long_read'+str(sn)+'.txt','w')
    for bin_id in bin_lr.keys():
        f.write(str(bin_id)+'\t'+str(bin_lr[bin_id])+'\n')
        try:
            f1=open(str(bin_id)+'_lr.fq','a')
        except:
            f1=open(str(bin_id)+'_lr.fq','w')
        f1.close()
    f.close()

    f=open('Long_read_bin'+str(sn)+'.txt','w')
    for lr in lr_bin.keys():
        f.write(str(lr)+'\t'+str(lr_bin[lr])+'\n')
    f.close()

    n, record_bin_line, record_bin_line2 = 0, {}, {}
    for line in open(long_read,'r'):
        n+=1
        record_bin_line[n]=[]

    print('Splitting long reads to different bins '+str(long_read))
    n = 0
    for line in open(long_read,'r'):
        n+=1
        m=n-1
        if m % 4 == 0:
            seq_id=str(line).strip().split(' ')[0].split('@')[1]
            if seq_id in lr_bin.keys():
                for bin_id in lr_bin[seq_id].keys():
                    record_bin_line[n].append(bin_id)
                    record_bin_line[n+1].append(bin_id)
                    record_bin_line[n+2].append(bin_id)
                    record_bin_line[n+3].append(bin_id)

    record_bin_line2=copy.deepcopy(record_bin_line)
    for line in record_bin_line2.keys():
        if len(record_bin_line2[line]) == 0:
            del record_bin_line[line]
    seq_num=int(len(record_bin_line)/4)
    print(str(seq_num)+' reads from '+str(long_read)+' will be splitted into different bins')

    n1=0
    for line in open(long_read,'r'):
        n1+=1
        if n1 in record_bin_line.keys():
            for bin_id in record_bin_line[n1]:
                f1=open(str(bin_id)+'_lr.fq','a')
                f1.write(line)
                f1.close()
    if n1 % 1000000 == 0:
        print('Parse', n1,'lines')
    print('Long reads '+str(long_read)+' splitting done!')
    return bin_lr

def mapping_sr(total_fa, datasets_list, fq, pair, mapping_tool, num_threads, batch, mp_run):
    f=open('Not_mapped_reads_'+str(batch)+'_'+str(mp_run)+'.txt','w')
    f.close()

    ### bowtie2
    if mapping_tool == 'bw2':
        os.system('bowtie2-build '+str(total_fa)+' '+str(total_fa))
        n = 0
        for item in datasets_list.keys():
            print('Mapping '+str(datasets_list[item]))
            n+=1
            os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(datasets_list[item][0])+' -2 '+str(datasets_list[item][1])+' -S '+str(item)+'.sam -q --no-unal') ### bowtie2
            parse_sam(str(item)+'.sam', fq, pair, n, 1, mp_run)
            os.system('rm '+str(item)+'.sam')
    elif mapping_tool == 'bwa':
        ### bwa
        os.system('bwa index '+str(total_fa))

        n = 0
        for item in datasets_list.keys():
            print('Mapping '+str(datasets_list[item]))
            n+=1
            os.system('bwa mem -t '+str(num_threads)+' '+str(total_fa)+' '+str(datasets_list[item][0])+' '+str(datasets_list[item][1])+' > '+str(item)+'.sam') ### bwa
            parse_sam_bwa(str(item)+'.sam', fq, pair, n, 2, mp_run)
            os.system('rm '+str(item)+'.sam')

def sr_polishing(bin_id, bin_seq, target_folder, sr_folder, num_threads, pwd, binset_folder, type, initial_bin_name, mapping_tool, mp_run, major_run, batch):
    r1=bin_seq[bin_id][0]
    r2=bin_seq[bin_id][1]
    xyz, prun=1, 0
    if type == 'lr_polishing':
        bin_name=initial_bin_name
        output_name_r=str(bin_id)+'_lr_polished'
        output_final=str(bin_id)+'_gf_lr_polished.fa'
        output_folder=pwd+'/'+binset_folder+'_gf_lr_polished'
        ss='_lr_polished'
    elif type == 'MAG_polishing':
        bin_name=str(initial_bin_name)
        output_name_r=str(bin_id)+'_mag_polished'
        output_final=str(bin_id)+'_mag_polished.fa'
        output_folder=pwd+'/'+binset_folder+'_MAGs_polished'
        ss='_mag_polished'
    
    if batch == 1:
        if mapping_tool == 'bw2':
            if mp_run == major_run:
                schedual_polish_iteration = 2
            elif mp_run < major_run:
                schedual_polish_iteration = 10
        elif mapping_tool == 'bwa':
            if mp_run == major_run:
                schedual_polish_iteration = 30
            elif mp_run < major_run:
                schedual_polish_iteration = 10

    else: ### batch == 2 or batch == 3
        if mp_run == major_run:
            schedual_polish_iteration = 20
        else:
            schedual_polish_iteration = 10

    try:
        iteration_num=int(str(initial_bin_name).split('_polished')[1].split('.fa')[0])
        prun=iteration_num
        print(bin_id+' start at iteration '+str(iteration_num))
        if iteration_num >= schedual_polish_iteration:
            print('Trying to move '+str(bin_id)+' to the final binset folder')
            os.system('cat '+str(initial_bin_name)+' '+str(output_final)+' > '+str(output_final)+'_m')
            # os.system('cat '+str(initial_bin_name)+'_bk '+str(output_final)+'_bk > '+str(output_final)+'_m')
            os.system('mv '+str(initial_bin_name)+' '+str(initial_bin_name)+'_bk')
            os.system('mv '+str(output_final)+' '+str(output_final)+'_bk')
            os.system('mv '+str(output_final)+'_m '+str(output_final))
            os.system('mv '+str(output_final)+' '+str(output_folder))
            f=open('Accompleted_bins_'+str(batch)+'_'+str(mp_run)+'.txt','a')
            f.write(str(output_final)+'\n')
            f.close()
            f=open('BASALT_log.txt','a')
            f.write('Finished '+str(bin_id)+' in run'+str(iteration_num)+'\n')
            f.close()
            xyz=0
    except:
        xyz=1

    while xyz == 1:
        prun+=1
        if prun <= schedual_polish_iteration:
            print('Processing '+str(bin_id)+' in polishing '+str(prun))
            output_name=output_name_r+str(prun)
            if prun == 1:
                try:
                    os.system('cp '+str(pwd)+'/'+str(target_folder)+'/'+str(bin_name)+' '+str(bin_name))
                except:
                    try:
                        prun=int(str(bin_name).split(str(ss))[1].split('.fa')[0])
                        output_name=output_name_r+str(prun)
                    except:
                        print('Error in polishing '+str(bin_id))
                        f=open('BASALT_log.txt','a')
                        f.write('Error in polishing '+str(bin_id)+'\n')
                        f.close()

            ###
            bin_contigs={}
            try:
                for record in SeqIO.parse(str(bin_name),'fasta'):
                    bin_contigs[record.seq]=0
            except:
                xyzzz=0
            # num=len(bin_contigs)
            ### bwa
            os.system('bwa index '+str(bin_name))
            os.system('bwa mem -t '+str(num_threads)+' '+str(bin_name)+' '+str(pwd)+'/'+str(sr_folder)+'/'+str(r1)+' '+str(pwd)+'/'+str(sr_folder)+'/'+str(r2)+' | samtools sort -@ '+str(num_threads)+' -O bam -o '+str(bin_id)+'.bam')
            ### bowtie2
            # os.system('samtools sort -@ 10 -O bam -o ONTmin_IT3_old.bam')
            os.system('samtools view -@ '+str(num_threads)+' -q 30 -b '+str(bin_id)+'.bam > '+str(bin_id)+'_sort.bam')
            os.system('samtools index -@  '+str(num_threads)+' '+str(bin_id)+'_sort.bam')
            # short read consensus call
            os.system('pilon --genome '+str(bin_name)+' --frags '+str(bin_id)+'_sort.bam --fix all --output '+str(output_name))

            redo, ix = 0, 0
            while ix < 3 and redo == 0:
                ix += 1
                try:
                    xyzzzz=0
                    for record in SeqIO.parse(str(output_name)+'.fasta','fasta'):
                        xyzzzz+=1
                        if xyzzzz == 4:
                            break
                    redo=1
                except:
                    print('Pilon polishing failed. Going to redo')
                    os.system('pilon --genome '+str(bin_name)+' --frags '+str(bin_id)+'_sort.bam --fix all --output '+str(output_name))

            ###
            remain_contigs, remain_contigs2, dp = {}, {}, 1
            try:
                xyzzzz=0
                for record in SeqIO.parse(str(output_name)+'.fasta','fasta'):
                    xyzzzz+=1
                    if xyzzzz == 4:
                        break
                
                fp=open(output_final,'a')
                for record in SeqIO.parse(str(output_name)+'.fasta','fasta'):
                    try:
                        bin_contigs[record.seq]+=1
                        fp.write('>'+str(str(record.id).replace('_pilon',''))+'\n'+str(record.seq)+'\n')
                        del bin_contigs[record.seq]
                    except:
                        remain_contigs[record.seq]=0
                        remain_contigs2[record.seq]=str(record.id)
                fp.close()
            except:
                print('Pilon polishing failed. Skip this bin '+str(str(output_name_r)))
                dp = 0
                fxyzz=open('BASALT_log.txt','a')
                fxyzz.write('Pilon polishing failed. Skip this bin '+str(str(output_name_r))+'\n')
                fxyzz.close()

            # if prun >= 4:
            #     prun2=int(prun)-2
            #     prun3=int(prun)-3
            #     output_name_pre2=output_name_r+str(prun2)+'.fasta'
            #     try:
            #         for record in SeqIO.parse(output_name_pre2,'fasta'):
            #             try:
            #                 remain_contigs[record.seq]+=1
            #                 # fp.write('>'+str(str(remain_contigs2[record.seq]).replace('_pilon',''))+'\n'+str(record.seq)+'\n')
            #                 # del remain_contigs[record.seq]
            #                 # del remain_contigs2[record.seq]
            #             except:
            #                 xyzzz=0
            #     except:
            #         xyzzz=0
                
            #     output_name_pre3=output_name_r+str(prun3)+'.fasta'
            #     try:
            #         for record in SeqIO.parse(output_name_pre3,'fasta'):
            #             try:
            #                 remain_contigs[record.seq]+=1
            #                 if remain_contigs[record.seq] == 2:
            #                     fp.write('>'+str(str(remain_contigs2[record.seq]).replace('_pilon',''))+'\n'+str(record.seq)+'\n')
            #                     del remain_contigs[record.seq]
            #                     del remain_contigs2[record.seq]
            #             except:
            #                 xyzzz=0
            #     except:
            #         xyzzz=0
            # fp.close()

            # if prun == 25 or prun == 10 or prun == 40 or prun == 50:
            #     os.system('cp '+str(output_name)+'.fasta '+str(output_name)+'.bk')
            # if prun >= 4:    
            #     os.system('rm '+str(output_name_pre2)+' '+str(output_name_pre3))
            os.system('rm '+str(bin_name)+'.sa '+str(bin_name)+'.pac '+str(bin_name)+'.ann '+str(bin_name)+'.amb '+str(bin_name)+'.bwt '+str(bin_id)+'_sort.bam.bai '+str(bin_id)+'_sort.bam '+str(bin_id)+'.bam')
            ###
            drep={}
            if len(remain_contigs) != 0:
                xyz=1
                bin_name=str(output_name)+'.fasta'
                if prun == schedual_polish_iteration:
                    print('Finished '+str(bin_id)+' in run'+str(prun))

                    fo=open(str(output_final),'a')
                    for seq in remain_contigs2.keys():
                        fo.write('>'+str(str(remain_contigs2[seq]).replace('_pilon',''))+'\n'+str(seq)+'\n')
                    fo.close()
                    
                    os.system('mv '+str(output_final)+' '+str(output_final)+'_org')
                    fo=open(str(output_final),'w')
                    for record in SeqIO.parse(str(output_final)+'_org','fasta'):
                        try:
                            drep[record.id]+=1
                        except:
                            drep[record.id]=0
                            fo.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    fo.close()
                    os.system('rm '+str(output_final)+'_org')
                    os.system('mv '+str(output_final)+' '+str(output_folder))

                    f=open('Accompleted_bins_'+str(batch)+'_'+str(mp_run)+'.txt','a')
                    f.write(str(output_final)+'\n')
                    f.close()
                    f=open('BASALT_log.txt','a')
                    f.write('Finished '+str(bin_id)+' in run'+str(prun)+'\n')
                    f.close()
                else:
                    fpo=open(str(output_name)+'.fasta','w')
                    for seq in remain_contigs2.keys():
                        fpo.write('>'+str(remain_contigs2[seq])+'\n'+str(seq)+'\n')
                    fpo.close()
            else:
                if dp == 1:
                    xyz=0
                    print('Finished '+str(bin_id)+' in run'+str(prun))

                    os.system('mv '+str(output_final)+' '+str(output_final)+'_org')
                    fo=open(str(output_final),'w')
                    for record in SeqIO.parse(str(output_final)+'_org','fasta'):
                        try:
                            drep[record.id]+=1
                        except:
                            drep[record.id]=0
                            fo.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    fo.close()
                    os.system('rm '+str(output_final)+'_org')
                    os.system('mv '+str(output_final)+' '+str(output_folder))

                    f=open('BASALT_log.txt','a')
                    f.write('Finished '+str(bin_id)+' in run'+str(prun)+'\n')
                    f.close()
                    f=open('Accompleted_bins_'+str(batch)+'_'+str(mp_run)+'.txt','a')
                    f.write(str(output_final)+'\n')
                    f.close()
                elif dp == 0:
                    f=open('BASALT_log.txt','a')
                    f.write('Unable to polish '+str(bin_id)+'\n')
                    f.close()
                    f=open('Accompleted_bins_'+str(batch)+'_'+str(mp_run)+'.txt','a')
                    f.write(str(output_final)+'\n')
                    f.close()
                    os.system('mv '+str(bin_id)+'.fa '+str(output_folder)+'/'+str(bin_id)+'_mag_polished.fa')

        else:
            xyz=0
            
def gf_lr_direct_polishing(gf_lr_fa, datasets_list, num_threads, pwd):
    fa_name=gf_lr_fa
    xyz, tp = 1, 0
    while xyz == 1:
        tp+=1
        if tp <= 20:
            print('Processing polishing iteration '+str(tp))
            record_seqs={}
            print('Recording reads from '+str(fa_name))
            for record in SeqIO.parse(fa_name,'fasta'):
                record_seqs[record.seq]=record.id
            x=len(record_seqs)
            print('There is '+str(x)+' seqs in '+str(fa_name))

            r=0
            for ds in datasets_list.keys():
                r+=1
                if r == 1:
                    output_name=fa_name+'_'+str(r)
                else:
                    fa_name_list=fa_name.split('_')
                    fa_name_list.remove(fa_name_list[-1])
                    fa_name2='_'.join(fa_name_list)
                    output_name=fa_name2+'_'+str(r)

                print('Processed polishing iteration '+str(tp)+' run '+str(r)+'. Input file is '+str(fa_name))
                os.system('bwa index '+str(fa_name))
                os.system('bwa mem -t '+str(num_threads)+' '+str(fa_name)+' '+str(datasets_list[ds][0])+' '+str(datasets_list[ds][1])+' | samtools sort -@ '+str(num_threads)+' -O bam -o '+str(fa_name)+'.bam')
                ### bowtie2
                # os.system('samtools sort -@ 10 -O bam -o ONTmin_IT3_old.bam')
                os.system('samtools view -@ '+str(num_threads)+' -q 30 -b '+str(fa_name)+'.bam > '+str(fa_name)+'_sort.bam')
                os.system('samtools index -@  '+str(num_threads)+' '+str(fa_name)+'_sort.bam')
                # short read consensus call
                os.system('pilon --genome '+str(fa_name)+' --frags '+str(fa_name)+'_sort.bam --fix all --output '+str(output_name))
                os.system('rm '+str(fa_name)+'.sa '+str(fa_name)+'.pac '+str(fa_name)+'.bwt '+str(fa_name)+'.ann '+str(fa_name)+'.amb '+str(fa_name)+'_sort.bam.bai '+str(fa_name)+'_sort.bam '+str(fa_name)+'.bam')
                try:
                    xyz=0
                    for record in SeqIO.parse(str(output_name)+'.fasta','fasta'):
                        xyz+=1
                    # if r != 1 and tp != 1:
                    #     os.system('rm '+str(fa_name))
                    fa_name=str(output_name)+'.fasta'
                except:
                    fa_name=fa_name

                print('Processed polishing iteration '+str(tp)+' run '+str(r)+'. Output file is '+str(fa_name))

            print('Recording reads from '+str(fa_name))
            for record in SeqIO.parse(fa_name,'fasta'):
                record_seqs[record.seq]=record.id
                y=len(record_seqs)
            print('There is '+str(y)+' seqs in '+str(fa_name))

            if x != y:
                xyz = 1
                f=open(gf_lr_fa+'_tp'+str(tp)+'.fa','w')
                for record in SeqIO.parse(fa_name,'fasta'):
                    try:
                        record_id2=str(record.id).replace('_pilon','')
                    except:
                        record_id2=str(record.id)
                    f.write('>'+record_id2+'\n'+str(record.seq)+'\n')
                f.close()
                os.system('rm '+fa_name)
                fa_name=gf_lr_fa+'_tp'+str(tp)+'.fa'
            else:
                xyz=0
                os.system('mv '+fa_name+' '+gf_lr_fa+'_polished.fa')
                print('Rename '+fa_name+' to '+gf_lr_fa+'_polished.fa')
        else:
            xyz=0
            os.system('mv '+fa_name+' '+gf_lr_fa+'_polished.fa')
            print('Rename '+fa_name+' to '+gf_lr_fa+'_polished.fa')
    
    f=open('BASALT_log.txt','a')
    f.write('Selected and Renamed '+fa_name+' to '+gf_lr_fa+'_polished.fa'+'\n')
    f.close()  
            
def fq_2_fa(bin_fq):
    fa=bin_fq.split('_lr.fq')[0]+'_lr.fa'
    f_name=open(fa,'w')
    n=0
    for line in open(bin_fq,'r'):
        n+=1
        m=n%4
        if m == 1:
            f_name.write('>'+str(line[1:]))
        elif m == 2:
            f_name.write(line)
    f_name.close()
    return fa

def filtration(binset_folder, bin_name, lr_fa, num_threads, pwd):
    # os.system('cp '+str(pwd)+'/'+str(binset_folder)+'_mod/'+str(bin_name)+' '+str(bin_name))
    # os.system('makeblastdb -dbtype nucl -in '+str(bin_name))
    # os.system('blastn -query '+str(lr_fa)+' -db '+str(bin_name)+' -outfmt 6 -evalue 1e-20 -max_target_seqs 10 -num_threads '+str(num_threads)+' -out '+str(bin_name)+'_blasted.txt')

    contigs_len={}
    for record in SeqIO.parse(bin_name, 'fasta'):
        contigs_len[record.id]=len(record.seq)
    for record in SeqIO.parse(lr_fa, 'fasta'):
        contigs_len[record.id]=len(record.seq)

    ct, f_b, e_lr = {}, {}, {}
    for line in open(str(bin_name)+'_blasted.txt', 'r'):
        que=str(line).strip().split('\t')[0]
        sub=str(line).strip().split('\t')[1]
        simi=float(str(line).strip().split('\t')[2])
        leng=int(str(line).strip().split('\t')[3])
        q_s=int(str(line).strip().split('\t')[6])
        q_e=int(str(line).strip().split('\t')[7])
        s_s=int(str(line).strip().split('\t')[8])
        s_e=int(str(line).strip().split('\t')[9])
        
        if s_e > s_s:
            s_e2=s_e
            s_s2=s_s
        else:
            s_s2=s_e
            s_e2=s_s
        delta_q=s_e2-s_s2
        q_perc=leng/contigs_len[sub]
        
        if simi >= 85 and leng >= 200:
            q_e_n=contigs_len[que]-50
            s_e_n=contigs_len[sub]-50
            if q_s <= 50 or q_e >= q_e_n:
                if s_s2 <= 50 or s_e2 >= s_e_n:
                    try:
                        ct[que][sub]=0
                    except:
                        ct[que]={}
                        ct[que][sub]=0
                    f_b[line]=0

            if q_perc >= 0.8:
                e_lr[line]=0

    f=open(str(bin_name)+'_gap_filling_lr.txt', 'w')
    f_que={}
    for que in ct.keys():
        if len(ct[que]) >= 2:
            f_que[que]=0
    
    for item in f_b.keys():
        que=item.split('\t')[0]
        if que in f_que.keys():
            f.write(item)
    f.close()

    f=open(str(bin_name)+'_elongation_lr.txt', 'w')
    for item in e_lr.keys():
        f.write(item)
    f.close()
    # os.system('rm '+bin_name)

def blast(binset_folder, bin_name, lr_fa, pwd):
    # os.system('cp '+str(pwd)+'/'+str(binset_folder)+'/'+str(bin_name)+' '+str(bin_name))
    os.system('cp '+str(pwd)+'/'+str(binset_folder)+'_mod/'+str(bin_name)+' '+str(bin_name))
    os.system('makeblastdb -dbtype nucl -in '+str(bin_name))
    os.system('blastn -query '+str(lr_fa)+' -db '+str(bin_name)+' -outfmt 6 -evalue 1e-20 -max_target_seqs 10 -num_threads 1 -out '+str(bin_name)+'_blasted.txt')
    os.system('rm '+str(bin_name)+'.nin '+str(bin_name)+'.nhr '+str(bin_name)+'.nsq ')

def checkm(bin_folder, pwd):
    print('Parsing '+bin_folder+' checkm output')
    refined_checkm={}
    os.chdir(pwd+'/'+str(bin_folder))
    for root, dirs, files in os.walk(pwd+'/'+str(bin_folder)):
        for file in files:
            if 'bin_stats_ext.tsv' in file:
                for line in open(file,'r'):
                    binID=str(line).strip().split('{\'')[0].strip()
                    genome_size=str(line).strip().split('Genome size\':')[1].split(',')[0].replace('}','').strip()
                    taxon=str(line).strip().split('lineage')[1].split('\'')[2].strip()
                    completeness=str(line).strip().split('Completeness')[1].split(':')[1].split(',')[0].strip()
                    contamination=str(line).strip().split('Contamination')[1].split(':')[1].split(',')[0].replace('}','').strip()
                    # GC=round(float(str(line).strip().split('\'GC\':')[1].split(', \'GCN4\'')[0].strip())*100, 1)
                    refined_checkm[str(binID)]={}
                    refined_checkm[str(binID)]['marker lineage']=str(taxon)
                    refined_checkm[str(binID)]['Completeness']=float(completeness)
                    refined_checkm[str(binID)]['Genome size']=float(eval(genome_size))
                    refined_checkm[str(binID)]['Contamination']=float(contamination)
    os.chdir(pwd)
    return refined_checkm

def bin_comparison(binset_folder, target_folder, origin_checkm, refined_checkm, pwd, type):
    removed_bins={}
    os.chdir(pwd+'/'+target_folder)
    f=open('bin_stats_ext.tsv', 'w')
    f2=open('Record_log.txt','w')
    os.chdir(pwd)
    for bin_name in refined_checkm.keys():
        bin_name_o=str(bin_name).split('_'+str(type)+'_lr')[0]
        bin_name_o_list=bin_name_o.split('.')
        bin_name_o_list.remove(bin_name_o_list[-1])
        bin_name_o2='.'.join(bin_name_o_list)
        if float(origin_checkm[bin_name_o2]['Contamination']) >= float(refined_checkm[bin_name]['Contamination']):
            removed_bins[bin_name_o2]=bin_name_o
            os.chdir(pwd+'/'+target_folder)
            f.write(bin_name_o2+'\t'+str(refined_checkm[bin_name])+'\n')
            os.system('mv '+str(bin_name)+'.fa '+bin_name_o)
            f2.write(bin_name+' was selected. Rename '+bin_name+'.fa to '+bin_name_o+'\n')
            os.chdir(pwd)
        else:
            # f.write(bin_name_o2+'\t'+str(refined_checkm[bin_name])+'\n')
            os.chdir(pwd+'/'+target_folder)
            os.system('rm '+str(bin_name)+'.fa')
            f2.write('Removed '+bin_name+'.fa. Copy the original bin into the folder.'+'\n')
            os.chdir(pwd)
            os.system('cp '+pwd+'/'+binset_folder+'/'+bin_name_o+' '+pwd+'/'+target_folder+'/'+bin_name_o)
    print('Moved '+str(len(removed_bins))+' original bins into gf folder')

    os.chdir(pwd+'/'+binset_folder)
    n=0
    for bin_name in origin_checkm.keys():
        if bin_name not in removed_bins.keys():
            n+=1
            os.system('cp '+pwd+'/'+binset_folder+'/'+bin_name+'.fa '+pwd+'/'+target_folder+'/'+bin_name+'.fa')    
            os.system('cp '+pwd+'/'+binset_folder+'/'+bin_name+'.fasta '+pwd+'/'+target_folder+'/'+bin_name+'.fasta')
            f.write(bin_name+'\t'+str(origin_checkm[bin_name])+'\n')
    os.chdir(pwd)
    print('Moved '+str(n)+' original bins into gf folder')
    f.close()
    f2.close()

def finding_extract_contigs(extra_fa, asb, mod_bin_contigs, contig_len, type, num_threads):
    gf_bin_extract_contigs={}
    os.system('makeblastdb -dbtype nucl -in '+str(asb))
    os.system('blastn -query '+str(extra_fa)+' -db '+str(asb)+' -outfmt 6 -evalue 1e-20 -max_target_seqs 10 -num_threads '+str(num_threads)+' -out '+str(asb)+'_'+str(type)+'_lr_blasted.txt')
        
    for line in open(str(asb)+'_'+str(type)+'_lr_blasted.txt','r'):
        simi=float(str(line).strip().split('\t')[2])
        if simi >= 85:
            contig=str(line).strip().split('\t')[1].strip()
            bin_name_list=str(line).strip().split('\t')[0].strip().split('||')[1].split('|')
            ali_len=int(str(line).strip().split('\t')[3].strip())
            try:
                perc=ali_len/contig_len[contig]
            except:
                perc=0
            
            if perc >= 0.9:
                for item in bin_name_list:
                    try:
                        try:
                            contig2=mod_bin_contigs[item][contig]
                            gf_bin_extract_contigs[item][contig2]=0
                        except:
                            gf_bin_extract_contigs[item][contig]=0
                    except:
                        gf_bin_extract_contigs[item]={}
                        try:
                            contig2=mod_bin_contigs[item][contig]
                            gf_bin_extract_contigs[item][contig2]=0
                        except:
                            gf_bin_extract_contigs[item][contig]=0
    return gf_bin_extract_contigs

def blast2(bin_id, mag_folder, lr_folder, pwd):
    bin_lr=bin_id+'_gf_lr_polished.fa'
    bin_name=bin_id+'_mag_polished.fa'
    os.system('cp '+str(mag_folder)+'/'+str(bin_name)+' '+pwd)
    os.system('makeblastdb -dbtype nucl -hash_index -in '+str(bin_name))
    os.system('blastn -db '+str(bin_name)+' -query '+str(lr_folder)+'/'+str(bin_lr)+' -evalue 1e-20 -outfmt 6 -num_threads 1 -max_target_seqs 10 -out '+str(bin_id)+'_lr_polished_VS_'+str(bin_id)+'_mag_polished.txt')
    os.system('rm '+str(bin_name)+'.nhr '+str(bin_name)+'.nin '+str(bin_name)+'.nog '+str(bin_name)+'.nsd '+str(bin_name)+'.nsi '+str(bin_name)+'.nsq '+str(bin_name)+'.nhi '+str(bin_name)+'.nhd '+str(bin_name))

# def fq_seq_record(binset_folder, pwd):
#     seq_record={}
#     for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_gf_lr_sr_bins_seq'):
#         for file in files:        
#             if '_seq_R1.fq' in file:
#                 print('Reading '+str(file))
#                 n1=0
#                 for line in open(pwd+'/'+binset_folder+'_gf_lr_sr_bins_seq/'+file,'r'):
#                     n1+=1
#                     if n1 % 4 == 1:
#                         seq_record[str(line.strip())[3:-2]]=0
#                         # print(str(line.strip())[3:-2])
#     return seq_record

def lr_fil_record(file, binset_folder, pwd):
    bin_name=file.split('_elongation_lr.txt')[0]
    bin_id=bin_name.split('.')[0]
    ffq=open(bin_id+'_fil_lr.fq','w')
    lr_seq={}
    for line in open(file,'r'):
        lr_id=str(line).strip().split('\t')[0]
        lr_seq[lr_id]=0
    fq=str(bin_id)+'_lr.fq'

    n, record_line = 0, {}
    for line in open(fq,'r'):
        n+=1
        if n % 4 == 1:
            record_line[n]=0
            record_line[n+1]=0
            record_line[n+2]=0
            record_line[n+3]=0

    n, n1 = 0, 0
    for line in open(fq,'r'):
        n+=1
        try:
            record_line[n]+=1
            ffq.write(line)
        except:
            n1=1
    ffq.close()
    os.system('mv '+str(bin_id)+'_fil_lr.fq '+str(pwd)+'/'+str(binset_folder)+'_filtrated_long_read')

def polishing_main(binset_folder, datasets_list, assembly_list, long_read, batch, mapping_tool, major_run, ram, num_threads):
    pwd=os.getcwd()
    try:
        f=open('BASALT_log.txt','a')
    except:
        f=open('BASALT_log.txt','w')
    f.close()

    polish_status_file=str(batch)+'_polishing_status.txt'
    try:
        f=open(polish_status_file,'a')
    except:
        f=open(polish_status_file,'w')
    f.close()

    if len(long_read) != 0:
        contig_len={}
        for asb in assembly_list:
            for record in SeqIO.parse(asb,'fasta'):
                contig_len[record.id]=len(record.seq)

        # x=0
        # for line in open(polish_status_file,'r'):
        #     if 'Short-read mapping done!' in line:
        #         x=1

        # if x == 0 and batch == 1:
            # f_not_mapped_reads=open('Polishing_not_mapped_reads.txt','w')
            # f_not_mapped_reads.close()
            
            # if batch == 2:
            #     sr_folder='BestBinset_outlier_refined_filtrated_retrieved_polished_sr_bins_seq'
            #     os.chdir(pwd+'/'+sr_folder)
            #     for root, dirs, files in os.walk(pwd+'/'+sr_folder):
            #         for file in files:
            #             bin_id=str(file).split('_seq_R')[0]
            #             hz=str(file).split('_seq_R')[1]
            #             new_bin_id=bin_name_change[bin_id]
            #             new_name=new_bin_id+'_seq_R'+hz
            #             os.system('mv '+str(file)+' '+str(pwd)+'/'+str(new_name))
            #     os.chdir(pwd)
            # else:
            # os.system('mkdir '+binset_folder+'_polished_sr_bins_seq')
            # mapping_sr(total_fa, datasets_list, fq, pair, num_threads, 1, batch)
            # for bin_id in bin_seq.keys():
            #     try:
            #         os.system('mv '+str(bin_seq[str(bin_id)][0])+' '+str(bin_seq[str(bin_id)][1])+' '+pwd+'/'+binset_folder+'_polished_sr_bins_seq/')
            #     except:
            #         y=0

            # f=open(polish_status_file,'a')
            # f.write('Short-read mapping done!'+'\n')
            # f.close()
            # f=open('BASALT_log.txt','a')
            # f.write('Short-read mapping done!'+'\n')
            # f.close()

        if batch == 1:
            print('Recording seqs')
            A=mod_bin(binset_folder)
            mod_bin_folder=A[0]
            total_fa=A[1]
            mod_bin_list=A[2]
            original_bins_checkm=A[3]
            bin_name_change=A[4]
            fq, pair, bin_checkm, bin_seq={}, {}, {}, {}
            bin_checkm=original_bins_checkm.copy()
            for bin_id in mod_bin_list:
                fq[str(bin_id)]={}
                pair[str(bin_id)]={}
                if bin_id not in bin_seq.keys():
                    bin_seq[str(bin_id)]=[]
                    bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R1.fq')
                    bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R2.fq')

            if len(datasets_list) != 0:
                x=0
                for line in open(polish_status_file,'r'):
                    if 'Long-read mapping done!' in line:
                        x=1

                if x == 0:
                    n, bin_lr=0, {}
                    for lrs in long_read:
                        n+=1
                        os.system('minimap2 -t '+str(num_threads)+' -ax map-ont Total_bins.fa '+str(lrs)+' > lr'+str(n)+'.sam')
                        print('Splitting long reads '+str(n))
                        bin_lr.update(parse_lr_sam('lr'+str(n)+'.sam', lrs, n))
                        os.system('rm lr'+str(n)+'.sam')
                    
                    f=open('Bin_id_with_lr.txt','w')
                    for item in bin_lr.keys():
                        f.write(str(item)+'\n')
                    f.close()
                    f=open(polish_status_file,'a')
                    f.write('Long-read mapping done!'+'\n')
                    f.close()
                    f=open('BASALT_log.txt','a')
                    f.write('Long-read mapping done!'+'\n')
                    f.close()

                x=0
                for line in open(polish_status_file,'r'):
                    if 'Long-read blast done!' in line:
                        x=1

                if x == 0:
                    bin_lr = {}
                    for line in open('Bin_id_with_lr.txt','r'):
                        bin_lr[line.strip()]=0

                    pool=Pool(processes=num_threads)
                    for lr in bin_lr.keys():
                        fq=lr+'_lr.fq'
                        lr_fa=fq_2_fa(fq)
                        bin_name=lr_fa.split('_lr')[0]+'.fa'
                        pool.apply_async(blast, args=(binset_folder, bin_name, lr_fa, pwd))
                    pool.close()
                    pool.join()

                    f=open(polish_status_file,'a')
                    f.write('Long-read blast done!'+'\n')
                    f.close()
                    f=open('BASALT_log.txt','a')
                    f.write('Long-read blast done!'+'\n')
                    f.close()            

                x=0
                for line in open(polish_status_file,'r'):
                    if 'Long-read filtration done!' in line:
                        x=1

                if x == 0:
                    bin_lr = {}
                    for line in open('Bin_id_with_lr.txt','r'):
                        bin_lr[line.strip()]=0

                    pool=Pool(processes=num_threads)
                    for lr in bin_lr.keys():
                        fq=lr+'_lr.fq'
                        lr_fa=fq_2_fa(fq)
                        bin_name=lr_fa.split('_lr')[0]+'.fa'
                        pool.apply_async(filtration, args=(binset_folder, bin_name, lr_fa, num_threads, pwd))
                    pool.close()
                    pool.join()

                    # for lr in bin_lr.keys():
                    #     fq=lr+'_lr.fq'
                    #     lr_fa=fq_2_fa(fq)
                    #     bin_name=lr_fa.split('_lr')[0]+'.fa'
                    #     filtration(binset_folder, bin_name, lr_fa, num_threads, pwd)

                    f=open(polish_status_file,'a')
                    f.write('Long-read filtration done!'+'\n')
                    f.close()
                    f=open('BASALT_log.txt','a')
                    f.write('Long-read filtration done!'+'\n')
                    f.close()            

                x=0
                for line in open(polish_status_file,'r'):
                    if 'Long-read extraction done!' in line:
                        x=1

                if x == 0:
                    gap_filling_bin_lr, elongation_bin_lr, gf = {}, {}, {}
                    for root, dirs, files in os.walk(pwd):
                        for file in files:      
                            if '_gap_filling_lr.txt' in file:
                                bin_name=file.split('_gap_filling_lr.txt')[0]
                                bin_lr_fa=bin_name.split('.fa')[0]+'_lr.fa'
                                bin_lr_fa2=bin_name.split('.fa')[0]+'_gf_lr.fa'
                                gap_filling_bin_lr[bin_name]={}
                                for line in open(file,'r'):
                                    lr=str(line).strip().split('\t')[0]
                                    gf[bin_name]={}
                                    gf[bin_name][lr]=0
                                    try:
                                        gap_filling_bin_lr[lr][bin_name]=0
                                    except:
                                        gap_filling_bin_lr[lr]={}
                                        gap_filling_bin_lr[lr][bin_name]=0

                                if bin_name in gf.keys():
                                    if len(gf[bin_name]) != 0:
                                        f=open(bin_lr_fa2, 'w')
                                        for record in SeqIO.parse(bin_lr_fa,'fasta'):
                                            if record.id in gf[bin_name].keys():
                                                f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                                        f.close()

                    os.system('mkdir '+str(binset_folder)+'_filtrated_long_read')
                    # pool=Pool(processes=num_threads)
                    # for root, dirs, files in os.walk(pwd):
                    #     for file in files:     
                    #         if '_elongation_lr.txt' in file:
                    #             pool.apply_async(lr_fil_record, args=(file))
                    # pool.close()
                    # pool.join()

                    for root, dirs, files in os.walk(pwd):
                        for file in files:     
                            if '_elongation_lr.txt' in file:
                                print('Processing '+str(file))
                                lr_fil_record(file, binset_folder, pwd)

                    ### Closure of elongation contigs since this process would increase uncertain problem
                    # for root, dirs, files in os.walk(pwd):
                    #     for file in files:     
                    #         if '_elongation_lr.txt' in file:
                    #             bin_name=file.split('_elongation_lr.txt')[0]
                    #             bin_lr_fa=bin_name.split('.fa')[0]+'_lr.fa'
                    #             bin_lr_fa2=bin_name.split('.fa')[0]+'_el_lr.fa'
                    #             elongation_bin_lr[bin_name], el={}, {}
                    #             for line in open(file,'r'):
                    #                 lr=str(line).strip().split('\t')[0]
                    #                 try:
                    #                     if lr not in gf[bin_name].keys():
                    #                         el[lr]=0
                    #                 except:
                    #                     el[lr]=0

                    #                 try:
                    #                     elongation_bin_lr[lr][bin_name]=0
                    #                 except:
                    #                     elongation_bin_lr[lr]={}
                    #                     elongation_bin_lr[lr][bin_name]=0

                    #             if len(el) != 0:
                    #                 f=open(bin_lr_fa2, 'w')
                    #                 for record in SeqIO.parse(bin_lr_fa,'fasta'):
                    #                     if record.id in el.keys():
                    #                         f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    #                 f.close()
                    ### Closure of elongation contigs since this process would increase uncertain problem

                    f=open(polish_status_file,'a')
                    f.write('Long-read extraction done!'+'\n')
                    f.close()
                    f=open('BASALT_log.txt','a')
                    f.write('Long-read extraction done!'+'\n')
                    f.close() 

                x=0
                for line in open(polish_status_file,'r'):
                    if 'Finding extract gap filling contigs done!' in line:
                        x=1

                if x == 0:
                    n, seqs, seqs2 = 0, {}, {}
                    for root, dirs, files in os.walk(pwd):
                        for file in files:     
                            if '_gf_lr.fa' in file:
                                try:
                                    for record in SeqIO.parse(file,'fasta'):
                                        n+=1
                                        bin_name=str(file).split('_gf_lr.fa')[0]
                                        try:
                                            seqs[record.seq]+=1
                                            seqs2[record.seq]+='|'+bin_name
                                        except:
                                            seqs[record.seq]=0
                                            seqs2[record.seq]=str(record.id)+'||'+bin_name
                                except:
                                    xyznn=0
                    print('Total '+str(n)+' seqs. Non-redundant seqs: '+str(len(seqs2)))
                    
                    f=open('Total_gf_lr.fa','w')
                    for seq in seqs2.keys():
                        f.write('>'+str(seqs2[seq])+'\n'+str(seq)+'\n')
                    f.close()

                    mod_contigs, mod_bin_contigs = {}, {}
                    for line in open(pwd+'/'+binset_folder+'_mod/Mod_contig_id.txt','r'):
                        mod_c=str(line).strip().split('\t')[0]
                        bin_name=mod_c.split('_')[0]
                        org_c=str(line).strip().split('\t')[1]
                        mod_contigs[org_c]=mod_c
                        try:
                            mod_bin_contigs[bin_name][org_c]=mod_c
                        except:
                            mod_bin_contigs[bin_name]={}
                            mod_bin_contigs[bin_name][org_c]=mod_c

                    gf_bin_extract_contigs = {}
                    for asb in assembly_list: #### These may have problems.
                        gf_bin_extract_contigs.update(finding_extract_contigs('Total_gf_lr.fa', asb, mod_bin_contigs, contig_len, 'gf', num_threads))
                    
                    f=open('GF_lr_connecting_contigs.txt', 'w')
                    f2=open('GF_lr_extract_contigs.txt', 'w')
                    for bin_name in gf_bin_extract_contigs.keys():
                        for contig in gf_bin_extract_contigs[bin_name].keys():
                            f.write(str(bin_name)+'\t'+str(contig)+'\n')
                            if str(bin_name) not in str(contig) and bin_name != 'Total':
                                f2.write(str(bin_name)+'\t'+str(contig)+'\n')
                    f.close()
                    f2.close()

                    f=open(polish_status_file,'a')
                    f.write('Finding extract gap filling contigs done!'+'\n')
                    f.close()
                    f=open('BASALT_log.txt','a')
                    f.write('Finding extract gap filling contigs done!'+'\n')
                    f.close()
                
                x=0
                for line in open(polish_status_file,'r'):
                    if 'Evaluation gap filling extract contigs done!' in line:
                        x=1

                if x == 0:
                    bin_extra_contigs, extra_contigs = {}, {}
                    for line in open('GF_lr_extract_contigs.txt', 'r'):
                        bin_name=str(line).strip().split('\t')[0]
                        contig=str(line).strip().split('\t')[1]
                        try:
                            bin_extra_contigs[bin_name][contig]=0
                            extra_contigs[contig]=0
                        except:
                            bin_extra_contigs[bin_name]={}
                            bin_extra_contigs[bin_name][contig]=0
                            extra_contigs[contig]=0

                    contig_seq={}
                    for asb in assembly_list:
                        for record in SeqIO.parse(asb,'fasta'):
                            try:
                                extra_contigs[record.id]+=1
                                contig_seq[record.id]=record.seq
                            except:
                                xyz=1
                    
                    if len(contig_seq) >= 1:
                        os.system('mkdir '+binset_folder+'_gf_lr')
                        mod_bin_id={}
                        for line in open(pwd+'/'+binset_folder+'_mod/Bin_name_mod.txt','r'):
                            mod_bin_id[str(line).strip().split('\t')[1].split('.fa')[0]]=str(line).strip().split('\t')[0]
                        
                        for bin_name in bin_extra_contigs.keys():
                            bin_org_name=mod_bin_id[bin_name]
                            os.system('cp /'+pwd+'/'+binset_folder+'/'+bin_org_name+' /'+pwd+'/'+binset_folder+'_gf_lr/'+bin_org_name+'_gf_lr.fa')
                            os.chdir(pwd+'/'+binset_folder+'_gf_lr')
                            f=open(bin_org_name+'_gf_lr.fa','a')
                            for contig in bin_extra_contigs[bin_name].keys():
                                f.write('>'+str(contig)+'\n'+str(contig_seq[contig])+'\n')
                            f.close()
                            os.chdir(pwd)
                        
                        os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset_folder)+'_gf_lr '+str(binset_folder)+'_gf_lr_checkm')
                        refined_checkm=checkm(str(binset_folder)+'_gf_lr_checkm/storage/', pwd)
                        origin_checkm=checkm(str(binset_folder), pwd)
                        bin_comparison(binset_folder, str(binset_folder)+'_gf_lr', origin_checkm, refined_checkm, pwd, 'gf')
                        f=open('BASALT_log.txt','a')
                        f.write(binset_folder+'_gf_lr folder create'+'\n')
                        f.close() 
                    else:
                        f=open('BASALT_log.txt','a')
                        f.write(binset_folder+'_gf_lr did not create since there is not qualified long read'+'\n')
                        f.close() 

                    os.system('rm Total_gf_lr.fa')
                    os.system('mkdir Gap_filling_long_read')
                    os.system('mv *_gf_lr.fa *_gap_filling_lr.txt Gap_filling_long_read')
                    os.system('mkdir BestBinset_long_read')
                    os.system('mv *_lr.fq BestBinset_long_read')
                    f=open(polish_status_file,'a')
                    f.write('Evaluation gap filling extract contigs done!'+'\n')
                    f.close()
                    f=open('BASALT_log.txt','a')
                    f.write('Evaluation gap filling extract contigs done!'+'\n')
                    f.close() 

        ################### Close the finding extract elongation contigs module ###################
        # x=0
        # for line in open(polish_status_file,'r'):
        #     if 'Finding extract elongation contigs done!' in line:
        #         x=1

        # if x == 0:
        #     f=open('Total_el_lr.fa','w')
        #     n, seqs, seqs2 = 0, {}, {}
        #     for root, dirs, files in os.walk(pwd):
        #         for file in files:     
        #             if '_el_lr.fa' in file:
        #                 for record in SeqIO.parse(file,'fasta'):
        #                     n+=1
        #                     bin_name=str(file).split('_el_lr.fa')[0]
        #                     try:
        #                         seqs[record.seq]+=1
        #                         seqs2[record.seq]+='|'+bin_name
        #                     except:
        #                         seqs[record.seq]=0
        #                         seqs2[record.seq]=str(record.id)+'||'+bin_name
        #     print('Total '+str(n)+' seqs. Non-redundant seqs: '+str(len(seqs2)))
        #     for seq in seqs2.keys():
        #         f.write('>'+str(seqs2[seq])+'\n'+str(seq)+'\n')
        #     f.close()

        #     mod_contigs, mod_bin_contigs = {}, {}
        #     for line in open(pwd+'/'+binset_folder+'_mod/Mod_contig_id.txt','r'):
        #         mod_c=str(line).strip().split('\t')[0]
        #         bin_name=mod_c.split('_')[0]
        #         org_c=str(line).strip().split('\t')[1]
        #         mod_contigs[org_c]=mod_c
        #         try:
        #             mod_bin_contigs[bin_name][org_c]=mod_c
        #         except:
        #             mod_bin_contigs[bin_name]={}
        #             mod_bin_contigs[bin_name][org_c]=mod_c

        #     el_bin_extract_contigs = {}
        #     for asb in assembly_list: #### These may have problems.
        #         el_bin_extract_contigs.update(finding_extract_contigs('Total_el_lr.fa', asb, mod_bin_contigs, contig_len, 'el', num_threads))
            
        #     f=open('EL_lr_connecting_contigs.txt', 'w')
        #     f2=open('EL_lr_extract_contigs.txt', 'w')
        #     for bin_name in el_bin_extract_contigs.keys():
        #         for contig in el_bin_extract_contigs[bin_name].keys():
        #             f.write(str(bin_name)+'\t'+str(contig)+'\n')
        #             if str(bin_name) not in str(contig):
        #                 f2.write(str(bin_name)+'\t'+str(contig)+'\n')
        #     f.close()
        #     f2.close()

        #     f=open(polish_status_file,'a')
        #     f.write('Finding extract elongation contigs done!'+'\n')
        #     f.close()
        #     f=open('BASALT_log.txt','a')
        #     f.write('Finding extract elongation contigs done!'+'\n')
        #     f.close()

        # x=0
        # for line in open(polish_status_file,'r'):
        #     if 'Evaluation elongation extract contigs done!' in line:
        #         x=1

        # if x == 0:
        #     bin_extra_contigs, extra_contigs = {}, {}
        #     for line in open('EL_lr_extract_contigs.txt', 'r'):
        #         bin_name=str(line).strip().split('\t')[0]
        #         contig=str(line).strip().split('\t')[1]
        #         try:
        #             bin_extra_contigs[bin_name][contig]=0
        #             extra_contigs[contig]=0
        #         except:
        #             bin_extra_contigs[bin_name]={}
        #             bin_extra_contigs[bin_name][contig]=0
        #             extra_contigs[contig]=0

        #     contig_seq={}
        #     for asb in assembly_list:
        #         for record in SeqIO.parse(asb,'fasta'):
        #             try:
        #                 extra_contigs[record.id]+=1
        #                 contig_seq[record.id]=record.seq
        #             except:
        #                 xyz=1
            
        #     if len(contig_seq) >= 1:
        #         os.system('mkdir '+binset_folder+'_el_lr')
        #         mod_bin_id={}
        #         for line in open(pwd+'/'+binset_folder+'_mod/Bin_name_mod.txt','r'):
        #             mod_bin_id[str(line).strip().split('\t')[1].split('.fa')[0]]=str(line).strip().split('\t')[0]
                
        #         for bin_name in bin_extra_contigs.keys():
        #             bin_org_name=mod_bin_id[bin_name]
        #             os.system('cp /'+pwd+'/'+binset_folder+'_gf_lr/'+bin_org_name+' /'+pwd+'/'+binset_folder+'_el_lr/'+bin_org_name+'_el_lr.fa')
        #             os.chdir(pwd+'/'+binset_folder+'_el_lr')
        #             f=open(bin_org_name+'_el_lr.fa','a')
        #             for contig in bin_extra_contigs[bin_name].keys():
        #                 f.write('>'+str(contig)+'\n'+str(contig_seq[contig])+'\n')
        #             f.close()
        #             os.chdir(pwd)
                
        #         os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset_folder)+'_el_lr '+str(binset_folder)+'_el_lr_checkm')
        #         refined_checkm=checkm(str(binset_folder)+'_el_lr_checkm/storage/', pwd)
        #         origin_checkm=checkm(str(binset_folder)+'_gf_lr', pwd)
        #         bin_comparison(binset_folder+'_gf_lr', binset_folder+'_el_lr', origin_checkm, refined_checkm, pwd, 'el')

        #     f=open(polish_status_file,'a')
        #     f.write('Evaluation elongation extract contigs done!'+'\n')
        #     f.close()
        #     f=open('BASALT_log.txt','a')
        #     f.write('Evaluation elongation extract contigs done!'+'\n')
        #     f.close()
        ################### Close the finding extract elongation contigs module ###################

    print('Batch '+str(batch)+' MAGs polishing start')
    f=open('BASALT_log.txt','a')
    f.write('Batch '+str(batch)+' MAGs polishing start'+'\n')
    f.close()

    if len(datasets_list) != 0:
        y=0
        for line in open(polish_status_file,'r'):
            if 'MAGs sr polishing done!' in line:
                y=1

        if y == 0:
            rx, mp_run = 1, 0
            while rx == 1:
                mp_run+=1
                if mp_run <= major_run: ### Primary run
                    x=0
                    for line in open(polish_status_file,'r'):
                        if 'Run'+str(mp_run)+' short-reads mapping done!' in line:
                            x=1

                    if x == 0:
                        print('Mapping and Polishing start at run'+str(mp_run))
                        y=0
                        for line in open(polish_status_file,'r'):
                            if 'Reads extraction at run '+str(mp_run)+' done!' in line:
                                y=1

                        if mp_run >= 2:
                            if y == 0:
                                prev_unmapped_reasds, map_await_ds ={}, {}
                                mp_run2=mp_run-1
                                for line in open('Not_mapped_reads_'+str(batch)+'_'+str(mp_run2)+'.txt','r'):
                                    prev_unmapped_reasds[str(line).strip()]=0
                                
                                # try:
                                #     for line in open('Not_mapped_reads'+str(mp_run)+'.txt','r'):
                                #         prev_unmapped_reasds2[str(line).strip()]=0
                                # except:
                                #     xyz=0
                                
                                # if len(prev_unmapped_reasds) == len(prev_unmapped_reasds2):
                                #     rx = 0

                                print('Re-writing previously un-mapped fq seqs')
                                try:
                                    fe1=open('Remained_seq1.fq','a')
                                    fe1.close()
                                    os.system('mv Remained_seq1.fq Remained_pre_seq1.fq')
                                    fe1=open('Remained_seq2.fq','a')
                                    fe1.close()
                                    os.system('mv Remained_seq2.fq Remained_pre_seq2.fq')
                                except:
                                    xx=0
                                fe1=open('Remained_seq1.fq','w')
                                fe2=open('Remained_seq2.fq','w')
                                fe1.close()
                                fe2.close()
                                unmapped_total = {}
                                for ds in datasets_list.keys():
                                    print('Scanning dataset: '+str(datasets_list[ds]))
                                    r1=datasets_list[ds][0]
                                    r2=datasets_list[ds][1]
                                    unmapped = {}
                                    # fxyz=open('IDs_'+str(ds)+'.txt','w')
                                    n1, record_num, record_num2 = 0, {}, {}
                                    for line in open(r1,'r'):
                                        n1 += 1
                                        if n1 % 4 == 1:
                                            ids=str(ds)+'_'+str(line.strip())[1:-2]
                                            # fxyz.write(str(ids)+'\n')
                                            
                                            # if ids in prev_unmapped_reasds.keys():
                                            try:
                                                # print(str(ids))
                                                prev_unmapped_reasds[str(ids)]+=1
                                                unmapped[str(ids)]=0
                                                unmapped_total[str(ids)]=0
                                                record_num[n1]=0
                                                record_num2[n1+1]=0
                                                record_num2[n1+2]=0
                                                record_num2[n1+3]=0
                                            except:
                                            # else:
                                                xyzzz=0

                                    # fxyz.close()
                                    n2=n1/4
                                    n3=n2-len(unmapped)
                                    print('Total '+str(n2)+' in dataset. '+str(n3)+' seqs have been mapped.')
                                    print('Total '+str(len(unmapped))+' unmapped reads in dataset. Writing the unmapped reads')
                                    if len(unmapped) != 0:
                                        fe1=open('Remained_seq1.fq','a')
                                        fe2=open('Remained_seq2.fq','a') 
                                        n1=0
                                        for line in open(r1,'r'):
                                            n1+=1
                                            try:
                                                record_num[n1]+=1
                                                fe1.write(str(str(line.strip())[0:-2])+' '+str(str(line.strip())[-1])+'\n')
                                            except:
                                                xyz=0

                                            try:
                                                record_num2[n1]+=1
                                                fe1.write(line)
                                            except:
                                                xyz=0


                                        n1=0
                                        for line in open(r2,'r'):
                                            n1+=1
                                            try:
                                                record_num[n1]+=1
                                                fe2.write(str(str(line.strip())[0:-2])+' '+str(str(line.strip())[-1])+'\n')
                                            except:
                                                xyz=0
                                                
                                            try:
                                                record_num2[n1]+=1
                                                fe2.write(line)
                                            except:
                                                xyz=0

                                        fe1.close()
                                        fe2.close()
                                if len(unmapped_total) != 0:
                                    map_await_ds={1:['Remained_seq1.fq','Remained_seq2.fq']}
                                    fzyz=open(polish_status_file,'a')
                                    fzyz.write('Reads extraction at run '+str(mp_run)+' done!'+'\n')
                                    fzyz.close()

                                else:
                                    os.system('rm Remained_seq1.fq Remained_seq2.fq')
                                    try:
                                        os.system('mv Remained_pre_seq1.fq Remained_seq1.fq ')
                                        os.system('mv Remained_pre_seq2.fq Remained_seq2.fq ')
                                    except:
                                        xx=0
                                    fx=open(polish_status_file,'a')
                                    fx.write('MAGs sr polishing done!'+'\n')
                                    fx.close()
                                    fx=open('BASALT_log.txt','a')
                                    fx.write('MAGs sr polishing done!'+'\n')
                                    fx.close()
                            else:
                                try:
                                    nx=0
                                    for line in open('Remained_seq1.fq','r'):
                                        nx+=1
                                        if nx >= 1:
                                            break
                                    map_await_ds={1:['Remained_seq1.fq','Remained_seq2.fq']}
                                except:
                                    map_await_ds={}

                                prev_unmapped_reasds = {}
                                mp_run2=mp_run-1
                                for line in open('Not_mapped_reads_'+str(batch)+'_'+str(mp_run2)+'.txt','r'):
                                    prev_unmapped_reasds[str(line).strip()]=0

                        # print(str(map_await_ds))
                        bin_seq = {}
                        if mp_run == 1:
                            print('Recording seqs')
                            if batch == 1:
                                try:
                                    A=mod_bin(binset_folder+'_gf_lr')
                                except:
                                    A=mod_bin(binset_folder) ### did not have binset_folder+'_gf_lr'
                            else: ### batch == 2 or batch == 3
                                A=mod_bin(binset_folder) ### 'BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC'

                            mod_bin_folder=A[0]
                            total_fa=A[1]
                            mod_bin_list=A[2]
                            original_bins_checkm=A[3]
                            bin_name_change=A[4]

                            fq, pair, bin_checkm, bin_seq={}, {}, {}, {}
                            bin_checkm=original_bins_checkm.copy()
                            for bin_id in mod_bin_list:
                                fq[str(bin_id)]={}
                                pair[str(bin_id)]={}
                                if bin_id not in bin_seq.keys():
                                    bin_seq[str(bin_id)]=[]
                                    bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R1.fq')
                                    bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R2.fq')
                                    try:
                                        f1=open(str(bin_id)+'_seq_R1.fq','a')
                                    except:
                                        f1=open(str(bin_id)+'_seq_R1.fq','w')

                                    try:
                                        f2=open(str(bin_id)+'_seq_R2.fq','a')
                                    except:
                                        f2=open(str(bin_id)+'_seq_R2.fq','w')
                                    f1.close()
                                    f2.close()

                            try:
                                sr_folder='BestBinset_sr_bins_seq'
                                os.chdir(pwd+'/'+sr_folder)
                                for root, dirs, files in os.walk(pwd+'/'+sr_folder):
                                    for file in files:
                                        bin_id=str(file).split('_seq_R')[0]
                                        hz=str(file).split('_seq_R')[1]
                                        new_bin_id=bin_name_change[bin_id]
                                        new_name=new_bin_id+'_seq_R'+hz
                                        os.system('mv '+str(file)+' '+str(pwd)+'/'+str(new_name))
                                os.chdir(pwd)

                            except:
                                os.system('mkdir BestBinset_sr_bins_seq')

                            mapping_sr(total_fa, datasets_list, fq, pair, mapping_tool, num_threads, batch, mp_run)
                        else: ### mp_run >= 2
                            if len(map_await_ds) != 0:
                                fq, pair ={}, {}
                                if batch == 1:
                                    fd=binset_folder+'_gf_lr_mod'
                                else: ### batch == 2 and batch == 3
                                    fd=binset_folder+'_mod' ### 'BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC_mod'

                                for line in open(pwd+'/'+fd+'/Bin_name_mod.txt','r'):
                                    bin_id=str(line).strip().split('\t')[1].split('.fa')[0]
                                    fq[str(bin_id)]={}
                                    pair[str(bin_id)]={}
                                    bin_seq[str(bin_id)]=[]
                                    bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R1.fq')
                                    bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R2.fq')
                                    f1=open(str(bin_id)+'_seq_R1.fq','w')
                                    f2=open(str(bin_id)+'_seq_R2.fq','w')
                                    f1.close()
                                    f2.close()

                                f1=open('Total_bins.fa','w')
                                os.chdir(binset_folder+'_MAGs_polished')
                                for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_MAGs_polished'):
                                    for file in files:
                                        if '_mag_polished.fa' in file:
                                            bin_id=file.split('_mag_polished.fa')[0]
                                            for record in SeqIO.parse(file, 'fasta'):
                                                # record_seq[bin_id][record.id]=str(record.seq)
                                                f1.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                                os.chdir(pwd)
                                f1.close()

                                mapping_sr('Total_bins.fa', map_await_ds, fq, pair, mapping_tool, num_threads, batch, mp_run)
                                for item in map_await_ds.keys():
                                    os.system('mv '+str(map_await_ds[item][0])+' '+str(map_await_ds[item][0])+'_'+str(mp_run))
                                    os.system('mv '+str(map_await_ds[item][1])+' '+str(map_await_ds[item][1])+'_'+str(mp_run))

                        bin_waiting_polish={}
                        for bin_id in bin_seq.keys():
                            r1=str(bin_seq[str(bin_id)][0])
                            r2=str(bin_seq[str(bin_id)][1])

                            if mp_run == 1 and batch != 2:
                                try:
                                    # os.system('mv '+str(r1)+' '+str(r2)+' '+pwd+'/'+binset_folder+'_gf_lr_sr_bins_seq')
                                    os.system('mv '+str(r1)+' '+str(r2)+' '+pwd+'/BestBinset_sr_bins_seq')
                                except:
                                    y=0

                            # elif len(map_await_ds) != 0: ### Merged old R1 reads seq with new mapped reads
                            else:
                                record_seq = []
                                try:
                                    for line in open(r1,'r'):
                                        record_seq.append(line)
                                except:
                                    xyz=0

                                if len(record_seq) < 4:
                                    os.system('rm '+r1+' '+r2)
                                else:
                                    bin_waiting_polish[bin_id]=1
                                    os.chdir(pwd+'/BestBinset_sr_bins_seq')

                                    f=open(r1, 'a')
                                    for ids in record_seq:
                                        f.write(ids)
                                    f.close()

                                    os.chdir(pwd)
                                    record_seq = []
                                    for line in open(r2,'r'):
                                        record_seq.append(line)

                                    os.chdir(pwd+'/BestBinset_sr_bins_seq')
                                        
                                    f=open(r2, 'a')
                                    for ids in record_seq:
                                        f.write(ids)
                                    f.close()
                                    os.chdir(pwd)
                                    os.system('rm '+r1+' '+r2)
                        
                        if len(bin_waiting_polish) != 0:
                            f=open('Bin_waiting_for_polish_batch.txt','w')
                            for bin_id in bin_waiting_polish.keys():
                                f.write(bin_id+'\n')
                            f.close()
                        else:
                            if mp_run != 1:
                                rx = 0

                        f=open(polish_status_file,'a')
                        f.write('Run'+str(mp_run)+' short-reads mapping done!'+'\n')
                        f.close()
                        f=open('BASALT_log.txt','a')
                        f.write('Run'+str(mp_run)+' short-reads mapping done!'+'\n')
                        f.close()

                        prev_unmapped_reasds2 = {}
                        try:
                            for line in open('Not_mapped_reads_'+str(batch)+'_'+str(mp_run)+'.txt','r'):
                                prev_unmapped_reasds2[str(line).strip()]=0
                        except:
                            print('There is not unmapped read in Not_mapped_reads'+str(mp_run)+'.txt')
                        try:
                            if len(prev_unmapped_reasds) == len(prev_unmapped_reasds2) or len(prev_unmapped_reasds2) == 0:
                                rx=0
                        except:
                            rx=1

                    ### Polishing bin before aligning long-read against MAG
                    if rx == 1:
                        x=0
                        for line in open(polish_status_file,'r'):
                            if 'Run'+str(mp_run)+' MAG sr polishing done!' in line:
                                x=1

                        if x == 0:
                            print('Processing run'+str(mp_run)+' polishing')
                            bin_waiting_polish={}
                            try:
                                for line in open('Bin_waiting_for_polish_batch.txt','r'):
                                    bin_waiting_polish[str(line).strip()]=0
                            except:
                                xyz=0

                            accomplish_bins={}
                            try:
                                f=open('Accompleted_bins_'+str(batch)+'_'+str(mp_run)+'.txt','a')
                                f.close()

                                for line in open('Accompleted_bins_'+str(batch)+'_'+str(mp_run)+'.txt','r'):
                                    accomplish_bins[str(line).strip()]=0
                                f=open('BASALT_log.txt','a')
                                f.write('Run'+str(mp_run)+' short-reads mapping '+str(len(accomplish_bins))+' bins accomplished'+'\n')
                                f.close()

                            except:
                                f=open('Accompleted_bins_'+str(batch)+'_'+str(mp_run)+'.txt','w')
                                f.close()

                            # sr_folder=binset_folder+'_gf_lr_sr_bins_seq/'
                            sr_folder='BestBinset_sr_bins_seq'
                            if mp_run == 1:
                                # input_bin_folder=binset_folder+'_mod' ### Be careful of input folder
                                if batch == 1:
                                    input_bin_folder=binset_folder+'_gf_lr_mod' ### Be careful of input folder
                                elif batch == 2:
                                    input_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC_mod'
                                elif batch == 3:
                                    input_bin_folder=binset_folder+'_mod'

                                os.system('mkdir '+binset_folder+'_MAGs_polished')
                            else:
                                mp_run2=mp_run-1
                                input_bin_folder=binset_folder+'_MP_'+str(mp_run2)
                                if len(accomplish_bins) == 0:
                                    os.system('mv '+binset_folder+'_MAGs_polished '+binset_folder+'_MP_'+str(mp_run2))
                                    os.system('mkdir '+binset_folder+'_MAGs_polished')
                                    os.chdir(pwd+'/'+binset_folder+'_MP_'+str(mp_run2))
                                    for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_MP_'+str(mp_run2)):
                                        for file in files:
                                            bin_id=file.split('_mag_polished.fa')[0]
                                            if bin_id in bin_waiting_polish.keys():
                                                os.system('mv '+file+' '+bin_id+'.fa')
                                            else:
                                                os.system('mv '+file+' '+pwd+'/'+binset_folder+'_MAGs_polished')
                                    os.chdir(pwd)
                            # for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_MP_'+str(mp_run2)):
                            #     for file in files:
                            #         if '_mag_polished.fa' in file:
                            #             os.system('mv '+pwd+'/'+binset_folder+'_'+str(mp_run2)+'/'+str(file)+' '+pwd+'/'+binset_folder+'_MAGs_polished/'+str(file))
                            
                            bin_seq={}
                            for root, dirs, files in os.walk(pwd+'/'+sr_folder):
                                for file in files:
                                    if '_seq_R1.fq' in file:
                                        bin_name=str(file).split('_seq_R1.fq')[0]
                                        bin_seq[bin_name]=[]
                                        bin_seq[bin_name].append(file)
                                        bin_seq[bin_name].append(bin_name+'_seq_R2.fq')

                            mod_bin_list=[]
                            if batch == 1:
                                input_bin_folder=binset_folder+'_gf_lr_mod' ### Be careful of input folder
                                input_bin_folder2=binset_folder+'_mod'
                            elif batch == 2:
                                input_bin_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished_re-assembly_OLC_mod'
                            elif batch == 3:
                                input_bin_folder=binset_folder+'_mod'

                            # for line in open(pwd+'/'+binset_folder+'_mod/Bin_name_mod.txt','r'):
                            try:
                                for line in open(pwd+'/'+input_bin_folder+'/Bin_name_mod.txt','r'): ### Needed to be checked
                                    bin_id=str(line).strip().split('\t')[1].strip().split('.fa')[0]
                                    if mp_run == 1:
                                        mod_bin_list.append(bin_id)
                                    else:
                                        if bin_id in bin_waiting_polish.keys():
                                            mod_bin_list.append(bin_id)
                            except:
                                for line in open(pwd+'/'+input_bin_folder2+'/Bin_name_mod.txt','r'):
                                    bin_id=str(line).strip().split('\t')[1].strip().split('.fa')[0]
                                    if mp_run == 1:
                                        mod_bin_list.append(bin_id)
                                    else:
                                        if bin_id in bin_waiting_polish.keys():
                                            mod_bin_list.append(bin_id)

                            fx=open('BASALT_log.txt','a')
                            fx.write(str(len(mod_bin_list))+' bin(s) in modified bin list: '+str(mod_bin_list)+'\n')
                            fx.close()

                            fold=int(num_threads/6)-1
                            folder_index=[1,1,1]
                            for i in range(1,fold):
                                folder_index.append(1.2*i)

                            i = 0
                            while i != -1:
                                i+=1
                                tpp=int(5*folder_index[i-1])
                                remain_bin, complete_bin = {}, {}
                                try:
                                    for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_MAGs_polished'):
                                        for file in files:
                                            if '_mag_polished.fa' in file:
                                                bin_id=file.split('_mag_polished.fa')[0]
                                                complete_bin[bin_id]=0
                                except:
                                    xyz=0

                                # print(str(len(complete_bin))+' bin(s). Completed bins: '+str(complete_bin))
                                fx=open('BASALT_log.txt','a')
                                fx.write(str(len(complete_bin))+' bin(s) polishing accomplish: '+str(complete_bin)+'\n')
                                fx.close()

                                for bin_id in mod_bin_list:
                                    if bin_id not in complete_bin.keys():
                                        remain_bin[bin_id]=0

                                print('There is '+str(len(remain_bin))+' remained in the run '+str(i))
                                fx=open('BASALT_log.txt','a')
                                fx.write('Mapping and polishing run'+str(mp_run)+'. There is '+str(len(remain_bin))+' remained in the run '+str(i)+'\n')
                                fx.close()

                                if len(remain_bin) == 0:
                                    i=-1
                                else:
                                    os.system('rm *.sa *.ann *.amb *.pac *.bwt *.bam *.bai')
                                    for root, dirs, files in os.walk(pwd):
                                        for file in files:
                                            if '_mag_polished' in file and '.fasta' in file:
                                                print('Recording '+str(file))
                                                bin_id=str(file).split('_mag_polished')[0]
                                                try:
                                                    bin_id_r=int(str(file).split('_mag_polished')[1].split('.fasta')[0])
                                                    try:
                                                        num=remain_bin[bin_id]
                                                        if bin_id_r > num:
                                                            remain_bin[bin_id]=bin_id_r
                                                    except:
                                                        xyz=0
                                                except:
                                                    xyz=0

                                # print(str(remain_bin))            
                                for bin_id in remain_bin.keys():
                                    num=remain_bin[bin_id]
                                    if num != 0:
                                        remain_bin[bin_id]=str(bin_id)+'_mag_polished'+str(remain_bin[bin_id])+'.fasta'
                                    else:
                                        remain_bin[bin_id]=str(bin_id)+'.fa'

                                # print(str(remain_bin))

                                num_projects_p=int(num_threads/2)
                                if tpp <= num_projects_p:
                                    num_projects=int(num_threads/tpp)
                                    pool=Pool(processes=num_projects)
                                    for bin_id in remain_bin.keys():
                                    # for bin_name in bin_lr:
                                        print('Processing '+str(bin_id))
                                        initial_bin_name=remain_bin[bin_id]
                                        try:
                                            fp=open(str(bin_id)+'_mag_polished.fa','a')
                                            fp.close()
                                        except:
                                            fp=open(str(bin_id)+'_mag_polished.fa','w')
                                            fp.close()
                                        pool.apply_async(sr_polishing, args=(bin_id, bin_seq, input_bin_folder, sr_folder, tpp, pwd, binset_folder, 'MAG_polishing', initial_bin_name, mapping_tool, mp_run, major_run, batch))
                                        # sr_polishing(bin_id, bin_seq, input_bin_folder, sr_folder, tpp, pwd, binset_folder, 'MAG_polishing', initial_bin_name, mp_run, batch)
                                    pool.close()
                                    pool.join()
                                    
                                    remain_bin = {}
                                    for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_MAGs_polished'):
                                        for file in files:
                                            if '_mag_polished.fa' in file:
                                                bin_id=file.split('_mag_polished.fa')[0]
                                                # os.system('mv '+file+' '+binset_folder+'_MAGs_polished')
                                                complete_bin[bin_id]=0
                                    
                                    for bin_id in mod_bin_list:
                                        if bin_id not in complete_bin.keys():
                                            remain_bin[bin_id]=0
                                    
                                    if len(remain_bin) == 0:
                                        i=-1

                                else:
                                    for bin_id in remain_bin.keys():
                                        print('Processing '+str(bin_id))
                                        initial_bin_name=remain_bin[bin_id]
                                        try:
                                            fp=open(str(bin_id)+'_mag_polished.fa','a')
                                            fp.close()
                                        except:
                                            fp=open(str(bin_id)+'_mag_polished.fa','w')
                                            fp.close()
                                        sr_polishing(bin_id, bin_seq, input_bin_folder, sr_folder, num_threads, pwd, binset_folder, 'MAG_polishing', initial_bin_name, mapping_tool, mp_run, major_run, batch)
                                    i=-1

                            os.system('rm hs_err_pid*')
                            # if mp_run != 1:
                            #     mp_run2=mp_run-1
                            #     os.system('mv '+binset_folder+'_temp '+binset_folder+'_'+str(mp_run2))
                            f=open(polish_status_file,'a')
                            f.write('Run'+str(mp_run)+' MAG sr polishing done!'+'\n')
                            f.close()
                            f=open('BASALT_log.txt','a')
                            f.write('Run'+str(mp_run)+' MAG sr polishing done!'+'\n')
                            f.close()
                            for i in range(1, 51):
                                os.system('rm *_mag_polished'+str(i)+'.fasta')
                            os.system('rm *.amb *.ann *.sa *.bwt *.pac *.bam *.bai *_mag_polished.fa')
                else:
                    rx = 0
            
            ####### ####### Below is process of polishing long-read. But I close the function since it is unecessary. The code is functional. It could be used in later ####### ####### 
            # x=0
            # for line in open(polish_status_file,'r'):
            #     if 'Long-reads polished done!' in line:
            #         x=1

            # if x == 0 and batch == 1:
            #     print('Long-reads polishing')
            #     bin_lr, bin_seq = {}, {}
            #     os.chdir(pwd+'/Gap_filling_long_read')
            #     for root, dirs, files in os.walk(pwd+'/Gap_filling_long_read'):
            #         for file in files:
            #             if '_gap_filling_lr.txt' in file:
            #                 n=0
            #                 for line in open(file,'r'):
            #                     n+=1
            #                 if n >= 1:
            #                     bin_name=str(file).split('_gap_filling_lr.txt')[0]
            #                     bin_lr[bin_name]={}
            #                     for line in open(file, 'r'):
            #                         contig=str(line).strip().split('\t')[0]
            #                         bin_lr[bin_name][contig]=0
            #                 else:
            #                     os.system('rm '+file)
            #     os.chdir(pwd)

            #     for bin_name in bin_lr.keys():
            #         bin_lr_name=bin_name.split('.fa')[0]+'_lr.fa'
            #         print('Recording seqs from '+str(bin_lr_name))
            #         try:
            #             x=0
            #             for line in open(bin_lr_name,'r'):
            #                 x+=1
            #                 if x == 1:
            #                     break
                        
            #             bin_gf_lr_total=bin_name.split('.fa')[0]+'_gf_lr_total.fa'
            #             f=open(bin_gf_lr_total,'w')
            #             for record in SeqIO.parse(bin_lr_name,'fasta'):
            #                 if record.id in bin_lr[bin_name].keys():
            #                     f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            #             f.close()
            #             os.system('mv '+str(bin_gf_lr_total)+' Gap_filling_long_read')
            #         except:
            #             print('There is not '+str(bin_lr_name)+' file in the folder')

            #     # sr_folder=binset_folder+'_gf_lr_sr_bins_seq/'
            #     sr_folder='BestBinset_outlier_refined_filtrated_retrieved_polished_sr_bins_seq'
            #     os.chdir(pwd+'/'+sr_folder)
            #     for root, dirs, files in os.walk(pwd+'/'+sr_folder):
            #         for file in files:
            #             if '_seq_R1.fq' in file:
            #                 bin_name=str(file).split('_seq_R1.fq')[0]
            #                 bin_seq[bin_name]=[]
            #                 bin_seq[bin_name].append(file)
            #                 bin_seq[bin_name].append(bin_name+'_seq_R2.fq')
            #     os.chdir(pwd)

            #     lr_folder='Gap_filling_long_read/'
            #     os.system('mkdir '+binset_folder+'_gf_lr_polished')

            #     fold=int(num_threads/6)
            #     folder_index=[1,1]
            #     for i in range(1,fold):
            #        folder_index.append(1.2*i)

            #     i = 0
            #     while i != -1:
            #         i+=1
            #         tpp=int(5*folder_index[i-1])
            #         remain_bin, complete_bin = {}, {}
            #         try:
            #             for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_gf_lr_polished'):
            #                 for file in files:
            #                     if '_gf_lr_polished.fa' in file:
            #                         bin_id=file.split('_gf_lr_polished.fa')[0]+'.fa'
            #                         complete_bin[bin_id]=0
            #         except:
            #             xyz=0

            #         for bin_name in bin_lr.keys():
            #             if bin_name not in complete_bin.keys():
            #                 remain_bin[bin_name]=0
            #         print('There is '+str(len(remain_bin))+' remained in the run '+str(i))
            #         fx=open('BASALT_log.txt','a')
            #         fx.write('There is '+str(len(remain_bin))+' remained in the run '+str(i)+'\n')
            #         fx.close()

            #         if len(remain_bin) == 0:
            #             i=-1
            #         else:
            #             for root, dirs, files in os.walk(pwd):
            #                 for file in files:
            #                     if '_lr_polished' in file and '.fasta' in file:
            #                         bin_name=str(file).split('_lr_polished')[0]+'.fa'
            #                         bin_id_r=int(str(file).split('_lr_polished')[1].split('.fasta')[0])
            #                         try:
            #                             num=remain_bin[bin_name]
            #                             if bin_id_r > num:
            #                                 remain_bin[bin_name]=bin_id_r
            #                         except:
            #                             xyz=0
                                    
            #         for bin_name in remain_bin.keys():
            #             bin_id=bin_name.split('.fa')[0]
            #             num=remain_bin[bin_name]
            #             if num != 0:
            #                 remain_bin[bin_name]=str(bin_id)+'_lr_polished'+str(remain_bin[bin_name])+'.fasta'
            #             else:
            #                 remain_bin[bin_name]=str(bin_id)+'_gf_lr_total.fa'

            #         num_projects_p=int(num_threads/2)
            #         if tpp <= num_projects_p:
            #             num_projects=int(num_threads/tpp)
            #             pool=Pool(processes=num_projects)
            #             for bin_name in remain_bin.keys():
            #             # for bin_name in bin_lr:
            #                 bin_id=bin_name.split('.fa')[0]
            #                 print('Processing '+str(bin_id))
            #                 initial_bin_name=remain_bin[bin_name]
            #                 try:
            #                     fp=open(str(bin_id)+'_gf_lr_polished.fa','a')
            #                     fp.close()
            #                 except:
            #                     fp=open(str(bin_id)+'_gf_lr_polished.fa','w')
            #                     fp.close()
            #                 pool.apply_async(sr_polishing, args=(bin_id, bin_seq, lr_folder, sr_folder, tpp, pwd, binset_folder, 'lr_polishing', initial_bin_name, mapping_tool, 1, major_run, batch))
            #                 # sr_polishing(bin_id, bin_seq, lr_folder, sr_folder, tpp, pwd, 'lr_polishing', initial_bin_name, batch)
            #             pool.close()
            #             pool.join()
                        
            #             remain_bin = {}
            #             for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_gf_lr_polished'):
            #                 for file in files:
            #                     if '_gf_lr_polished.fa' in file:
            #                         bin_name=file.split('_gf_lr_polished.fa')[0]+'.fa'
            #                         # os.system('mv '+file+' '+binset_folder+'_gf_lr_polished')
            #                         complete_bin[bin_name]=0
                        
            #             for bin_name in bin_lr.keys():
            #                 if bin_name not in complete_bin.keys():
            #                     remain_bin[bin_name]=0
                        
            #             if len(remain_bin) == 0:
            #                 i=-1

            #         else:
            #             for bin_name in remain_bin.keys():
            #                 bin_id=bin_name.split('.fa')[0]
            #                 print('Processing '+str(bin_id))
            #                 initial_bin_name=remain_bin[bin_name]
            #                 try:
            #                     fp=open(str(bin_id)+'_gf_lr_polished.fa','a')
            #                     fp.close()
            #                 except:
            #                     fp=open(str(bin_id)+'_gf_lr_polished.fa','w')
            #                     fp.close()
            #                 sr_polishing(bin_id, bin_seq, lr_folder, sr_folder, num_threads, pwd, binset_folder, 'lr_polishing', initial_bin_name, mapping_tool, 1, major_run, batch)
            #             i=-1

            #     os.system('rm hs_err_pid* *_blasted.txt *_gf_lr_total.fa *.fa_elongation_lr.txt *_lr.fa')
            #     for i in range(1,51):
            #         os.system('rm *_lr_polished'+str(i)+'.fasta')
                
            #     f=open(polish_status_file,'a')
            #     f.write('Long-reads polished done!'+'\n')
            #     f.close()
            #     f=open('BASALT_log.txt','a')
            #     f.write('Long-reads polished done!'+'\n')
            #     f.close()       
        #######  ####### ####### ####### ####### ####### #######------------------------------------------------------  ####### ####### ####### ####### #######

        os.system('mkdir Polish_'+str(batch)+'_file')
        os.system('mv Accompleted_bins* Not_mapped_reads_* GF_lr_* Bin_long_read* Bin_id_with_* Long_read_bin* *_polishing_status.txt Bin_waiting_for_polish_batch.txt Polish_'+str(batch)+'_file')
        os.system('rm *bt2 bin* Total_bins.fa')
    print('Polishing done!')

if __name__ == '__main__': 
    # binset_folder='BestBinset_outlier_refined_filtrated_retrieved_OLC_re-assembly_OLC'
    binset_folder='BestBinset_outlier_refined'
    # binset_folder='bin7_test'
    # binset_folder='1_Opera_unpolished_cat_contigs_MAGs_polished_1'
    # assembly_list=['1_RH_insert_cat_270_final.contigs.fa','2_RH_S001_insert_270_final.contigs.fa','3_RH_S002_insert_270_final.contigs.fa','4_RH_S003_insert_270_final.contigs.fa','5_RH_S004_insert_270_final.contigs.fa','6_RH_S005_insert_270_final.contigs.fa']
    assembly_list=['1_assembly.fa','4_assembly.fa']
    # datasets_list={'1':['PE_r1_RH_S001_insert_270_1.fastq','PE_r2_RH_S001_insert_270_2.fastq'], '2':['PE_r1_RH_S002_insert_270_1.fastq','PE_r2_RH_S002_insert_270_2.fastq'], '3':['PE_r1_RH_S003_insert_270_1.fastq','PE_r2_RH_S003_insert_270_2.fastq'],'4':['PE_r1_RH_S004_insert_270_1.fastq','PE_r2_RH_S004_insert_270_2.fastq'],'5':['PE_r1_RH_S005_insert_270_1.fastq','PE_r2_RH_S005_insert_270_2.fastq']}
    # datasets_list={'1':['RH_S001_insert_270_mate1.fq','RH_S001_insert_270_mate2.fq'], '2':['RH_S002_insert_270_mate1.fq','RH_S002_insert_270_mate2.fq'], '3':['RH_S003_insert_270_mate1.fq','RH_S003_insert_270_mate2.fq'],'4':['RH_S004_insert_270_mate1.fq','RH_S004_insert_270_mate2.fq'],'5':['RH_S005_insert_270_mate1.fq','RH_S005_insert_270_mate2.fq']}
    datasets_list={'1':['SRR12358675_sub_1.fastq','SRR12358675_sub_2.fastq'],'2':['SRR12358676_sub_1.fastq','SRR12358676_sub_2.fastq']}
    # long_read=['anonymous_reads1.fq','anonymous_reads2.fq','anonymous_reads3.fq','anonymous_reads4.fq','anonymous_reads5.fq'] ### Write the name of long read here, if there is not, just let it to be blank
    long_read=['SRR12358673_sub.fastq','SRR12358674_sub.fastq','SRR12358678_sub.fastq']
    major_run=3
    num_threads=60
    ram=250
    mapping_tool='bwa' ### bw2 or bwa
    batch=1 ### batch 1: 1st polishing for hybrid assembly, which polishs the assembly before reassembly;
    ### batch 2: 2nd polishing for hybrid assembly, which polishs the assembly after reassembly.
    ### batch 3: polishing for SR assembly
    polishing_main(binset_folder, datasets_list, assembly_list, long_read, batch, mapping_tool, major_run, ram, num_threads)