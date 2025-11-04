#!/usr/bin/env python
from Bio import SeqIO
import sys, os, threading, copy, math
from multiprocessing import Pool

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
    f2=open('Bin_name_mod_quality_report.tsv','w')
    f2.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    n, record_seq, mod_bin_list, mod_bin_dict, bin_contig_num = 0, {}, [], {}, {}
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
                f1.write(str(file)+'\t'+'bin'+str(n)+'.fa'+'\n')
                for record in SeqIO.parse(file, 'fasta'):
                    m+=1
                    f.write('bin'+str(n)+'_'+str(m)+'\t'+str(record.id)+'\n')
                    record_seq['bin'+str(n)]['bin'+str(n)+'_'+str(m)]=str(record.seq)
                bin_contig_num['bin'+str(n)]=m

    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            if 'quality_report.tsv' in str(file) and str(file) != 'Bin_name_mod_quality_report.tsv':
                n=0
                for line in open(str(file),'r'):
                    n+=1
                    if n >= 2:
                        checkm_name=str(line).strip().split('\t')[0]
                        mod_bin_checkm_name=str(mod_bin_dict[checkm_name])
                        bins_checkm[mod_bin_checkm_name]={}

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

                        bins_checkm[mod_bin_checkm_name]['N50']=int(N50)
                        bins_checkm[mod_bin_checkm_name]['Completeness']=float(completeness)
                        bins_checkm[mod_bin_checkm_name]['Genome size']=int(eval(genome_size))
                        bins_checkm[mod_bin_checkm_name]['Contamination']=float(contamination)
                        f2.write(mod_bin_checkm_name+'\t'+str(genome_size)+'\t'+str(completeness)+'\t'+str(contamination)+'\t'+str(N50)+'\n')
    f.close()
    f1.close()
    f2.close()

    os.system('mv Mod_contig_id.txt Bin_name_mod.txt Bin_name_mod_quality_report.tsv '+pwd+'/'+str(binset_folder)+'_mod')
    
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
    return str(binset_folder)+'_mod', 'Total_bins.fa', mod_bin_list, bins_checkm

def parse_sam(sam_file, fq, pair, n):
    print('Reading reads id')
    f_not_mapped_reads=open('Not_mapped_reads.txt','a')
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
            bin_id=flist[2].split('_')[0]
            read_id=flist[0]
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
            read_num=reads.split('_')[-1] #
            read_id_name=str(n)+'_'+read_id.split('_')[0] #
            fq_seq=flist[9]+'\n'+'+'+'\n'+flist[10]+'\n'
            try:
                fq[bin_id][read_id_name]+=1
                f1=open(str(bin_id)+'_seq_R'+str(read_num)+'.fq','a')
                f1.write('@'+str(reads)[:-2]+' '+str(read_num)+'\n'+str(fq_seq))
                f1.close()
                if fq[bin_id][read_id_name] == 2:
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

def mapping_sr(total_fa, datasets_list, fq, pair, num_threads):
    os.system('bowtie2-build '+str(total_fa)+' '+str(total_fa))
    n = 0
    for item in datasets_list.keys():
        n+=1
        os.system('bowtie2 -p '+str(num_threads)+' -x '+str(total_fa)+' -1 '+str(datasets_list[item][0])+' -2 '+str(datasets_list[item][1])+' -S '+str(item)+'.sam -q --no-unal')
        parse_sam(str(item)+'.sam', fq, pair, n)
        os.system('rm '+str(item)+'.sam')

def parse_checkm(checkm_containing_folder,pwd):
    #pwd=os.getcwd()
    bins_checkm={}

    os.chdir(pwd+'/'+checkm_containing_folder)
    for root, dirs, files in os.walk(pwd+'/'+checkm_containing_folder):
        for file in files:        
            if 'quality_report.tsv' in file:
                n=0
                for line in open(file,'r'):
                    n+=1
                    if n >= 2:
                        genome_ids=str(line).strip().split('\t')[0]
                        if '_genomes.0' not in genome_ids:
                            bins_checkm[genome_ids]={}
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

                            bins_checkm[genome_ids]['N50']=int(N50)
                            bins_checkm[genome_ids]['Completeness']=float(completeness)
                            bins_checkm[genome_ids]['Genome size']=int(eval(genome_size))
                            bins_checkm[genome_ids]['Contamination']=float(contamination)
    return bins_checkm

def reassembly(bin_seq, reassembly_bin_folder, num_threads, bins_seq_folder, long_read, ram, pwd):
    for item in bin_seq.keys():
        os.system('spades.py -1 '+str(bin_seq[item][0])+' -2 '+str(bin_seq[item][1])+' -o '+str(item)+'_spades_reassembly --careful -t '+str(num_threads)+' -m '+str(ram))
        os.chdir(str(pwd)+'/'+str(item)+'_spades_reassembly')
        xxxx=0
        if os.path.isfile("contigs.fasta") == True:
            f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
            for record in SeqIO.parse('contigs.fasta', 'fasta'):
                if len(record.seq) >= 1000:
                    xxxx+=1
                    f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            f.close()
        else:
            os.chdir(str(pwd))
            os.system('spades.py -1 '+str(bin_seq[item][0])+' -2 '+str(bin_seq[item][1])+' -o '+str(item)+'_spades_reassembly -t '+str(num_threads)+' -m '+str(ram))
            os.chdir(str(pwd)+'/'+str(item)+'_spades_reassembly')
            if os.path.isfile("contigs.fasta") == True:
                f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
                for record in SeqIO.parse('contigs.fasta', 'fasta'):
                    if len(record.seq) >= 1000:
                        xxxx+=1
                        f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                f.close()
        os.chdir(str(pwd))

        if xxxx >= 2:
            os.system('mv '+str(item)+'_spades_reassembly/'+str(item)+'_SPAdes_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_re-assembly_contigs.fa')
        ### os.system('cp '+str(item)+'_spades_reassembly/mismatch_corrector/contigs/corrected_contigs.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_re-assembly_corrected_contigs.fa')
        ### os.system('/home/emma/MEGAHIT-1.2.2-beta-Linux-static/bin/megahit -1 '+str(bin_seq[item][0])+' -2 '+str(bin_seq[item][1])+' -o '+str(item)+'_megahit_reassembly --min-contig-len 1000 -t '+str(num_threads))
	    ### os.system('mv '+str(item)+'_megahit_reassembly/final.contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_megahit_re-assembly_contigs.fa')

        os.system('fq2fa --merge --filter '+str(bin_seq[item][0])+' '+str(bin_seq[item][1])+' idba.fa')
        n=0
        for record in SeqIO.parse('idba.fa','fasta'):
            n+=1
            if n == 1:
                seq_len=len(record.seq)
        
        if seq_len <= 125:
            os.system('idba_ud -r idba.fa -o '+str(item)+'_idba_reassembly --num_threads '+str(num_threads)+' --min_contig 1000')
        else:
            os.system('idba_ud -l idba.fa -o '+str(item)+'_idba_reassembly --num_threads '+str(num_threads)+' --min_contig 1000')

        os.system('mv '+str(item)+'_idba_reassembly/contig.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_IDBA_re-assembly_contigs.fa')

        os.system('rm -rf '+str(item)+'_spades_reassembly')
        os.system('rm -rf '+str(item)+'_idba_reassembly')
        os.system('mv '+str(bin_seq[item][0])+' '+str(bin_seq[item][1])+' '+str(bins_seq_folder))

def unicycler_mul(item, pwd, sr_folder, bin_seq, bin_lr, t_p_p, reassembly_bin_folder, bins_seq_folder):
    os.system('unicycler -1 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][0])+' -2 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][1])+' -l '+str(bin_lr[item])+' -o '+str(item)+'_unicycler_reassembly -t '+str(t_p_p)+' --mode conservative --no_pilon')
    os.system('mv '+str(item)+'_unicycler_reassembly/assembly.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_UNICYCLER_re-assembly_contigs.fa')
    # os.system('mv '+str(bin_seq[item][0])+' '+str(bin_seq[item][1])+' '+str(pwd)+'/'+str(sr_folder))
    os.system('mv '+str(bin_lr[item])+' '+str(pwd)+'/'+str(bins_seq_folder))
    os.system('rm -rf '+str(item)+'_unicycler_reassembly')

def reassembly_lr(bin_seq, bin_lr, reassembly_bin_folder, num_threads, bins_seq_folder, sr_folder, ram, pwd):
    try:
        os.system('java -Xmx'+str(ram)+'G -jar pilon-1.23.jar')
    except:
        print('pilon running in default ram')

    num_project=1
    if num_threads >= 40:
        if num_threads < 60:
            num_project=2
        else:
            num_project=math.ceil(num_threads/30)
    if ram >= 64:
        num_project2=math.ceil(ram/55)
        if num_project2 < num_project:
            num_project=num_project2

    if num_project == 1:
        for item in bin_lr.keys():
            if item in bin_seq.keys():
                ### Re-assembly using Unicycler
                print('Reassembling '+str(item))
                # os.system('mv '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][0])+' '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][1])+' '+str(pwd))
                os.system('unicycler -1 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][0])+' -2 '+str(pwd)+'/'+str(sr_folder)+'/'+str(bin_seq[item][1])+' -l '+str(bin_lr[item])+' -o '+str(item)+'_unicycler_reassembly -t '+str(num_threads)+' --mode conservative --no_pilon')
                os.system('mv '+str(item)+'_unicycler_reassembly/assembly.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_UNICYCLER_re-assembly_contigs.fa')
                # os.system('mv '+str(bin_seq[item][0])+' '+str(bin_seq[item][1])+' '+str(pwd)+'/'+str(sr_folder))
                os.system('mv '+str(bin_lr[item])+' '+str(pwd)+'/'+str(bins_seq_folder))
                os.system('rm -rf '+str(item)+'_unicycler_reassembly')
                ### Re-assembly using Flye
                ### os.system('flye --nano-raw '+str(bin_lr[item])+' --threads '+str(num_threads)+' -o '+str(item)+'_flye_reassembly --iterations 5')
                ### os.system('mv '+str(item)+'_flye_reassembly/assembly.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_FLYe_re-assembly_contigs.fa')
            else:
                xyz=0
                ### os.system('flye --nano-raw '+str(bin_lr[item])+' --threads '+str(num_threads)+' -o '+str(item)+'_flye_reassembly --iterations 5')
                ### os.system('mv '+str(item)+'_flye_reassembly/assembly.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_FLYe_re-assembly_contigs.fa')
    else:
        print('Processing', str(num_project), 'simultaneously')
        pool=Pool(processes=num_project)
        t_p_p=math.ceil(num_threads/num_project)
        for item in bin_lr.keys():
            if item in bin_seq.keys():
                print('Reassembling '+str(item))
                pool.apply_async(unicycler_mul, args=(item, pwd, sr_folder, bin_seq, bin_lr, t_p_p, reassembly_bin_folder, bins_seq_folder))
        pool.close()
        pool.join()

        os.system('mv '+item+' '+str(bins_seq_folder))
    # os.chdir(str(pwd))
    os.system('rm -rf '+str(item)+'_unicycler_reassembly')
    # os.system('rm -rf '+str(item)+'_flye_reassembly')

def bin_comparison(paired_bins, bin_checkm):
    pwd=os.getcwd()
    f=open('Reassembled_bins_comparison.txt','w')
    best_bin, best_bin_checkm={}, {}
    for item in paired_bins.keys():
        best_bin_checkm_name_list=item.split('.')
        best_bin_checkm_name_list.remove(best_bin_checkm_name_list[-1])
        best_bin_checkm_name='.'.join(best_bin_checkm_name_list)
        f.write(str(item)+'\t'+str(bin_checkm[best_bin_checkm_name])+'\n')
        for item2 in paired_bins[item]:
            reass_bin_checkm_name_list=item2.split('.')
            reass_bin_checkm_name_list.remove(reass_bin_checkm_name_list[-1])
            reass_bin_checkm_name='.'.join(reass_bin_checkm_name_list)
            f.write(str(item2)+'\t'+str(bin_checkm[reass_bin_checkm_name])+'\n')
            best_bin_cpn=bin_checkm[best_bin_checkm_name]['Completeness']
            best_bin_ctn=bin_checkm[best_bin_checkm_name]['Contamination']
            best_bin_ml=bin_checkm[best_bin_checkm_name]['N50']
            reass_bin_cpn=bin_checkm[reass_bin_checkm_name]['Completeness']
            reass_bin_ctn=bin_checkm[reass_bin_checkm_name]['Contamination']
            reass_bin_ml=bin_checkm[reass_bin_checkm_name]['N50']
        
            delta_cpn_ctn_bestbin=float(best_bin_cpn)-5*float(best_bin_ctn)
            delta_cpn_ctn_reass_bin=float(reass_bin_cpn)-float(5*reass_bin_ctn)

            if '_SPAdes_' in reass_bin_checkm_name or '_UNICYCLER_' in reass_bin_checkm_name:
                if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
                    best_bin_checkm_name=best_bin_checkm_name
                elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
                    best_bin_checkm_name=reass_bin_checkm_name
                elif delta_cpn_ctn_bestbin == delta_cpn_ctn_reass_bin:
                    if reass_bin_ml > best_bin_ml:
                        best_bin_checkm_name=reass_bin_checkm_name
                    else:
                        best_bin_checkm_name=best_bin_checkm_name
                else:
                    continue
            elif '_IDBA_' in reass_bin_checkm_name:
                delta_idba = delta_cpn_ctn_reass_bin-delta_cpn_ctn_bestbin
                if delta_idba < 3:
                    best_bin_checkm_name=best_bin_checkm_name
                elif delta_idba >= 3:
                    best_bin_checkm_name=reass_bin_checkm_name
                else:
                    continue
                # if best_bin_ml > reass_bin_ml:
                #    best_bin_checkm_name=best_bin_checkm_name
                # elif best_bin_ml < reass_bin_ml:
                #    best_bin_checkm_name=reass_bin_checkm_name
                # else:
                #    best_bin_checkm_name=reass_bin_checkm_name
            # elif '_UNICYCLER_' in reass_bin_checkm_name:
            #     if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=best_bin_checkm_name
            #     elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=reass_bin_checkm_name
            #     elif delta_cpn_ctn_bestbin == delta_cpn_ctn_reass_bin:
            #         if reass_bin_ml > best_bin_ml:
            #             best_bin_checkm_name=reass_bin_checkm_name
            #         else:
            #             best_bin_checkm_name=best_bin_checkm_name
            #     else:
            #         continue
            # elif '_FLYe_' in reass_bin_checkm_name:
            #     if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=best_bin_checkm_name
            #     elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
            #         best_bin_checkm_name=reass_bin_checkm_name
            #     else:
            #         continue

        best_bin[best_bin_checkm_name+'.fa']=best_bin_checkm_name
        best_bin_checkm[best_bin_checkm_name]=bin_checkm[best_bin_checkm_name].copy()
    f.close()
    return best_bin, best_bin_checkm

def re_assembly_main(binset_folder, datasets_list, long_read, hybri_reassembly, ram, num_threads):
    pwd=os.getcwd()
    try:
        f=open('Assembly_status.txt','a')
    except:
        f=open('Assembly_status.txt','w')
    f.close()

    try:
        os.mkdir(binset_folder+'_re-assembly_binset')
    except:
        print(binset_folder+'_re-assembly_binset exists')
    
    assembled_bins={}
    try:
        os.mkdir(binset_folder+'_sr_bins_seq')
    except:
        print(binset_folder+'_sr_bins_seq exists')
        for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_sr_bins_seq'):
            for file in files:        
                if '.fq' in file and '_seq_' in file and 'bin' in file:
                    assembled_bins[str(file).split('_seq_')[0].strip()]=0

    A=mod_bin(binset_folder)
    mod_bin_folder=A[0]
    total_fa=A[1]
    mod_bin_list=A[2]
    original_bins_checkm=A[3]
    fq, pair, bin_seq, bin_checkm={}, {}, {}, {}
    bin_checkm=original_bins_checkm.copy()
    for bin_id in mod_bin_list:
        fq[str(bin_id)]={}
        pair[str(bin_id)]={}
        bin_seq[str(bin_id)]=[]
        bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R1.fq')
        bin_seq[str(bin_id)].append(str(bin_id)+'_seq_R2.fq')
        if len(assembled_bins) == 0:
            f1=open(str(bin_id)+'_seq_R1.fq','w')
            f2=open(str(bin_id)+'_seq_R2.fq','w')
            f1.close()
            f2.close()

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Short-read mapping done!' in line:
            x=1

    if x == 0:
        f_not_mapped_reads=open('Not_mapped_reads.txt','w')
        f_not_mapped_reads.close()

        mapping_sr(total_fa, datasets_list, fq, pair, num_threads)
        f=open('Assembly_status.txt','a')
        f.write('Short-read mapping done!'+'\n')
        f.close()

    bin_seq2={}
    for item in bin_seq:
        if item not in assembled_bins.keys():
            bin_seq2[item]=bin_seq[item]

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Short-read assembly done!' in line:
            x=1

    if x == 0:
        reassembly(bin_seq2, binset_folder+'_re-assembly_binset', num_threads, binset_folder+'_sr_bins_seq', long_read, ram, pwd)
        f=open('Assembly_status.txt','a')
        f.write('Short-read assembly done!'+'\n')
        f.close()

    if len(long_read) != 0 and hybri_reassembly == 'y':
        print('Reassemblying long reads')
        assembled_bins={}
        try:
            os.mkdir(binset_folder+'_lr_bins_seq')
        except:
            print(binset_folder+'_lr_bins_seq exists')
            for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_lr_bins_seq'):
                for file in files:        
                    if '_lr.fq' in file and 'bin' in file:
                        assembled_bins[str(file).split('_lr')[0].strip()]=0
        os.chdir(pwd)

        if len(assembled_bins) == 0:
            print('Mapping long reads')
            x=0
            for line in open('Assembly_status.txt','r'):
                if 'Long-read mapping done!' in line:
                    x=1

            if x == 0:
                n, bin_lr=0, {}
                for lrs in long_read:
                    n+=1
                    os.system('minimap2 -t '+str(num_threads)+' -ax map-ont Total_bins.fa '+str(lrs)+' > lr'+str(n)+'.sam')
                    print('Splitting long reads '+str(n))
                    bin_lr.update(parse_lr_sam('lr'+str(n)+'.sam', lrs, n))

                f=open('Assembly_status.txt','a')
                f.write('Long-read mapping done!'+'\n')
                f.close()

                for i in range(1, n+1):
                    os.system('rm lr'+str(n)+'.sam')

        x=0
        for line in open('Assembly_status.txt','r'):
            if 'Long-read assembly done!' in line:
                x=1

        if x == 0:
            n, bin_lr=0, {}
            for lrs in long_read:
                n+=1
                for line in open('Bin_long_read'+str(n)+'.txt','r'):
                    bin_id=str(line).strip().split('\t')[0]
                    bin_lr[bin_id]=0
 
            bin_lr2={}
            for item in bin_lr.keys():
                if item not in assembled_bins.keys():
                    bin_lr2[item]=str(item)+'_lr.fq'
            
            reassembly_lr(bin_seq, bin_lr2, binset_folder+'_re-assembly_binset', num_threads, binset_folder+'_lr_bins_seq', binset_folder+'_sr_bins_seq', ram, pwd)
            f=open('Assembly_status.txt','a')
            f.write('Long-read assembly done!'+'\n')
            f.close()

    x=0
    for line in open('Assembly_status.txt','r'):
        if 'Assembly quality evaluation done!' in line:
            x=1
    
    if x == 0:
        print('Checking reassembled bins')
        os.chdir(str(binset_folder)+'_re-assembly_binset')
        for root, dirs, files in os.walk(pwd+'/'+str(binset_folder)+'_re-assembly_binset'):
            for file in files:
                try:
                    ff=open(file+'_filtrated.fa','w')
                    for record in SeqIO.parse(file, 'fasta'):
                        if len(record.seq) != 0:
                            ff.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
                    ff.close()
                    os.system('mv '+str(file)+'_filtrated.fa '+str(file))
                except:
                    xyzzzz=0
        os.chdir(pwd)
        
        # os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset_folder)+'_re-assembly_binset '+str(binset_folder)+'_re-assembly_binset_checkm')
        os.system('checkm2 predict -t '+str(num_threads)+' -i '+str(binset_folder)+'_re-assembly_binset -x fa -o '+str(binset_folder)+'_re-assembly_binset_checkm')
        f=open('Assembly_status.txt','a')
        f.write('Assembly quality evaluation done!'+'\n')
        f.close()

    os.chdir(pwd)
    reassembly_bin_checkm=parse_checkm(binset_folder+'_re-assembly_binset_checkm', pwd)
    bin_checkm.update(reassembly_bin_checkm)

    # checkm_t=open('bin_checkm_t.txt','w')
    # for item in bin_checkm.keys():
    #     checkm_t.write(str(item)+'\t'+str(bin_checkm[item])+'\n')
    # checkm_t.close()

    paired_bins, best_bin, best_bin_checkm={}, {}, {}
    os.chdir(pwd+'/'+binset_folder+'_re-assembly_binset')
    for root, dirs, files in os.walk(pwd+'/'+binset_folder+'_re-assembly_binset'):
        for file in files:
            if '_SPAdes_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_SPAdes_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
            elif '_IDBA_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_IDBA_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
            elif '_UNICYCLER_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_UNICYCLER_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
            elif '_FLYe_re-assembly_contigs.fa' in file:
                n_re=0
                for line in open(file,'r'):
                    n_re+=1
                    if n_re == 2:
                        break
                if n_re == 2:
                    o_bin=str(file).split('_FLYe_re-assembly_contigs.fa')[0]+'.fa'
                    if o_bin not in paired_bins.keys():
                        paired_bins[o_bin]=[]
                        paired_bins[o_bin].append(str(file))
                    else:
                        paired_bins[o_bin].append(str(file))
    os.chdir(pwd)

    # paired_bin_t=open('Paired_bins_t.txt','w')
    # for item in paired_bins.keys():
    #     paired_bin_t.write(str(item)+'\t'+str(paired_bins[item])+'\n')
    # paired_bin_t.close()

    best_bin_after_c=bin_comparison(paired_bins, bin_checkm)
    best_bin=best_bin_after_c[0].copy()
    best_bin_checkm=best_bin_after_c[1].copy()

    try:
        os.mkdir(binset_folder+'_re-assembly')
    except:
        print(binset_folder+'_re-assembly exists')
    
    selected_bins={}
    for item in best_bin.keys():
        if '_re-assembly_contigs.fa' in item:
            s_name=item.split('_')[0]+'.fa'
            selected_bins[s_name]=1
            os.system('cp '+pwd+'/'+binset_folder+'_re-assembly_binset/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')
        else:
            selected_bins[item]=1
            os.system('cp '+pwd+'/'+str(binset_folder)+'_mod/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')

    os.chdir(pwd+'/'+str(binset_folder)+'_mod')
    for root, dirs, files in os.walk(pwd+'/'+str(binset_folder)+'_mod'):
        for file in files:
            if '.fa' in file:
                if file not in selected_bins.keys():
                    os.system('cp '+pwd+'/'+str(binset_folder)+'_mod/'+file+' '+pwd+'/'+binset_folder+'_re-assembly')
                    item_checkm_name=file.split('.fa')[0]
                    best_bin_checkm[item_checkm_name]=bin_checkm[item_checkm_name]
    os.chdir(pwd)

    os.chdir(pwd+'/'+binset_folder+'_re-assembly')
    f=open('Best_binset_after_re-assembly_quality_report.tsv','w')
    f.write('Bin_ID'+'\t'+'Genome_size'+'\t'+'Completeness'+'\t'+'Contamination'+'\t'+'N50'+'\n')
    for item in best_bin_checkm.keys():
        # f.write(str(item)+'\t'+str(best_bin_checkm[item])+'\n')
        f.write(str(item)+'\t'+str(best_bin_checkm[item]['Genome size'])+'\t'+str(best_bin_checkm[item]['Completeness'])+'\t'+str(best_bin_checkm[item]['Contamination'])+'\t'+str(best_bin_checkm[item]['N50'])+'\n')
    f.close()
    os.chdir(pwd)

    # try:
    #     os.system('rm -rf '+binset_folder+'_sr_bins_seq')
    # except:
    #     x=1

    # try:
    #     os.system('rm -rf '+binset_folder+'_lr_bins_seq')
    # except:
    #     x=1
    print('Re-assembly done!')

if __name__ == '__main__': 
    binset_folder='1_Opera_unpolished_cat_contigs.fasta_BestBinsSet_outlier_refined_filtrated_retrieved'
    datasets_list={'1':['PE_r1_RH_S001_insert_270_mate1.fq','PE_r2_RH_S001_insert_270_mate2.fq'], '2':['PE_r1_RH_S002_insert_270_mate1.fq','PE_r2_RH_S002_insert_270_mate2.fq'], '3':['PE_r1_RH_S003_insert_270_mate1.fq','PE_r2_RH_S003_insert_270_mate2.fq'],'4':['PE_r1_RH_S004_insert_270_mate1.fq','PE_r2_RH_S004_insert_270_mate2.fq'],'5':['PE_r1_RH_S005_insert_270_mate1.fq','PE_r2_RH_S005_insert_270_mate2.fq']}
    long_read=['anonymous_reads1.fq','anonymous_reads2.fq','anonymous_reads3.fq','anonymous_reads4.fq','anonymous_reads5.fq'] ### Write the name of long read here, if there is not, just let it to be blank
    hybri_reassembly='n' ### Use Unicycler to re-assembly. e.g. --hybrid y / --hybrid n; defalt no
    num_threads=20
    ram=250
    re_assembly_main(binset_folder, datasets_list, long_read, hybri_reassembly, ram, num_threads)