#!/usr/bin/env python
from Bio import SeqIO
import sys, os, threading, copy, math
from multiprocessing import Pool
# from glob import glob


def hybrid_parse_checkm(checkm_containing_folder,pwd):
    #pwd=os.getcwd()
    bins_checkm={}

    os.chdir(pwd+'/'+checkm_containing_folder)
    for root, dirs, files in os.walk(pwd+'/'+checkm_containing_folder):
        for file in files:        
            if 'bin_stats_ext.tsv' in file: 
                for line in open(file,'r'):
                    genome_ids=str(line).strip().split('\t')[0]
                    if '_genomes.0' not in genome_ids:
                        bins_checkm[genome_ids]={}

                        try:
                            marker_lineage=str(line).strip().split('\'marker lineage\': \'')[1].strip().split('\'')[0]
                            bins_checkm[genome_ids]['marker lineage']=marker_lineage
                        except:
                            print('marker lineage error')
                            bins_checkm[genome_ids]['marker lineage']=root

                        try:    
                            completeness=str(line).strip().split('\'Completeness\': ')[1].split(', ')[0]
                            bins_checkm[genome_ids]['Completeness']=float(completeness)
                        except:
                            bins_checkm[genome_ids]['Completeness']=0

                        try:
                            genome_size=str(line).strip().split('\'Genome size\':')[1].strip().split(', ')[0]
                            bins_checkm[genome_ids]['Genome size']=int(genome_size)
                        except:
                            bins_checkm[genome_ids]['Genome size']=0

                        bins_checkm[genome_ids]['Contamination']=float(str(line).strip().split('Contamination\': ')[1].split('}')[0].split(',')[0])

                        try:
                            bins_checkm[genome_ids]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split(',')[0].strip())
                        except:
                            bins_checkm[genome_ids]['Mean scaffold length']=float(str(line).strip().split('Mean scaffold length\':')[1].split('}')[0].strip())
    os.chdir(pwd)
    return bins_checkm

def assembly_mul(bins_seq_folder, bin_seq, item, reassembly_bin_folder, pwd, num_threads, ram):
    # try:
    #     os.system('mv '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq.gz '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq.gz '+pwd)
    #     os.system('gzip -d '+str(item)+'_seq_R1.fq.gz')
    #     os.system('gzip -d '+str(item)+'_seq_R2.fq.gz')
    # except:
    #     print(str(item)+' corrected reads does not exist')

    try:
        fx=open('Basalt_log.txt','a')
    except:
        fx=open('Basalt_log.txt','w')
    fx.write('Assembling '+str(item)+' using SPAdes'+'/n')
    os.system('spades.py -1 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq -2 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq -o '+str(item)+'_spades_reassembly --careful -t '+str(num_threads)+' -m '+str(ram))
    os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/contigs.fasta '+str(item)+'_contigs.fasta')
    
    xxxx=0
    if os.path.isfile(str(item)+'_contigs.fasta') == True:
        f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
        for record in SeqIO.parse(str(item)+'_contigs.fasta', 'fasta'):
            if len(record.seq) >= 1000:
                xxxx+=1
                f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
        f.close()

        os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R1* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
        os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R2* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')
    else:
        os.system('spades.py -1 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq -2 '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq -o '+str(item)+'_spades_reassembly -t '+str(num_threads)+' -m '+str(ram))
        # os.chdir(str(pwd)+'/'+str(item)+'_spades_reassembly')
        os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/contigs.fasta '+str(item)+'_contigs.fasta')
        if os.path.isfile(str(item)+'_contigs.fasta') == True:
            f=open(str(item)+'_SPAdes_re-assembly_contigs.fa','w')
            for record in SeqIO.parse(str(item)+'_contigs.fasta', 'fasta'):
                if len(record.seq) >= 1000:
                    xxxx+=1
                    f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
            f.close()

        os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R1* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
        os.system('mv '+pwd+'/'+str(item)+'_spades_reassembly/corrected/'+str(item)+'_seq_R2* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')
    
    os.system('rm '+str(item)+'_contigs.fasta')

    fx.close()
    if xxxx >= 2:
        os.system('mv '+str(item)+'_SPAdes_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_re-assembly_contigs.fa')
    ### os.system('cp '+str(item)+'_spades_reassembly/mismatch_corrector/contigs/corrected_contigs.fasta '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_re-assembly_corrected_contigs.fa')
    ### os.system('/home/emma/MEGAHIT-1.2.2-beta-Linux-static/bin/megahit -1 '+str(bin_seq[item][0])+' -2 '+str(bin_seq[item][1])+' -o '+str(item)+'_megahit_reassembly --min-contig-len 1000 -t '+str(num_threads))
    ### os.system('mv '+str(item)+'_megahit_reassembly/final.contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_megahit_re-assembly_contigs.fa')

    os.system('fq2fa --merge --filter '+str(pwd)+'/'+str(bins_seq_folder)+'/'+str(bin_seq[item][0])+' '+str(pwd)+'/'+str(bins_seq_folder)+'/'+str(bin_seq[item][1])+' '+str(item)+'_idba.fa')
    n=0
    for record in SeqIO.parse(str(item)+'_idba.fa','fasta'):
        n+=1
        if n == 1:
            seq_len=len(record.seq)
    
    if seq_len <= 125:
        os.system('idba_ud -r '+str(item)+'_idba.fa -o '+str(item)+'_idba_reassembly --num_threads '+str(num_threads)+' --min_contig 1000')
    else:
        os.system('idba_ud -l '+str(item)+'_idba.fa -o '+str(item)+'_idba_reassembly --num_threads '+str(num_threads)+' --min_contig 1000')

    os.system('mv '+pwd+'/'+str(item)+'_idba_reassembly/contig.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_IDBA_re-assembly_contigs.fa')

    os.system('rm -rf '+str(item)+'_spades_reassembly')
    os.system('rm -rf '+str(item)+'_idba_reassembly')
    # os.system('mv '+str(bin_seq[item][0])+' '+str(bin_seq[item][1])+' '+str(bins_seq_folder))
    os.system('rm '+str(item)+'_idba.fa')

def SR_reassembly(bin_seq, reassembly_bin_folder, num_threads, bins_seq_folder, long_read, ram, pwd):
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
    
    os.system('mkdir SPAdes_corrected_reads')
    # for item in bin_seq.keys():
    #     os.system('spades.py -1 '+str(pwd)+'/'+str(bins_seq_folder)+'/'+str(bin_seq[item][0])+' -2 '+str(pwd)+'/'+str(bins_seq_folder)+'/'+str(bin_seq[item][1])+' -o '+str(item)+'_SPAdes_corrected_reads --only-error-correction -t '+str(num_threads)+' -m '+str(ram))
    #     os.system('mv '+pwd+'/'+str(item)+'_SPAdes_corrected_reads/corrected/'+str(item)+'_seq_R1* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
    #     os.system('mv '+pwd+'/'+str(item)+'_SPAdes_corrected_reads/corrected/'+str(item)+'_seq_R2* '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')
    #     os.system('rm -rf '+str(item)+'_SPAdes_corrected_reads')

    print('Processing', str(num_project), 'simultaneously')
    pool=Pool(processes=num_project)
    t_p_p=math.ceil(num_threads/num_project)
    for item in bin_seq.keys():
        print('Reassembling '+str(item))
        pool.apply_async(assembly_mul, args=(bins_seq_folder, bin_seq, item, reassembly_bin_folder, pwd, t_p_p, ram))
        # assembly_mul('SPAdes_corrected_reads', bins_seq_folder, bin_seq, item, reassembly_bin_folder, pwd, t_p_p, ram)
    pool.close()
    pool.join()

def hybrid_bin_comparison(paired_bins, bin_checkm):
    pwd=os.getcwd()
    f=open('Reassembled_bins_comparison.txt','w')
    best_bin, best_bin_checkm={}, {}
    for item in paired_bins.keys():
        for item2 in paired_bins[item]:
            if '_polished' in item2:
                best_bin_checkm_name_list=item2.split('.')
                best_bin_checkm_name_list.remove(best_bin_checkm_name_list[-1])
                best_bin_checkm_name='.'.join(best_bin_checkm_name_list)
                # best_bin_cpn=bin_checkm[best_bin_checkm_name]['Completeness']
                # best_bin_ctn=bin_checkm[best_bin_checkm_name]['Contamination']
                # best_bin_ml=bin_checkm[best_bin_checkm_name]['Mean scaffold length']
                f.write(str(item2)+'\t'+str(bin_checkm[best_bin_checkm_name])+'\n')

        for item2 in paired_bins[item]:
            if '_polished' not in item2:   
                reass_bin_checkm_name_list=item2.split('.')
                reass_bin_checkm_name_list.remove(reass_bin_checkm_name_list[-1])
                reass_bin_checkm_name='.'.join(reass_bin_checkm_name_list)
                f.write(str(item2)+'\t'+str(bin_checkm[reass_bin_checkm_name])+'\n')
                best_bin_cpn=bin_checkm[best_bin_checkm_name]['Completeness']
                best_bin_ctn=bin_checkm[best_bin_checkm_name]['Contamination']
                best_bin_ml=bin_checkm[best_bin_checkm_name]['Mean scaffold length']
                reass_bin_cpn=bin_checkm[reass_bin_checkm_name]['Completeness']
                reass_bin_ctn=bin_checkm[reass_bin_checkm_name]['Contamination']
                reass_bin_ml=bin_checkm[reass_bin_checkm_name]['Mean scaffold length']
        
                delta_cpn_ctn_bestbin=float(best_bin_cpn)-5*float(best_bin_ctn)
                delta_cpn_ctn_reass_bin=float(reass_bin_cpn)-float(5*reass_bin_ctn)

                delta_cpn_ctn_bestbin_x=float(best_bin_cpn)-float(best_bin_ctn)
                delta_cpn_ctn_reass_bin_x=float(reass_bin_cpn)-float(reass_bin_ctn)

                if '_hybird_' in reass_bin_checkm_name:
                    if delta_cpn_ctn_bestbin > delta_cpn_ctn_reass_bin:
                        if delta_cpn_ctn_bestbin >= 40:
                            best_bin_checkm_name=best_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name
                    elif delta_cpn_ctn_bestbin < delta_cpn_ctn_reass_bin:
                        if delta_cpn_ctn_reass_bin >= 40:
                            best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name
                    elif delta_cpn_ctn_bestbin == delta_cpn_ctn_reass_bin:
                        if reass_bin_ml > best_bin_ml:
                            best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            best_bin_checkm_name=best_bin_checkm_name
                    else:
                        continue
                else:
                    delta_idba = delta_cpn_ctn_reass_bin-delta_cpn_ctn_bestbin
                    if delta_idba < 3:
                        if delta_cpn_ctn_bestbin >= 40:
                            best_bin_checkm_name=best_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name
                    elif delta_idba >= 3 and delta_idba < 6:
                        ratio=reass_bin_ml/best_bin_ml
                        if ratio >= 0.5:
                            if delta_cpn_ctn_reass_bin >= 40:
                                best_bin_checkm_name=reass_bin_checkm_name
                            else:
                                if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                    best_bin_checkm_name=best_bin_checkm_name
                                elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                    best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            best_bin_checkm_name=best_bin_checkm_name
                    elif delta_idba >= 6:
                        if delta_cpn_ctn_reass_bin >= 40:
                            best_bin_checkm_name=reass_bin_checkm_name
                        else:
                            if delta_cpn_ctn_bestbin_x >= delta_cpn_ctn_reass_bin_x:
                                best_bin_checkm_name=best_bin_checkm_name
                            elif delta_cpn_ctn_bestbin_x < delta_cpn_ctn_reass_bin_x: 
                                best_bin_checkm_name=reass_bin_checkm_name 
                    else:
                        continue

        best_bin[best_bin_checkm_name+'.fa']=best_bin_checkm_name
        best_bin_checkm[best_bin_checkm_name]=bin_checkm[best_bin_checkm_name].copy()
    f.close()
    return best_bin, best_bin_checkm

def hybrid_assembly_mul(sr_folder, bin_seq, item, bin_lr_reads, lr_folder, reassembly_bin_folder, pwd, num_threads, ram):
    try:
        os.system('gzip -d '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq.gz')
        os.system('gzip -d '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq.gz')
    except:
        print(str(item)+' corrected reads does not exist')

    # try:
    #     os.system('mv '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R1.fq.gz '+pwd+'/'+bins_seq_folder+'/'+str(item)+'_seq_R2.fq.gz '+pwd)
    #     os.system('gzip -d '+str(item)+'_seq_R1.fq.gz')
    #     os.system('gzip -d '+str(item)+'_seq_R2.fq.gz')
    # except:
    #     print(str(item)+' corrected reads does not exist')


    if os.path.isfile(pwd+'/'+sr_folder+'/'+str(item)+'_seq_R1.fq') == True:
        os.system('spades.py -1 '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R1.fq -2 '+pwd+'/SPAdes_corrected_reads/'+str(item)+'_seq_R2.fq --nanopore '+pwd+'/'+lr_folder+'/'+str(bin_lr_reads[item])+' -o '+str(item)+'_spades_hybrid_reassembly --careful --only-assembler -t '+str(num_threads)+' -m '+str(ram))
        # os.system('rm '+str(item)+'_seq_R1.fq '+str(item)+'_seq_R2.fq')
    else:
        os.system('spades.py -1 '+pwd+'/'+sr_folder+'/'+str(bin_seq[item][0])+' -2 '+pwd+'/'+sr_folder+'/'+str(bin_seq[item][1])+' --nanopore '+pwd+'/'+lr_folder+'/'+str(bin_lr_reads[item])+' -o '+str(item)+'_spades_hybrid_reassembly --careful -t '+str(num_threads)+' -m '+str(ram))
    os.system('mv '+str(pwd)+'/'+str(item)+'_spades_hybrid_reassembly/contigs.fasta '+str(item)+'_contigs.fasta')
    # os.chdir(str(pwd)+'/'+str(item)+'_spades_hybrid_reassembly')
    xxxx=0
    if os.path.isfile(str(item)+'_contigs.fasta') == True:
        f=open(str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa','w')
        for record in SeqIO.parse(str(item)+'_contigs.fasta', 'fasta'):
            if len(record.seq) >= 1000:
                xxxx+=1
                f.write('>'+str(record.id)+'\n'+str(record.seq)+'\n')
        f.close()
        os.system('rm '+str(item)+'_contigs.fasta')
    # os.chdir(str(pwd))
    
    if xxxx >= 1:
        os.system('mv '+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa')
        # os.system('mv '+str(item)+'_spades_hybrid_reassembly/'+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa '+str(reassembly_bin_folder)+'/'+str(item)+'_SPAdes_hybrid_re-assembly_contigs.fa')
    os.system('rm -rf '+str(item)+'_spades_hybrid_reassembly')

def hybrid_re_assembly_main(binset_folder, sr_folder, lr_folder, ram, num_threads):
    pwd=os.getcwd()
    bin_checkm={}
    try:
        f=open('Hybrid_re-assembly_status.txt','a')
    except:
        f=open('Hybrid_re-assembly_status.txt','w')
    f.close()

    reassembly_bin_folder=binset_folder+'_re-assembly_binset'
    try:
        os.mkdir(reassembly_bin_folder)
    except:
        print(reassembly_bin_folder+' exists')

    bin_seq={}
    os.chdir(pwd+'/'+binset_folder)
    for root, dirs, files in os.walk(pwd+'/'+binset_folder):
        for file in files:
            if '_mag_polished.fa' in file:
                bin_id=str(file).split('_mag_polished.fa')[0]
                bin_seq[bin_id]=[]
                os.system('cp '+str(file)+' '+pwd+'/'+reassembly_bin_folder)

    os.chdir(pwd+'/'+sr_folder)
    for root, dirs, files in os.walk(pwd+'/'+sr_folder):
        for file in files:
            if '_seq_R1.fq' in file:
                bin_id=str(file).split('_seq_R1.fq')[0]
                bin_seq[bin_id].append(file)
            elif '_seq_R2.fq' in file:
                bin_id=str(file).split('_seq_R2.fq')[0]
                bin_seq[bin_id].append(file)
    os.chdir(pwd)
    
    x=0
    for line in open('Hybrid_re-assembly_status.txt','r'):
        if 'Short-read assembly done!' in line:
            x=1

    if x == 0:
        SR_reassembly(bin_seq, reassembly_bin_folder, num_threads, sr_folder, lr_folder, ram, pwd)
        f=open('Hybrid_re-assembly_status.txt','a')
        f.write('Short-read assembly done!'+'\n')
        f.close()

    if lr_folder != '':
        x=0
        for line in open('Hybrid_re-assembly_status.txt','r'):
            if 'Hybrid-assembly done!' in line:
                x=1

        if x == 0:
            bin_lr_reads={}
            os.chdir(pwd+'/'+lr_folder)
            for root, dirs, files in os.walk(pwd+'/'+lr_folder):
                for file in files:
                    if '_lr.fq' in file:
                        bin_id=str(file).split('_lr.fq')[0]
                        bin_lr_reads[bin_id]=file
            os.chdir(pwd)

            ####
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
            
            print('Processing', str(num_project), 'simultaneously')
            pool=Pool(processes=num_project)
            t_p_p=math.ceil(num_threads/num_project)
            for item in bin_lr_reads.keys():
                print('Reassembling '+str(item))
                pool.apply_async(hybrid_assembly_mul, args=(sr_folder, bin_seq, item, bin_lr_reads, lr_folder, reassembly_bin_folder, pwd, t_p_p, ram))
            pool.close()
            pool.join()

            f=open('Hybrid_re-assembly_status.txt','a')
            f.write('Hybrid-assembly done!'+'\n')
            f.close()

    x=0
    for line in open('Hybrid_re-assembly_status.txt','r'):
        if 'Re-assembly quality evaluation done!' in line:
            x=1
    
    if x == 0:
        os.system('checkm lineage_wf -t '+str(num_threads)+' -x fa '+str(binset_folder)+'_re-assembly_binset '+str(binset_folder)+'_re-assembly_binset_checkm')
        f=open('Hybrid_re-assembly_status.txt','a')
        f.write('Re-assembly quality evaluation done!'+'\n')
        f.close()

    reassembly_bin_checkm=hybrid_parse_checkm(binset_folder+'_re-assembly_binset_checkm/storage',pwd)
    bin_checkm.update(reassembly_bin_checkm)

    x=0
    for line in open('Hybrid_re-assembly_status.txt','r'):
        if 'Bin selecting done!' in line:
            x=1
    
    if x == 0:
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
                elif '_SPAdes_hybrid_re-assembly_contigs.fa' in file:
                    n_re=0
                    for line in open(file,'r'):
                        n_re+=1
                        if n_re == 2:
                            break
                    if n_re == 2:
                        o_bin=str(file).split('_SPAdes_hybrid_re-assembly_contigs.fa')[0]+'.fa'
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
                elif '_mag_polished.fa' in file:
                    n_re=0
                    for line in open(file,'r'):
                        n_re+=1
                        if n_re == 2:
                            break
                    if n_re == 2:
                        o_bin=str(file).split('_mag_polished.fa')[0]+'.fa'
                        if o_bin not in paired_bins.keys():
                            paired_bins[o_bin]=[]
                            paired_bins[o_bin].append(str(file))
                        else:
                            paired_bins[o_bin].append(str(file))
        os.chdir(pwd)

        best_bin_after_c=hybrid_bin_comparison(paired_bins, bin_checkm)
        best_bin=best_bin_after_c[0].copy()
        best_bin_checkm=best_bin_after_c[1].copy()

        try:
            os.mkdir(binset_folder+'_re-assembly')
        except:
            print(binset_folder+'_re-assembly exists')
        
        selected_bins={}
        for item in best_bin.keys():
            os.system('cp '+pwd+'/'+reassembly_bin_folder+'/'+item+' '+pwd+'/'+binset_folder+'_re-assembly')

        os.chdir(pwd+'/'+binset_folder+'_re-assembly')
        f=open('Best_binset_after_re-assembly_bin_stats_ext.tsv','w')
        for item in best_bin_checkm.keys():
            f.write(str(item)+'\t'+str(best_bin_checkm[item])+'\n')
        f.close()
        os.chdir(pwd)

        f=open('Hybrid_re-assembly_status.txt','a')
        f.write('Bin selecting done!'+'\n')
        f.close()
    print('Re-assembly done!')

if __name__ == '__main__': 
    binset_folder='BestBinset_outlier_refined_filtrated_retrieved_MAGs_polished'
    # datasets_list={'1':['PE_r1_RH_S001_insert_270_mate1.fq','PE_r2_RH_S001_insert_270_mate2.fq'], '2':['PE_r1_RH_S002_insert_270_mate1.fq','PE_r2_RH_S002_insert_270_mate2.fq'], '3':['PE_r1_RH_S003_insert_270_mate1.fq','PE_r2_RH_S003_insert_270_mate2.fq'],'4':['PE_r1_RH_S004_insert_270_mate1.fq','PE_r2_RH_S004_insert_270_mate2.fq'],'5':['PE_r1_RH_S005_insert_270_mate1.fq','PE_r2_RH_S005_insert_270_mate2.fq']}
    # long_read=['anonymous_reads1.fq','anonymous_reads2.fq','anonymous_reads3.fq','anonymous_reads4.fq','anonymous_reads5.fq'] ### Write the name of long read here, if there is not, just let it to be blank
    sr_folder='BestBinset_outlier_refined_filtrated_retrieved_polished_sr_bins_seq'
    lr_folder='BestBinset_outlier_refined_filtrated_retrieved_long_read'
    num_threads=60
    ram=250
    hybrid_re_assembly_main(binset_folder, sr_folder, lr_folder, ram, num_threads)