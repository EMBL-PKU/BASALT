#!/usr/bin/env python
# -*- coding: UTF-8 -*-
#coding=utf-8

import pandas as pd
# import seaborn as sns
# import matplotlib.pyplot as plt
from Bio import SeqIO
import os, sys, threading
from sklearn.decomposition import PCA
import numpy as np
from multiprocessing import Pool
# from Outlier_remover import *

def checkm_tax_wf(marker_lineage_level_name, marker_lineage_taxon, num, folder_name):
    os.system('checkm taxonomy_wf -t '+str(num)+' -x fa '+str(marker_lineage_level_name)+' '+str(marker_lineage_taxon)+' '+str(folder_name)+' '+str(folder_name)+'_checkm')

if __name__ == '__main__': 
    marker_lineage_level_name_num, threads_num, tax_num, folder_name_num = 0, 0, 0, 0
    for i in range(1, len(sys.argv)):
        if '-l' in str(sys.argv[i]):
            marker_lineage_level_name_num=int(i)+1
        elif '-t' in str(sys.argv[i]):
            threads_num=int(i)+1
        elif '-c' in str(sys.argv[i]):
            tax_num=int(i)+1
        elif '-f' in str(sys.argv[i]):
            folder_name_num=int(i)+1
        else:
            continue

    marker_lineage_level_name=str(sys.argv[marker_lineage_level_name_num])
    num=int(sys.argv[threads_num])
    marker_lineage_taxon=str(sys.argv[tax_num])
    folder_name=str(sys.argv[folder_name_num])
    checkm_tax_wf(marker_lineage_level_name, marker_lineage_taxon, num, folder_name)