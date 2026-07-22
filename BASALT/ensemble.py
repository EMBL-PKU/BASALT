#!/usr/bin/env python

"""
Inference-time ensemble for BASALT DL contamination classifiers.

This script loads multiple trained MLP checkpoints, performs ensemble
prediction on the test dataset, and writes the predicted labels to
an output file. It is used internally during BASALT outlier removal.
"""

import math
import argparse
import os
from pathlib import Path

from collections import Counter
from glob import glob

import scipy.io as scio
import torch
import numpy as np
import torch.optim as optim
from tensorboardX import SummaryWriter
from torch.utils.data import random_split
from tqdm import tqdm

from model import MLP
import torch.optim.lr_scheduler as lr_scheduler
from os.path import join
from my_dataset import MyDataSet_test
from basalt_runtime import require_model_directory
from utils import train_one_epoch, evaluate_ensemble, get_log_dir, del_best_ckpt, save_confusion_mat


def main(args):
    op=args.output
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    # device = torch.device("cpu")

    dataset = MyDataSet_test(type=args.norm_type, use_256=args.use_256, fea=args.fea, split='test', dec=args.dec)

    batch_size = args.batch_size
    nw = 2
    print(f'Using {nw} dataloader workers every process')

    val_loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=False, drop_last=False,
                                             pin_memory=True, num_workers=nw, persistent_workers=True)

    model = MLP(input_size=dataset.datas.shape[-1], num_classes=args.num_classes).to(device)
    prs = []
    # 0.944 0.792 # 0.942 0.742 # 0.933 0.746 # 0.924 0.748 # 0.922 0.740
    ensemble_dict = {5: '5_92_76_ensemble', 4: '4_90_71_ensemble', 3: '3_89_70_ensemble',
                     2: '2_89_70_ensemble', 1: '1_89_69_ensemble'}
    if args.ckpt_dir:
        model_dir = Path(args.ckpt_dir).expanduser().resolve()
    else:
        model_dir = require_model_directory()

    descriptor = model_dir / (ensemble_dict[dataset.n_col] + '.csv')
    if not descriptor.is_file():
        raise FileNotFoundError(
            "BASALT model descriptor not found: {}. Re-run the canonical "
            "BASALT_models_download.py installer.".format(descriptor)
        )
    ckpts = np.atleast_1d(np.loadtxt(descriptor, dtype=str))

    for ckpt in tqdm(ckpts):
        # print(str(ckpt))
        checkpoint_root = (model_dir / ensemble_dict[dataset.n_col]).resolve()
        checkpoint = (checkpoint_root / str(ckpt)).resolve()
        if checkpoint_root not in checkpoint.parents:
            raise ValueError(
                "Unsafe checkpoint path in BASALT model descriptor: {}".format(ckpt)
            )
        if not checkpoint.is_file():
            raise FileNotFoundError(
                "BASALT model checkpoint not found: {}. Re-run the canonical "
                "BASALT_models_download.py installer.".format(checkpoint)
            )
        model.load_state_dict(
            torch.load(checkpoint, map_location='cpu', weights_only=True)
        )
        pr = evaluate_ensemble(model=model,
                               data_loader=val_loader,
                               device=device)
        prs.append(pr)
    prs = np.array(prs)
    prs = (prs.sum(0) > (len(ckpts)//2 - args.mode)).astype(int)

    pr = np.array(prs)

    j = 0
    class_dict = {v: k for k, v in dataset.class_dict.items()}
    for i, line in tqdm(enumerate(dataset.lines)):
        l = line.split('\t')
        if len(l) != 1:
            dataset.lines[i] = '\t'.join([class_dict[pr[j]], *l])
            j += 1
    # with open('Predicted_potential_outlier.txt', 'w') as f:
    with open(op, 'w') as f:
        for l in dataset.lines:
            f.write(l)

    # class_dict = {'Real': 0, 'Contaminated': 1}
    # genes = []
    # labels, datas = [], []
    # with open(r"D:\新建文件夹\new_cls_pair\5-Opera_unpolishedCAT_BASALT_BestBinset_bin_contigs_type_depth.txt", 'r') as f:
    #     lines = f.readlines()
    # for i, line in tqdm(enumerate(lines)):
    #     l = line.split('\t')
    #     if len(l) == 1:
    #         if i != 0:
    #             labels.append(np.vstack(temp_labels))
    #         genes.append(l[0].strip())
    #         temp_labels, temp_datas = [], []
    #     else:
    #         label = l[0]
    #         if label in list(class_dict.keys()):
    #             label = class_dict[label]
    #             temp_labels.append(label)
    #         else:
    #             label = 0
    #             temp_labels.append(label)
    # labels.append(np.vstack(temp_labels))
    # labels = np.vstack(labels).flatten()
    # print(Counter(list(labels)))
    # tr = labels.flatten()
    # mask_0 = tr == 0
    # tr_0, pr_0 = tr[mask_0], pr[mask_0]
    # acc_0 = np.sum(tr_0==pr_0)/tr_0.shape[0]
    # mask_1 = tr == 1
    # tr_1, pr_1 = tr[mask_1], pr[mask_1]
    # acc_1 = np.sum(tr_1==pr_1)/tr_1.shape[0]
    # print(acc_0, acc_1)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('--norm_type', default='absmmn', choices=['mmn', 'absmmn', 'none'])
    parser.add_argument('--use_256', type=bool, default=True)
    parser.add_argument('--fea', type=bool, default=True)
    parser.add_argument('--dec', type=int, default=3)
    parser.add_argument('--mode', type=int, default=2)
    parser.add_argument('--ckpt_dir', type=str, default='', help='path to model checkpoint')
    parser.add_argument('--num_classes', type=int, default=2)
    parser.add_argument('--batch_size', type=int, default=512)
    parser.add_argument('-o', '--output', type=str, dest='output', default='Predicted_potential_outlier.txt')

    parser.add_argument('--device', default='cpu', help='device id (i.e. 0 or 0,1 or cpu)')

    opt = parser.parse_args()

    main(opt)
