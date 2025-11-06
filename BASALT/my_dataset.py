#!/usr/bin/env python
import os.path
import pickle
from collections import Counter
from glob import glob
from itertools import combinations

import scipy.io as scio
import numpy as np
from PIL import Image
import torch
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from torch.utils.data import Dataset
from tqdm import tqdm
import os
from utils import norm


def read_raw_data(split):
    lines = []
    if split == 'train':
        for p in glob(r'./train/*.txt'):
            with open(p, 'r') as f:
                lines += f.readlines()
    elif split == 'val':
        for p in glob(r'./val/*.txt'):
            with open(p, 'r') as f:
                lines += f.readlines()
    elif split == 'test':
        for p in glob(r'./BestBinset_outlier_test/*.txt'):
        # for p in glob(r'./1_assembly_sample1.fa_200_metabat_genomes.8.fa_coverage_TNFs_test_folder/*.txt'):
        # for p in glob(r'./test/*.txt'):
            with open(p, 'r') as f:
                lines += f.readlines()
    else:
        raise NameError
    return lines

class MyDataSet(Dataset):
    def __init__(self, type='mmn', use_256=False, fea=False, use_col=3, split='train'):
        class_dict = {'Real': 0, 'Contaminated': 1}
        self.class_dict = class_dict
        u = 't' if use_256 else 'f'
        f = 't' if fea else 'f'
        cache = type+'-'+u+'-'+f+'-'+split+'-'+str(use_col)+'.mat'
        if os.path.exists(cache):
            mat = scio.loadmat(cache)
            self.datas = mat['data']
            self.labels = mat['label'].flatten()
            print(Counter(list(self.labels)))
            self.datas = torch.tensor(self.datas).to(torch.float)
            self.labels = torch.tensor(self.labels).to(torch.long)


        else:
            genes = []
            labels, datas = [], []
            lines = read_raw_data(split)
            self.lines = lines
            for i, line in tqdm(enumerate(lines)):
                l = line.split('\t')
                if len(l) == 1:
                    if i != 0:
                        datas.append(norm(np.vstack(temp_datas), type))
                        labels.append(np.vstack(temp_labels))
                    genes.append(l[0].strip())
                    temp_labels, temp_datas = [], []
                else:
                    label = l[0]
                    if label in list(class_dict.keys()):
                        label = class_dict[label]
                        ind = list(eval(l[3]).keys())
                        cs = [c for c in combinations(ind, use_col)]
                        for c in cs:
                            data = [eval(l[3])[j] for j in c]
                            if fea and len(data) > 1:
                                data += get_time_domain_features(np.array(data))
                            if use_256:
                                data += list(eval(l[5]))
                            temp_labels.append(label)
                            temp_datas.append(data)
            datas.append(norm(np.vstack(temp_datas), type))
            labels.append(np.vstack(temp_labels))

            self.datas = np.vstack(datas)
            self.labels = np.vstack(labels).flatten()
            print(Counter(list(self.labels)))
            self.datas = torch.tensor(self.datas).to(torch.float)
            self.labels = torch.tensor(self.labels).to(torch.long)

    def __len__(self):
        return len(self.datas)

    def __getitem__(self, idx):
        return self.datas[idx], self.labels[idx]
class MyDataSet_(Dataset):
    def __init__(self, type='mmn', use_256=False, fea=False, split='train'):
        class_dict = {'Real': 0, 'Contaminated': 1}
        self.class_dict = class_dict

        genes = []
        labels, datas = [], []
        lines = read_raw_data(split)
        self.lines = lines
        for i, line in tqdm(enumerate(lines)):
            l = line.split('\t')
            if len(l) == 1:
                if i != 0:
                    datas.append(norm(np.vstack(temp_datas), type))
                    labels.append(np.vstack(temp_labels))
                genes.append(l[0].strip())
                temp_labels, temp_datas = [], []
            else:
                label = l[0]
                if label in list(class_dict.keys()):
                    label = class_dict[label]
                    ind = list(eval(l[3]).keys())
                    self.n_col = len(ind)
                    data = [eval(l[3])[j] for j in ind]
                    if fea and len(data) > 1:
                        data += get_time_domain_features(np.array(data))
                    if use_256:
                        data += list(eval(l[5]))
                    temp_labels.append(label)
                    temp_datas.append(data)
        datas.append(norm(np.vstack(temp_datas), type))
        labels.append(np.vstack(temp_labels))

        self.datas = np.vstack(datas)
        self.labels = np.vstack(labels).flatten()
        print(Counter(list(self.labels)))
        self.datas = torch.tensor(self.datas).to(torch.float)
        self.labels = torch.tensor(self.labels).to(torch.long)

    def __len__(self):
        return len(self.datas)

    def __getitem__(self, idx):
        return self.datas[idx], self.labels[idx]

class MyDataSet_test(Dataset):
    def __init__(self, type='mmn', use_256=False, fea=False, split='train', dec=3):
        self.dec = dec
        class_dict = {'Real': 0, 'Contaminated': 1}
        self.class_dict = class_dict
        genes = []
        datas = []
        lines = read_raw_data(split)
        self.lines = lines
        for i, line in tqdm(enumerate(lines)):
            l = line.split('\t')
            if len(l) == 1:
                if i != 0:
                    datas.append(norm(np.vstack(temp_datas), type))
                genes.append(l[0].strip())
                temp_labels, temp_datas = [], []
            else:
                data = eval(l[2])
                index = list(data.keys())
                self.n_col = len(index)
                data = [data[j] for j in index]
                if fea and len(data) > 1:
                    data += get_time_domain_features(np.array(data))
                if use_256:
                    data += list(eval(l[4]))
                temp_datas.append(data)
        datas.append(norm(np.vstack(temp_datas), type))

        self.datas = np.vstack(datas)
        if self.n_col > 5:
            self.pca()
        self.datas = torch.tensor(self.datas).to(torch.float)

    def __len__(self):
        return len(self.datas)

    def __getitem__(self, idx):
        return self.datas[idx]

    def pca(self):
        x = self.datas[:, :self.n_col]

        pca = PCA(n_components=self.dec)
        X_pca = pca.fit_transform(x)
        self.datas = np.hstack((X_pca, self.datas[:, self.n_col:]))

        # tsne = TSNE(n_components=3, random_state=42)
        # X_tsne = tsne.fit_transform(x)
        # self.datas = np.hstack((X_tsne, self.datas[:, self.n_col:]))

        self.n_col = self.dec

def get_time_domain_features(data):
    mean_ = np.mean(data)  # 1.均值
    var_ = np.var(data)  # 2.方差
    std_ = np.std(data)  # 3.标准差
    max_ = np.max(data)  # 4.最大值
    min_ = np.min(data)  # 5.最小值

    len_ = len(data)

    absXbar = np.sum(abs(data))
    x_r = np.sum(np.sqrt(abs(data)))
    S = np.sum((data - mean_) ** 3)
    K = np.sum((data - mean_) ** 4)
    x_rms = np.sum(data ** 2)

    x_p = np.max(abs(data))  # 6.峰值
    x_rms = np.sqrt(x_rms / (len_ + 1e-7))  # 7.均方根值
    absXbar = absXbar / (len_ + 1e-7)  # 8.绝对平均值
    x_r = (x_r / (len_ + 1e-7)) ** 2  # 9.方根幅值
    W = x_rms / (mean_ + 1e-7)  # 10.波形指标
    C = x_p / (x_rms + 1e-7)  # 11.峰值指标
    I = x_p / (mean_ + 1e-7)  # 12.脉冲指标
    L = x_p / (x_r + 1e-7)  # 13.裕度指标
    S = S / ((len_ - 1) * std_ ** 3 + 1e-7)  # 14.偏斜度
    K = K / ((len_ - 1) * std_ ** 4 + 1e-7)  # 15.峭度

    fea = [mean_, var_, std_, max_, min_, x_p, x_rms, absXbar, x_r, W, C, I, L, S, K]
    return fea

if __name__ == '__main__':
    dataset = MyDataSet(type='absmmn', use_256=True, fea=True, use_col=2, split='val')
    scio.savemat('./absmmn-t-t-val-1.mat', {'data': dataset.datas.numpy(), 'label': dataset.labels.numpy()})
    # print(len(dataset))
    # batch_size = 256
    # train_loader = torch.utils.data.DataLoader(dataset, batch_size=batch_size, shuffle=True,
    #                                            pin_memory=True, num_workers=2, persistent_workers=True)
    # for x, y in train_loader:
    #     print(x.shape, y.shape)