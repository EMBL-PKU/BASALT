#!/usr/bin/env python
import os
import sys
from glob import glob

from matplotlib import pyplot as plt
from sklearn.metrics import confusion_matrix, accuracy_score, ConfusionMatrixDisplay
from torch import nn
import torch
from torch.nn import functional as F
from tqdm import tqdm
import numpy as np

import os
import requests
from tqdm import tqdm

def download_model(url, local_dir=None):
    if local_dir is None:
        user_dir = os.path.expanduser('~')
        local_dir = f"{user_dir}/.cache/BASALT"
    local_path = f"{local_dir}/{os.path.basename(url)}"  # 指定本地保存路径和文件名

    if os.path.exists(local_path):
        print(f"File already exists in {local_path}.")
    else:
        os.makedirs(local_dir, exist_ok=True)
        print(f"File will be saved in {local_path}.")
        response = requests.get(url, stream=True)  # 向URL发送请求
        total_size = int(response.headers.get('content-length', 0))  # 获取文件总大小
        block_size = 1024  # 每次下载的块大小
        progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True)  # 创建进度条对象

        with open(local_path, 'wb') as f:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data))  # 更新进度条
                f.write(data)  # 将下载的数据写入到本地文件中

        progress_bar.close()  # 关闭进度条


    return local_dir

def del_best_ckpt(log_dir):
    for i in glob(os.path.join(log_dir, '*best*')):
        os.remove(i)

# def norm(t, type='mmn'):
#     if type == 'mmn':
#         return (t - np.min(t, axis=0)) / (np.max(t, axis=0) - np.min(t, axis=0))
#     elif type == 'absmmn':
#         t = np.abs(t-np.mean(t, axis=0))
#         return t/np.max(t, axis=0)
#     else:
#         return t

def norm(t, type='mmn'):
    if type == 'mmn':
        return (t - np.min(t, axis=0)) / (np.max(t, axis=0) - np.min(t, axis=0) + 1e-8)
    elif type == 'absmmn':
        t = np.abs(t-np.mean(t, axis=0))
        return t/(np.max(t, axis=0) + 1e-8)
    else:
        return t


def get_log_dir(set='train', comment=''):
    log_dir = os.path.join('runs', set, comment)
    if os.path.exists(log_dir):
        raise FileExistsError
    os.makedirs(log_dir, exist_ok=True)
    return log_dir


def train_one_epoch(model, optimizer, data_loader, device, epoch, weight, loss_fun):
    model.train()
    if loss_fun == 'ce':
        loss_function = torch.nn.CrossEntropyLoss(weight=torch.tensor([weight, 1-weight], dtype=torch.float).to(device))
        # loss_function = torch.nn.CrossEntropyLoss()
    elif loss_fun == 'fl':
        loss_function = focal_loss(alpha=weight, num_classes=2)
    mean_loss = torch.zeros(1).to(device)
    optimizer.zero_grad()

    data_loader = tqdm(data_loader, file=sys.stdout)

    for step, data in enumerate(data_loader):
        x, labels = data

        pred = model(x.to(device))

        loss = loss_function(pred, labels.to(device))
        loss.backward()
        mean_loss = (mean_loss * step + loss.detach()) / (step + 1)  # update mean losses

        data_loader.desc = "[train epoch {}] mean loss {}".format(epoch, round(mean_loss.item(), 3))

        if not torch.isfinite(loss):
            print('WARNING: non-finite loss, ending training ', loss)
            sys.exit(1)

        optimizer.step()
        optimizer.zero_grad()

    return mean_loss.item()


@torch.no_grad()
def evaluate(model, data_loader, device, epoch, weight, loss_fun):
    model.eval()

    if loss_fun == 'ce':
        loss_function = torch.nn.CrossEntropyLoss(weight=torch.tensor([weight, 1-weight], dtype=torch.float).to(device))
        # loss_function = torch.nn.CrossEntropyLoss()
    elif loss_fun == 'fl':
        loss_function = focal_loss(alpha=weight, num_classes=2)

    tr, pr = [], []

    mean_loss = torch.zeros(1).to(device)
    data_loader = tqdm(data_loader, file=sys.stdout)

    for step, data in enumerate(data_loader):
        x, labels = data
        pred = model(x.to(device))

        loss = loss_function(pred, labels.to(device))
        mean_loss = (mean_loss * step + loss.detach()) / (step + 1)  # update mean losses
        data_loader.desc = "[val epoch {}] mean loss {}".format(epoch, round(mean_loss.item(), 3))

        pred = torch.max(pred, dim=1)[1]
        # sum_num += torch.eq(pred, labels.to(device)).sum()
        tr += list(labels.detach().cpu().numpy().flatten())
        pr += list(pred.detach().cpu().numpy().flatten())

    tr = np.array(tr)
    pr = np.array(pr)
    # print(accuracy_score(tr, pr))
    mask_0 = tr == 0
    tr_0, pr_0 = tr[mask_0], pr[mask_0]
    acc_0 = np.sum(tr_0==pr_0)/tr_0.shape[0]
    mask_1 = tr == 1
    tr_1, pr_1 = tr[mask_1], pr[mask_1]
    acc_1 = np.sum(tr_1==pr_1)/tr_1.shape[0]

    return acc_0, acc_1

def save_confusion_mat(log_dir, matrix, classes=None, title=None):
    plt.figure()
    disp = ConfusionMatrixDisplay(confusion_matrix=matrix, display_labels=classes)
    disp.plot()
    if title:
        plt.title(title)
    name = title or 'confusion_mat'
    plt.savefig(f'{log_dir}/{name}.png')


@torch.no_grad()
def evaluate_ensemble(model, data_loader, device):
    model.eval()


    tr, pr = [], []

    for step, data in enumerate(data_loader):
        x = data
        pred = model(x.to(device))
        pred = torch.max(pred, dim=1)[1]
        # sum_num += torch.eq(pred, labels.to(device)).sum()
        # tr += list(labels.detach().cpu().numpy().flatten())
        pr += list(pred.detach().cpu().numpy().flatten())

    # tr = np.array(tr)
    pr = np.array(pr)

    # print(accuracy_score(tr, pr))

    # mask_0 = tr == 0
    # tr_0, pr_0 = tr[mask_0], pr[mask_0]
    # acc_0 = np.sum(tr_0==pr_0)/tr_0.shape[0]
    # mask_1 = tr == 1
    # tr_1, pr_1 = tr[mask_1], pr[mask_1]
    # acc_1 = np.sum(tr_1==pr_1)/tr_1.shape[0]
    #
    # print(f'Contaminated total: {np.sum(mask_1)}\n'
    #       f'Contaminated correctly removed: {np.sum(tr_1==pr_1)}\n'
    #       f'Contaminated wrong remained: {np.sum(mask_1)-np.sum(tr_1==pr_1)}\n'
    #       f'Real total: {np.sum(mask_0)}\n'
    #       f'Real wrong removed: {np.sum(mask_0)-np.sum(tr_0==pr_0)}\n'
    #       f'Real correctly remained: {np.sum(tr_0==pr_0)}\n')

    return pr


class focal_loss(nn.Module):
    def __init__(self, alpha=0.25, gamma=2, num_classes=3, size_average=True):
        super(focal_loss, self).__init__()
        self.size_average = size_average
        if isinstance(alpha, list):
            assert len(alpha) == num_classes
            self.alpha = torch.Tensor(alpha)
        else:
            assert alpha < 1
            self.alpha = torch.zeros(num_classes)
            self.alpha[0] += alpha
            self.alpha[1:] += (1 - alpha)

        self.gamma = gamma

    def forward(self, preds, labels):
        preds = preds.view(-1, preds.size(-1))
        self.alpha = self.alpha.to(preds.device)
        preds_logsoft = F.log_softmax(preds, dim=1)  # log_softmax
        preds_softmax = torch.exp(preds_logsoft)  # softmax

        preds_softmax = preds_softmax.gather(1, labels.view(-1, 1))  # 这部分实现nll_loss ( crossempty = log_softmax + nll )
        preds_logsoft = preds_logsoft.gather(1, labels.view(-1, 1))
        self.alpha = self.alpha.gather(0, labels.view(-1))
        loss = -torch.mul(torch.pow((1 - preds_softmax), self.gamma),
                          preds_logsoft)  # torch.pow((1-preds_softmax), self.gamma) 为focal loss中 (1-pt)**γ

        loss = torch.mul(self.alpha, loss.t())
        if self.size_average:
            loss = loss.mean()
        else:
            loss = loss.sum()
        return loss
