#!/usr/bin/env python
import torch.nn as nn

class LBR(nn.Module):
    def __init__(self, hidden_size):
        super(LBR, self).__init__()
        self.fc1 = nn.Linear(hidden_size, 4*hidden_size)
        self.bn1 = nn.BatchNorm1d(4*hidden_size)
        self.at1 = nn.ReLU()
        self.fc2 = nn.Linear(4*hidden_size, 4*hidden_size)
        self.bn2 = nn.BatchNorm1d(4*hidden_size)
        self.at2 = nn.ReLU()
        self.fc3 = nn.Linear(4*hidden_size, hidden_size)
    def forward(self, x):
        res = x
        x = self.fc1(x)
        x = self.bn1(x)
        x = self.at1(x)
        x = self.fc2(x)
        x = self.bn2(x)
        x = self.at2(x)
        x = self.fc3(x)
        return x + res



class MLP(nn.Module):
    def __init__(self, num_classes, input_size, hidden_size=256):
        super(MLP, self).__init__()
        self.embed = nn.Linear(input_size, hidden_size)
        self.m1 = LBR(hidden_size)
        self.m2 = LBR(hidden_size)
        self.cls = nn.Linear(hidden_size, num_classes)

    def forward(self, x):
        x = self.embed(x)
        x = self.m1(x) + x
        x = self.m2(x) + x
        x = self.cls(x)
        return x


