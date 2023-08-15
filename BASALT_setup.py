#!/usr/bin/env python
import os

pwd=os.getcwd()
os.system('git clone https://github.com/EMBL-PKU/BASALT.git')

os.chdir('BASALT')
print('Creating BASALT environment')
print('This will take a while...')
os.system('conda env create -n BASALT --file basalt_env.yml')

print('Downloading models')
try:
    import requests
    from tqdm import tqdm
    os.system('python BASALT_models_download.py')
except:
    os.system('pip install tqdm')
    os.system('pip install requests')
    os.system('python BASALT_models_download.py')

print('BASALT installation done!')
