#!/usr/bin/env python
import os

pwd=os.getcwd()
os.system('git clone https://github.com/EMBL-PKU/BASALT.git')

os.chdir('BASALT')
print('Creating BASALT environment')
os.system('conda env create -n BASALT --file basalt_env.yml')

os.system('unzip BASALT_script.zip')
os.system('chmod -R 777 BASALT_script')
condaenv=os.popen('conda info --envs').read()
a=condaenv.split('\n')
for item in a:
    a=str(item).split('/')
    name=a[0].strip()
    if name == 'BASALT':
        path='/'.join(a[1:len(a)])
        path='/'+path+'/bin'

os.system('mv BASALT_script/* '+path)

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
