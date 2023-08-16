#!/usr/bin/env python
import os

pwd=os.getcwd()

def git_clone(repo_url, max_attempts=3):
    for attempt in range(1, max_attempts + 1):
        result = os.system(f'git clone {repo_url}')
        if result == 0:
            return
        else:
            print(f"Attempt {attempt} failed.")
            if attempt == max_attempts:
                raise Exception("There was an issue with the network connection. Please attempt again.")

repo_url = "https://github.com/EMBL-PKU/BASALT.git"

git_clone(repo_url)

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
