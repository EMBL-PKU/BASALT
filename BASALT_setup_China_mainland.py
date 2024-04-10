#!/usr/bin/env python
import os
try:
    import requests
    from tqdm import tqdm
except:
    os.system('pip install tqdm')
    os.system('pip install requests')
    import requests
    from tqdm import tqdm

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

def download_model(local_dir=None):
    if local_dir is None:
        user_dir = os.path.expanduser('~')
        # url = "https://github.com/LinB203/test/releases/download/v1/pytorch_model.bin" 
        url = "https://figshare.com/ndownloader/files/41093033"
        local_dir = f"{user_dir}/.cache"
    local_path = f"{local_dir}/BASALT.zip"

    if os.path.exists(local_path):
        print("File already exists.")
    else:
        os.makedirs(local_dir, exist_ok=True)
        print(f"File will be saved in {local_path}.")
        response = requests.get(url, stream=True)
        total_size = int(response.headers.get('content-length', 0))
        block_size = 1024
        progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True)

        with open(local_path, 'wb') as f:
            for data in response.iter_content(block_size):
                progress_bar.update(len(data)) 
                f.write(data)

        progress_bar.close()

        print("File downloaded successfully.")

if __name__ =='__main__':
    pwd=os.getcwd()
    
    print('Creating BASALT environment')
    print('This will take a while...')
    # os.system('site=https://mirrors.tuna.tsinghua.edu.cn/anaconda')
    # os.system('conda config --add channels $/{site/}/pkgs/free/')
    # os.system('conda config --add channels $/{site/}/pkgs/main/')
    # os.system('conda config --add channels $/{site/}/cloud/conda-forge/')
    # os.system('conda config --add channels $/{site/}/cloud/bioconda/')
    
    os.system('conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/free/')
    os.system('conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/pkgs/main/')
    os.system('conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge/')
    os.system('conda config --add channels https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda/')
    
    git_clone(repo_url)
    os.chdir('BASALT')
    os.system('conda env create -n BASALT --file basalt_env.yml')

    condaenv=os.popen('conda info --envs').read()
    a=condaenv.split('\n')
    for item in a:
        a=str(item).split('/')
        name=a[0].strip()
        if name == 'BASALT':
            path='/'.join(a[1:len(a)])
            path='/'+path+'/bin/*'
            os.system('chmod -R 777 '+path)

    # print('Downloading models')
    # print('This will take a while...')
    # download_model()
    # user_dir = os.path.expanduser('~')
    # destination_folder = f"{user_dir}/.cache"
    # source_file = destination_folder+"/BASALT.zip"

    # # os.system unzip 
    # exit_code = os.system(f"unzip -o {source_file} -d {destination_folder}")

    # if exit_code == 0:
    #     print(f"successfully unzip")
    # else:
    #     print(f"unzip failed: please check your unzip configuration")
        
    print('BASALT scripts and environment installation done!')
    print('You need to manually download the BASALT models following the instruction.')