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
    git_clone(repo_url)
    try:
        flog=open('BASALT_installation_log.txt','a')
    except:
        flog=open('BASALT_installation_log.txt','w')

    stats=[]
    for line in open('BASALT_installation_log.txt','r'):
        stats.append(line.strip())

    if 'BASALT environment installation done!' not in stats:
        os.chdir('BASALT')
        print('Creating BASALT environment')
        print('This will take a while...')
        os.system('conda env create -n BASALT --file basalt_env.yml')

        # os.system('unzip BASALT_script.zip')
        # os.system('chmod -R 777 BASALT_script')
        condaenv=os.popen('conda info --envs').read()
        a=condaenv.split('\n')
        for item in a:
            a=str(item).split('/')
            name=a[0].strip()
            if name == 'BASALT':
                path='/'.join(a[1:len(a)])
                path='/'+path+'/bin'
                os.system('chmod -R 777 '+path+'/*')

        # os.system('mv BASALT_script/* '+path)
        flog('BASALT environment installation done!')
    flog.close()

    if 'BASALT model download done!' not in stats:
        print('Downloading models')
        print('This will take a while...')
        download_model()
        user_dir = os.path.expanduser('~')
        destination_folder = f"{user_dir}/.cache"
        source_file = destination_folder+"/BASALT.zip"

        # os.system unzip 
        exit_code = os.system(f"unzip -o {source_file} -d {destination_folder}")

        if exit_code == 0:
            print(f"successfully unzip")
        else:
            print(f"unzip failed: please check your unzip configuration")

        flog=open('BASALT_installation_log.txt','a')    
        flog('BASALT model download done!')
        flog.close()

    print('BASALT installation done!')
    flog=open('BASALT_installation_log.txt','a')    
    flog('BASALT installation done!')
    flog.close()
