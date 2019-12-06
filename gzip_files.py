import os
from multiprocessing import Pool
import subprocess

input_dir = "./gwas_results"
def gzip_file(fname):
    subprocess.check_output(["gzip", fname])

files_to_zip = []

for root, dirs, files in os.walk(input_dir):
    for file in files:
        if ".gz" not in file and ".DS_Store" not in file:
            files_to_zip.append(os.path.join(root, file))

pool = Pool(6)
pool.map(gzip_file, files_to_zip)
pool.close()
pool.join()
