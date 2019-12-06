import pandas as pd
import numpy as np
import subprocess
import os
import sys
import errno
import glob
from multiprocessing import Pool


# ---------------------------------

# configs:

gwas_methods = ['stratified', 'mega']
n_causal_snps = 12
shared_snps = 'ss_0_75'

# inputs:

input_dir = "./simulations/gcta_output/%d_causal_snps/%s/" % (n_causal_snps, shared_snps)

bed_file = "/Users/szabad/data/1000G/affy_6_biallelic_snps_maf005_thinned_aut"
cluster_assignment_file = "./inputs/hdbscan_labels_min10_1000G_UMAP_PC3_NC2_NN15_MD0.5_20184421291.txt"
cluster_files_dir = "./simulations/keep_files"
clust_covar_file = "./inputs/affy_6_biallelic_snps_maf005_thinned_aut_pcs_10.eigenvec"

# output:

output_dir = input_dir.replace("simulations", "gwas_results").replace("gcta", "plink")

# ---------------------------------


def makedir(cdir):
    try:
        os.makedirs(cdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def combine_phenotypes(fs, h, plink_outdir):

    dfs = []

    for c_dir in glob.glob(os.path.join(input_dir, fs, "*/")):
        df = pd.read_csv(os.path.join(c_dir, h, "sim_out.phen"), sep="\s",
                         engine="python", header=None)
        dfs.append(df)

    dfs = pd.concat(dfs)

    pheno_file = os.path.join(plink_outdir, "combined.pheno")

    dfs.to_csv(pheno_file, header=False, index=False, sep=" ")

    return pheno_file


def run_gwas(fs, h, method):

    pheno_outdir = os.path.join(output_dir, fs, h)
    makedir(pheno_outdir)

    pheno_file = combine_phenotypes(fs, h, pheno_outdir)

    plink_outdir = os.path.join(pheno_outdir, method)
    makedir(plink_outdir)

    plink_cmds = ["/Users/szabad/software/plink2/plink2",
                  "--allow-no-sex",
                  "--bfile", bed_file,
                  "--linear",
                  "--chr", "22",
                  "--pheno", pheno_file,
                  "--covar", clust_covar_file,
                  "--covar-variance-standardize",
                  "--adjust"]

    if method == 'stratified':
        for c in clust_ids:
            clust_dir = os.path.join(plink_outdir, c)
            makedir(clust_dir)
            plink_cmds += ["--out", os.path.join(clust_dir, "gwas_out"),
                           "--keep", os.path.join(cluster_files_dir, c + ".csv")]

            subprocess.check_output(plink_cmds)
            subprocess.check_output(["/Users/szabad/software/homebrew/bin/pigz", "-r", clust_dir])
            plink_cmds = plink_cmds[:-4]
    else:

        plink_cmds += ["--out", os.path.join(plink_outdir, "gwas_out")]

        subprocess.check_output(plink_cmds)
        subprocess.check_output(["/Users/szabad/software/homebrew/bin/pigz", "-r", plink_outdir])


if __name__ == '__main__':
    # Grab the simulation configurations from directory structure:
    try:
        frac_shared = next(os.walk(input_dir))[1]
        clust_ids = next(os.walk(os.path.join(input_dir, frac_shared[0])))[1]
        hsq = next(os.walk(os.path.join(input_dir, frac_shared[0], clust_ids[0])))[1]
    except Exception as e:
        print(str(e))
        sys.exit()

    args = []

    for fs in frac_shared:
        for h in hsq:
            for gwm in gwas_methods:
                args.append((fs, h, gwm))

    pool = Pool(2)

    pool.starmap(run_gwas, args)

    pool.close()
    pool.join()
