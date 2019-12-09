import pandas as pd
import numpy as np
import os
import subprocess
from multiprocessing import Pool
import errno


def makedir(cdir):
    try:
        os.makedirs(cdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def select_snps_random(n, snp_df=None):
    return np.random.choice(snp_df['SNP'].values, size=n, replace=False)


def generate_simulation(outdir, c_loci_file, h, c):

    h_outdir = os.path.join(outdir, "hsq_" + str(h).replace(".", "_"))

    cmds = ["/Users/szabad/software/gcta_1.92.4beta2_mac/bin/gcta",
            "--bfile", bed_file,
            "--simu-qt",
            "--simu-causal-loci", c_loci_file,
            "--simu-hsq", str(h),
            "--simu-rep", str(n_replicates),
            "--keep", os.path.join(keep_file_dir, str(c) + ".csv"),
            "--out", os.path.join(h_outdir, "sim_out")]

    subprocess.check_output(cmds)


# ---------------------------------

# inputs:

bed_file = "/Users/szabad/data/1000G/affy_6_biallelic_snps_maf005_thinned_aut"
snp_file = "~/data/1000G/affy_6_biallelic_snps_maf005_thinned_aut.bim"
cluster_assignment_file = "./inputs/hdbscan_labels_min10_1000G_UMAP_PC3_NC2_NN15_MD0.5_20184421291.txt"

# configs:

chr_n = 22
n_replicates = 100
n_causal_snps = 500
shared_snps = [1.0, 0.9, 0.75]
shared_effect_sizes = [1.0, 0.75, 0.5, 0.25]
hsq = [0.1, 0.3, 0.6, 0.9]
snp_select_func = select_snps_random

# output:

keep_file_dir = "./simulations/keep_files/"
simulation_output_dir = "./simulations/gcta_output/%d_causal_snps" % n_causal_snps

# ---------------------------------

makedir(keep_file_dir)

# Generate the keep files:

clust_df = pd.read_csv(cluster_assignment_file, header=None,
                       names=['FID', 'IID', 'CLUSTER'], sep="\s",
                       engine='python')

cluster_ids = np.unique(clust_df['CLUSTER'])

# Generate keep files for GCTA:
for c in cluster_ids:
    ddf = clust_df.loc[clust_df['CLUSTER'] == c, ]
    ddf.to_csv(os.path.join(keep_file_dir, str(c) + ".csv"), index=False,
               header=False, sep=' ')


# ---------------------------------

# Read SNP dataframe:

snp_df = pd.read_csv(snp_file, header=None,
                     names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], sep="\s",
                     engine='python')

if chr_n is not None:
    snp_df = snp_df.loc[snp_df['CHR'] == chr_n]

# ---------------------------------

if __name__ == '__main__':
    args = []

    for ss in shared_snps:

        n_shared_causal = int(np.round(ss*n_causal_snps))
        shared_causal_snps = snp_select_func(n_shared_causal)

        c_snp_df = snp_df.loc[snp_df['SNP'].isin(shared_causal_snps), ['SNP']].reset_index(drop=True)
        c_snp_df['ES'] = 0.

        if ss < 1.:
            rem_snp_df = snp_df.loc[~snp_df['SNP'].isin(shared_causal_snps), ['SNP']]

        for fs in shared_effect_sizes:

            n_shared_effect = int(np.round(fs * len(c_snp_df)))
            c_snp_df.iloc[np.arange(n_shared_effect), -1] = np.random.normal(size=n_shared_effect)

            for c in cluster_ids:
                if fs < 1.:
                    c_snp_df.iloc[np.arange(n_shared_effect, len(c_snp_df)), -1] = np.random.normal(
                        size=len(c_snp_df) - n_shared_effect)

                if ss < 1.:
                    nonshared_causal_snps = snp_select_func(n_causal_snps - n_shared_causal,
                                                            rem_snp_df)

                    rem_snp_df = rem_snp_df.loc[~rem_snp_df['SNP'].isin(nonshared_causal_snps), ['SNP']]

                    nonshared_causal_snps = snp_df.loc[snp_df['SNP'].isin(nonshared_causal_snps),
                                                       ['SNP']].reset_index(drop=True)

                    nonshared_causal_snps['ES'] = np.random.normal(size=len(nonshared_causal_snps))
                    final_c_snp_df = c_snp_df.append(nonshared_causal_snps).reset_index(drop=True)

                else:
                    final_c_snp_df = c_snp_df

                outdir = os.path.join(simulation_output_dir, "ss_" + str(ss).replace(".", "_"),
                                      "fs_" + str(fs).replace(".", "_"), str(c))
                makedir(outdir)

                c_loci_file = os.path.join(outdir, "causal_loci.csv")
                final_c_snp_df.to_csv(c_loci_file, index=False, header=False, sep="\t")

                for h in hsq:
                    h_outdir = os.path.join(outdir, "hsq_" + str(h).replace(".", "_"))
                    makedir(h_outdir)

                    args.append((outdir, c_loci_file, h, c))

    pool = Pool(6)

    pool.starmap(generate_simulation, args)

    pool.close()
    pool.join()
