import os
import sys
import pandas as pd

n_causal_snps = 12
input_dir = "./simulations/gcta_output/%d_causal_snps" % n_causal_snps

try:
    snp_shared = next(os.walk(input_dir))[1]
    frac_shared = next(os.walk(os.path.join(input_dir, snp_shared[0])))[1]
    clust_ids = next(os.walk(os.path.join(input_dir, snp_shared[0], frac_shared[0])))[1]
except Exception as e:
    print(str(e))
    sys.exit()

for ss in snp_shared:
    print("=========================")
    print(ss)
    for fs in frac_shared:
        print("-----------")
        print(fs)
        clust_dfs = []

        for c in clust_ids:
            df = pd.read_csv(os.path.join(input_dir, ss, fs, c, "causal_loci.csv"), header=None,
                             sep="\t", engine='python', names=['SNP', 'ES'])
            clust_dfs.append(df)

        clust_snps = [set(cdf['SNP'].values) for cdf in clust_dfs]
        clust_es = [set(cdf['ES'].values) for cdf in clust_dfs]

        shared_causal = clust_snps[0].intersection(*clust_snps[1:])
        print("# Shared Causal:", len(shared_causal))
        print(float(len(shared_causal)) / len(clust_dfs[0]))

        # ----------------------

        shared_es = clust_es[0].intersection(*clust_es[1:])
        print("# Shared ES:", len(shared_es))
        print(float(len(shared_es)) / len(shared_causal))
