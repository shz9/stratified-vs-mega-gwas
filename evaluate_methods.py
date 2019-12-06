import pandas as pd
import numpy as np
import os
import sys
import glob
import errno
from multiprocessing import Pool
import matplotlib.pylab as plt
import seaborn as sns


pd.set_option('display.max_columns', 10)

# ---------------------------------

# configs:
clust_size_threshold = 100
n_causal_snps = 12
shared_snps = 'ss_0_75'
pvalue_correction = 'FDR_BH'
scores = ['Recall', 'Precision', 'Non-shared Recall', 'Non-shared Precision']

# inputs:

cluster_assignment_file = "./inputs/hdbscan_labels_min10_1000G_UMAP_PC3_NC2_NN15_MD0.5_20184421291.txt"
sim_input = "./simulations/gcta_output/%d_causal_snps/%s/" % (n_causal_snps, shared_snps)
gwas_input = "./gwas_results/plink_output/%d_causal_snps/%s/" % (n_causal_snps, shared_snps)

# output:

output_dir = gwas_input.replace("gwas_results/plink_output", "evaluation")

# ---------------------------------


def makedir(cdir):
    try:
        os.makedirs(cdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def score_method(inf_file, true_snps, fs):

    inferred_causal = pd.read_csv(inf_file, sep="\s+", engine='python')
    inferred_causal = inferred_causal.loc[inferred_causal[pvalue_correction] < .05, ]

    if len(inferred_causal) < 1:
        return None
    else:

        recall = float(len(inferred_causal.loc[inferred_causal['ID'].isin(true_snps)])) / len(true_snps)
        precision = float(len(inferred_causal.loc[inferred_causal['ID'].isin(true_snps)])) / len(inferred_causal)

        if fs != 1:
            fs_true_snps = true_snps[int(fs*n_causal_snps):]
        else:
            fs_true_snps = true_snps

        fs_recall = float(len(inferred_causal.loc[inferred_causal['ID'].isin(fs_true_snps)])) / len(fs_true_snps)
        fs_precision = float(len(inferred_causal.loc[inferred_causal['ID'].isin(fs_true_snps)])) / len(inferred_causal)

        return recall, precision, fs_recall, fs_precision


def evaluate_overall_fishing_performance():

    res_df = pd.DataFrame(columns=['Heritability', 'GWAS Strategy', 'Fraction Shared'] + scores)

    for fs in frac_shared:

        fsf = float(fs.replace("fs_", "").replace("_", "."))

        for h in hsq:

            print(fs, h)

            hf = float(h.replace("hsq_", "").replace("_", "."))

            for m in gwas_methods:

                args = []

                if m == 'stratified':
                    for c in clust_ids:
                        true_causal = pd.read_csv(os.path.join(sim_input, fs, str(c), "causal_loci.csv"),
                                                  sep="\t", engine='python',
                                                  header=None, names=['SNP', 'ES'])
                        args += [(fname, true_causal['SNP'].values, fsf)
                                 for fname in glob.glob(os.path.join(gwas_input,
                                                                     fs, h, m,
                                                                     str(c), "*.linear.adjusted"))]
                else:
                    causal_set = []

                    for c in clust_ids:
                        true_causal = pd.read_csv(os.path.join(sim_input, fs, str(c), "causal_loci.csv"),
                                                  sep="\t", engine='python',
                                                  header=None, names=['SNP', 'ES'])
                        causal_set += list(true_causal['SNP'].values)

                    causal_set = np.array(list(set(causal_set)))

                    args = [(fname, causal_set, fsf)
                            for fname in glob.glob(os.path.join(gwas_input,
                                                                fs, h,
                                                                m, "*.linear.adjusted"))]

                pool = Pool(5)

                res = pool.starmap(score_method, args)

                pool.close()
                pool.join()

                for score in res:
                    if score is not None:
                        res_df.loc[len(res_df)] = [hf, m, fsf] + list(score)

    print(res_df.head())
    print(res_df.describe())

    for sc in scores:
        sns.catplot(x="Heritability", y=sc,
                    hue="GWAS Strategy", col="Fraction Shared",
                    data=res_df,
                    kind="violin")

        plt.savefig(os.path.join(output_dir, sc + ".pdf"))
        plt.close()


if __name__ == '__main__':

    print("# Causal SNPs: %d" % n_causal_snps)

    clust_df = pd.read_csv(cluster_assignment_file, header=None,
                           names=['FID', 'IID', 'CLUSTER'], sep="\s",
                           engine='python')

    clust_size = clust_df.groupby('CLUSTER').size()

    print(clust_size)

    makedir(output_dir)

    # Grab the simulation configurations from directory structure:
    try:
        frac_shared = next(os.walk(sim_input))[1]
        clust_ids = next(os.walk(os.path.join(sim_input, frac_shared[0])))[1]
        hsq = next(os.walk(os.path.join(sim_input, frac_shared[0], clust_ids[0])))[1]
        gwas_methods = next(os.walk(os.path.join(gwas_input, frac_shared[0], hsq[0])))[1]
    except Exception as e:
        print(str(e))
        sys.exit()

    clust_ids = [c for c in clust_ids if c in clust_size[clust_size > clust_size_threshold].index.astype(str)]

    print(clust_ids)

    evaluate_overall_fishing_performance()
