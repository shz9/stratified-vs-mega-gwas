import pandas as pd
import numpy as np
import os
import sys
import glob
import errno
from multiprocessing import Pool

pd.set_option('display.max_columns', 10)

# ---------------------------------

# configs:
clust_size_threshold = 100

# inputs:

snp_bim_file = "~/data/1000G/affy_6_biallelic_snps_maf005_thinned_aut.bim"
cluster_assignment_file = "./inputs/hdbscan_labels_min10_1000G_UMAP_PC3_NC2_NN15_MD0.5_20184421291.txt"
sim_input = "./simulations/gcta_output/"
gwas_input = "./gwas_results/plink_output/"

# output:

output_dir = gwas_input.replace("gwas_results/plink_output", "evaluation")

# ---------------------------------


def makedir(cdir):
    try:
        os.makedirs(cdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def calculate_power(snp, sighits_set, ld_window=1e5):
    """
    LD window length selected according to recommendation by Fang et al. 2019 (AJHG)
    """

    snp_loc = snp_meta_df.loc[snp_meta_df['SNP'] == snp, 'BP']
    valid_snps = snp_meta_df.loc[np.abs(snp_meta_df['BP'] - snp_loc) <= ld_window, 'SNP']

    if np.any(np.isin(valid_snps, list(sighits_set))):
        return 1.
    else:
        return 0.


def score_analysis_methods(input_dir, true_snps, params, pval_threshold=0.05,
                           pval_correction='BONF'):

    print(input_dir)
    fs, h = params

    # -----------------------------------------------------------
    # Generate sets of true causal snps for this configuration:

    # Sets based on SNP ID
    causal_sets = [set(ci['SNP']) for ci in true_snps]
    causal_union = set.union(*causal_sets)
    causal_intersect = set.intersection(*causal_sets)

    causal_nonshared_snp = [c.difference(causal_intersect) for c in causal_sets]
    causal_nonshared_snp_union = causal_union.difference(causal_intersect)

    # Sets based on SNP ID and effect size
    causal_es_sets = [set(ci.itertuples(index=False, name=None)) for ci in true_snps]
    causal_es_union = set.union(*causal_es_sets)
    causal_es_intersect = set.intersection(*causal_es_sets)
    causal_es_intersect_list = set([ci[0] for ci in causal_es_intersect])

    causal_nonshared_es = [set([ci[0] for ci in list(c.difference(causal_es_intersect))])
                           for c in causal_es_sets]
    causal_nonshared_es_union = set([ci[0] for ci in list(causal_es_union.difference(causal_es_intersect))])

    # -----------------------------------------------------------
    # Create a dataframe to compare statistical power of both methods:

    causal_snp_df = pd.DataFrame({'SNP': list(causal_union)})
    causal_snp_df['Heritability'] = h
    causal_snp_df['Fraction Shared'] = fs
    causal_snp_df['Presence'] = causal_snp_df['SNP'].apply(
        lambda x: ['Exclusive', 'Shared'][x in causal_intersect])
    causal_snp_df['Effect Size'] = causal_snp_df['SNP'].apply(
        lambda x: ['Exclusive', 'Shared'][x in causal_es_intersect_list])

    # -----------------------------------------------------------
    scores = []

    for method in gwas_methods:

        causal_snp_df[method] = 0.
        causal_snp_df[method + "_LD"] = 0.

        if method == 'stratified':
            rep_files = glob.glob(os.path.join(input_dir, method, clust_ids[0], "*.linear.adjusted.gz"))

            # Loop over replicates:
            for rf in rep_files:
                base_rf = os.path.basename(rf)
                sighits = []
                for c in clust_ids:
                    df = pd.read_csv(os.path.join(input_dir, method, c, base_rf), sep="\s+", engine='python')
                    sighits.append(set(df.loc[df[pval_correction] < pval_threshold, 'ID'].values))

                sighits_union = set.union(*sighits)

                causal_snp_df[method] += causal_snp_df['SNP'].apply(
                    lambda x: calculate_power(x, sighits_union, ld_window=0.) / len(rep_files))
                causal_snp_df[method + "_LD"] += causal_snp_df['SNP'].apply(
                    lambda x: calculate_power(x, sighits_union) / len(rep_files))

                try:
                    nonshared_snp_recall = float(len(set.union(*[sighits[i].intersection(causal_nonshared_snp[i])
                                                                 for i in range(len(clust_ids))]))) / \
                                           len(causal_nonshared_snp_union)
                except Exception as e:
                    nonshared_snp_recall = None

                try:
                    nonshared_es_recall = float(len(set.union(*[sighits[i].intersection(causal_nonshared_es[i])
                                                                for i in range(len(clust_ids))]))) / \
                                          len(causal_nonshared_es_union)
                except Exception as e:
                    nonshared_es_recall = None

                try:
                    overall_precision = float(len(set.union(*[sighits[i].intersection(causal_sets[i])
                                                               for i in range(len(clust_ids))]))) / len(sighits_union)
                except Exception as e:
                    overall_precision = None

                scores.append({
                    'GWAS Strategy': method,
                    'Heritability': h,
                    'Fraction Shared ES': fs,
                    'Overall Recall': float(len(set.union(*[sighits[i].intersection(causal_sets[i])
                                                            for i in range(len(clust_ids))]))) / len(causal_union),
                    'Overall Precision': overall_precision,
                    'Non-shared SNP Recall': nonshared_snp_recall,
                    'Non-shared ES Recall': nonshared_es_recall
                })

        else:
            rep_files = glob.glob(os.path.join(input_dir, method, "*.linear.adjusted.gz"))

            for rf in rep_files:
                df = pd.read_csv(rf, sep="\s+", engine='python')

                sighits = set(df.loc[df[pval_correction] < pval_threshold, 'ID'].values)

                causal_snp_df[method] += causal_snp_df['SNP'].apply(
                    lambda x: calculate_power(x, sighits, ld_window=0.) / len(rep_files))
                causal_snp_df[method + "_LD"] += causal_snp_df['SNP'].apply(
                    lambda x: calculate_power(x, sighits) / len(rep_files))

                try:
                    nonshared_snp_recall = float(len(sighits.intersection(causal_nonshared_snp_union))) / \
                                           len(causal_nonshared_snp_union)
                except Exception as e:
                    nonshared_snp_recall = None

                try:
                    nonshared_es_recall = float(len(sighits.intersection(causal_nonshared_es_union))) / \
                                          len(causal_nonshared_es_union)
                except Exception as e:
                    nonshared_es_recall = None

                try:
                    overall_precision = float(len(sighits.intersection(causal_union))) / len(sighits)
                except Exception as e:
                    overall_precision = None

                scores.append({
                    'GWAS Strategy': method,
                    'Heritability': h,
                    'Fraction Shared ES': fs,
                    'Overall Recall': float(len(sighits.intersection(causal_union))) / len(causal_union),
                    'Overall Precision': overall_precision,
                    'Non-shared SNP Recall': nonshared_snp_recall,
                    'Non-shared ES Recall': nonshared_es_recall
                })

    return pd.DataFrame(scores), causal_snp_df


def evaluate_overall_performance():

    for nc in n_causal_snps:
        for ss in snp_shared:

            config_output_dir = os.path.join(output_dir, nc, ss)
            makedir(config_output_dir)

            args = []

            for fs in frac_shared:
                fsf = float(fs.replace("fs_", "").replace("_", "."))

                for h in hsq:
                    hf = float(h.replace("hsq_", "").replace("_", "."))

                    true_causal = [pd.read_csv(os.path.join(sim_input, nc, ss, fs, str(c), "causal_loci.csv"),
                                               sep="\t", engine='python',
                                               header=None, names=['SNP', 'ES'])
                                   for c in clust_ids]

                    args.append([os.path.join(gwas_input, nc, ss, fs, h), true_causal, (fsf, hf)])

            pool = Pool(6)
            res = pool.starmap(score_analysis_methods, args)
            pool.close()
            pool.join()

            score_df = pd.concat([r[0] for r in res])
            score_df.to_csv(os.path.join(config_output_dir, "score_df.csv"))
            power_df = pd.concat([r[1] for r in res])
            power_df.to_csv(os.path.join(config_output_dir, "power_df.csv"))


if __name__ == '__main__':

    snp_meta_df = pd.read_csv(snp_bim_file, header=None,
                              names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'], sep="\s",
                              engine='python')
    snp_meta_df = snp_meta_df.loc[snp_meta_df['CHR'] == 22, ]
    print(len(snp_meta_df))

    clust_df = pd.read_csv(cluster_assignment_file, header=None,
                           names=['FID', 'IID', 'CLUSTER'], sep="\s",
                           engine='python')

    clust_size = clust_df.groupby('CLUSTER').size()
    makedir(output_dir)

    # Grab the simulation configurations from directory structure:
    try:
        n_causal_snps = next(os.walk(sim_input))[1]
        snp_shared = next(os.walk(os.path.join(sim_input, n_causal_snps[0])))[1]
        frac_shared = next(os.walk(os.path.join(sim_input, n_causal_snps[0], snp_shared[0])))[1]
        clust_ids = next(os.walk(os.path.join(sim_input, n_causal_snps[0], snp_shared[0], frac_shared[0])))[1]
        hsq = next(os.walk(os.path.join(sim_input, n_causal_snps[0], snp_shared[0], frac_shared[0], clust_ids[0])))[1]
        gwas_methods = next(os.walk(os.path.join(gwas_input, n_causal_snps[0], snp_shared[0], frac_shared[0], hsq[0])))[1]
    except Exception as e:
        print(str(e))
        sys.exit()

    clust_ids = [c for c in clust_ids if c in clust_size[clust_size > clust_size_threshold].index.astype(str)]

    evaluate_overall_performance()
