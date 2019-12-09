import os
import errno
import sys
import pandas as pd
import numpy as np
import matplotlib.pylab as plt
import seaborn as sns

# --------------------------------------
# Config:
score_names = ['Overall Recall', 'Overall Precision', 'Non-shared SNP Recall', 'Non-shared ES Recall']

# Inputs:
evaluation_dir = "./evaluation/"

# --------------------------------------


def makedir(cdir):
    try:
        os.makedirs(cdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise


def plot_scores(nc, ss):

    nci = int(nc.split("_")[0])
    ssf = float(ss.replace("ss_", "").replace("_", "."))

    input_dir = os.path.join(evaluation_dir, nc, ss)

    makedir(os.path.join(input_dir, "plots"))
    score_df = pd.read_csv(os.path.join(input_dir, "score_df.csv"))

    for sc in score_names:
        g = sns.catplot(x="Heritability", y=sc,
                        hue="GWAS Strategy", col="Fraction Shared ES",
                        data=score_df,
                        kind="violin")
        g.fig.subplots_adjust(top=.9)
        plt.suptitle("%d Causal SNPs, %d%% shared across clusters" % (nci, int(100 * ssf)))
        plt.savefig(os.path.join(input_dir, "plots", sc + ".pdf"))
        plt.close()

    score_df_m = score_df.loc[score_df["Fraction Shared ES"].isin([0.5, 1.0]), ]

    for sc in score_names:
        g = sns.catplot(x="Heritability", y=sc,
                        hue="GWAS Strategy", col="Fraction Shared ES",
                        data=score_df_m,
                        kind="violin")
        g.fig.subplots_adjust(top=.9)
        plt.suptitle("%d Causal SNPs, %d%% shared across clusters" % (nci, int(100 * ssf)))
        plt.savefig(os.path.join(input_dir, "plots", sc + ".png"))
        plt.close()


def plot_power(nc, ss, m_suffix="", filter_0_1=True):

    nci = int(nc.split("_")[0])
    ssf = float(ss.replace("ss_", "").replace("_", "."))

    input_dir = os.path.join(evaluation_dir, nc, ss)

    makedir(os.path.join(input_dir, "plots"))
    power_df = pd.read_csv(os.path.join(input_dir, "power_df.csv"))

    power_df['mega' + m_suffix] = power_df['mega' + m_suffix].clip(upper=1.0, lower=0.0)
    power_df['stratified' + m_suffix] = power_df['stratified' + m_suffix].clip(upper=1.0, lower=0.0)

    if filter_0_1:
        power_df = power_df.loc[~((power_df['mega' + m_suffix] == 0.) & (power_df['stratified' + m_suffix] == 0.)) |
                                ((power_df['mega' + m_suffix] == 1.) & (power_df['stratified' + m_suffix] == 1.)), ]

    power_df = power_df.sort_values('Presence', ascending=False)
    x = np.linspace(0.0, 1., 100)

    power_df.rename(columns={'mega' + m_suffix: 'Mega Analysis Power',
                             'stratified' + m_suffix: 'Stratified Analysis Power'}, inplace=True)

    g = sns.FacetGrid(power_df, col="Heritability", row='Fraction Shared', margin_titles=True)
    g = g.map(sns.scatterplot, "Mega Analysis Power", "Stratified Analysis Power",
              hue='Presence', palette={"Shared": "#95a5a6", "Exclusive": "#9b59b6"},
              style='Effect Size', markers={"Shared": "o", "Exclusive": "X"},
              data=power_df,
              legend='full')

    for ax in g.axes.flat:
        ax.plot(x, x, color='dimgray', linestyle='--')

    g.fig.subplots_adjust(top=0.92, right=0.85)
    if m_suffix == "_LD":
        plt.suptitle("Power (w/ LD): %d Causal SNPs, %d%% shared across clusters" % (nci, int(100 * ssf)))
    else:
        plt.suptitle("Power: %d Causal SNPs, %d%% shared across clusters" % (nci, int(100 * ssf)))
    plt.legend(loc='center left', bbox_to_anchor=(1.15, 0.7))
    plt.savefig(os.path.join(input_dir, "plots", "power" + m_suffix + ".pdf"))
    plt.close()

    power_df_m = power_df.loc[power_df['Fraction Shared'].isin([0.5, 1.0]), ]

    g = sns.FacetGrid(power_df_m, col="Heritability", row='Fraction Shared', margin_titles=True)
    g = g.map(sns.scatterplot, "Mega Analysis Power", "Stratified Analysis Power",
              hue='Presence', palette={"Shared": "#95a5a6", "Exclusive": "#9b59b6"},
              style='Effect Size', markers={"Shared": "o", "Exclusive": "X"},
              data=power_df_m,
              legend='full')

    for ax in g.axes.flat:
        ax.plot(x, x, color='dimgray', linestyle='--')

    g.fig.subplots_adjust(top=0.92, right=0.85)
    if m_suffix == "_LD":
        plt.suptitle("Power (w/ LD): %d Causal SNPs, %d%% shared across clusters" % (nci, int(100 * ssf)))
    else:
        plt.suptitle("Power: %d Causal SNPs, %d%% shared across clusters" % (nci, int(100 * ssf)))
    plt.legend(loc='center left', bbox_to_anchor=(1.15, 0.7))
    plt.savefig(os.path.join(input_dir, "plots", "power" + m_suffix + ".png"))
    plt.close()

    """
    for h in np.unique(power_df['Heritability']):
        for fs in np.unique(power_df['Fraction Shared']):
            df = power_df.loc[(power_df['Heritability'] == h) & (power_df['Fraction Shared'] == fs), ]
            sns.scatterplot("mega", "stratified", hue="Presence", style="Effect Size",
                            palette={"Shared": "#95a5a6", "Exclusive": "#9b59b6"},
                            markers={"Shared": "o", "Exclusive": "X"},
                            data=df)
            plt.plot(x, x, color='black', linestyle='--')

            plt.xlabel("Mega Analysis Power")
            plt.ylabel("Stratified Analysis Power")

            plt.savefig(os.path.join(input_dir, "plots", "h_%.1f_fs_%.2f_power_filter.pdf" % (h, fs)))
            plt.close()
    """


if __name__ == '__main__':
    try:
        n_causal_snps = next(os.walk(evaluation_dir))[1]
        snp_shared = next(os.walk(os.path.join(evaluation_dir, n_causal_snps[0])))[1]
    except Exception as e:
        print(str(e))
        sys.exit()

    for nc in n_causal_snps:
        for ss in snp_shared:
            plot_scores(nc, ss)
            plot_power(nc, ss)
            plot_power(nc, ss, m_suffix="_LD")
