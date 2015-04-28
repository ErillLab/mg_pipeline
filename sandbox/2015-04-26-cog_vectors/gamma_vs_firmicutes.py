# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 22:02:57 2015

@author: Talmo
"""
import sys
sys.path.append("../..")
from igc_pipeline import *
from tqdm import tqdm
from glob import glob
import pandas as pd
import numpy as np

#%% Configuration
#samples = get_all_samples(HMP=False).index.tolist()
samples = get_all_samples(HMP=False).index.tolist()[0:100]
#samples = get_MetaHit()[0:85]

# From Cornish et al. (2014)
putative_cogs = ["COG0389", "COG0468", "COG0556", "COG1974", "COG0210", "COG2001", "COG0732", "COG0178", "COG0187", "COG4974", "COG0632", "COG0497", "COG0399", "COG0653", "COG1136", "COG0463", "COG1609", "COG0745", "COG0582"]
SOS_cogs = ["COG0174", "COG0177", "COG0178", "COG0187", "COG0210", "COG0214", "COG0270", "COG0305", "COG0322", "COG0330", "COG0380", "COG0389", "COG0417", "COG0420", "COG0463", "COG0468", "COG0477", "COG0491", "COG0497", "COG0531", "COG0550", "COG0556", "COG0582", "COG0587", "COG0591", "COG0596", "COG0606", "COG0629", "COG0632", "COG0642", "COG0749", "COG0776", "COG0784", "COG0790", "COG0817", "COG0850", "COG0863", "COG1028", "COG1066", "COG1125", "COG1132", "COG1199", "COG1200", "COG1201", "COG1219", "COG1253", "COG1322", "COG1349", "COG1372", "COG1388", "COG1404", "COG1414", "COG1443", "COG1479", "COG1533", "COG1573", "COG1585", "COG1609", "COG1653", "COG1674", "COG1877", "COG1957", "COG1961", "COG1971", "COG1974", "COG1982", "COG1988", "COG2132", "COG2137", "COG2200", "COG2255", "COG2318", "COG2372", "COG2764", "COG2818", "COG2827", "COG3066", "COG3141", "COG3226", "COG3279", "COG3324", "COG3449", "COG3600", "COG3636", "COG3657", "COG3668", "COG3798", "COG3857", "COG3905", "COG4188", "COG4277", "COG4335", "COG4535", "COG4544", "COG4799", "COG4948", "COG5321", "COG5404", "COG5615"]

#%% Compute
def get_cog_probs(samples, pssm):
    LL_ratios = pd.Series()
    n = pd.Series()
    cog_probs = pd.DataFrame(columns=["n", "post_prob"])
    for sample in tqdm(samples):
        try:
            genes = load_sample_genes(sample)
            if "class" not in genes:
                continue
            operons = get_operons(sample)
            scores = get_sample_scores(sample, pssm)
            genes2operon = get_genes2operon(genes, operons)
        except:
            continue
        
        # Filter by head completeness
        genes = genes.loc[np.unique(np.hstack(operons.genes[operons.head_completeness != "Lack 5'-end"]))].dropna(how="all")
        genes.index.name = "gene_name"
        
        # Filter by taxonomy
        #genes = genes[genes.phylum == "Firmicutes"]
        genes = genes[genes["class"] == "Gammaproteobacteria"]
        
        # Use old assignments
        genes["eggNOG"] = genes["eggNOG_old"]
        assert (genes["eggNOG"] == genes["eggNOG_old"]).all()
        
        # Group genes by COG
        grouped = genes.reset_index().groupby("eggNOG")
        
        # Compute log likelihoods
        get_LL_ratio = lambda gene_names: pssm.LL_ratio(np.hstack(scores.loc[genes2operon[gene_names]]))
        sample_ratios = grouped["gene_name"].agg(get_LL_ratio)
        
        # Aggregate
        idx = LL_ratios.index.union(sample_ratios.index)
        LL_ratios = LL_ratios.get(idx).fillna(1) * sample_ratios.get(idx).fillna(1)
        n = n.get(idx).fillna(0)
        n.loc[sample_ratios.index] += 1
    
    # Compute posteriors
    cog_probs.post_prob = LL_ratios.map(lambda LLR: 1 / (1 + LLR * pssm.pb / pssm.pf))
    cog_probs.sort("post_prob", ascending=False, inplace=True)
    cog_probs.n = n
    return cog_probs

def plot_cogs(cogs):
    k = len(cogs)
    plt.figure()
    plt.bar(range(k), cogs)
    plt.xticks(np.arange(k) + 0.5, cogs.index, rotation="vertical")
    plt.ylim(0, 1)
    plt.subplots_adjust(bottom=0.20)
    plt.ylabel("Posterior probability of being regulated")

#%% Compare
Firmicutes_probs = get_cog_probs(samples, Firmicutes_LexA)
Firmicutes_probs.to_csv("Firmicutes_probs_0-100_tax-filtered_head-filtered_old.csv")

#Gamma_probs = get_cog_probs(samples[0:100], GammaProteobacteria_LexA)
#Gamma_probs.to_csv("Gamma_probs_0-100_tax-filtered_head-filtered.csv")

#%% Merge
#cog_probs = Firmicutes_probs.copy().rename(columns={"post_prob": "post_prob_Firmicutes"})
#cog_probs["post_prob_Gamma"] = Gamma_probs.post_prob
#cog_probs.to_csv("LexA-Firmicutes_vs_Gammaproteobacteria.csv")

#%% Load and merge
#cog_probs = pd.merge(pd.read_csv("Firmicutes_probs_0-100.csv", index_col="eggNOG"), pd.read_csv("Gamma_probs_0-100.csv", index_col="eggNOG", usecols=["eggNOG", "post_prob"]), left_index=True, right_index=True, suffixes=("_Firmicutes", "_Gamma"))
#cog_probs.to_csv("LexA-Firmicutes_vs_Gammaproteobacteria_0-100.csv")

#%% Plot
plot_figs = False
if plot_figs:
    # Top Firmicutes
    plot_cogs(Firmicutes_probs.post_prob.sort(ascending=False, inplace=False).head(25))
    plt.title(Firmicutes_LexA)
    
    # Top Gamma
    plot_cogs(cog_probs.post_prob_Gamma.sort(ascending=False, inplace=False).head(25))
    plt.title(GammaProteobacteria_LexA)
    
    # Putative Firmicutes
    plot_cogs(cog_probs.post_prob_Firmicutes[putative_cogs])
    
    # Putative Gamma
    plot_cogs(cog_probs.post_prob_Gamma[putative_cogs])
    
    # SOS Firmicutes
    plot_cogs(cog_probs.post_prob_Firmicutes[SOS_cogs])
    
    # SOS Gamma
    plot_cogs(cog_probs.post_prob_Gamma[SOS_cogs])