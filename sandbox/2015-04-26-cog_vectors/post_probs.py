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

#%% Configuration
pssm = Firmicutes_LexA
samples = get_all_with_scores(pssm)

# From Cornish et al. (2014)
putative_cogs = ["COG0389", "COG0468", "COG0556", "COG1974", "COG0210", "COG2001", "COG0732", "COG0178", "COG0187", "COG4974", "COG0632", "COG0497", "COG0399", "COG0653", "COG1136", "COG0463", "COG1609", "COG0745", "COG0582"]
SOS_cogs = ["COG0174", "COG0177", "COG0178", "COG0187", "COG0210", "COG0214", "COG0270", "COG0305", "COG0322", "COG0330", "COG0380", "COG0389", "COG0417", "COG0420", "COG0463", "COG0468", "COG0477", "COG0491", "COG0497", "COG0531", "COG0550", "COG0556", "COG0582", "COG0587", "COG0591", "COG0596", "COG0606", "COG0629", "COG0632", "COG0642", "COG0749", "COG0776", "COG0784", "COG0790", "COG0817", "COG0850", "COG0863", "COG1028", "COG1066", "COG1125", "COG1132", "COG1199", "COG1200", "COG1201", "COG1219", "COG1253", "COG1322", "COG1349", "COG1372", "COG1388", "COG1404", "COG1414", "COG1443", "COG1479", "COG1533", "COG1573", "COG1585", "COG1609", "COG1653", "COG1674", "COG1877", "COG1957", "COG1961", "COG1971", "COG1974", "COG1982", "COG1988", "COG2132", "COG2137", "COG2200", "COG2255", "COG2318", "COG2372", "COG2764", "COG2818", "COG2827", "COG3066", "COG3141", "COG3226", "COG3279", "COG3324", "COG3449", "COG3600", "COG3636", "COG3657", "COG3668", "COG3798", "COG3857", "COG3905", "COG4188", "COG4277", "COG4335", "COG4535", "COG4544", "COG4799", "COG4948", "COG5321", "COG5404", "COG5615"]

#%% One sample
sample = samples[105]

# Load data
genes = load_sample_genes(sample)
operons = get_operons(sample)
scores = get_sample_scores(sample, pssm)
genes2operon = get_genes2operon(genes, operons)

# Group genes by COG
grouped = genes.reset_index().groupby("eggNOG")

# Compute posterior for each COG group
get_prob = lambda gene_names: pssm.post_prob(np.hstack(scores.loc[genes2operon[gene_names]]))
cog_probs = grouped["gene_name"].agg({"n": len, "post_prob": get_prob})
cog_probs.sort("post_prob", ascending=False, inplace=True)

#%% Plot top COGs
def plot_cogs(cogs):
    k = len(cogs)
    plt.figure()
    plt.bar(range(k), cogs.post_prob)
    plt.xticks(np.arange(k) + 0.5, cogs.index, rotation="vertical")
    plt.subplots_adjust(bottom=0.20)
    plt.ylabel("Posterior probability of being regulated")

plot_cogs(cog_probs.head(20))
plot_cogs(cog_probs.loc[SOS_cogs].sort("post_prob", ascending=False).dropna())
plot_cogs(cog_probs.loc[putative_cogs])

#%% Combinations
score_sets = [[10, 11, 12], [13, 14, 15]]

# Concatenated scores
conc_prob = pssm.post_prob(np.hstack(score_sets))

# Scores computed separately by LL ratio
sep_prob = 1 / (1 + np.prod(map(pssm.LL_ratio, score_sets)) * pssm.pb / pssm.pf)

assert(conc_prob == sep_prob)

#%% Multiple samples
LL_ratios = pd.Series()
for sample in tqdm(samples[0:400]):
    genes = load_sample_genes(sample)
    if "eggNOG_old" not in genes:
        continue
    operons = get_operons(sample)
    scores = get_sample_scores(sample, pssm)
    genes2operon = get_genes2operon(genes, operons)

    #TODO: Filter by head completeness
    #scores.get(operons.head_completeness != "Lack 5'-end").dropna()
    # Group genes by COG
    grouped = genes.reset_index().groupby("eggNOG")
    
    # Compute log likelihoods
    get_LL_ratio = lambda gene_names: pssm.LL_ratio(np.hstack(scores.loc[genes2operon[gene_names]]))
    sample_ratios = grouped["gene_name"].agg(get_LL_ratio)
    
    # Aggregate
    idx = LL_ratios.index.union(sample_ratios.index)
    LL_ratios = LL_ratios.get(idx).fillna(1) * sample_ratios.get(idx).fillna(1)

# Compute posteriors
cog_probs = LL_ratios.map(lambda LLR: 1 / (1 + LLR * pssm.pb / pssm.pf))
cog_probs.sort(ascending=False)

#%% Plot top COGs
def plot_cogs2(cogs):
    k = len(cogs)
    plt.figure()
    plt.bar(range(k), cogs)
    plt.xticks(np.arange(k) + 0.5, cogs.index, rotation="vertical")
    plt.subplots_adjust(bottom=0.20)
    plt.ylabel("Posterior probability of being regulated")

plot_cogs2(cog_probs.head(20))
plot_cogs2(cog_probs.loc[SOS_cogs].dropna().sort(ascending=False, inplace=False))
plot_cogs2(cog_probs.loc[putative_cogs])