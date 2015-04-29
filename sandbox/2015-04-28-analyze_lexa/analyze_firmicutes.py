# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 14:32:51 2015

@author: Talmo
"""

import sys
sys.path.append("../..")
from igc_pipeline import *
from tqdm import tqdm
from glob import glob
import pandas as pd
import numpy as np

# From Cornish et al. (2014)
putative_cogs = ["COG0389", "COG0468", "COG0556", "COG1974", "COG0210", "COG2001", "COG0732", "COG0178", "COG0187", "COG4974", "COG0632", "COG0497", "COG0399", "COG0653", "COG1136", "COG0463", "COG1609", "COG0745", "COG0582"]
SOS_cogs = ["COG0174", "COG0177", "COG0178", "COG0187", "COG0210", "COG0214", "COG0270", "COG0305", "COG0322", "COG0330", "COG0380", "COG0389", "COG0417", "COG0420", "COG0463", "COG0468", "COG0477", "COG0491", "COG0497", "COG0531", "COG0550", "COG0556", "COG0582", "COG0587", "COG0591", "COG0596", "COG0606", "COG0629", "COG0632", "COG0642", "COG0749", "COG0776", "COG0784", "COG0790", "COG0817", "COG0850", "COG0863", "COG1028", "COG1066", "COG1125", "COG1132", "COG1199", "COG1200", "COG1201", "COG1219", "COG1253", "COG1322", "COG1349", "COG1372", "COG1388", "COG1404", "COG1414", "COG1443", "COG1479", "COG1533", "COG1573", "COG1585", "COG1609", "COG1653", "COG1674", "COG1877", "COG1957", "COG1961", "COG1971", "COG1974", "COG1982", "COG1988", "COG2132", "COG2137", "COG2200", "COG2255", "COG2318", "COG2372", "COG2764", "COG2818", "COG2827", "COG3066", "COG3141", "COG3226", "COG3279", "COG3324", "COG3449", "COG3600", "COG3636", "COG3657", "COG3668", "COG3798", "COG3857", "COG3905", "COG4188", "COG4277", "COG4335", "COG4535", "COG4544", "COG4799", "COG4948", "COG5321", "COG5404", "COG5615"]

# Configuration
pssm = Firmicutes_LexA
lexa_cog = "COG1974"
use_old = False

#%%
sample = "MH0011"
print "Sample:", sample

# Load data
ORFs = get_ORFs(sample)
ORF_seqs = SeqIO.index(get_sample_paths(sample)["ORFs_fasta"], "fasta")
genes = load_sample_genes(sample)
operons = get_operons(sample)
genes2operon = get_genes2operon(genes, operons)
scores = get_sample_scores(sample, pssm)
probs = pd.read_csv("Firmicutes_probs_0-1000_tax-filtered_head-filtered.csv", index_col="eggNOG")
LL_ratios = pd.Series()

all_genes = genes.copy()

# Use old assignments
if use_old:
    genes["eggNOG"] = genes["eggNOG_old"]
    assert (genes["eggNOG"] == genes["eggNOG_old"]).all()

print "Total LexA genes:", sum(genes.eggNOG == lexa_cog)

# Filter by head completeness
genes = genes.loc[np.unique(np.hstack(operons.genes[operons.head_completeness != "Lack 5'-end"]))].dropna(how="all")
genes.index.name = "gene_name"
print "LexA genes (head complete):", sum(genes.eggNOG == lexa_cog)

# Filter by taxonomy
genes = genes[(genes.phylum == "Firmicutes") | (genes.phylum == "Actinobacteria")]
#genes = genes[genes["class"] == "Gammaproteobacteria"]
print "LexA genes (head complete + Firmicutes + Actino):", sum(genes.eggNOG == lexa_cog)

# Group genes by COG
grouped = genes.reset_index().groupby("eggNOG")

# Compute log-likelihood ratios and posteriors
get_LL_ratio = lambda gene_names: pssm.LL_ratio(np.hstack(scores.loc[genes2operon[gene_names]]))
sample_ratios = grouped["gene_name"].agg(get_LL_ratio)
idx = LL_ratios.index.union(sample_ratios.index)
LL_ratios = LL_ratios.get(idx).fillna(1) * sample_ratios.get(idx).fillna(1)
cog_probs = LL_ratios.map(lambda LLR: 1 / (1 + LLR * pssm.pb / pssm.pf))

lexa_genes = grouped.get_group(lexa_cog)
lexa_gene_names = lexa_genes["gene_name"]
lexa_orfs = ORFs.loc[lexa_gene_names]
lexa_operons = operons.loc[genes2operon[lexa_gene_names]]
lexa_scores = scores.loc[genes2operon[lexa_gene_names]]
lexa_scores_flat = np.hstack(lexa_scores)
lexa_ratios = sample_ratios.loc[lexa_cog]
lexa_prob = cog_probs.loc[lexa_cog]

#print lexa_operons.promoter_seq.map(len)
print "LexA operons:", len(lexa_operons)
#print
print "PSSM:", pssm
print "LexA scores (sample): mean = %.2f, std = %.2f, max = %.2f" % (np.mean(lexa_scores_flat), np.std(lexa_scores_flat), np.max(lexa_scores_flat))
print "LexA LL (sample):", lexa_ratios
print "LexA prob (sample):", lexa_prob
print "LexA prob (all):", probs.post_prob.loc[lexa_cog]

#%% Visualize
# Scores
def plot_scores(scores):
    plt.figure()
    for x in scores:
        plt.plot(x)
    plt.ylim(ymin=0)
    plt.legend(scores.index)
    plt.title("Scores for promoters of operons with LexA in sample %s" % sample)
    plt.ylabel("Score (soft-max)"); plt.xlabel("Promoter position")

def plot_operons(ops):
    from matplotlib import collections
    fig, axes = plt.subplots(len(ops))
    #plt.title("Operons with LexA in sample %s" % sample)
    for i in xrange(len(ops)):
        operon = ops.iloc[i]
        ax = axes[i]
        starts = ORFs.loc[operon.genes].start
        ends = ORFs.loc[operon.genes].end
        
        y = 0
        pts = zip([(x, y) for x in starts], [(x, y) for x in ends])
    
        c = sns.color_palette()
        lc = collections.LineCollection(pts, colors=c, linewidths=20)
        
        ax.add_collection(lc)
        ax.autoscale()
        ax.yaxis.set_major_locator(plt.NullLocator())
        
    plt.xlabel("Position in scaftig (bp)")
    
    
plot_scores(lexa_scores)
#plot_operons(lexa_operons)

#%% Pull out sequences
#lexa_seqs = {gene_name: str(ORF_seqs[gene_name].seq) for gene_name in lexa_gene_names}
#fasta_file = open("LexA_seqs-Firmicutes-%s.fa" % sample, "w")
#for gene_name, seq in lexa_seqs.items():
#    print >>fasta_file, ">" + gene_name
#    print >>fasta_file, seq
#fasta_file.close()