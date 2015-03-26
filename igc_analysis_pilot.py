# -*- coding: utf-8 -*-
"""
Preliminary analysis routines

Created on Wed Mar 25 16:22:26 2015

@author: Talmo
"""

from igc_pipeline import *
import seaborn as sns

#%%
sample = "MH0001"
PSSM = Firmicutes_LexA
sample = get_unique_sample(sample)
    
# Load sample data
genes = load_sample_genes(sample)
operons = predict_operons(sample)
scores = get_sample_scores(sample, PSSM)

# Filter by genes with taxonomy and COGs
valid_genes = genes[(genes.phylum != "unknown") & (genes.eggNOG != "unknown")]
valid_operons = operons[operons.genes.apply(lambda x: len(valid_genes.index.intersection(x)) > 0)]
valid_scores = scores.loc[valid_operons.index]

# Compute soft max
soft_max = valid_scores.applymap(np.exp).sum(1).map(np.log)
soft_max1d = np.hstack(soft_max)

#%% Break down by taxonomy
tax = valid_genes.groupby(['phylum', 'genus'])
