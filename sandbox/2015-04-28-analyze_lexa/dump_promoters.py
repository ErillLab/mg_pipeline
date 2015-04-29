# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 17:46:07 2015

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
reca_cog = "COG0468"

samples = get_all_samples(HMP=False).index.tolist()[0:100]

#%%
def get_promoters(sample, cog=lexa_cog):
    # Load data
    genes = load_sample_genes(sample)
    operons = get_operons(sample)
    genes2operon = get_genes2operon(genes, operons)
    
    # Filter by head completeness
    genes = genes.loc[np.unique(np.hstack(operons.genes[operons.head_completeness != "Lack 5'-end"]))].dropna(how="all")
    
    # Filter by taxonomy
    genes = genes[genes.phylum == "Firmicutes"]
    #genes = genes[genes["class"] == "Gammaproteobacteria"]
    
    # Get COG genes
    cog_genes = genes[genes["eggNOG"] == cog]
    cog_gene_names = cog_genes.index.tolist()
    
    # Get COG operons
    cog_operons = operons.loc[genes2operon[cog_gene_names]]
    
    return cog_operons.promoter_seq
    
#%%
cog = reca_cog

promoter_ids = []
promoters_fasta = ""
for sample in tqdm(samples):
    promoter_seqs = get_promoters(sample, cog)
    for operon, seq in promoter_seqs.iteritems():
        promoter_id = "%s|%d" % (sample, operon)
        if promoter_id not in promoter_ids:        
            promoters_fasta += ">%s\n%s\n" % (promoter_id, seq)
            promoter_ids.append(promoter_id)
open("RecA_promoters-Firmicutes-0-100.fa", "w").write(promoters_fasta)