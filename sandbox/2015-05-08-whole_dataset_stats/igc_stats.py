# -*- coding: utf-8 -*-
"""
Created on Fri May  8 12:02:51 2015

@author: Talmo
"""

import sys
import os
sys.path.append("../..")
from igc_pipeline import *
from tqdm import tqdm
from glob import glob
import numpy as np
import pandas as pd
import seaborn as sns

samples = [sample for sample in get_all_samples(HMP=False) if has_data(sample).all()]

stats = pd.DataFrame(index=samples,
                     columns=["ORFs", "genes", "genes+has_cog", "operons", "operons+head_complete", 
                     "operons+full_promoter", "operons+head_complete+full_promoter", 
                     "ORFs+head_complete+full_promoter", "genes+head_complete+full_promoter",
                     "genes+has_cog+head_complete+full_promoter"])

for sample in tqdm(samples):
    
    ORFs = get_operons(sample)
    genes = get_genes(sample)
    operons = get_operons(sample)
    
    stats.at[sample, "ORFs"] = len(ORFs)
    stats.at[sample, "genes"] = len(genes)
    stats.at[sample, "genes+has_cog"] = sum(genes.eggNOG.notnull())
    stats.at[sample, "operons"] = len(operons)
    stats.at[sample, "operons+head_complete"] = sum(operons.head_completeness != "Lack 5'-end")
    stats.at[sample, "operons+full_promoter"] = sum(operons.promoter_seq.apply(len) >= 300)
    stats.at[sample, "operons+head_complete+full_promoter"] = sum((operons.head_completeness != "Lack 5'-end") & (operons.promoter_seq.apply(len) >= 300))
    stats.at[sample, "ORFs+head_complete+full_promoter"] = np.unique(np.hstack(operons[(operons.head_completeness != "Lack 5'-end") & (operons.promoter_seq.apply(len) >= 300)].genes)).shape[0]
    stats.at[sample, "genes+head_complete+full_promoter"] = sum(genes.loc[np.unique(np.hstack(operons[(operons.head_completeness != "Lack 5'-end") & (operons.promoter_seq.apply(len) >= 300)].genes))].gene_id.notnull())
    stats.at[sample, "genes+has_cog+head_complete+full_promoter"] = sum(genes.loc[np.unique(np.hstack(operons[(operons.head_completeness != "Lack 5'-end") & (operons.promoter_seq.apply(len) >= 300)].genes))].eggNOG.notnull())
    
stats.to_csv("2015-05-08-igc_stats.csv")