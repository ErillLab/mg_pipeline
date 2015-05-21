# -*- coding: utf-8 -*-
"""
Created on Tue May 19 13:49:39 2015

@author: Talmo
"""

import sys
sys.path.append("../..")
from igc_pipeline import *
from tqdm import tqdm
import pandas as pd
import numpy as np

probs_file = "Gamma_probs_0-1000_tax-filtered_head-filtered_annotated.csv"
probs = pd.read_csv(probs_file)
samples = get_all_samples()


sample = samples[0]

genes = get_genes(sample)
operons = get_operons(sample)
genes2operons = get_genes2operon(genes, operons)
ORFs = get_ORFs(sample)

COG = probs.eggNOG.iloc[2]
COG_genes = genes[genes.eggNOG == COG]
COG_operons_id = genes2operons[COG_genes.index]
COG_operons = operons.loc[COG_operons_id]
COG_operon_genes = genes.loc[np.hstack(COG_operons.genes)]