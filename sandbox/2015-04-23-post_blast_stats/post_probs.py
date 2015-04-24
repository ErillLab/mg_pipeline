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

#%% Find BLASTed
blast_out_path = IGC_path + "Taxonomy_BLAST/has_cog-eggNOG/"
samples = get_all_samples()
blasted = [p.replace(blast_out_path, "").replace(".tbl", "") for p in glob(blast_out_path + "*.tbl")]
not_blasted = [sample for sample in samples if sample not in blasted]

