# -*- coding: utf-8 -*-
"""
This script queries the Integrated Gene Catalog for some basic stats.

Created on Wed Feb 18 14:21:39 2015

@author: Talmo
"""

import pandas as pd

#%% Load data
p = "E:\\metagenomics\\Li et al. (2014)\\3.IGC.AnnotationInfo\\IGC.annotation.summary.v2.head"
t = pd.read_table(p, names=['gene_id', 'gene_name', 'gene_length', 'completeness', 'cohort_origin', 'phylum', 'genus', 'kegg', 'eggNOG', 'sample_freq', 'individual_freq', 'eggNOG_funccat', 'kegg_funccat', 'cohort_assembled'])

lexA = t[t.eggNOG == 'COG1974']
recA = t[t.eggNOG == 'COG0468']

#%% Get stats for lexA
nog = lexA
stats = pd.Series()

stats['n_total'] = nog.shape[0]
stats['n_complete'] = sum(nog['completeness'] == 'Complete')
stats['n_has_phylum'] = sum(nog['phylum'] != 'unknown')
print stats
print nog.phylum.value_counts()
lexA.to_csv('IGC.lexA.csv')

#%% Get stats for recA
nog = recA
stats = pd.Series()

stats['n_total'] = nog.shape[0]
stats['n_complete'] = sum(nog['completeness'] == 'Complete')
stats['n_has_phylum'] = sum(nog['phylum'] != 'unknown')
print stats
print nog.phylum.value_counts()
lexA.to_csv('IGC.recA.csv')
