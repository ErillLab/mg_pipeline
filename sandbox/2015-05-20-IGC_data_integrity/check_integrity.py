# -*- coding: utf-8 -*-
"""
Created on Wed May 20 16:19:51 2015

@author: Talmo
"""

import sys
sys.path.append("../..")
from igc_pipeline import *
import pandas as pd
import os

#%% Configuration
samples = get_all_samples()
sample = samples[0]

#%% Load data
scaftigs = get_scaftigs(sample, index_only=False, convert_to_str=False)
ORFs = get_ORFs(sample)
ORF_seqs = get_ORF_seqs(sample, index_only=False, convert_to_str=False)
genes = get_genes(sample)
operons = get_operons(sample)
genes2operon = get_genes2operon(genes, operons)

# Raw string versions
scaftigs_str = {scaf: str(seq.seq).upper() for scaf, seq in scaftigs.items()}
ORF_seqs_str = {orf: str(seq.seq).upper() for orf, seq in ORF_seqs.items()}

#%% Stats
print "sample:", sample
print "scaftigs:", len(scaftigs)
print "ORFs:", len(ORFs)
print "ORF_seqs:", len(ORF_seqs)
print "genes:", len(genes)
print "operons:", len(operons)
print

#%% ORFs in ORF table can be found in ORF sequences and vice-versa
# 1. ORFs == ORF_seqs
orfs_seqs = set(ORF_seqs.keys())
orfs_table = set(ORFs.index.tolist())
orfs_all = orfs_seqs.union(orfs_table)
orfs_both = orfs_seqs.intersection(orfs_table)
orfs_seqs_only = orfs_seqs.difference(orfs_table)
orfs_table_only = orfs_table.difference(orfs_seqs)

print "ORFs ∪ ORF_seqs:", len(orfs_all)
print "ORFs ∩ ORF_seqs:", len(orfs_both)
print "ORFs - ORF_seqs:", len(orfs_table_only)
print "ORF_seqs - ORFs:", len(orfs_seqs_only)
print

# 2. ORFs in scaftigs
# All ORF sequences can be found in scaftigs
#   - Get all scaftigs referenced by ORFs and check that they are all in scaftigs
#   - Check if actual sequences can be found in scaftigs
#   - (Optional) Find scaftigs without ORFs?

# 3. genes in ORFs
# All genes are found in ORFs
#   - Get all ORFs referenced by genes and check that they are all in ORFs

# 4. genes in operons.genes

# 5. operons.genes in ORFs

# 6. operons.promoter_seq in scaftigs

