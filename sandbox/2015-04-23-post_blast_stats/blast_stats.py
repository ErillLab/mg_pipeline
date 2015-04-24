# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 21:57:52 2015

@author: Talmo
"""
import sys
import os
sys.path.append("../..")
from igc_pipeline import *
from tqdm import tqdm
from glob import glob

#%% Find BLASTed
blast_out_path = IGC_path + "Taxonomy_BLAST/has_cog-eggNOG/"
samples = get_all_samples()
blasted = [p.replace(blast_out_path, "").replace(".tbl", "") for p in glob(blast_out_path + "*.tbl") if os.stat(p).st_size > 0]
not_blasted = [sample for sample in samples if sample not in blasted]

#%% Aggregate stats
samples = blasted
num_samples = len(samples)

COGs = pd.Series()
COG_stats = pd.Series(index=["updated", "no_change", "dropped", "missing", "old_total", "total"]).fillna(0)
#no_change_cogs = 0; updated_cogs = 0; missing_cogs = 0; total_cogs = 0

levels = ["phylum", "class", "order", "genus"]
taxonomy = {level: pd.Series() for level in levels}
tax_stats = pd.DataFrame(index=levels, columns=["new", "updated", "no_change", "dropped", "missing", "old_total", "total"]).fillna(0)

for sample in tqdm(samples):
    genes = load_sample_genes(sample)
    
    # COGs
    COGs = COGs.add(genes.groupby("eggNOG").size(), fill_value=0)
    idx = genes.eggNOG_old != "unknown"
    COG_stats.no_change += sum(genes.eggNOG_old[idx] == genes.eggNOG[idx])
    COG_stats.updated += sum((genes.eggNOG[idx] != genes.eggNOG_old[idx]) & ~genes.eggNOG[idx].isnull())
    COG_stats.dropped += sum(genes.eggNOG.isnull() & idx)
    COG_stats.missing += sum(genes.eggNOG.isnull())
    COG_stats.total += sum(~genes.eggNOG.isnull()) # total = updated + no_change
    COG_stats.old_total += sum(idx) # old_total = updated + no_change + dropped
    
    # Taxonomy
    for level in levels:
        taxonomy[level] = taxonomy[level].add(genes.groupby(level).size(), fill_value=0)
        tax_stats["missing"][level] += sum(genes[level].isnull())
        tax_stats["total"][level] += sum(~genes[level].isnull())
        if level + "_old" in genes:
            idx = genes[level + "_old"] != "unknown"
            tax_stats["new"][level] += sum(~idx & ~genes[level].isnull())
            tax_stats["updated"][level] += sum(idx & ~genes[level].isnull())
            tax_stats["dropped"][level] += sum(idx & genes[level].isnull())
            tax_stats["no_change"][level] += sum(genes[level][idx] == genes[level + "_old"][idx])
            tax_stats["old_total"][level] += sum(idx)
        else:
            tax_stats["new"][level] += sum(~genes[level].isnull())

#%% Save
COGs.to_csv("2014-04-23-COGs.csv")
for level in levels:
    taxonomy[level].to_csv("2014-04-23-Taxonomy_%s.csv" % level)

#%% Analyze
print "Samples: %d" % num_samples
print

print "COGs:"
print COG_stats.to_string()
LexA_COG = "COG1974"
print "LexA [%s]: %d" % (LexA_COG, COGs[LexA_COG])
print

print "Taxonomy:"
print tax_stats.T.to_string()
