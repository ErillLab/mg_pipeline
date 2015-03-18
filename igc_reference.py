# -*- coding: utf-8 -*-
"""
This script queries the Integrated Gene Catalog for some basic stats.

Created on Wed Feb 18 14:21:39 2015

@author: Talmo
"""
import sys
from time import time, clock
import pandas as pd
import numpy as np
from Bio import SeqIO
from PSSMScorer import PSSMScorer

#base_path = "E:\\metagenomics\\"
base_path = "..\\" # metagenomics folder
IGC_path = base_path + "Li et al. (2014)\\"

#%% Load data
p = IGC_path + "3.IGC.AnnotationInfo\\IGC.annotation.summary.v2"
t = pd.read_table(p, names=['gene_id', 'gene_name', 'gene_length', 'completeness', 'cohort_origin', 'phylum', 'genus', 'kegg', 'eggNOG', 'sample_freq', 'individual_freq', 'eggNOG_funccat', 'kegg_funccat', 'cohort_assembled'], index_col='gene_id')

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
nog.to_csv('IGC.lexA.csv')

#%% Get stats for recA
nog = recA
stats = pd.Series()

stats['n_total'] = nog.shape[0]
stats['n_complete'] = sum(nog['completeness'] == 'Complete')
stats['n_has_phylum'] = sum(nog['phylum'] != 'unknown')
print stats
print nog.phylum.value_counts()
nog.to_csv('IGC.recA.csv')

#%% Get scaftig
gene_name = lexA.iloc[0].gene_name
sample, gene = gene_name.split('_')
sample_path = IGC_path + "4.IndividualAssmeblies\\V1.CD44-0.scaftig.more500.corrected.break.fa"
orf_path = IGC_path + "5.IndividualORFs\\V1.CD44-0.scaftig.more500.corrected.break.fa.more100.fa"

# Load sample scaffolds/scafftigs
scaffolds = SeqIO.index(sample_path, 'fasta')

# Load FASTA headers for ORFs
with open(orf_path) as f:
    lines = [line[1:-2] for line in f if line[0] == '>']
ORF_pattern = '(?P<gene_name>[^ ]+)[ ]+\\[(?P<orf_type>[^\]]+)\\][ ]+locus=(?P<scaffold>[^:]+):(?P<start>[^:]+):(?P<end>[^:]+):(?P<strand>[^:\[]+)\\[(?P<completeness>[^\[\]]+)'
ORFs = pd.Series(lines).str.extract(ORF_pattern).convert_objects(convert_numeric=True)

#%% Operon prediction
threshold_IGI = 50 # max intergenic interval (bp)
promoter_region = (-250, +50) # bp rel to gene start
operon_tables = []
#operons = pd.DataFrame()

# Group by 'scaffold' column and do groupwise operations
i = 0; t = clock()
for scaffold_name, scaf_ORFs in ORFs.groupby('scaffold'):
#    if scaf_ORFs.shape[0] > 10:
        #break
    
    # Sort by start values
    
    # Compute intergenic intervals
    IGIs = scaf_ORFs.start[1:].values - scaf_ORFs.end[:-1].values
    
    # Find operon gene indices
    head_genes = np.where(IGIs > threshold_IGI)[0] + 1
    operons_start = np.hstack(([0], head_genes))
    operons_end = np.hstack((head_genes, [scaf_ORFs.shape[0]]))
    head_gene_start = scaf_ORFs.start.iloc[operons_start]
    
    # Get promoter regions
    promoter_start = (head_gene_start + promoter_region[0]).clip(1).reset_index(drop=True)
    promoter_end = (head_gene_start + promoter_region[1]).reset_index(drop=True)
    promoter_seq = [scaffolds[scaffold_name][s-1:e].seq.tostring() for s, e in zip(promoter_start, promoter_end)]
    
    # Build operon table
    scaffold_operons = pd.DataFrame()
    scaffold_operons['genes'] = [scaf_ORFs.index[op_start:op_end].tolist() for op_start, op_end in zip(operons_start, operons_end)]
    scaffold_operons['promoter_start'] = promoter_start
    scaffold_operons['promoter_end'] = promoter_end
    scaffold_operons['promoter_seq'] = promoter_seq
    
    
    # Append to list of operon tables
    operon_tables.append(scaffold_operons)
    
    i = i + 1
    if i % 1000 == 0:
        print i
    
# Merge operon tables
operons = pd.concat(operon_tables, ignore_index=True)
operons.index.name = 'operon'
print "%.2fs" % (clock() - t)
    


#%% Extract sequence data for specific gene
# Get gene information
gene_orf = ORFs[ORFs['gene_name'] == gene_name]

# Get scaffold
scaffold_name = gene_orf.scaffold.iat[0]
scaffold = scaffolds[scaffold_name]
gene_start = gene_orf.start.iat[0] - 1
gene_end = gene_orf.end.iat[0]
gene = scaffold[gene_start:gene_end]

# Get promoter
promoter_region = (-250, +50)
promoter_start = max(gene_start + promoter_region[0], 0)
promoter_end = min(gene_start + promoter_region[1], gene_end)
promoter = scaffold[promoter_start:promoter_end]

#%% Score promoter region
lexA_sites = base_path + "binding_sites\\lexA.Alphaproteobacteria.txt"
pssm = PSSMScorer(lexA_sites)
print pssm.search(promoter)