# -*- coding: utf-8 -*-
"""
Checking against the gene summary table.

Created on Sat Mar 21 08:20:48 2015

@author: Talmo
"""
import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from igc_pipeline import *

#%% Pre-processing
# Load big gene table
genes = load_gene_summary()

# Match against regex
t = time()
sample_matches = genes.gene_name.str.extract(samples_regex)
log('Extracted sample regex.', t)

# Get unique names for samples
t = time()
samples = sample_matches.apply(get_unique_sample)
log('Got unique samples.', t)

# Find genes with no sample match
no_match = genes.gene_name[samples.isnull()]

# Find which samples are unmatched
matched_samples = samples[~samples.isnull()].unique()
unmatched_samples = samples_index[~samples_index.Assemblies_filenames.isnull()].index.diff(matched_samples)

#%% Search for unmatched gene names in unmatched ORF files
no_match_idx = pd.Index(no_match.values)
search_matches = {}
for i, sample in enumerate(unmatched_samples):
    print "[%d/%d]" % (i+1, len(unmatched_samples)),
    ORFs = parse_ORFs(sample)
    search_matches[sample] = ORFs.index.intersection(no_match_idx)
    if len(search_matches[sample]) > 0:
        print "Hit:", sample

#%% Search for unmatched gene names in matched ORF files
no_match_idx = pd.Index(no_match.values)
for i, sample in enumerate(matched_samples[653:]):
    print "[%d/%d]" % (i+1, len(matched_samples)),
    ORFs = parse_ORFs(sample)
    search_matches[sample] = ORFs.index.intersection(no_match_idx)
    if len(search_matches[sample]) > 0:
        print "Hit:", sample
