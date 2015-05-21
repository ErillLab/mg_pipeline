# -*- coding: utf-8 -*-
"""
This script checks whether each sample's ORF FASTA files or scaftig assembly
files exist and cross-references with the sample index.

Created on Wed May 20 12:09:53 2015

@author: Talmo
"""

import sys
sys.path.append("../..")
from igc_pipeline import *
import pandas as pd
import os
from glob import glob

#%% Completeness of samples index (samples_index.csv)
samples_idx = pd.read_csv(samples_index_path)
all_samples = samples_idx.sample
samples_missing_ORFs_path = all_samples[samples_idx.ORFs_filenames.isnull()]
samples_missing_Assembly_path = all_samples[samples_idx.Assemblies_filenames.isnull()]
samples_missing_any = all_samples[samples_idx.isnull().any(axis=1)]
samples_with_paths = all_samples[~samples_idx.isnull().any(axis=1)]

print "All samples (reported in samples_index.csv):", len(all_samples)
print "Missing ORF (FASTA) path:", len(samples_missing_ORFs_path)
print "Missing Assembly (scaftig) path:", len(samples_missing_Assembly_path)
print "Missing anything:", len(samples_missing_any)
print "Not missing anything:", len(samples_with_paths)
print
print "Missing ORF <=> Missing assembly:", (samples_missing_ORFs_path == samples_missing_Assembly_path).all()
print "Samples with missing paths:"
print ", ".join(samples_missing_any)
print

#%% ORFs files
ORF_filenames = pd.DataFrame()
ORF_filenames["original_path"] = samples_idx.set_index("sample").ORFs_filenames[samples_with_paths]
ORF_filenames["standard_path"] = ORF_filenames.index + ".fna"
ORF_filenames = orfs_fasta_path + ORF_filenames
ORFs_exist = ORF_filenames.applymap(os.path.exists)

other_files = [p for p in glob(orfs_fasta_path + "*") if ~(ORF_filenames == p).any().any() and not os.path.isdir(p)]

print "ORFs (FASTA) folder:", orfs_fasta_path
print "Samples with original path: %d/%d" % (ORFs_exist.original_path.sum(), len(ORFs_exist))
print "Samples with standard path: %d/%d" % (ORFs_exist.standard_path.sum(), len(ORFs_exist))
print "Other files in folder (%d): %s" % (len(other_files), ", ".join(other_files))
print

#%% Assembly files
Assembly_filenames = pd.DataFrame()
Assembly_filenames["original_path"] = samples_idx.set_index("sample").Assemblies_filenames[samples_with_paths]
Assembly_filenames["standard_path"] = Assembly_filenames.index + ".fna"
Assembly_filenames = scaftigs_path + Assembly_filenames
Assembly_exist = Assembly_filenames.applymap(os.path.exists)

other_files = [p for p in glob(scaftigs_path + "*") if ~(Assembly_filenames == p).any().any() and not os.path.isdir(p)]

print "Assembly (scaftigs) folder:", scaftigs_path
print "Samples with original path: %d/%d" % (ORFs_exist.original_path.sum(), len(ORFs_exist))
print "Samples with standard path: %d/%d" % (ORFs_exist.standard_path.sum(), len(ORFs_exist))
print "Other files in folder (%d): %s" % (len(other_files), ", ".join(other_files))
print