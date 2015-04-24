# -*- coding: utf-8 -*-
"""
Some validation for sample <=> filename detection.

Created on Sat Mar 21 00:38:32 2015

@author: Talmo
"""
import os
import pandas as pd

#%% Configuration
# This script's directory
base_path = os.path.dirname(os.path.abspath(__file__))

# Modified gut_microbiome_samples_20140613.xls with inferred filenames
xls_path = os.path.join(base_path, 'gut_microbiome_samples_20140613.xls')

# Output file containing unified samples index table
samples_index_path = os.path.join(base_path, 'samples_index.csv')

#%% Parse
# Read Excel file
xls = pd.read_excel(xls_path)

# Sample data
samples = xls.loc[:, ['Name', 'Alias', 'SRASample', 'BioSample', 'SRAStudy', 'Study reference/citation', 'NCBI_taxon_ID', 'NCBI_taxon_name']]
samples.rename(columns={'Name': 'sample'}, inplace=True)
samples.dropna(how='all', inplace=True)

# ORFs
ORFs = xls.loc[:, ['ORFs_filenames', 'Inferred_name']]
ORFs.rename(columns={'Inferred_name': 'sample'}, inplace=True)
ORFs.dropna(how='all', inplace=True)

# Assemblies
assemblies = xls.loc[:, ['Assemblies_filenames', 'Inferred_name.1']]
assemblies.rename(columns={'Inferred_name.1': 'sample'}, inplace=True)
assemblies.dropna(how='all', inplace=True)

#%% Validation
print "=== ORFs.sample not in samples.sample:"
print ORFs[~ORFs.sample.isin(samples.sample)]
print "=== assemblies.sample not in samples.sample:", assemblies.sample.isin(samples.sample).all()
print assemblies[~assemblies.sample.isin(samples.sample)]
print

#%% Merge ORFs
samples = pd.merge(samples, ORFs, on='sample', how='outer')

#%% Missing assemblies
print "=== Null assemblies:"
print assemblies[assemblies.sample.isnull()]
print "=== Samples not in assemblies:"
print samples.loc[~samples.sample.isin(assemblies.sample[~assemblies.sample.isnull()]), ["sample", "Alias", "ORFs_filenames"]]

#%% Merge assemblies
# Use Alias to predict missing samples
missing = assemblies[assemblies.sample.isnull()].copy()
missing['Alias'] = missing.Assemblies_filenames.str.extract('(.*)[.]scaf')
missing['sample'] = samples.reset_index().set_index('Alias').loc[missing.Alias, 'sample'].tolist()
assemblies.loc[assemblies.sample.isnull(), 'sample'] = missing['sample']

# Merge to rest of samples
samples = pd.merge(samples, assemblies, how='outer', on='sample')

# Validate
print 'All assemblies in samples:', assemblies.Assemblies_filenames.isin(samples.Assemblies_filenames).all()

#%% Save
samples['Study reference/citation'] = samples['Study reference/citation'].str.encode('utf-8')
samples.rename(columns={'Alias':'alias', 'Study reference/citation': 'study'}, inplace=True)
samples.to_csv(samples_index_path, index=False)

#%% Test
loaded_samples = pd.read_csv(samples_index_path)
print 'Loaded samples == samples:'
print ((loaded_samples == samples) | (loaded_samples.isnull() & samples.isnull())).all()