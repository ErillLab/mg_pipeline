# -*- coding: utf-8 -*-
"""
This script automates processing of the IGC dataset.

Created on Tue Mar 17 19:13:08 2015

@author: Talmo
"""

import os
import sys
import time
from glob import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from PSSMScorer import PSSMScorer

# From the timeit module:
if sys.platform == 'win32':
    time = time.clock # On Windows, the best timer is time.clock
else:
    time = time.time # On most other platforms the best timer is time.time
    
#%% Configuration
# Output verbosity
verbosity = 2 # 0 = no output, 1 = minimal, 2 = debugging

# Base paths
base_path = "/home/cuda/2TB/metagenomics/" # metagenomics folder
IGC_path = base_path + "IGC/"

# Data paths
gene_summary_path = IGC_path + "3.IGC.AnnotationInfo/IGC.annotation.summary.v2"
scaftigs_path = IGC_path + "4.IndividualAssmeblies/"
orfs_fasta_path = IGC_path + "5.IndividualORFs/"

# Processed data paths
genes_path = IGC_path + "Genes/"
orfs_path = IGC_path + "ORFs/"
operons_path = IGC_path + "Operons/"

# PSSM scoring
binding_sites_path = base_path + 'binding_sites/'
Firmicutes_LexA = PSSMScorer(binding_sites_path + 'Firmicutes_LexA.txt')
GammaProteobacteria_LexA = PSSMScorer(binding_sites_path + 'GammaProteobacteria_LexA.txt')

# Operon prediction
threshold_IGI = 50 # max intergenic interval (bp)
promoter_region = (-250, +50) # bp relative to gene start


def log(msg, start_time=None, verbosity_level=1):
    """ Convenience function for printing pipeline messages.
    
    Args:
        msg:
            the message string
        start_time:
            if provided, displays time elapsed since this time
        verbosity_level:
            if the verbosity is less than this, nothing will be printed
        
    Example:
        >>> start_time = time()
        >>> log("Done processing.", start_time) 
        Done processing. [10.94s] """
        
    if verbosity_level <= verbosity:
        if start_time is not None:
            msg = msg + " [%.2fs]" % (time() - start_time)
        print msg

def get_sample_paths(sample):
    """ Returns paths to sample data files. """
    paths = {}
    
    # TODO: Better path detection
    try:
        paths["scaftigs"] = glob(scaftigs_path + "*" + sample + "*")[0]
    except:
        log("Could not find sample scaftigs.")
    try:
        paths["ORFs_fasta"] = glob(orfs_fasta_path + "*" + sample + "*")[0]
    except:
        log("Could not find sample ORF FASTA.")
    
    # Processed data
    paths["genes"] = genes_path + sample + ".csv"
    paths["ORFs"] = orfs_path + sample + ".csv"
    paths["operons"] = operons_path + sample + ".csv"
    
    return paths

def get_MetaHit():
    """ Returns the names of the samples from the MetaHit database. """
    # Search for MetaHit scaftigs
    MH_paths = pd.Series(glob(scaftigs_path + "MH*"))
    
    # Extract sample names from paths
    MH_samples = MH_paths.str.extract("(MH\\d+)").values
    
    # Sort by sample number
    MH_samples = sorted(MH_samples, key=lambda x: int(x[2:]))
    return MH_samples
    

#%% Genes processing
def load_gene_summary(limit=None, sample_col=True):
    """ Loads the integrated gene summary.
    
    Warning: This may take a while as there are 9879896 rows in the table.
    
    Args:
        limit: int, default None
            The number of rows to return. If None, returns all.
        sample_col: bool, default True
            If True, adds a column containing the sample name extracted from 
            the gene_name.
    
    Returns:
        gene_summary: DataFrame
            The IGC gene summaries table indexed by gene_id.
        
    Example:
        >>> print load_gene_summary(1)
    ================  ===========  ============  =============  =======  =======  ======  =======  ==============  ===============  ================  =======================================================================  ================
    gene_name         gene_length  completeness  cohort_origin  phylum   genus    kegg    eggNOG   sample_freq     individual_freq  eggNOG_funccat    kegg_funccat                                                             cohort_assembled
    ================  ===========  ============  =============  =======  =======  ======  =======  ==============  ===============  ================  =======================================================================  ================
    T2D-6A_GL0083352  88230        Complete      CHN            unknown  unknown  K01824  COG5184  0.224151539069  0.236448598131   Lipid Metabolism  Cell cycle control, cell division, chromosome partitioning;Cytoskeleton  EUR;CHN;USA
    ================  ===========  ============  =============  =======  =======  ======  =======  ==============  ===============  ================  =======================================================================  ================
    """
    
    t = time()
    gene_summary_cols = ['gene_id', 'gene_name', 'gene_length', 'completeness', \
        'cohort_origin', 'phylum', 'genus', 'kegg', 'eggNOG', 'sample_freq', \
        'individual_freq', 'eggNOG_funccat', 'kegg_funccat', 'cohort_assembled']
    gene_summary = pd.read_table(gene_summary_path, names=gene_summary_cols, index_col='gene_id', nrows=limit)
    
    # Add sample column to gene_summary table
    idx = gene_summary.gene_name.map(lambda x: x[0].isdigit())
    samples = pd.Series(index=gene_summary.index)
    samples[idx] = gene_summary[idx].gene_name.str.split('.').str[0]
    samples[~idx] = gene_summary[~idx].gene_name.str.split('_').str[0]
    #idx = samples.map(lambda x: len(x) == 1)
    #samples[idx] = samples[idx].map(lambda x: x[0].split('.'))
    gene_summary['sample'] = samples
    
    log("Loaded %d genes from summary." % (gene_summary.shape[0]), t)
    return gene_summary

def save_individual_genes(overwrite=False, gene_summary=None):
    """ Splits the integrated gene summary into individual tables for each 
    sample.
    
    Saves data to: genes_path/[sample].csv. 
    
    """
    
    # Create parent folder if it does not exist
    if not os.path.exists(genes_path):
        os.makedirs(genes_path)
    
    # Load the combined gene summary table
    if gene_summary is None:
        gene_summary = load_gene_summary()
    
    # Save the list of samples
    #gene_summary['sample'].unique().tofile(IGC_path + 'samples.txt', '\n')
    
    # Group by sample and save data
    t = time()
    gene_summary = gene_summary.groupby('sample', sort=False)
    for sample, genes in gene_summary:
        sample_genes_path = genes_path + sample + ".csv"
        
        if overwrite or not os.path.exists(sample_genes_path):
            genes.to_csv(sample_genes_path)
            
        #log('%s: %d' % (sample, genes.shape[0]), None, 2)
        log('%s: %s' % (sample, genes_path + sample + ".csv"), None, 2)
    
    log("Saved genes for %d samples individually." % gene_summary.ngroups, t)

    
#%% ORF and Operon processing
def parse_ORFs(sample, overwrite=False):
    """ Parses the FASTA headers in the ORF FASTA file and saves the data.
    
    Args:
        sample: str
            The name of the sample
        overwrite: boolean, default False
            If True, will overwrite saved ORF table, otherwise will attempt to
            load the table from cache.
            
    Returns:
        ORFs: DataFrame
            Table with data on the open reading frames and the columns:
                gene_name, orf_type, scaffold, start, end, completeness
                
    Example:
        >>> ORFs = parse_ORFs('MH0001')
        Parsed 21452 ORFs from MH0001. [0.44s]
        >>> ORFs.head(1)
        ...
        """
    t = time()
    
    # Get paths
    sample_paths = get_sample_paths(sample)
    
    # Check for cache
    if not overwrite and os.path.exists(sample_paths["ORFs"]):
        ORFs = pd.read_csv(sample_paths["ORFs"], index_col=0)
        log("Loaded %d cached ORFs from %s." % (ORFs.shape[0], sample), t)
        return ORFs
    
    # Read in FASTA headers
    with open(sample_paths["ORFs_fasta"]) as f:
        lines = [line[1:-2] for line in f if line[0] == '>']
        
    # Parse FASTA headers for ORFs
    ORF_pattern = '(?P<gene_name>[^ ]+)[ ]+\\[(?P<orf_type>[^\]]+)\\][ ]+locus=(?P<scaffold>[^:]+):(?P<start>[^:]+):(?P<end>[^:]+):(?P<strand>[^:\[]+)\\[(?P<completeness>[^\[\]]+)'
    ORFs = pd.Series(lines).str.extract(ORF_pattern).convert_objects(convert_numeric=True)
    ORFs.set_index("gene_name", inplace=True)
    
    # Create parent folder if it does not exist
    if not os.path.exists(orfs_path):
        os.makedirs(orfs_path)
    
    # Save
    ORFs.to_csv(sample_paths["ORFs"])
    log("Parsed %d ORFs from %s." % (ORFs.shape[0], sample), t)
    
    return ORFs
    
    
def predict_operons(sample, overwrite=False):
    """ Predicts operons for the given sample. """
    t = time()
    
    # Get paths
    sample_paths = get_sample_paths(sample)
    
    # Check for cache
    if not overwrite and os.path.exists(sample_paths["operons"]):
        operons = pd.read_csv(sample_paths["operons"], index_col=0)
        log("Loaded %d cached operons from %s." % (operons.shape[0], sample), t)
        return operons
    
    # Load sample scaffolds (scafftigs)
    scaffolds = SeqIO.index(sample_paths["scaftigs"], 'fasta')
    
    # Load the predicted open reading frames
    ORFs = parse_ORFs(sample)
    
    # Group by 'scaffold' column and do groupwise operations
    ORFs = ORFs.groupby('scaffold', sort=False)
    n_ORFs = ORFs.ngroups
    operon_tables = []; i = 0
    for scaffold_name, scaf_ORFs in ORFs:
        for strand, strand_ORFs in scaf_ORFs.groupby('strand', sort=False):
            # Sort by start values
            strand_ORFs.sort('start', inplace=True)
        
            # Compute intergenic intervals
            IGIs = strand_ORFs.start[1:].values - strand_ORFs.end[:-1].values
        
            # Find operon gene indices
            head_genes = np.where(IGIs > threshold_IGI)[0] + 1
            operons_start = np.hstack(([0], head_genes))
            operons_end = np.hstack((head_genes, [scaf_ORFs.shape[0]]))
            head_gene_start = scaf_ORFs.start.iloc[operons_start]
        
            # Get promoter regions
            promoter_start = (head_gene_start + promoter_region[0]).clip(1).reset_index(drop=True)
            promoter_end = (head_gene_start + promoter_region[1]).reset_index(drop=True)
            promoter_seq = [str(scaffolds[scaffold_name][s-1:e].seq) for s, e in zip(promoter_start, promoter_end)]
        
            # Build operon table
            strand_operons = pd.DataFrame()
            strand_operons['genes'] = [strand_ORFs.index[op_start:op_end].tolist() for op_start, op_end in zip(operons_start, operons_end)]
            strand_operons['strand'] = strand
            strand_operons['promoter_start'] = promoter_start
            strand_operons['promoter_end'] = promoter_end
            strand_operons['promoter_seq'] = promoter_seq
            strand_operons['head_completeness'] = strand_ORFs.ix[op_start, 'completeness']
        
            # Append to list of operon tables
            operon_tables.append(strand_operons)
        
        i = i + 1
        if verbosity >= 2 and i % (n_ORFs / 5) == 0:
            log("%.2f%%: %d/%d" % (float(i) / n_ORFs * 100, i, n_ORFs), None, 2)
        
    # Merge operon tables
    operons = pd.concat(operon_tables, ignore_index=True)
    operons.index.name = 'operon'
    
    # Create parent folder if it does not exist
    if not os.path.exists(operons_path):
        os.makedirs(operons_path)
        
    # Save
    operons.to_csv(sample_paths["operons"])
    
    log("Predicted %d operons for %s." % (operons.shape[0], sample), t)
    return operons


#%% PSSM scoring
def score_sample(sample, PSSM):
    """ Scores the promoters of all the operons in a sample using the PSSM. 
    
    Args:
        sample: str
            Name of the sample to score
        PSSM: PSSMScorer
            PSSM to use for scoring
        """
        
    # Load operons for sample
    operons = predict_operons(sample)
    
    # Score all promoters
    t = time()
    scores = operons.promoter_seq.apply(PSSM.score_all)
    
    log("Scored %d promoters for %s." % (scores.shape[0], sample), t)
    return scores
    