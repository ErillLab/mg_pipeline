# -*- coding: utf-8 -*-
"""
This script automates processing of the IGC dataset.

Created on Tue Mar 17 19:13:08 2015

@author: Talmo
"""

import os
import sys
import time
import ast
from glob import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
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
script_dir_path = os.path.dirname(os.path.abspath(__file__)) + os.path.sep
IGC_path = base_path + "IGC/"

# Samples index
samples_index_path = script_dir_path + "samples_index.csv"

# IGC Data paths
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

#%% General
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

# Load pre-processed samples index
samples_index = pd.read_csv(samples_index_path, index_col="sample")
samples_aliases = samples_index["alias"].reset_index().set_index("alias")
data_samples = samples_index[~(samples_index.Assemblies_filenames.isnull() | samples_index.ORFs_filenames.isnull())]

def get_unique_sample(query_sample, error=False):
    """ Returns a unique sample name by resolving ambiguity with aliases.
    
    Args:
        query_sample: str
            Sample name to look up
            
    Returns:
        sample: str
            Unique sample name that indexes into samples_index table.
            Returns None if query_sample is not found or raises LookupError.
            
    Example:
        >>> get_unique_sample("MH192")
        'MH0192'
        >>> samples_index.loc[get_unique_sample("MH192")]
        ...
    """
    
    if query_sample in samples_index.index:
        # Sample name is already unique
        return query_sample
    elif query_sample in samples_aliases.index:
        # Get sample name from alias
        return samples_aliases.at[query_sample, "sample"]
    else:
        # Sample name wasn't found
        if error:
            raise LookupError("Sample not found: %s" % query_sample)
        return None

def get_sample_paths(sample):
    """ Returns paths to sample data files. 
    
    Args:
        sample: str
            Sample to query.
            This name will be disambiguated by get_unique_sample().
            Raises an exception if not found.
            
    Returns:
        paths: Series
            Absolute paths to the files associated with the sample. Keys:
            
            :sample: Unique sample name (see ``get_unique_sample()``)
            :scaftigs: Assembled scaftigs (in ``scaftigs_path``)
            :ORFs_fasta: ORF FASTA files (in ``orfs_fasta_path``)
            :genes: Summary of genes in sample (see ``save_sample_genes()``)
            :ORFs: ORF table without sequences (see ``parse_ORFs()``)
            :operons: Predicted operons (see ``predict_operons()``)
            
    Example:
        >>> get_sample_paths('MH192')
        scaftigs      /home/cuda/2TB/metagenomics/IGC/4.IndividualAs...
        ORFs_fasta    /home/cuda/2TB/metagenomics/IGC/4.IndividualAs...
        genes          /home/cuda/2TB/metagenomics/IGC/Genes/MH0192.csv
        ORFs            /home/cuda/2TB/metagenomics/IGC/ORFs/MH0192.csv
        operons       /home/cuda/2TB/metagenomics/IGC/Operons/MH0192...
        dtype: object
    """
    paths = pd.Series()
    
    # Get unique sample name
    sample = get_unique_sample(sample, True)
    
    # IGC Data
    paths["scaftigs"] = scaftigs_path + samples_index.at[sample, "Assemblies_filenames"]
    paths["ORFs_fasta"] = orfs_fasta_path + samples_index.at[sample, "ORFs_filenames"]
    
    # Processed data
    paths["genes"] = genes_path + sample + ".csv"
    paths["ORFs"] = orfs_path + sample + ".csv"
    paths["operons"] = operons_path + sample + ".csv"
    
    # Make sure paths are absolute
    paths = paths.apply(os.path.abspath)
    
    return paths

def has_data(sample):
    """ Checks if the sample has its data files. 
    
    Args:
        sample: str
            The sample to be queried
            
    Returns:
        data_exists: Series
            Same keys as returned by ``get_sample_paths()`` but will be True
            for the ones for which the file exists.
    """
    return get_sample_paths(sample).apply(os.path.exists)

def check_samples_data():
    """ Returns a table of all samples and which data files are present. """
    return data_samples.reset_index().sample.apply(has_data).set_index(data_samples.index)

def get_MetaHit(study=2010):
    """ Returns the names of the samples from the MetaHit database. 
    
    Args:
        study: int or str, default 2010
            Year of the study to return. Studies:
            
            :2010: Qin et al. (2010) | doi:10.1038/nature08821
            :2012: Qin et al. (2012) | doi:10.1038/nature11450
            :2010: Le Chatelier et al. (2013) | doi:10.1038/nature12506
            :"all": Returns all studies.
    
    Returns:
        MetaHit_samples: list
            The samples corresponding to the study selected
    """
    
    queries = {"2010": "A human gut microbial gene catalogue", \
               "2012": "A metagenome-wide association", \
               "2013": "Richness of human gut"}
    
    if study == 'all':
        return [item for inner_list in [get_MetaHit(k) for k in queries.keys()] for item in inner_list]
    
    return samples_index[samples_index.study.str.contains(queries[str(study)])].index.tolist()
    

#%% Genes processing
# Matches all sample names and aliases (use get_unique_sample() to disambiguate)
samples_regex = "(" + "|".join(set(samples_index.index.tolist() + samples_index.alias.tolist())) + ")"
def load_gene_summary(limit=None, extract_samples=False):
    """ Loads the integrated gene summary.
    
    Warning: This may take a while as there are 9879896 rows in the table.
    
    Args:
        limit: int, default None
            The number of rows to return. If None, returns all.
        extract_sample: bool, default False
            Tries to extract the sample that each gene comes from and appends 
            it to the "sample" column (see ``extract_sample()``).
    
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
    log("Loaded %d genes from summary." % (gene_summary.shape[0]), t)
    
    if extract_samples:
        gene_summary["sample"] = extract_sample(gene_summary["gene_name"])
    
    return gene_summary

def extract_sample(gene_names):
    """ Returns a Series of samples given a series of gene_names. """
    
    # Match against regex
    t = time()
    samples = gene_names.str.extract(samples_regex)
    log("Extracted sample regex.", t)
    
    # Get unique names for samples
    t = time()
    samples = samples.apply(get_unique_sample)
    log("Got unique samples.", t)
    
    return samples


def save_sample_genes(overwrite=False, gene_summary=None):
    """ Splits the integrated gene summary into individual tables for each 
    sample.
    
    Saves data to: [genes_path]/[sample].csv
    
    """
    
    # Create parent folder if it does not exist
    if not os.path.exists(genes_path):
        os.makedirs(genes_path)
    
    # Load the combined gene summary table
    if gene_summary is None:
        gene_summary = load_gene_summary()

    # Extract samples from gene names
    if "sample" not in gene_summary.columns:
        gene_summary["sample"] = extract_sample(gene_summary["gene_name"])

    # Drop genes without sample
    n_total = len(gene_summary)
    gene_summary = gene_summary[~gene_summary.sample.isnull()]
    n_after = len(gene_summary)
    log("%d/%d genes without sample (%.2f%%)." % (n_total - n_after, n_total, 100.0*(n_total - n_after)/n_total), None, 2)
    
    empty_samples = samples_index[~(samples_index.Assemblies_filenames.isnull() | samples_index.ORFs_filenames.isnull())].index.diff(gene_summary.sample.unique())
    log("%d/%d samples without genes (%.2f%%)." % (len(empty_samples), len(samples_index), 100.0*len(empty_samples)/len(samples_index)), None, 2)
    
    # Group by sample and save data
    t = time()
    gene_summary = gene_summary.groupby("sample", sort=False)
    for sample, sample_genes in gene_summary:
        sample_genes_path = get_sample_paths(sample).genes
        
        if overwrite or not os.path.exists(sample_genes_path):
            sample_genes.to_csv(sample_genes_path)
            
    log("Saved genes for %d samples individually." % gene_summary.ngroups, t)

def load_sample_genes(sample):
    """ Loads the genes for a sample. """
    # Get path
    sample_genes_path = get_sample_paths(sample).genes
    
    # Check if it exists
    if not os.path.exists(sample_genes_path):
        raise EnvironmentError("Genes file for %s does not exist: %s" % (sample, sample_genes_path))
        
    # Load and return
    sample_genes = pd.read_csv(get_sample_paths(sample).genes, index_col="gene_name")
    sample_genes.sort_index(inplace=True)
    return sample_genes
        

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
    
    # Get unique sample and paths
    sample = get_unique_sample(sample)
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
        operons = pd.read_csv(sample_paths["operons"], index_col=0, converters={'genes':ast.literal_eval})
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
    
def plot_scores(scores):
    # Analysis
    all_scores = np.hstack(scores.values)
    plt.figure()
    plt.hist(all_scores, bins=20, range=(10, 21), normed=True)
    plt.title('Firmicutes_LexA (Sample: MH0001)')
    plt.xlabel('PSSM Score')
    plt.ylabel('Number of sites')
    
    high_scores = all_scores[all_scores > 16.0]
    high_scores.sort()
    plt.figure()
    plt.bar(range(len(high_scores)), high_scores)
    plt.title('Firmicutes_LexA (Sample: MH0001) > 16.0 bits')
    plt.xlabel('Top Sites')
    plt.ylabel('PSSM Score')
    
def batch_score(samples):
    for sample in samples:
        try:
            scores = score_sample(sample, Firmicutes_LexA)
            scores.to_csv(IGC_path + 'Scores/Firmicutes_LexA/'+sample+'.csv')
        except:
            print 'Failed:', sample