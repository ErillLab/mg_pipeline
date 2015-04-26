# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 18:26:07 2015

@author: Talmo
"""

from igc_pipeline import *
from Bio import Entrez
from tqdm import tqdm
from glob import glob
import os

#%% Configuration
# Paths
eggnog_path = "/home/cuda/2TB/eggNOG/"
members_path = eggnog_path + "NOG.members.txt"
funccats_path = eggnog_path + "NOG.funccat.txt"
blast_out_path = IGC_path + "Taxonomy_BLAST/has_cog-eggNOG/"

# NCBI
Entrez.email = "erill@umbc.edu"

#%% Cache data files
t = time()
# Load eggNOG files
eggnog_members = pd.read_table(members_path, index_col="protein name")
eggnog_funccat = pd.read_table(funccats_path, names=["cog", "funccat"], index_col="cog")

# Load taxonomy lookup table
eggnog_tax = pd.read_csv(eggnog_tax_path, index_col="tax_id")

log("Loaded eggNOG data files.", t)

#%% Find not BLASTed
samples = get_all_samples()
blasted = [p.replace(blast_out_path, "").replace(".tbl", "") for p in glob(blast_out_path + "*.tbl") if os.stat(p).st_size > 0]
not_blasted = [sample for sample in samples if sample not in blasted]

#%% Processing functions
def extract_ORFs_with_COGs(sample):
    """ Extracts ORFs for genes that have COG assignment and saves their
    FASTA sequences. """
    sample = get_unique_sample(sample)
    sample_paths = get_sample_paths(sample)
    
    # Load sample data
    genes = load_sample_genes(sample)
    
    # Get ORFs with COGs
    t = time()
    has_cog = genes[genes.eggNOG != "unknown"].index.tolist()
    ORF_seqs = SeqIO.index(sample_paths["ORFs_fasta"], "fasta")
    valid_orfs = [ORF_seqs[orf] for orf in has_cog]

    # Create containing folder if doesn't exist
    if not os.path.exists(orfs_fasta_path + "has_cog"):
        os.makedirs(orfs_fasta_path + "has_cog")
    
    # Save to file
    valid_orfs_path = orfs_fasta_path + "has_cog/" + sample + ".fna"
    SeqIO.write(valid_orfs, valid_orfs_path, "fasta")
    log("Saved ORFs with COGs for %s." % sample, t)

def query_taxonomy(tax_ids, use_first_result=True):
    """ Returns the taxonomy for by querying the NCBI Taxonomy database. """
    
    # Make sure tax_ids is a list of strings
    if type(tax_ids) != list:
        tax_ids = [tax_ids]
    tax_ids = [str(int(tax_id)) for tax_id in tax_ids]
    
    # Query the taxonomy database via Entrez
    results = pd.DataFrame(Entrez.parse(Entrez.efetch("taxonomy", id=tax_ids)))
    
    # Resolve multiple results per query
    if use_first_result:
        results = results.groupby("TaxId").first().reset_index()
    
    # Return taxonomy tables
    taxonomy = [pd.DataFrame(lineage) for lineage in results.LineageEx]
    if len(taxonomy) == 1:
        taxonomy = taxonomy[0]
    return taxonomy

def update_genes_with_eggnog(sample):
    """ Updates gene table for a given sample using BLAST output from querying
    the eggNOG database. """
        
    # Load BLAST output
    blast_out = pd.read_table(blast_out_path + sample + ".tbl", names=blast_columns)
    
    # Get hits by lowest e-value
    keep_cols = ["hit", "e_value"]
    blast_hits = blast_out.groupby("query").first()[keep_cols]
    
    # Get COGs and functinal categories
    blast_hits = blast_hits.join(eggnog_members["#nog name"], on="hit")
    blast_hits = blast_hits.join(eggnog_funccat, on="#nog name")
    
    # Parse out eggNOG taxonomy ID
    blast_hits = blast_hits.join(blast_hits["hit"].str.extract("(?P<tax_id>\d+)\.(?P<gene_id>.*)").convert_objects(convert_numeric=True))
    
    # Get taxonomy for hits
    hits_tax = blast_hits.join(eggnog_tax, on="tax_id")
    hits_tax.rename(columns={"hit": "eggNOG_id", "e_value": "eggNOG_e_value", "#nog name": "eggNOG", "funccat": "eggNOG_funccat"}, inplace=True)
    hits_tax.drop(["tax_id", "gene_id"], axis=1, inplace=True)
    hits_tax.index.name = "gene_name"
    
    # Load current genes
    genes = load_sample_genes(sample)
    
    # Update with new taxonomy
    genes = genes.join(hits_tax, lsuffix="_old")
    
    # Save
    sample_genes_path = get_sample_paths(sample).genes
    genes.to_csv(sample_genes_path)
    
#%% Batch process samples
samples = get_all_samples()

failed = []; done = []
for sample in tqdm(samples):
    try:
        update_genes_with_eggnog(sample)
        done.append(sample)
    except:
        failed.append(sample)
        continue
    