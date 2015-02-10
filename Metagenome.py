# -*- coding: utf-8 -*-
"""
This class contains the methods for processing a metagenome.

Created on Mon Feb 09 23:19:42 2015

@author: Talmo
"""
import re
import os
#import sys
from time import clock
import json

#import numpy as np
import pandas as pd
#import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline

import mg_rast
#import pssm
#import database as db
#import taxonomy as tax

""" Configuration """
# Folder where metagenomes are stored
metagenomes_path = "E:\\metagenomics\\metagenomes"

# Path to BLAST database for eggNOG
eggnog_db_path = "E:\\eggNOG\\blast_db\\eggnogv4"

# Path to eggNOG reference tables
eggnog_NOG_members_path = "E:\\eggNOG\\NOG.members.txt"
eggnog_NOG_funccat_path = "E:\\eggNOG\\NOG.funccat.txt"
eggnog_bactNOG_members_path = "E:\\eggNOG\\all.members\\bactNOG.members.txt"
eggnog_bactNOG_funccat_path = "E:\\eggNOG\\all.funccat\\bactNOG.funccat.txt"


""" Static methods """
def get_paths(mg_id):
    """ Returns the paths to files for the specified metagenome. """
    
    # Paths
    paths = {}
    paths["base"] = metagenomes_path + "\\" + mg_id
    paths["metadata"] = paths["base"] + "\\" + mg_id + "_metadata.json"
    paths["299"] = paths["base"] + "\\" + mg_id + ".299.screen.passed.fna"
    paths["350"] = paths["base"] + "\\" + mg_id + ".350.genecalling.faa"
    paths["650"] = paths["base"] + "\\" + mg_id + ".650.protein.sims"
    paths["operon_hits"] = paths["base"] + "\\" + mg_id + "_operon-hits.faa"
    paths["eggnog_blast"] = paths["base"] + "\\" + mg_id + "_eggnogv4-BLAST.tbl"
    paths["db"] = paths["base"] + "\\" + mg_id + ".db"
    
    paths["hits"] = paths["base"] + "\\hits.csv"
    paths["operons"] = paths["base"] + "\\operons.csv"
    
    # Make paths absolute
    for key, value in paths.items():
        paths[key] = os.path.abspath(paths[key])
    
    return paths


class Metagenome:
    """ This class can be used to process a metagenome from MG-RAST. """
    
    def __init__(self, mg_id, verbosity=1):
        """ Initializes the class. """
        self.mg_id = mg_rast.validate_mg_id(mg_id)
        self.id = self.mg_id
        self.verbosity = verbosity
        self.paths = get_paths(mg_id)

        self.scaffolds = None # 299 -> load_scaffolds()
        self.genes = None # 350 -> load_genes()
        self.hits = None # 650 -> load_hits()
        self.operons = None # predict_operons()
        
        self.n_scaffolds = None
        self.n_hits = None
        self.n_operons = None
        
        print "Initialized metagenome: " + self.mg_id
        self.load_metadata()
        
    def check_path(self, what):
        """ Checks if the specified metagenome file exists. """
        return os.path.exists(self.paths[what])
        
    def load_metadata(self, force=False):
        """ Loads MG-RAST metadata associated with metagenome. """
        # Create base directory if it doesn't exist
        if not self.check_path('base'):
            os.makedirs(self.paths['base'])
            
        # Query API for metadata and save
        if force or not self.check_path('metadata'):
            metadata = mg_rast.get_metadata(self.mg_id, parse=False)
            open(self.paths['metadata'], 'w').write(metadata)
            if self.verbosity > 0:
                print "Saved metadata for metagenome."
        else:
            metadata = open(self.paths['metadata']).read()
            if self.verbosity > 0:
                print "Loaded cached metadata for metagenome."
        
        # Deserialize data
        self.metadata = json.loads(metadata)

        # Parse out basic info about metagenome
        self.name = self.metadata["name"]
        self.project = self.metadata["project"][0]
        
        if self.verbosity > 0:
            print " -- Name: " + self.name
            print " -- Project: " + self.project
        
    def download_data(self, force=False):
        """ Downloads the relevant files from MG-RAST. """
        
        # Create base directory if it doesn't exist
        if not self.check_path('base'):
            os.makedirs(self.paths['base'])
            
        # Download data files
        stages = ['299.1', '350.1', '650.1']
        for stage in stages:
            if force or not self.check_path(stage[:-2]):
                if self.verbosity > 0:
                    t = clock(); print 'Downloading stage %s...' % stage,
                mg_rast.download(self.mg_id, stage, folder_path=metagenomes_path)
                if self.verbosity > 0:
                    print "Done. [%.2fs]" % (clock() - t)
            else:
                if self.verbosity > 0:
                    print 'Skipping download for existing stage %s.' % stage
                    
    def load_scaffolds(self):
        """ Loads the scaffolds containing raw sequences from the FASTA data. 
        Sequences come from stage 299 in the MG-RAST pipeline.
        
        Updates scaffolds object after processing. """
        
        # Check if the data file exists
        if not self.check_path('299'):
            self.download_data()
            
        t = clock()
        self.scaffolds = SeqIO.index(self.paths['299'], 'fasta')
        if self.verbosity > 0:
            print 'Scaffolds indexed. [%.2fs]' % (clock() - t)
        
        return self.scaffolds
        
    def load_genes(self):
        """ Loads the predicted proteins from the FASTA file.
        Gene predictions come from stage 350 in the MG-RAST pipeline.
        
        Updates the genes object after processing. """
        
        # Check if the data file exists
        if not self.check_path('350'):
            self.download_data()
        
        t = clock()
        self.genes = SeqIO.index(self.paths['350'], 'fasta')
        if self.verbosity > 0:
            print 'Genes indexed. [%.2fs]' % (clock() - t)
                
    def parse_blat(self):
        """ Loads and parses gene predictions (BLAT hits) from file. 
        BLAT hits come from stage 650 in the MG-RAST pipeline.
        
        Updates the hits table after processing. """
        
        # Check if the data file exists
        if not self.check_path('650'):
            self.download_data()
        
        t = clock()
        blat_cols = ['query_id', 'hit_id', 'percentage_identity', 'alignment_length', 'num_mismatches', 'num_gap_openings', 'q.start', 'q.end', 's.start', 's.end', 'e_value', 'score_bits']
        blat = pd.read_table(self.paths['650'], names=blat_cols)
        #blat = blat[blat.query_id.str.contains('scaffold')] # only rows containing scaffold
        blat = blat[~blat.query_id.str.contains('aa')] # only rows not starting with aa
        
        if self.verbosity > 0:        
            print 'BLAT loaded. [%.2fs]' % (clock() - t)
            print 'Parsing BLAT hits...'; t = clock()
            
        n_hits = blat.query_id.nunique()
        hits = pd.DataFrame(index=range(n_hits), columns=['query_id', 'scaffold', 'hit_start', 'hit_end', 'strand', 'md5', 'md5_e_val', 'operon', 'operon_pos'])
        i = 0
        for query_id, query_hits in blat.groupby('query_id'):
            # Parse query ID
            scaffold, hit_start, hit_end, strand = re.match("([\w\d_.]+)_(\d+)_(\d+)_([+-])$", query_id).groups()
            
            # Get the best hit
            idx = query_hits.e_value.argmin()
            md5 = query_hits.hit_id[idx]
            md5_e_val = query_hits.e_value[idx]
            
            # Save (fast scalar indexing)
            hits.at[i, 'query_id'] = query_id
            hits.at[i, 'scaffold'] = scaffold
            hits.at[i, 'hit_start'] = int(hit_start)
            hits.at[i, 'hit_end'] = int(hit_end)
            hits.at[i, 'strand'] = strand
            hits.at[i, 'md5'] = md5
            hits.at[i, 'md5_e_val'] = float(md5_e_val)
            
            # Slow saving (has lots of overhead):
            #hits.ix[i] = [scaffold, int(hit_start), int(hit_end), strand, md5, float(md5_e_val)]
            
            # Report progress
            i = i + 1
            if i % 2000 == 0 and self.verbosity > 1:
                print i, (1.0 * i / n_hits * 100)
                
        if self.verbosity > 0:
            print 'Parsed BLAT hits. [%.2fs]' % (clock() - t)
                
        #TODO: Update database with hits
        hits.to_csv(self.paths['hits'])
        
        # Update instance variable
        self.hits = hits
        self.update_hit_stats()
        
        return hits, blat

    def load_hits(self):
        """ Loads predicted genes from file or parses the BLAT file to create it. """
        in_db = False
        if in_db:
            #TODO: Load from database
            pass
        elif self.check_path('hits'):
            # Load from file
            self.hits = pd.read_csv(self.paths['hits'])
            self.update_hit_stats()
            
            if self.verbosity > 0:
                print 'Loaded hits from file.'
        else:
            # Create the hits table from the BLAT file
            self.parse_blat()
            
    def update_hit_stats(self):
        """ Updates statistics on hits. """
        
        # Make sure hits are loaded
        if self.hits is None:
            self.load_hits()
        
        # Update hit statistics
        self.n_scaffolds = self.hits.scaffold.nunique()
        self.n_hits = self.hits.shape[0]
        
        # Update hit statistics dependent on operons
        if self.operons is not None:
            self.n_scaffolds_in_operons = self.hits[~self.hits.operon.isnull()].scaffold.nunique()
            self.n_hits_in_operons = self.hits[~self.hits.operon.isnull()].shape[0]
            
    def predict_operons(self, min_intergenic_distance=50, promoter_region=(-250, +50)):
        """ Predicts operons based on intergenic distance.
        
        Updates the operons table after processing.
        
        Parameters:
            min_intergenic_distance = 50 # bp
            promoter_region = (-250, +50) # relative to gene start
        """

        # Make sure scaffolds are loaded
        if self.scaffolds is None:
            self.load_scaffolds()
        
        # Make sure hits are loaded
        if self.hits is None:
            self.load_hits()
        
        if self.verbosity > 0:
            print 'Predicting operons...'; t = clock()
        
        # Pre-allocate 2x the number of scaffolds
        operons = pd.DataFrame(index=range(self.n_scaffolds * 2), \
            columns=['scaffold', 'hits', 'strand', 'promoter_start', 'promoter_end', 'promoter_seq'])
        
        # Loop through hits grouped by scaffolds
        i = -1; g = 0
        for scaffold, scaffold_hits in self.hits.groupby('scaffold'):
            # Get sequence
            scaffold_seq = self.scaffolds[scaffold].seq
            
            # Loop through scaffold hits grouped by strand
            for strand, strand_hits in scaffold_hits.groupby('strand'):
                # Sort by start position
                strand_hits = strand_hits.sort('hit_start')
            
                # Build operons
                n_strand_hits = strand_hits.shape[0]
                for j in range(n_strand_hits):
                    idx = strand_hits.index[j]
                    hit_start = strand_hits.hit_start[idx]
                    #hit_end = strand_hits.hit_end[idx]
                    last_hit_end = strand_hits.hit_end.iat[j - 1]
                    
                    # Start new operon
                    if j == 0 or (hit_start - last_hit_end) > min_intergenic_distance:
                        i = i + 1
                        promoter_start = max(0, hit_start + promoter_region[0])
                        promoter_end = min(hit_start + promoter_region[1], len(scaffold_seq))
                        
                        # Save (fast scalar indexing)
                        operons.at[i, 'scaffold'] = scaffold
                        operons.at[i, 'strand'] = strand
                        operons.at[i, 'promoter_start'] = promoter_start
                        operons.at[i, 'promoter_end'] = promoter_end
                        operons.at[i, 'promoter_seq'] = scaffold_seq[promoter_start:promoter_end]
                        
                        # Reset operon
                        operon = []
                        
                    # Add current hit to end of the operon
                    operon.append(idx)
                    
                    # Update operon
                    operons.at[i, 'hits'] = operon
                    g = g + 1
            
            if g % 2000 == 0 and self.verbosity > 1:
                print g, (1.0 * g / self.n_hits * 100)
            
        # Drop empty rows
        operons = operons.ix[:i]
        
        # Drop operons without upstream promoter (fragmented genes)
        operons = operons[operons.promoter_start > 0]
        
        # Update index
        n_operons = len(operons.index)
        operons.index = range(n_operons)
        
        # Update instance variable
        self.operons = operons
        self.update_hits_in_operons()
        self.update_operon_stats()

        # Save operons
        #TODO: Save to database
        operons.to_csv(self.paths['operons'])
        
        if self.verbosity > 0:
            print 'Predicted operons and extracted promoter regions. [%.2fs]' % (clock() - t)
            
    def update_hits_in_operons(self):
        """ Adds information from the operon table to the hits table. """
        # Make sure hits are loaded
        if self.hits is None:
            self.load_hits()
        
        # Make sure operons are loaded
        if self.operons is None:
            self.load_operons()
            
        # Update hits
        for i, operon_hits in zip(self.operons.index, self.operons.hits.values):
            for pos, hit_idx in enumerate(operon_hits):
                self.hits.at[hit_idx, 'operon'] = i
                self.hits.at[hit_idx, 'operon_pos'] = pos
        
        # Update stats
        #self.update_operon_stats()
        
        # Save
        #TODO: Database
        self.hits.to_csv(self.paths['hits'])
        if self.verbosity > 0:
            print "Updated hits table with operons."
        
    def hits_in_operons(self):
        """ Returns the rows from the hits table that are contained in an operon. """
        # Make sure hits are loaded
        if self.hits is None:
            self.load_hits()
        
        # Make sure operons are loaded
        if self.operons is None:
            self.load_operons()
            
        return self.hits[~self.hits.operon.isnull()]
            
    def load_operons(self):
        """ Loads or creates operons table. """
        in_db = False
        if in_db:
            #TODO: Load from database
            pass
        elif self.check_path('operons'):
            # Load from file
            self.operons = pd.read_csv(self.paths['operons'])
            self.update_operon_stats()
            #self.update_hits_in_operons()
            
            if self.verbosity > 0:
                print 'Loaded operons from file.'
        else:
            # Create operons table
            self.predict_operons()
        
    def update_operon_stats(self):
        """ Updates stats based on operon table. """
        # Make sure operons are loaded
        if self.operons is None:
            self.load_operons()
        
        # Update operon statistics
        self.n_operons = self.operons.shape[0]
        self.operon_lengths = self.operons.hits.apply(len)
        
        # Update hits in operons
        self.n_scaffolds_in_operons = self.hits[~self.hits.operon.isnull()].scaffold.nunique()
        self.n_hits_in_operons = self.hits[~self.hits.operon.isnull()].shape[0]
        
    def dump_proteins(self, force=False):
        """ Saves protein sequences of hits with operons to FASTA file.
        This can be used as input for BLASTing.
        
        Format of operon_hits.faa is:
            >operon_id
            AA_SEQ
            ...
        """
        if not force and self.check_path('operon_hits'):
            print 'Skipping overwriting existing protein sequences.'
            return
        
        # Make sure genes are loaded
        if self.genes is None:
            self.load_genes()
            
        # Make sure operons have been predicted
        if self.operons is None:
            self.load_operons()
        
        if self.verbosity > 0:
            print 'Saving protein sequences for BLASTing...'; t = clock()
            
        # Save predicted amino acid sequences for genes in operons
        with open(self.paths['operon_hits'], 'w') as f:
            for idx in self.hits[~self.hits.operon.isnull()].index:
                f.write('>%d\n%s\n' % (idx, self.genes[self.hits.ix[idx].query_id].seq))
                
        if self.verbosity > 0:
            print 'Finished saving protein sequences. [%.2fs]' % (clock() - t)
            
    def blast_eggnog(self):
        """ Runs BLAST against the eggNOG database. """
        
        # Make sure operon hits have been dumped
        if not self.check_path('operon_hits'):
            self.dump_proteins()
        
        # Set up command
        quote = lambda x: "\"" + x + "\""
        blastp_cmd = NcbiblastpCommandline(query=quote(self.paths["operon_hits"]), \
                                            db=quote(eggnog_db_path), \
                                            out=quote(self.paths["eggnog_blast"]), \
                                            outfmt=6, \
                                            num_threads=8, \
                                            evalue=0.001)
                                            
        # Start timing the run
        if self.verbosity > 0:
            start = clock(); print "BLASTing %s (this may take a while)..." % self.mg_id
        
        # Execute
        blastp_cmd()
        
        # Finished
        time_elapsed = clock() - start
        if self.verbosity > 0:
            print "Done BLASTing %s. [%.2f min]" % (time_elapsed / 60.0)
            
    def parse_eggnog_blast(self):
        """ Parses the output of BLAST against the eggNOG database. """
        
        # Check if BLAST output exists
        if not self.check_path('eggnog_blast'):
            print 'BLAST file does not exist: ' + self.paths['eggnog_blast']
            return
        
        # Read in BLAST output into dataframe
        blast_cols = ['query_id', 'hit_id', 'percentage_identity', 'alignment_length', 'num_mismatches', 'num_gap_openings', 'q.start', 'q.end', 's.start', 's.end', 'e_value', 'score_bits']
        blast = pd.read_table(self.paths['eggnog_blast'], names=blast_cols)
        
        # Load membership tables
        NOG_members = pd.read_table(eggnog_NOG_members_path)
            #names=['cog', 'id', 'start_pos', 'end_pos'])
        bactNOG_members = pd.read_table(eggnog_bactNOG_members_path)
            #names=['cog', 'id', 'start_pos', 'end_pos'])
            
        # Load functional category tables
        NOG_funccat = pd.read_table(eggnog_NOG_funccat_path, \
            names=['cog', 'funccat'])
        bactNOG_funccat = pd.read_table(eggnog_bactNOG_funccat_path, \
            names=['cog', 'funccat'])
        
        # Create hit_nogs table
        hit_nogs = pd.DataFrame(index=blast.query_id.unique(), \
            columns=['tax_id', 'protein_id', 'e_val', 'NOG', 'NOG_funccat', 'bactNOG', 'bactNOG_funccat'])
        for hit in hit_nogs.index:
            identifier = blast.at[hit, 'hit_id']
            tax_id, protein_id = identifier.split('.')
            e_val = blast.at[hit, 'e_value']

            # BLAST hit
            hit_nogs.at[hit, 'tax_id'] = tax_id
            hit_nogs.at[hit, 'protein_id'] = protein_id
            hit_nogs.at[hit, 'e_val'] = e_val
            
            # NOG Membership
            # How do we do this fast?
            #hit_nogs.at[hit, 'NOG'] = NOG_members['#nog name'][NOG_members['protein name'] == mg.blast.ix[1, 'hit_id']]

            # Functional category

        # Update stats

        # Save
            
    def process(self):
        """ Runs all processing steps of the pipeline. """
        self.download_data()
        self.load_scaffolds()
        self.load_genes()
        self.load_hits()
        self.load_operons()
        self.dump_proteins()
        