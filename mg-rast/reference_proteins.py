# -*- coding: utf-8 -*-
"""
This script is used to BLAST reference genes to eggNOG.

Created on Tue Feb 10 10:37:07 2015

@author: Talmo
"""

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastpCommandline
from time import clock
import pandas as pd

#%% Blast reference genes
# Path to BLAST database for eggNOG
eggnog_db_path = "E:\\eggNOG\\blast_db\\eggnogv4"

# Amino acid FASTA file with reference proteins
ref_path = "E:\\metagenomics\\proteins\\lexA.Alphaproteobacteria.faa"
ref_faa = list(SeqIO.parse(ref_path, 'fasta'))

# Build BLAST command
output_path = ref_path[:-3] + "tbl"
blastp_cmd = NcbiblastpCommandline(query=ref_path, \
                                    db=eggnog_db_path, \
                                    out=output_path, \
                                    outfmt=6, \
                                    num_threads=8, \
                                    evalue=0.0001)
# Run the command
t = clock()
blastp_cmd()
print "Finished BLASTing. [%.2fs]" % (clock() - t)

#%% Parse output
# Parse BLAST
blast_cols = ['query_id', 'hit_id', 'percentage_identity', 'alignment_length', 'num_mismatches', 'num_gap_openings', 'q.start', 'q.end', 's.start', 's.end', 'e_value', 'score_bits']
blast = pd.read_table(output_path, names=blast_cols)
top_hits = blast.groupby('query_id').first()

# Load eggNOG reference tables
members_path = "E:\\eggNOG\\all.members\\aproNOG.members.txt"
funccat_path = "E:\\eggNOG\\all.funccat\\aproNOG.funccat.txt"
members = pd.read_table(members_path)
funccat = pd.read_table(funccat_path)

#%% Search for protein name
for accession, protein_name in top_hits.hit_id.iteritems():
    NOG = members[members['protein name'] == protein_name]
    if NOG.shape[0] > 0:
        print accession, "==>", protein_name, "==>", NOG.reset_index().ix[0, 1]


#%% BLAST LexA repressor [Agrobacterium fabrum str. C58] (NP_354402.1)
ref_path = "E:\\metagenomics\\proteins\\NP_354402.1.faa"

# Build BLAST command
eggnog_db_path = "E:\\eggNOG\\blast_db\\eggnogv4"
output_path = ref_path[:-3] + "tbl"
blastp_cmd = NcbiblastpCommandline(query=ref_path, \
                                    db=eggnog_db_path, \
                                    out=output_path, \
                                    outfmt=6, \
                                    num_threads=8, \
                                    evalue=0.0001)
# Run the command
t = clock()
blastp_cmd()

# Load BLAST results
blast_cols = ['gene', 'protein_name', 'percentage_identity', 'alignment_length', 'num_mismatches', 'num_gap_openings', 'q.start', 'q.end', 's.start', 's.end', 'e_value', 'score_bits']
blast = pd.read_table(output_path, names=blast_cols, usecols=['gene', 'protein_name', 'e_value'])

print "Finished BLASTing. [%.2fs]" % (clock() - t)

#%% RecA protein [Agrobacterium fabrum str. C58] (NP_354855.2)
ref_path = "E:\\metagenomics\\proteins\\NP_354855.2.faa"

# Build BLAST command
eggnog_db_path = "E:\\eggNOG\\blast_db\\eggnogv4"
output_path = ref_path[:-3] + "tbl"
blastp_cmd = NcbiblastpCommandline(query=ref_path, \
                                    db=eggnog_db_path, \
                                    out=output_path, \
                                    outfmt=6, \
                                    num_threads=8, \
                                    evalue=0.0001)
# Run the command
t = clock()
blastp_cmd()

# Load BLAST results
blast_cols = ['gene', 'protein_name', 'percentage_identity', 'alignment_length', 'num_mismatches', 'num_gap_openings', 'q.start', 'q.end', 's.start', 's.end', 'e_value', 'score_bits']
blast = pd.read_table(output_path, names=blast_cols, usecols=['gene', 'protein_name', 'e_value'])

print "Finished BLASTing. [%.2fs]" % (clock() - t)

#%% Parse output
# Reference paths
members_path = "E:\\eggNOG\\all.members\\bactNOG.members.txt"
funccat_path = "E:\\eggNOG\\all.funccat\\bactNOG.funccat.txt"
#members_path = "E:\\eggNOG\\all.members\\aproNOG.members.txt"
#funccat_path = "E:\\eggNOG\\all.funccat\\aproNOG.funccat.txt"
#members_path = "E:\\eggNOG\\NOG.members.txt"
#funccat_path = "E:\\eggNOG\\NOG.funccat.txt"

# Load eggNOG reference tables
members = pd.read_table(members_path, names=['nog', 'protein_name', 'start_pos', 'end_pos'], usecols=['nog', 'protein_name'], skiprows=1)
funccat = pd.read_table(funccat_path, names=['nog', 'funccat'])

# Add NOG and functional categories
NOGs = pd.merge(blast, members, how='inner', left_on='protein_name', right_on='protein_name')
        
# Get the functional categories
NOGs = pd.merge(NOGs, funccat, how='left', left_on='nog', right_on='nog')

# Output NOG distribution
#print NOGs.groupby('nog').describe()
print NOGs.groupby('nog').count()

# LexA
#bactNOG
#bactNOG04582     287
#bactNOG108766      2
#bactNOG247172      2
#aproNOG
#aproNOG01160    158
#NOG
#COG1974    297

# RecA
#bactNOG
#nog                                               
#bactNOG00450   305
#aproNOG                                               
#aproNOG01143   172
#COG                                          
#COG0468   305




