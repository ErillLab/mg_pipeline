# -*- coding: utf-8 -*-
"""
This class provides a simple wrapper for Biopython's PSSM scoring functions.

Created on Wed Feb 11 11:05:08 2015

@author: Talmo
"""

from Bio import Seq, SeqRecord, Alphabet, motifs
import numpy as np

class PSSMScorer:
    def __init__(self, binding_sites, pseudocounts=1):
        """ Creates a PSSM scorer object. Accepts binding sites in the form of
        a path to a text file containing one site per line, or a list of 
        Biopython Bio.Seq objects."""
        
        self.alphabet = Alphabet.IUPAC.unambiguous_dna
        
        if type(binding_sites) == str:
            self.seqs = [Seq.Seq(site.strip(), self.alphabet) for site in open(binding_sites).readlines()]
        elif type(binding_sites) == list:
            self.seqs = binding_sites
            
        # Construct motif
        self.motif = motifs.create(self.seqs)
        self.motif.pseudocounts = pseudocounts
        
        # Construct PSSM and reverse PSSM
        self.pssm = self.motif.pssm
        self.pssm_r = self.pssm.reverse_complement()
        
    def convert_seq(self, seq):
        """ Converts sequence strings to Biopython Seq objects with appropriate
        alphabet. """
        
        if type(seq) == str:
            return Seq.Seq(seq.strip(), self.alphabet)
        elif type(seq) == SeqRecord.SeqRecord:
            return self.convert_seq(seq.seq)
        elif type(seq) == Seq.Seq and seq.alphabet != self.alphabet:
            seq.alphabet = self.alphabet
            return seq
        else:
            return seq
        
    def score(self, seq):
        """ Scores a sequence and returns the best score between the forward
        and reverse strand. """
        seq = self.convert_seq(seq)
        return max(self.pssm.calculate(seq), self.pssm_r.calculate(seq))
        
    def score_all(self, seq, strand=False):
        """ Scores all sites in a sequence and returns an array of scores. """
        if strand:
            return self.search(seq, -np.inf)
        else:
            return self.search(seq, -np.inf)[:, 1]
        
    def search(self, seq, threshold=0.0):
        """ Search for the sites in the sequence with a score above a threshold.
        Searches on both strands."""
        return np.array(list(self.pssm.search(self.convert_seq(seq), both=True, threshold=threshold)))

    
