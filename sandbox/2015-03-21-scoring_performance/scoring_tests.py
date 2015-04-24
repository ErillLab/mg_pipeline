# -*- coding: utf-8 -*-
"""
Scoring benchmarks.

Created on Sat Mar 21 13:24:38 2015

@author: Talmo
"""

import os
import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))
from igc_pipeline import *

#%% Load data
sample = "MH0001"
operons = predict_operons(sample)

#%% Baseline (Biopython)
def biopython():
    t = time()
    scores = operons.promoter_seq.apply(lambda x: Firmicutes_LexA.score_all(x, True))
    log("[biopython]  %d seqs, %d bp" % (len(seqs), sum(seqs_pd.apply(len))), t)
    return scores

#%% Pandas
# Format data
seqs = operons.promoter_seq.str.upper()
#pssm = pd.DataFrame(Firmicutes_LexA.pssm) #Fast: pssm.at[pos, base]
pssm = dict(Firmicutes_LexA.pssm)
pssm_r = dict(Firmicutes_LexA.pssm_r)

def calculate(seq):
    m = len(pssm)
    n = len(seq)
    scores = np.zeros(n-m+1)
    scores_r = np.zeros(n-m+1)
    for i in xrange(n-m+1):
        for pos in xrange(m):
            scores[i] += pssm[seq[i+pos]][pos]
            scores_r[i] += pssm_r[seq[i+pos]][pos]
    return scores, scores_r

t = time()
#scored = seqs.apply(calculate)
scored = seqs.apply(calculate).apply(pd.Series, args=([["+","-"]]))
#scored = seqs.apply(lambda x: pd.Series(calculate(x), index=["+","-"]))
log("[pandas] %d seqs, %d bp" % (len(seqs), sum(seqs_pd.apply(len))), t)

#%% NumPy
# Convert sequences to matrix
seqs = operons.promoter_seq.str.upper().as_matrix()

# NumPy
keys = np.array(Firmicutes_LexA.pssm.keys())
pssm = np.array(Firmicutes_LexA.pssm.values())
pssm_r = np.array(Firmicutes_LexA.pssm_r.values())