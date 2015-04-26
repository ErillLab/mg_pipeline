# -*- coding: utf-8 -*-
"""
Created on Sat Apr 25 20:45:56 2015

@author: Talmo
"""

import sys
sys.path.append("../..")
from igc_pipeline import *
from PSSMScorer import *
from tqdm import tqdm
import numpy as np
from multiprocessing import Pool

#%% Configuration
pssm = GammaProteobacteria_LexA
samples = get_all_with_scores(pssm)
n_processes = 8

#%% Mean
t = time()
def accum_mu(sample):
    operons = get_operons(sample)
    valid_operons = operons.index[operons.head_completeness != "Lack 5'-end"]
    scores = np.hstack(get_sample_scores(sample, pssm, soft_max=True).loc[valid_operons])
    total = np.sum(scores)
    n = np.float64(len(scores))
    return total, n
    
p = Pool(n_processes)
result = p.map(accum_mu, samples, 1)
totals = np.sum(result, axis=0)
p.close()
mu = totals[0] / totals[1]
print "mu = %f [%.2fs | %.2fs / sample]" % (mu, time() - t, (time() - t) / len(samples))

#%% Standard deviation
t = time()
def accum_sigma(sample):
    operons = get_operons(sample)
    valid_operons = operons.index[operons.head_completeness != "Lack 5'-end"]
    scores = np.hstack(get_sample_scores(sample, pssm, soft_max=True).loc[valid_operons])
    total = np.sum(abs(scores - mu) ** 2)
    n = np.float64(len(scores))
    return total, n

p = Pool(n_processes)
result = p.map(accum_sigma, samples, 1)
totals = np.sum(result, axis=0)
p.close()
sigma = np.sqrt(totals[0] / totals[1])
print "sigma = %f [%.2fs | %.2fs / sample]" % (sigma, time() - t, (time() - t) / len(samples))

#%% Results:
#Firmicutes_LexA [16 bp | 79 seqs]
#mu = -17.918493763638413
#sigma = 8.2211415419612841
#
# 