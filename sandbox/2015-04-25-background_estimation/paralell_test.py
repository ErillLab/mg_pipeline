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
pssm = Firmicutes_LexA
samples = get_all_with_scores(pssm)

#%% Parallel test

def f(sample):
    n = np.float64(0)
    total = np.float64(0)

    operons = get_operons(sample)
    valid_operons = operons.index[operons.head_completeness != "Lack 5'-end"]
    scores = np.hstack(get_sample_scores(sample, pssm, soft_max=True).loc[valid_operons])
    total += np.sum(scores)
    n += len(scores)

    return total, n

def naive(k):
    result = map(f, samples[0:k])
    totals = np.sum(result, axis=0)
    mu = totals[0] / totals[1]
    return mu
    
def par(k):
    p = Pool(8)
    result = p.map(f, samples[0:k], 1)
    totals = np.sum(result, axis=0)
    p.close()
    mu = totals[0] / totals[1]
    return mu

def par_async(k):
    p = Pool(8)
    result = p.map(f, samples[0:k], 1)
    totals = np.sum(result, axis=0)
    p.close()
    mu = totals[0] / totals[1]
    return mu

#%time naive(8)
#%time par(8)
#%time par_async(8)

