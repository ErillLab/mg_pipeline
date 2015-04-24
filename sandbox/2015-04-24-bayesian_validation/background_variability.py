# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 21:49:40 2015

@author: Talmo
"""
# Test variability of background
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.append("../..")
from PSSMScorer import PSSMScorer
from tqdm import *

pssm = PSSMScorer('/home/cuda/2TB/metagenomics/binding_sites/Firmicutes_LexA.txt')

scores = 15
#n = np.repeat([1e3, 1e4, 1e5], 5)
#n = np.repeat(np.linspace(1e4, 1e6, 10), 5)
n = np.repeat(np.logspace(4, 6, 100), 10)
p = np.zeros(n.shape)
for i in trange(len(n)):
    pssm.initialize_estimator(n[i])
    p[i] = pssm.post_prob(scores)
    
# Plot
plt.figure()
plt.plot(n, p, 'o')
plt.title("Background variability (%s)" % pssm.name)
plt.xlabel("# random seqs generated to estimate background")
plt.ylabel("Posterior probability of score = 15")
