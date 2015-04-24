# -*- coding: utf-8 -*-
"""
Created on Fri Apr 24 00:10:09 2015

@author: Talmo
"""

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import random
import sys
sys.path.append("../..")
from PSSMScorer import PSSMScorer, random_seq
from tqdm import *

pssm = PSSMScorer('/home/cuda/2TB/metagenomics/binding_sites/Firmicutes_LexA.txt')

#%% Single site curve
scores = np.arange(0, 25, 0.5)
probs = map(pssm.post_prob, scores)

plt.figure()
plt.plot(scores, probs)
plt.title("Single site scores\n%s" % pssm)
plt.xlabel("Site PSSM score")
plt.ylabel("Posterior probability of being bound")

#%% Motif sites
probs = map(pssm.post_prob, pssm.score_self(True))

plt.figure()
plt.title("PSSM motif sites\n%s" % pssm)
plt.boxplot(probs, labels=["Motif"])
plt.ylabel("Posterior probability")

#%% True versus random scores
n = 1000 # seqs (promoters)
l = 300 # promoter length

# Generate random background sequences
bg_seqs = [random_seq(l) for i in xrange(n)]

# Generate random promoters with sites sampled from the motif
make_segment = lambda x: random_seq(l / x - pssm.m) + str(random.choice(pssm.seqs))
make_fg = lambda n_fg_sites: "".join(make_segment(n_fg_sites) for i in range(n_fg_sites))
fg1_seqs = [make_fg(1) for i in xrange(n)]
fg2_seqs = [make_fg(2) for i in xrange(n)]
fg3_seqs = [make_fg(3) for i in xrange(n)]

# Compute scores
bg_scores = [pssm.score(seq, True) for seq in bg_seqs]
fg1_scores = [pssm.score(seq, True) for seq in fg1_seqs]
fg2_scores = [pssm.score(seq, True) for seq in fg2_seqs]
fg3_scores = [pssm.score(seq, True) for seq in fg3_seqs]

# Compute posteriors
bg_post = [pssm.post_prob(score) for score in bg_scores]
fg1_post = [pssm.post_prob(score) for score in fg1_scores]
fg2_post = [pssm.post_prob(score) for score in fg2_scores]
fg3_post = [pssm.post_prob(score) for score in fg3_scores]

plt.figure()
plt.boxplot([bg_post, fg1_post, fg2_post, fg3_post], labels=["Background only", "BG + 1 site", "BG + 2 sites", "BG + 3 sites"])
plt.title("Random promoters [%d bp]\n%s" % (l, pssm))
plt.xlabel("n = %d promoters per group" % n)
plt.ylabel("Posterior probability")
