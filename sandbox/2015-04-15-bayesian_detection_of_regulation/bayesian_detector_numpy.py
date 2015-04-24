# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 21:15:51 2015

@author: Talmo
"""
import numpy as np
import scipy
import random
import sys
sys.path.append("../..")
from PSSMScorer import PSSMScorer

pssm = PSSMScorer('/home/cuda/2TB/metagenomics/binding_sites/Firmicutes_LexA.txt')

# Params
alpha = 1.0/300 # frequency of lexA?
pf = 1/100.0 # foreground probability
pb = 1 - pf # background probability
num_random = 10000 # number of random sequences to score to estimate background

# Scores
pssm_scores = pssm.score_self(True)
random_site = lambda l: "".join([random.choice("ACGT") for i in xrange(l)])
background_scores = [pssm.score(random_site(pssm.length), True) for i in xrange(num_random)]

# Distributions
mu_y, sigma_y = np.mean(pssm_scores), np.std(pssm_scores)
mu_x, sigma_x = np.mean(background_scores), np.std(background_scores)
pdf_y = scipy.stats.distributions.norm(mu_y, sigma_y).pdf
pdf_x = scipy.stats.distributions.norm(mu_x, sigma_x).pdf

L_b = lambda scores: pdf_x(scores)
L_f = lambda scores: alpha * pdf_y(scores) + (1 - alpha) * pdf_x(scores)
Ls_b = lambda scores: np.exp(np.sum(np.log(L_b(scores))))
Ls_f = lambda scores: np.exp(np.sum(np.log(L_f(scores))))

LL_ratio = lambda scores: np.exp(np.sum(np.log(L_b(scores)) - np.log(L_f(scores))))

def post_prob(scores):
    """Compute posterior probability of score array"""
    return 1/(1 + LL_ratio(scores) * pb / pf)