from motifs import Escherichia_coli
from utils import make_pssm,score_seq,mean,sd,random_site,score_genome,wc
from collections import Counter
import scipy
from scipy import stats
import random
from math import exp,log

def prod(xs):
    return reduce(lambda x,y:x*y,xs)

alpha = 1.0/300
pf = 1/100.0
pb = 1 - pf
pssm = make_pssm(Escherichia_coli.LexA)
pssm_scores = [log(exp(score_seq(pssm,seq)) + exp(score_seq(pssm,seq))) for seq in Escherichia_coli.LexA]
background_scores = [(lambda seq:log(exp(score_seq(pssm,seq)) + exp(score_seq(pssm,seq))))(random_site(16))
                     for i in range(10000)]
mu_y,sigma_y = mean(pssm_scores),sd(pssm_scores)
mu_x,sigma_x = mean(background_scores),sd(background_scores) # can be read off pssm but I'm lazy
pdf_y = scipy.stats.distributions.norm(mu_y,sigma_y).pdf
pdf_x = scipy.stats.distributions.norm(mu_x,sigma_x).pdf
L_b = lambda xs:pdf_x(xs)
L_f = lambda z:alpha*pdf_y(z) + (1-alpha)*pdf_x(z)
Ls_b = lambda xs:exp(sum(log(L_b(x)) for x in xs))
Ls_f = lambda xs:exp(sum(log(L_f(x)) for x in xs))
LL_ratio = lambda xs:exp(sum(log(L_b(x))-log(L_f(x)) for x in xs))

def post_prob_ref(scores):
    """included only for reference purposes"""
    return Ls_f(scores)*pf/(Ls_f(scores)*pf + Ls_b(scores)*pb)
    

def post_prob(scores):
    """Compute posterior probability of score array"""
    return 1/(1+LL_ratio(scores)*pb/pf)
    
def sanity_check(background_num,site_num):
    background_scores = score_genome(pssm,random_site(background_num))
    site_scores = [score_seq(pssm,random.choice(Escherichia_coli.LexA)) for i in range(site_num)]
    return background_scores + site_scores
