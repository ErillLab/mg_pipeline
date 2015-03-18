# -*- coding: utf-8 -*-
"""


Created on Tue Feb 10 01:05:51 2015

@author: Talmo
"""
from Metagenome import Metagenome
from PSSMScorer import PSSMScorer
import pandas as pd

# Metagenomes
metagenomes = []

# UMC Utrecht, NL, ICU patient gut microbiome
metagenomes.append('mgm4508945.3')  # Day 14 ICU
metagenomes.append('mgm4508946.3')  # Day 16 ICU
metagenomes.append('mgm4508947.3')  # Day 28 neurology ward

# Gut microbiota in the Irish Elderly and its links to health and diet
metagenomes.append('mgm4491401.3') # EM338_DH_DG2_3A
metagenomes.append('mgm4491410.3') # EM283_C_DG2_3A #
metagenomes.append('mgm4491412.3') # EM272_DH_DG2_3A
metagenomes.append('mgm4491413.3') # EM268_C_DG2_1B
metagenomes.append('mgm4491417.3') # EM232_R_DG1_2A #
metagenomes.append('mgm4491421.3') # EM209_R_NA_NA
metagenomes.append('mgm4491423.3') # EM205_C_DG1_2A
metagenomes.append('mgm4491477.3') # EM203_C_DG1_3A
metagenomes.append('mgm4491482.3') # EM177_C_DG2_1B
metagenomes.append('mgm4491487.3') # EM172_C_DG2_1A
metagenomes.append('mgm4491562.3') # EM039_C_NA_NA

# Waseca farm
metagenomes.append('mgm4441091.3')

# Phosphorus Removing (EBPR) Sludge Community
metagenomes.append('mgm4441092.3')

# Acid Mine Drainage Metagenome
metagenomes.append('mgm4441137.3') # 5-Way (CG) Acid Mine Drainage Biofilm
metagenomes.append('mgm4441138.3') # UBA Acid Mine Drainage Biofilm


#%% Process metagenomes
for mg_id in metagenomes:
    mg = Metagenome(mg_id)
    mg.process()

#%% Gene and NOG counts
n_genes = {}; n_NOGs = {}
for mg_id in metagenomes:
    mg = Metagenome(mg_id, verbosity=0)
    if mg.check_path('eggnog_blast'):
        mg.load_NOGs()
        print mg
        print "  Genes: %d / NOGs: %d" % (mg.NOGs.gene.nunique(), mg.NOGs.nog.nunique())
        n_genes[mg_id] = mg.NOGs.gene.nunique()
        n_NOGs[mg_id] = mg.NOGs.nog.nunique()
n_genes = pd.DataFrame(n_genes.values(), index=n_genes.keys(), columns=['n_genes'])
n_NOGs = pd.DataFrame(n_NOGs.values(), index=n_NOGs.keys(), columns=['n_NOGs'])
stats = n_genes.join(n_NOGs)

#%% Count TF NOGs
NOGs = ["COG0468", "COG1974"]

NOG_counts = {NOG: {} for NOG in NOGs}
for mg_id in metagenomes:
    mg = Metagenome(mg_id, verbosity=0)
    if mg.check_path('eggnog_blast'):
        mg.load_NOGs()
        print mg
        for NOG in NOGs:
            n = sum(mg.NOGs.nog == NOG)
            print "  %s: %d" % (NOG, n)
            NOG_counts[NOG][mg.id] = n
NOG_counts = pd.DataFrame(NOG_counts)
stats = stats.join(NOG_counts)

#%% PSSM scores
lexA_COG = "COG1974"

# Create PSSM
base_path = "E:\\metagenomics\\binding_sites\\"
# "lexA.Actinobacteria.txt"
sites = {"lexA": base_path + "lexA.Alphaproteobacteria.txt"}
lexA = PSSMScorer(sites["lexA"])

# Score sites
mg = Metagenome("mgm4491401.3")
mg.process()

# Get COG genes
lexA_NOGs = mg.NOGs[mg.NOGs.nog == lexA_COG].gene
lexA_genes = mg.hits.ix[lexA_NOGs][["operon", "operon_pos"]]

# Get corresponding operons
lexA_operons = mg.operons.ix[lexA_genes.operon]

# Get promoter region for those operons
lexA_promoters = lexA_operons.promoter_seq

# Score
