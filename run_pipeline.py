# -*- coding: utf-8 -*-
"""


Created on Tue Feb 10 01:05:51 2015

@author: Talmo
"""
from Metagenome import Metagenome

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


# Process metagenome
#mg = Metagenome(metagenomes[3])
#mg.process()

for mg_id in metagenomes:
    mg = Metagenome(mg_id)
    mg.process()
    