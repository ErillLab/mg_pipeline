# -*- coding: utf-8 -*-
"""
Created on Wed May 20 15:00:44 2015

@author: Talmo
"""

import sys
sys.path.append("../..")
from igc_pipeline import *
import pandas as pd
import os
from glob import glob

samples = get_all_samples()
sample = samples[0]

original_paths = get_sample_paths(sample, original_paths=True)
new_paths = get_sample_paths(sample, original_paths=False)

if not os.path.exists(new_paths["scaftigs"]):
    if os.path.exists(original_paths["scaftigs"]):
        try:
            os.rename(original_paths["scaftigs"], new_paths["scaftigs"])
            log("Renamed assembled scaftig for sample %s." % sample, verbosity_level=2)
        except:
            log("Could not rename original assembled scaftig for sample %s." % sample, verbosity_level=1)
    else:
        log("Could not find original assembled scaftig for sample %s." % sample, verbosity_level=1)

if not os.path.exists(new_paths["ORFs_fasta"]):
    if os.path.exists(original_paths["ORFs_fasta"]):
        try:
            os.rename(original_paths["ORFs_fasta"], new_paths["ORFs_fasta"])
            log("Renamed ORF sequences for sample %s." % sample, verbosity_level=2)
        except:
            log("Could not rename original ORF sequences for sample %s." % sample, verbosity_level=1)
    else:
        log("Could not find original ORF sequences for sample %s." % sample, verbosity_level=1)
