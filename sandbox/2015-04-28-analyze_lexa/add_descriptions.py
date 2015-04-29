# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 12:21:01 2015

@author: Talmo
"""
import pandas as pd

#%% Paths
eggnog_path = "/home/cuda/2TB/eggNOG/"
members_path = eggnog_path + "NOG.members.txt"
funccats_path = eggnog_path + "NOG.funccat.txt"
descriptions_path = eggnog_path + "NOG.description.txt"

#%% Load data
descriptions = pd.read_csv(descriptions_path, sep="\t", names=["eggNOG", "description"], index_col="eggNOG")
funccats = pd.read_csv(funccats_path, sep="\t", names=["eggNOG", "funccat"], index_col="eggNOG")

Firmicutes_probs = pd.read_csv("Firmicutes_probs_0-1000_tax-filtered_head-filtered.csv", index_col="eggNOG")
Gamma_probs = pd.read_csv("Gamma_probs_0-1000_tax-filtered_head-filtered.csv", index_col="eggNOG")

#%% Add annotations
Firmicutes_probs = Firmicutes_probs.join(funccats)
Gamma_probs = Gamma_probs.join(funccats)

Firmicutes_probs = Firmicutes_probs.join(descriptions)
Gamma_probs = Gamma_probs.join(descriptions)

#%% Save
Firmicutes_probs.to_csv("Firmicutes_probs_0-1000_tax-filtered_head-filtered_annotated.csv")
Gamma_probs.to_csv("Gamma_probs_0-1000_tax-filtered_head-filtered_annotated.csv")
