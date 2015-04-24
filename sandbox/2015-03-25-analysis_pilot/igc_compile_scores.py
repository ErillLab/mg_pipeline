# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 11:26:48 2015

@author: Talmo
"""

from igc_pipeline import *

def compile_scores(n0=0, n1=10):
    samples = get_MetaHit()[n0:n1]
    PSSM = Firmicutes_LexA
    
    data_cols = ["COG", "sample", "operons", "sm_scores"]
    data = []
    for sample in samples:
        t = time()
        sample = get_unique_sample(sample)
        
        # Load sample data
        genes = load_sample_genes(sample)
        operons = predict_operons(sample)
        scores = get_sample_scores(sample, PSSM)
        SM = lambda o: np.log(np.sum(np.exp(np.hstack(scores.loc[o].values))))
        
        # Filter by genes with taxonomy and COGs
        valid_genes = genes[genes.eggNOG != "unknown"]
        #valid_genes = genes[(genes.phylum != "unknown") & (genes.eggNOG != "unknown")]
        valid_operons = operons[operons.genes.apply(lambda x: len(valid_genes.index.intersection(x)) > 0)]
        
        # Group by COGs and pre-allocate data
        COGs = valid_genes.groupby("eggNOG")
        n_COGs = COGs.ngroups + sum([k.count(";") for k in COGs.groups.keys()])
        sample_data = pd.DataFrame(index=range(n_COGs), columns=data_cols)
        sample_data["sample"] = sample
        i = 0
        
        # Filter by COG
        for COG, COG_genes in COGs:
            # Get operons containing genes
            cog_operons = valid_operons[valid_operons.genes.apply(lambda x: len(COG_genes.index.intersection(x)) > 0)]
            #log("%s: %d genes in %d operons" % (COG, len(COG_genes), len(cog_operons)))
            
            # Compute soft max of scores (one per operon)
            sm_scores = map(SM, cog_operons.index)
            
            # Handle multiple COGs, ex: "COG0142;NOG136421"
            for _COG in COG.split(";"):
                # Save data
                sample_data.at[i, "COG"] = _COG
                sample_data.at[i, "operons"] = cog_operons.index.tolist()
                sample_data.at[i, "sm_scores"] = sm_scores
                i += 1
        
        # Append data to other samples
        data.append(sample_data)
        log("Processed scores for %s." % sample, t)
        
    data = pd.concat(data, ignore_index=True)
    
    data.set_index("COG", inplace=True)

    data["n_genes"] = data.sm_scores.apply(len)
    data["mean_sm"] = data.sm_scores.apply(np.mean)
    data["min_sm"] = data.sm_scores.apply(np.min)
    data["max_sm"] = data.sm_scores.apply(np.max)
    data.sort("max_sm", ascending=False, inplace=True)

    return data

#%%
data = compile_scores(15, 20)
data.to_hdf("sm_scores-samples15-20-has_cog.h5", "table")