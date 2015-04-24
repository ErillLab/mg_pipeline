# -*- coding: utf-8 -*-
"""
Preliminary analysis routines

Created on Wed Mar 25 16:22:26 2015

@author: Talmo
"""

from igc_pipeline import *
import seaborn as sns

#%% Processing
#sample = "MH0001"
sample = "O2.UC44-2"
PSSM = Firmicutes_LexA
sample = get_unique_sample(sample)
    
# Load sample data
genes = load_sample_genes(sample)
ORFs = parse_ORFs(sample)
operons = predict_operons(sample)
scores = get_sample_scores(sample, PSSM)

# Filter by genes with taxonomy and COGs
valid_genes = genes[(genes.phylum != "unknown") & (genes.eggNOG != "unknown")]
valid_operons = operons[operons.genes.apply(lambda x: len(valid_genes.index.intersection(x)) > 0)]
valid_scores = scores.loc[valid_operons.index]

# Compute soft max
soft_max = valid_scores.applymap(np.exp).sum(1).map(np.log)
soft_max1d = np.hstack(soft_max)

#%% Top hits
# Get max score of each operon
max_scores = soft_max.apply(max)
max_scores.sort(ascending=False)

# Get top hits
top_hits = max_scores.head().index
top_operons = operons.loc[top_hits]

# Get genes from the operon with the best score
best_operon = top_operons.iloc[0]
operon_genes = genes.loc[best_operon.genes].dropna()

#%% Unpack hits
thresh = 10.0
o = 11453 # operon

hit_scores = soft_max.loc[o][soft_max.loc[o] > thresh]
has_hits = hit_scores.any()

hit_genes = operons.loc[o, "genes"]

cols = "gene", "score", "operon", "operon_pos", "sample", "COG", "phylum", "genus"
hits = pd.DataFrame(columns=cols)

hits["gene"] = hit_genes
hits["score"] = hit_scores.max()
hits["operon"] = o
hits["operon_pos"] = range(len(hit_genes))
hits["sample"] = sample
gene_data = genes.loc[hit_genes, ["eggNOG", "phylum", "genus"]]
gene_data[gene_data == "unknown"] = np.NaN
hits[["COG", "phylum", "genus"]] = gene_data.values

print hits

#%% Visualize operons
ops = top_operons

from matplotlib import collections
fig, axes = plt.subplots(len(ops))
for i in xrange(len(ops)):
    operon = ops.iloc[i]
    ax = axes[i]
    starts = ORFs.loc[operon.genes].start
    ends = ORFs.loc[operon.genes].end
    
    y = 0
    pts = zip([(x, y) for x in starts], [(x, y) for x in ends])

    c = sns.color_palette()
    lc = collections.LineCollection(pts, colors=c, linewidths=20)
    
    ax.add_collection(lc)
    ax.autoscale()
    ax.yaxis.set_major_locator(plt.NullLocator())
    
plt.xlabel("Position in scaftig (bp)")


#%% Full soft max idea
hit = max_scores.index[0]
control = max_scores.index[-1]

o = hit
both_strands = np.hstack(scores.loc[o].values)
soft_maxed = np.log(np.sum(np.exp(both_strands)))

SM = lambda o: np.log(np.sum(np.exp(np.hstack(scores.loc[o].values))))

# All operons sorted by their max scores
sm_scores_sorted = map(SM, max_scores.index)

plt.plot(max_scores, b)
plt.xlabel("Max operon score"); plt.ylabel("Soft max of all operon scores")

#%% Full soft max on more data
#%%prun -T "test"
lexA_COG = "COG1974"
pcna_COG = "COG0592" # DNA polymerase sliding clamp subunit (PCNA homolog)
control_COG = "COG1803" # Methylglyoxal synthase

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
#_data = data

#%%
data = compile_scores(0, 5)
data.to_hdf("sm_scores-samples0-5-has_cog.h5", "table")

#%%
filenames = ["sm_scores-samples0-10.h5", "sm_scores-samples30-40-has_cog.h5", "sm_scores-samples40-50-has_cog.h5"]

data = pd.concat((pd.read_hdf(f, "table") for f in filenames))
data.reset_index(inplace=True)

#%%
g = data.groupby("COG")
threshold = 12.0
num_scores_per_cog = g.apply(lambda x0: sum(x0.sm_scores.apply(lambda x2: sum([x1 > threshold for x1 in x2]))))
num_scores_per_cog.sort(ascending=False)

# COGs of interest (Table S6)
COIs = ["COG0389", "COG0468", "COG0556", "COG1974", "COG0210", "COG2001", "COG0732", "COG0178", "COG0187", "COG4974", "COG0632", "COG0497", "COG0399", "COG0653", "COG1136", "COG0463", "COG1609", "COG0745", "COG0582"]
print num_scores_per_cog[COIs]