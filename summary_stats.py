import pandas as pd
from tqdm import tqdm
from glob import glob
table_summary = {}
LexA_COG = "COG1974"

def fill_table(table, keep_unknowns=False):
    """
    This method is in charge of filling the global variable called table_summary.
    Table summary is a dictionary that clusters the sample data into specific categories.
    The map is below: (K stands for K and V stands for variable)
                                  - <K> N (Number of Samples)
                                |
                                  - <K> Phylum (Pandas dataframe that holds all the phylum data)
                        <K>EUR  |
    <V>table_summary ->           - <K> Genus (Pandas dataframe that holds all the genus data)
                        <K>CHN  |
                                  - <K> eggNog (Pandas dataframe that holds all the cog data)

    params:
        table 
            A pandas data frame table that must have the following headers: 
            cohort_origin, genus, phylum, eggNOG (order does not matter)
    """
    if not keep_unknowns:
        table.drop(table.index[table.phylum == "unknown"], inplace=True)
        table.drop(table.index[table.genus == "unknown"], inplace=True)
        table.drop(table.index[table.eggNOG == "unknown"], inplace=True)
    for cohort in table["cohort_origin"].unique():
        cohort_table = table.loc[table["cohort_origin"] == cohort]

        if cohort not in table_summary:
            table_summary[cohort] = {}

        if "N" not in table_summary[cohort]:
            table_summary[cohort]["N"] = 1
        else:
            table_summary[cohort]["N"] += 1

        if "phylum" not in table_summary[cohort]:
            table_summary[cohort]["phylum"] = pd.DataFrame(columns=["phylum"])

        #count the number of unique values in the phylum column
        table_phylum = cohort_table.loc[table['genus'] != "unknown",['phylum']].apply(pd.Series.value_counts,axis=0)
        table_summary[cohort]["phylum"] = table_summary[cohort]["phylum"].add(table_phylum,fill_value=0)
    
        if "genus" not in table_summary[cohort]:
            table_summary[cohort]["genus"] = pd.DataFrame(columns=["genus"])

        #count the number of unique values in the genus column
        table_genus = cohort_table.loc[table['genus'] != "unknown",['genus']].apply(pd.Series.value_counts,axis=0)
        table_summary[cohort]["genus"] = table_summary[cohort]["genus"].add(table_genus,fill_value=0)
   
        if "eggNOG" not in table_summary[cohort]:
            table_summary[cohort]["eggNOG"] = pd.DataFrame(columns=["eggNOG"])

        index = 0
        cogs_table = {}

        #forced to do this by each row because one gene belongs to multiple cogs
        for row in cohort_table["eggNOG"]:
            if "gene_cog" not in table_summary[cohort]:
                table_summary[cohort]["gene_cog"] = 0
            table_summary[cohort]["gene_cog"] += 1
            if ";" in row:
                cogs = row.split(";")
                for cog in cogs:
                    if cog not in cogs_table:
                        cogs_table[cog] = 0
                    cogs_table[cog] += 1
            else:
                if row not in cogs_table:
                    cogs_table[row] = 0
                cogs_table[row] += 1

        #create data frame and create the table
        cog_pd = pd.DataFrame.from_dict(cogs_table, orient="index")
        cog_pd.columns = ["eggNOG"]
        table_summary[cohort]["eggNOG"] = table_summary[cohort]["eggNOG"].add(cog_pd,fill_value=0)
 
def summarize(data, verbal=False, using_files=True):
    """
    This function calls the function above which is titled fill_table and prints out
    the summary data. This function can either take in a list of files or it can take in a list of pandas dataframes

    params: 
        data 
            An ambigious variable that can either be a list of files or a list of panda dataframes
        verbal
            Print out the summary table True for print False for don't print
        using_files
            Boolean variable to tell the fucntion if you plan to use a list of file name or a list of panda dataframes
            True for file names False for list of pandas dataframes
    """

    if using_files:
        for file_name in tqdm(data):
            fill_table(pd.read_csv(file_name))
    else:
        for table in tqdm(data):
            fill_table(table)

    for cluster in table_summary:
        #total_genes = sum(table_summary[cluster]["phylum"].values) # number of genes
        #total_genes = table_summary[cluster]["N"] # number of samples
        total_genes = table_summary[cluster]["eggNOG"].eggNOG.sum() # number of genes in COGs with duplicates
        
        phylum_percent = table_summary[cluster]["phylum"].apply(lambda x: x/total_genes * 100)
        phylum_percent.columns = ["percent"]
        table_summary[cluster]["phylum"] = pd.concat([table_summary[cluster]["phylum"],phylum_percent],axis=1)

        #Read above for fix
        genus_percent = table_summary[cluster]["genus"].apply(lambda x: x/total_genes * 100)
        genus_percent.columns = ["percent"]
        table_summary[cluster]["genus"] = pd.concat([table_summary[cluster]["genus"],genus_percent],axis=1)

        #read above for fix
        cog_percent = table_summary[cluster]["eggNOG"].apply(lambda x: x/table_summary[cluster]["gene_cog"] * 100)
        cog_percent.columns = ["percent"]
        table_summary[cluster]["eggNOG"] = pd.concat([table_summary[cluster]["eggNOG"],cog_percent],axis=1)

        #Print the data
        if verbal:
            print "Cluster %s:\n" % cluster
            print "Number of Samples: %d\n" % table_summary[cluster]["N"]
            print "Taxonomy:"
            print table_summary[cluster]["phylum"].sort("percent", ascending=False)
            print "----------------------------------"
            print table_summary[cluster]["genus"].sort("percent", ascending=False)
            print "-----------------------------------"
            print "COGS:"
            print table_summary[cluster]["eggNOG"].sort("percent", ascending=False)
            print "------------------------------------"
            print "End Summary"

#summarize(glob.glob("~/2TB/metagenomics/IGC/Genes/*.csv"), True)
genes_path = '/home/cuda/2TB/metagenomics/IGC/Genes/'
summarize(glob(genes_path + "*.csv"), True)
print "LexA (%s):" % LexA_COG
for cluster in table_summary.keys():
    print "%s:" % cluster,
    lexA_counts = table_summary[cluster]["eggNOG"].loc[LexA_COG]
    print "%d (%f%%)" % (lexA_counts.eggNOG, lexA_counts.percent)