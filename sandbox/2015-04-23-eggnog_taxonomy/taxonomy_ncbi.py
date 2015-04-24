import pandas as pd
from Bio import Entrez
from tqdm import tqdm
Entrez.email = "erill@umbc.edu"

def parse_taxfile(file_name, column_number=0):
        """
            This function takes in a file name and returns a pandas dataframe for all taxonomy. Make sure the filename and directory is correct. 
            param: 
                    file_name   
                        The name of the file
                    column_number
                        The column number that contains the tax id
                        
            returns:
                pandas dataframe where the columns are: tax_id, kingdom, phylum, class, subclass, order, family, genus, species
    
            example:
                pandas_object = parse_taxfile("eggnogv4.species.txt",1)
        """
        columns = ["tax_id", "superkingdom", "kingdom", "phylum", "class", "subclass", "order", "family", "genus", "species"]
        taxonomy_db = pd.DataFrame(columns=columns)
        for line in tqdm(open(file_name).readlines()):
            try:
                if "#" not in line:
                    line_attr = line.split("\t")
                    tax_dict = get_taxonomy(line_attr[column_number])
                    taxonomy_db = taxonomy_db.append(pd.DataFrame(tax_dict,index=[0]), ignore_index=True)
            except e:
                print e
                return taxonomy_db
        return taxonomy_db

def get_taxonomy(tax_id):
    """
        This fuction takes in a tax id and returns a dictionary containing the taxonomic information.
        In this dictionary the keys are the taxonomic heirachy (kingdom, superkingdom, phylum, etc.)

        params: tax_id
            NCBI id for the taxonomy database

        returns:
            a dictionary containing the taxonomy heirachy

        example:
            tax_dict = get_taxonomy("8883")
    """
    tax_dict = {}
    handle = Entrez.efetch("taxonomy", id=tax_id)
    tax_dict.update({'tax_id':tax_id})
    records = Entrez.parse(handle)
    for record in records:
        for tax_class in record['LineageEx']:
            if tax_class['Rank'] != 'no rank':
                tax_dict.update({tax_class['Rank']:tax_class['ScientificName']})
    return tax_dict

eggnog_taxonomy = parse_taxfile("eggnogv4.species.txt",1)
eggnog_taxonomy.to_hdf("eggnogv4_taxonomy.h5", "table")
eggnog_taxonomy.to_csv("eggnogv4_taxonomy.csv")