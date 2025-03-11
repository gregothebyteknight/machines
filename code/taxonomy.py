
# DESCRIPTION
# This script reads a CSV file containing taxid data and a tree file in itol format.
# It filters the taxid data by kingdom (Bacteria) and translates the nodes names in the tree file to taxids.
# The new tree file is saved in the data folder with the prefix "new_" in the name.

# Also there is another mode if your file lacks of "Kingdom" information column
# Then you can import modules `get_lineage` and `tax_filt`

# IMPORT LIBRARIES
import warnings
warnings.filterwarnings("ignore")
import argparse

import pandas as pd
from modules import get_lineage, tax_filt, tree_tax

parser = argparse.ArgumentParser()
parser.add_argument('--taxid_file', type = argparse.FileType('r'), required = True,
                    help = 'Path to the CSV file containing taxid data')
parser.add_argument('--tree_file', type = argparse.FileType('r'), default = None,
                    help = 'Path to the itol file to change the nodes names to taxids')
parser.add_argument('--get_tax', type = str, default = None,
                    help = 'Get csv table with taxid and selected tax')
args = parser.parse_args()

# Read the CSV into a DataFrame using pandas
tax_tab = pd.read_csv(args.taxid_file, header = None, sep = "\t")[[0, 1]]

# FILTERING THE DATAFRAME
tax_tab = tax_tab.rename(columns = {0: "protid",1: "taxid"}) # naming the column with taxid
tax_tab = pd.DataFrame(tax_tab).dropna(subset = ["taxid"]) # remove rows with empty taxid

# CLEANING THE TAXID COLUMN
if tax_tab.iloc[-1,-1] == "Search has CONVERGED!":
    tax_tab = tax_tab.iloc[:-1,] # removing the last row if it is the message "Search has CONVERGED!"
# rest only the first taxid in row if there are more than one
tax_tab["taxid"] = tax_tab["taxid"].str.split(";").str[0].astype(int)

# INFERING THE TAXONOMY
tax_tab['taxonomy'] = tax_tab.taxid.apply(get_lineage)
tax_tab = tax_filt(tax_tab)
tax_tab = tax_tab.join(tax_tab.taxonomy.str.split(",", expand = True, n = 11),  
                       how = 'left', lsuffix = '_left', rsuffix = '_right')
tax_tab.columns = ['protid', 'taxid', 'taxonomy', 0, 1, "kingdom", "phylum", "class", 
                   "order", "family", "genus", "species", 9, 10, 11]

# SELECTING BY KINGDOM
tax_tab = tax_tab[tax_tab["kingdom"] == "Bacteria"] # filter by kingdom

# TRANSLATION OF NODES NAMES
# Create a dictionary with the taxid as key and the saccver as value
tax_dict = tax_tab.set_index('protid')['taxid'].to_dict()

# OVERWRITING ITOL TREE FILE
tree_path = args.tree_file 
if tree_path != None:
    tree_tax(tree_path, tax_dict)

# GETTING CSV WITH TAX INFORMATION
tax_type = args.get_tax
if tax_type != None:
    tax_tab[['taxid', tax_type]].to_csv('../data/tax_to_lineage.csv', index = False)