
# IMPORT LIBRARIES
import warnings
warnings.filterwarnings("ignore")
import argparse
import ete3
import os

import pandas as pd
from modules import tax_filt, get_lineage

parser = argparse.ArgumentParser()
parser.add_argument('--taxid_file', type = argparse.FileType('r'), required = True,
                    help='Path to the CSV file containing taxid data')
parser.add_argument('--tree_file', type = argparse.FileType('r'), required = True,
                    help='Path to the itol file to change the nodes names to taxids')
args = parser.parse_args()

# Read the CSV into a DataFrame using pandas
tax_tab = pd.read_csv(args.taxid_file, header = None, sep = "\t")

# FILTERING THE DATAFRAME
tax_tab = tax_tab.rename(columns={0: "prot_id", 1: "tax_id", 2: "kingdom"}) # naming the column with taxid
tax_tab = tax_tab.dropna(subset = ["tax_id"]) # remove rows with empty taxid

# CLEANING THE TAXID COLUMN
if tax_tab.iloc[-1,-1] == "Search has CONVERGED!":
    tax_tab = tax_tab.iloc[:-1,] # removing the last row if it is the message "Search has CONVERGED!"
# rest only the first taxid in row if there are more than one
tax_tab["tax_id"] = tax_tab["tax_id"].str.split(";").str[0].astype(int)

# SELECTING BY KINGDOM
tax_tab = tax_tab[tax_tab["kingdom"] == "Bacteria"] # filter by kingdom
# selection only prot_id and taxid columns
tax_tab = pd.DataFrame(tax_tab[["prot_id", "tax_id"]])

# TRANSLATION OF NODES NAMES
# Create a dictionary with the taxid as key and the saccver as value
tax_dict = tax_tab.set_index('prot_id')['tax_id'].to_dict()

# Read the itol file
t = ete3.Tree(args.tree_file.name, format = 3)
for leaf in t:
    new_id = tax_dict.get(leaf.name, leaf.name)  # If no matching taxid found, keep the original saccver
    leaf.name = new_id

# Write the new tree
name, ext = os.path.splitext(os.path.basename(args.tree_file.name))
t.write(outfile = f"../data/new_{name}{ext}", format = 3)
print(f"new_{name}{ext} created")