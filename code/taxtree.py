
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('--taxid_file', type = argparse.FileType('r'), required = True,
                    help='Path to the CSV file containing taxid data')
parser.add_argument('--tree_file', type = argparse.FileType('r'), required = True,
                    help='Path to the itol file to change the nodes names to taxids')
args = parser.parse_args()

# Read the CSV into a DataFrame using pandas
df = pd.read_csv(args.taxid_file)


