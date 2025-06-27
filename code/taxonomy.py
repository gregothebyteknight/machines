
# DESCRIPTION
# This script reads a CSV file containing taxid data and a tree file in itol format.
# It filters the taxid data by kingdom (Bacteria) and translates the nodes names in the tree file to taxids.
# The new tree file is saved in the data folder with the prefix "new_" in the name.

# Also there is another mode if your file lacks of "Kingdom" information column
# Then you can import modules `get_lineage` and `tax_filt`

# IMPORT LIBRARIES
import warnings
warnings.filterwarnings("ignore") # Consider managing warnings more selectively
import argparse
from typing import Dict, Any, IO, List
import pandas as pd
import os # For path manipulation

from modules import get_lineage, tax_filt, tree_tax, generate_itol_simple_annotation

# Standard taxonomic ranks often encountered. The order matters for mapping.
# NCBITaxa output might include higher levels like "cellular organisms" or "root" at the beginning.
STANDARD_RANKS = ["superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species"]


def main():
    parser = argparse.ArgumentParser(
        description="Process taxonomic information, update phylogenetic trees, and export taxonomic data."
    )
    parser.add_argument('--taxid_file', type=argparse.FileType('r'), required=True,
                        help='Path to the CSV file containing taxid data (tab-separated, no header, protid in col 0, taxid in col 1)')
    parser.add_argument('--tree_file', type=argparse.FileType('r'), default=None,
                        help='Path to the tree file (e.g., Newick) to update with taxids')
    parser.add_argument('--get_tax', type=str, default=None,
                        help='Specify a taxonomic rank (e.g., "phylum", "genus") to export with taxids to a CSV file.')
    parser.add_argument('--output_dir', type=str, default='../data',
                        help='Directory for output files (modified tree, taxonomy CSV). Default: ../data')
    parser.add_argument('--output_csv_name', type=str, default=None,
                        help='Optional: Filename for the exported taxonomy CSV. If not provided, a name is generated based on --get_tax.')
    parser.add_argument('--filter_kingdom', type=str, default="Bacteria",
                        help='Kingdom to filter by (e.g., "Bacteria", "Archaea", "Viruses"). Set to empty string "" to disable kingdom filtering. Default: "Bacteria"')
    parser.add_argument('--make_itol_annotation', type=str, default=None,
                        help='Column name from the taxonomy table to use for iTOL simple annotation (e.g., "phylum"). Requires --protid_column_name to be set if not "protid".')
    parser.add_argument('--itol_output_file', type=str, default="itol_simple_annotation.txt",
                        help='Filename for the iTOL annotation file. Will be saved in --output_dir. Default: itol_simple_annotation.txt')
    parser.add_argument('--protid_column_name', type=str, default="protid",
                        help='Name of the column containing protein IDs (for matching tree leaves in iTOL annotation). Default: "protid"')


    args = parser.parse_args()

    # Ensure output directory exists
    if not os.path.exists(args.output_dir):
        try:
            os.makedirs(args.output_dir)
            print(f"Created output directory: {args.output_dir}")
        except OSError as e:
            print(f"Error creating output directory {args.output_dir}: {e}")
            return


    # Read the CSV into a DataFrame using pandas
    try:
        tax_tab: pd.DataFrame = pd.read_csv(args.taxid_file, header=None, sep="\t", usecols=[0, 1])
    except pd.errors.EmptyDataError:
        print(f"Error: The taxid file {args.taxid_file.name} is empty or not formatted correctly.")
        return
    except Exception as e:
        print(f"Error reading taxid file {args.taxid_file.name}: {e}")
        return
    finally:
        args.taxid_file.close()


    # FILTERING THE DATAFRAME
    tax_tab = tax_tab.rename(columns={0: "protid", 1: "taxid"})
    tax_tab = tax_tab.dropna(subset=["taxid"])

    # CLEANING THE TAXID COLUMN
    if not tax_tab.empty and isinstance(tax_tab.iloc[-1, -1], str) and "Search has CONVERGED!" in tax_tab.iloc[-1, -1]:
        tax_tab = tax_tab.iloc[:-1, :]

    tax_tab["taxid"] = tax_tab["taxid"].astype(str).str.split(";").str[0]
    tax_tab = tax_tab[tax_tab["taxid"].str.match(r'^\d+$')]
    if tax_tab.empty:
        print("No valid numeric taxids found after cleaning.")
        return
    tax_tab["taxid"] = tax_tab["taxid"].astype(int)


    # INFERRING THE TAXONOMY
    print("Inferring taxonomy (this may take some time)...")
    tax_tab['full_lineage_str'] = tax_tab["taxid"].apply(get_lineage)
    # Apply tax_filt, which expects a 'taxonomy' column. We'll rename full_lineage_str temporarily.
    tax_tab_for_filt = tax_tab.rename(columns={'full_lineage_str': 'taxonomy'})
    tax_tab_for_filt = tax_filt(tax_tab_for_filt)
    tax_tab['cleaned_lineage_str'] = tax_tab_for_filt['taxonomy']

    # Expand cleaned lineage string into separate columns
    lineage_df = tax_tab['cleaned_lineage_str'].str.split(",", expand=True)

    # Name new columns temporarily as rank_0, rank_1, ...
    temp_rank_cols = [f"temp_rank_{i}" for i in range(lineage_df.shape[1])]
    lineage_df.columns = temp_rank_cols
    tax_tab = tax_tab.join(lineage_df)

    # Attempt to map temp_rank columns to STANDARD_RANKS
    # This assumes the lineage from get_lineage starts with high-level ranks
    # and roughly aligns with STANDARD_RANKS.
    # E.g., if lineage is "cellular organisms,Bacteria,Terrabacteria,...",
    # we want to map "Bacteria" to "kingdom" if possible.
    # We iterate through STANDARD_RANKS and see if any temp_rank column typically contains it.
    # This is heuristic. A more robust way would be if get_lineage returned ranks.

    # Simple approach: Assume the first few columns of lineage_df might correspond to STANDARD_RANKS
    # This renaming is a best-effort for clarity.
    rename_map: Dict[str, str] = {}
    # Offset can be adjusted if NCBITaxa typically returns e.g. "root" or "cellular organisms" first.
    # Let's try to find "Bacteria", "Archaea", "Eukaryota", "Viruses" in the first few temp_ranks
    # and map that column to "kingdom" or "superkingdom".

    # A simple heuristic: try to match standard rank names based on typical output.
    # This part is tricky because ete3's get_lineage just gives ordered names.
    # For now, let's assume the user will rely on args.get_tax using one of the
    # STANDARD_RANKS and we will try to find a column that seems to match.
    # The original code effectively picked 'rank_1' (second field) as kingdom.

    # Let's rename the first few temp_rank_cols to the first few STANDARD_RANKS
    # This is a simplification and might not always be accurate.
    max_ranks_to_map = min(len(temp_rank_cols), len(STANDARD_RANKS))
    for i in range(max_ranks_to_map):
        rename_map[temp_rank_cols[i]] = STANDARD_RANKS[i]
    tax_tab = tax_tab.rename(columns=rename_map)


    # Filter by Kingdom if specified
    if args.filter_kingdom and "kingdom" in tax_tab.columns:
        print(f"Filtering by kingdom: {args.filter_kingdom}")
        tax_tab_filtered = tax_tab[tax_tab["kingdom"] == args.filter_kingdom].copy()
        if tax_tab_filtered.empty:
            print(f"No entries found for Kingdom '{args.filter_kingdom}'.")
            # Decide: continue with empty dataframe or stop? For now, continue.
            # tax_tab = tax_tab_filtered # This would make tax_tab empty
        else:
            print(f"Filtered down to {len(tax_tab_filtered)} entries for Kingdom '{args.filter_kingdom}'.")
            tax_tab = tax_tab_filtered
    elif args.filter_kingdom:
        print(f"Warning: Kingdom column 'kingdom' not found after renaming. Skipping kingdom filtering for '{args.filter_kingdom}'.")


    # TRANSLATION OF NODES NAMES FOR THE TREE
    protid_to_taxid_map: Dict[str, Any] = tax_tab.set_index('protid')['taxid'].to_dict()

    if args.tree_file:
        # tree_tax will attempt to create output_dir if it doesn't exist via modules.py logic
        # It returns the path on success, or an error message string on failure.
        output_tree_message = tree_tax(args.tree_file, protid_to_taxid_map, output_dir=args.output_dir)

        if output_tree_message.startswith("ERROR:"):
            print(output_tree_message) # Print the error message from tree_tax
        else:
            # This was the success case path from tree_tax
            print(f"Modified tree processing complete. Path: {output_tree_message}")
            # NOTE: To visualize the tree directly from script (if ete3 GUI tools are available, or to save to image):
            # try:
            #     from ete3 import TreeStyle, NodeStyle
            #     t = ete3.Tree(output_tree_message, format=3) # Reload the modified tree
            #     ts = TreeStyle()
            #     ts.show_leaf_name = True
            #     # Example: Save to PNG
            #     # tree_image_path = os.path.join(args.output_dir, f"rendered_tree_{os.path.splitext(os.path.basename(args.tree_file.name))[0]}.png")
            #     # t.render(tree_image_path, tree_style=ts)
            #     # print(f"Rendered tree image saved to: {tree_image_path}")
            #     # Or show interactively (if environment supports it)
            #     # t.show(tree_style=ts)
            # except ImportError:
            #     print("ete3 plotting tools not fully available to render/show tree image directly.")
            # except Exception as e:
            #     print(f"Error during optional tree rendering/showing: {e}")

        # Close the tree_file regardless of success or failure of tree_tax, if it was opened.
        # argparse.FileType handles opening.
        try:
            args.tree_file.close()


    # GETTING CSV WITH TAX INFORMATION
    if args.get_tax:
        selected_rank_col_name = args.get_tax # User provides standard name e.g. "phylum"

        if selected_rank_col_name in tax_tab.columns and 'taxid' in tax_tab.columns:
            if args.output_csv_name:
                output_csv_filename = args.output_csv_name
            else:
                output_csv_filename = f"tax_to_{selected_rank_col_name}.csv"

            output_csv_full_path = os.path.join(args.output_dir, output_csv_filename)

            try:
                tax_tab[['taxid', selected_rank_col_name]].to_csv(output_csv_full_path, index=False)
                print(f"Taxonomy data for '{selected_rank_col_name}' saved to {output_csv_full_path}")
            except Exception as e:
                print(f"Error saving taxonomy data to CSV {output_csv_full_path}: {e}")
        else:
            print(f"Error: Column '{selected_rank_col_name}' or 'taxid' not found in the processed table.")
            print(f"Available columns for export: { [col for col in tax_tab.columns if col in STANDARD_RANKS or col == 'taxid'] }")
            print(f"All available columns: {tax_tab.columns.tolist()}")

    # GENERATE ITOL ANNOTATION FILE
    if args.make_itol_annotation:
        annotation_value_column = args.make_itol_annotation
        if args.protid_column_name not in tax_tab.columns:
            print(f"Error: Protein ID column '{args.protid_column_name}' not found in table for iTOL annotation.")
        elif annotation_value_column not in tax_tab.columns:
            print(f"Error: Annotation value column '{annotation_value_column}' not found in table for iTOL annotation.")
            print(f"Available columns: {tax_tab.columns.tolist()}")
        else:
            itol_file_path = os.path.join(args.output_dir, args.itol_output_file)
            print(f"\nGenerating iTOL simple annotation file for column '{annotation_value_column}'...")
            # Using a default color for now, can be made configurable
            # The generate_itol_simple_annotation function handles file writing and error messages.
            generate_itol_simple_annotation(
                data_df=tax_tab,
                id_column=args.protid_column_name,
                value_column=annotation_value_column,
                output_filepath=itol_file_path,
                dataset_label=f"{annotation_value_column}_annotation",
                color="#0000FF" # Example color: Blue
            )

if __name__ == '__main__':
    main()