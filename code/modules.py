
# IMPORT LIBRARIES
import re
import os
from typing import List, Dict, Any, IO
import ete3
import pandas as pd
from Bio.Align import MultipleSeqAlignment
from Bio.SeqRecord import SeqRecord

# ==============================================================================
# Taxonomy Utilities
# ==============================================================================

def tax_filt(tax_tab: pd.DataFrame) -> pd.DataFrame:
    """
    Filters a DataFrame containing taxonomy information by removing specified patterns
    from the 'taxonomy' column.

    Args:
        tax_tab (pd.DataFrame): DataFrame with taxonomic data, expecting a 'taxonomy' column.

    Returns:
        pd.DataFrame: The DataFrame with cleaned 'taxonomy' column.
    """
    patterns: List[str] = [
        'FCB group,',
        'TACK group,',
        'Bacteroidota/Chlorobiota group,',
        'Cyanobacteriota/Melainabacteria group,',
        'Terrabacteria group,',
        'Stenosarchaea group,',
        'Bacillus cereus group,',
        'Leptospirillum sp. Group III,',
        'Aneurinibacillus group,',
        'Rhizobium/Agrobacterium group,',
        'Caldisericota/Cryosericota group,',
        'PVC group,',
        'environmental samples,',
        'Chlorobium/Pelodictyon group,',
        'Leptolyngbya group,'
    ]

    pattern_regex: str = '|'.join(map(re.escape, patterns))
    # Ensure 'taxonomy' column is string type before replace
    tax_tab['taxonomy'] = tax_tab['taxonomy'].astype(str).str.replace(pattern_regex, '', regex=True)
    return tax_tab

def get_lineage(tax_id: int) -> str:
    """
    Retrieves the full taxonomic lineage for a given NCBI taxon ID.

    Args:
        tax_id (int): The NCBI taxon ID.

    Returns:
        str: A comma-separated string representing the lineage.
             Returns an empty string if lineage retrieval fails.
    """
    try:
        ncbi = ete3.NCBITaxa()
        lineage_ids: List[int] = ncbi.get_lineage(tax_id)
        if lineage_ids is None:
            return ""
        lineage_translator: Dict[int, str] = ncbi.get_taxid_translator(lineage_ids)
        lineage_names: List[str] = [lineage_translator[taxid] for taxid in lineage_ids if taxid in lineage_translator]
        return ",".join(lineage_names)
    except Exception as e:
        print(f"Error retrieving lineage for taxid {tax_id}: {e}")
        return ""

# ==============================================================================
# Tree Utilities
# ==============================================================================

def tree_tax(tree_file: IO[str], tax_dict: Dict[str, Any], output_dir: str = "../data") -> str:
    """
    Updates leaf names in a phylogenetic tree based on a taxonomy dictionary
    and writes the modified tree to a new file.

    Args:
        tree_file (IO[str]): An open file object for the input tree file (e.g., Newick format).
        tax_dict (Dict[str, Any]): Dictionary mapping original leaf names to new names (e.g., taxids).
        output_dir (str): Directory to save the modified tree. Defaults to "../data".

    Returns:
        str: The path to the written tree file.
    """
    tree = ete3.Tree(tree_file.name, format=3) # format=3 often means Newick with names & branch lengths
    for leaf in tree:
        # If no matching taxid found, keep the original name
        new_id = tax_dict.get(leaf.name, leaf.name)
        leaf.name = str(new_id) # Ensure name is a string

    # Ensure output directory exists
    try:
        if not os.path.exists(output_dir): # Check explicitly to print message only when creating
            os.makedirs(output_dir)
            # print(f"Created directory for tree output: {output_dir}") # Optional: for debugging
        # If it already exists, os.makedirs(output_dir, exist_ok=True) would also work silently.
    except OSError as e:
        print(f"Error creating directory {output_dir}: {e}")
        # Depending on desired behavior, could re-raise or return an error indicator
        return f"ERROR: Could not create directory {output_dir}"


    base_name = os.path.basename(tree_file.name)
    name, ext = os.path.splitext(base_name)
    output_tree_path = os.path.join(output_dir, f"tax_{name}{ext}")

    try:
        tree.write(outfile=output_tree_path, format=3)
        print(f"Modified tree created: {output_tree_path}")
    except Exception as e: # Catching general exception from tree.write
        print(f"Error writing tree to {output_tree_path}: {e}")
        return f"ERROR: Could not write tree to {output_tree_path}"

    return output_tree_path

# ==============================================================================
# Sequence & Motif Utilities
# ==============================================================================

def seq_with_res(alignment: MultipleSeqAlignment, residue_conditions: Dict[int, str]) -> None:
    """
    Identifies and prints sequences in an alignment that match specific residues at given positions.

    Args:
        alignment (MultipleSeqAlignment): The multiple sequence alignment object.
        residue_conditions (Dict[int, str]): A dictionary where keys are 0-based
                                             positions and values are the expected
                                             single-letter amino acid codes.
    """
    counter = 0
    for record in alignment:
        seq_str = str(record.seq)
        match = True
        for position, expected_residue in residue_conditions.items():
            if position >= len(seq_str) or seq_str[position] != expected_residue:
                match = False
                break

        if match:
            counter += 1
            print(f"Sequence ID: {record.id}")
            # print(f"Full Sequence: {seq_str}") # Uncomment if full sequence is needed
            print(f"Description: {record.description}")
            # Print matched residues for verification
            matched_residues_info = {pos: seq_str[pos] for pos in residue_conditions}
            print(f"Matched residues at positions: {matched_residues_info}")
            print()  # empty line for readability

    print(f"Number of sequences matching all conditions: {counter}")


def motive_table(alignment: MultipleSeqAlignment, motif_positions: Dict[str, int]) -> pd.DataFrame:
    """
    Creates a DataFrame summarizing specific amino acid residues at defined positions
    for each sequence in a multiple sequence alignment.

    Args:
        alignment (MultipleSeqAlignment): The multiple sequence alignment object.
        motif_positions (Dict[str, int]): A dictionary where keys are labels for the
                                           positions (e.g., 'H1_pos1') and values are
                                           the 0-based integer positions in the sequence.
                                           Example: {'17': 17, '60': 60, ...}

    Returns:
        pd.DataFrame: A DataFrame where each row corresponds to a sequence,
                      and columns include 'ProteinID' and the residues at the
                      specified motif positions.
    """
    result_list: List[Dict[str, Any]] = []

    for record in alignment:
        seq_str = str(record.seq)
        current_motif_data: Dict[str, Any] = {'ProteinID': record.id}

        for label, pos_0_indexed in motif_positions.items():
            if 0 <= pos_0_indexed < len(seq_str):
                current_motif_data[label] = seq_str[pos_0_indexed]
            else:
                current_motif_data[label] = '-' # Use '-' for out-of-bounds positions

        result_list.append(current_motif_data)

    motif_df = pd.DataFrame(result_list)
    return motif_df

def define_site(df_row: pd.Series, pos1_label: str = '28', pos2_label: str = '127') -> str:
    """
    Defines a site group based on amino acids at two specified positions in a DataFrame row.
    The input DataFrame row is expected to be a result from `motive_table`.

    Args:
        df_row (pd.Series): A row from the DataFrame generated by `motive_table`.
        pos1_label (str): The column label in df_row for the first position (e.g., '28').
        pos2_label (str): The column label in df_row for the second position (e.g., '127').

    Returns:
        str: The site group classification (e.g., '(K|R)M_group', 'other').
    """
    # Ensure the labels exist in the Series, otherwise default to a gap '-'
    res1 = df_row.get(pos1_label, '-')
    res2 = df_row.get(pos2_label, '-')
    seq_pair = f"{res1}{res2}"

    # Using a dictionary for cleaner mapping
    site_map: Dict[str, str] = {
        'KM': '(K|R)M_group', 'RM': '(K|R)M_group', # Consolidate (K|R)M
        'KS': '(K|R)S_group', 'RS': '(K|R)S_group', # Consolidate (K|R)S
        'K-': '(K|R)-_group', 'R-': '(K|R)-_group', # Consolidate (K|R)-
        'KG': 'KG_group',
        # 'RM': 'RM_group', # Already covered by (K|R)M_group
        'RA': 'RA_group',
        # 'RS': 'RS_group', # Already covered by (K|R)S_group
        'RG': 'RG_group',
        'AM': 'AM_group',
        'AA': 'AA_group',
        'AS': 'AS_group',
        'A-': 'A-_group',
        'AG': 'AG_group'
    }

    # Special handling for patterns like (K|R)M etc. needs careful ordering or direct check
    if res1 in ('K', 'R') and res2 == 'M':
        return '(K|R)M_group'
    if res1 in ('K', 'R') and res2 == 'S':
        return '(K|R)S_group'
    if res1 in ('K', 'R') and res2 == '-':
        return '(K|R)-_group'

    return site_map.get(seq_pair, 'other')