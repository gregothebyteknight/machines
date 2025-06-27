
"""
Analyzes a multiple sequence alignment to identify specific motifs and classify sites.

This script reads a multiple sequence alignment, extracts residues from predefined
positions to form a motif table, and then defines site groups based on the
residues at two specific positions (defaulting to 'pos28' and 'pos127' from the
motif table).

The script can print site group statistics to the console or save them to a file.
It also includes an unused dictionary `res_dict` and an example of how it could
be used with the `seq_with_res` function from `modules.py` for more specific
sequence filtering based on exact residue matches at multiple positions.
"""

# IMPORT MODULES AND VARIABLES
import argparse
from typing import Dict, IO
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
import pandas as pd

# Assuming modules.py is in the same directory or PYTHONPATH
from modules import motive_table, define_site
# from modules import seq_with_res # For the example usage

# This dictionary defines specific residues at particular 0-indexed positions.
# It is intended for use with the `seq_with_res` function (from modules.py)
# to filter sequences that match these exact residues at these exact positions.
# It is not directly used in the main workflow of this script otherwise.
res_dict: Dict[int, str] = {
    17: 'D',  # Example: position 17 (0-indexed) should be 'D'
    60: 'D',
    61: 'E',
    28: 'K',
    59: 'S',
    127: 'M',
    15: 'R',
    79: 'K',
    99: 'R',
    5: 'E',
    116: 'R'
}

# Define the 0-indexed positions to be extracted by motive_table.
# The keys (e.g., 'pos5') will be used as column names in the resulting DataFrame.
# These positions are critical for the subsequent site definition.
motif_positions_to_extract: Dict[str, int] = {
    'pos5': 5,
    'pos15': 15,
    'pos17': 17,
    'pos28': 28,  # Used by define_site default
    'pos59': 59,
    'pos60': 60,
    'pos61': 61,
    'pos79': 79,
    'pos99': 99,
    'pos116': 116,
    'pos127': 127  # Used by define_site default
}

def main():
    parser = argparse.ArgumentParser(description="Analyze motifs and define sites in a sequence alignment.")
    parser.add_argument('--align_file', type=argparse.FileType('r'), required=True,
                        help='Path to the multiple sequence alignment file (e.g., Clustal format)')
    parser.add_argument('--stats_outfile', type=argparse.FileType('w'), default=None,
                        help='Optional: Path to save the site group statistics (e.g., as CSV). Prints to console if not provided.')
    # TODO: Add arguments for pos1_label and pos2_label for define_site if customization is needed.

    args = parser.parse_args()

    try:
        alignment: MultipleSeqAlignment = AlignIO.read(args.align_file, "clustal")
    except FileNotFoundError: # Should be caught by argparse.FileType('r') but good practice
        print(f"Error: Alignment file not found at '{args.align_file.name}'.")
        return # exit() is not ideal in functions
    except ValueError as e:
        print(f"Error reading alignment file '{args.align_file.name}': {e}. Ensure it's a valid Clustal format.")
        return
    finally:
        if args.align_file: # Ensure it was successfully opened before trying to close
            args.align_file.close()

    print(f"Alignment length: {alignment.get_alignment_length()}")

    # Generate the table of residues at specified motif positions
    motif_df: pd.DataFrame = motive_table(alignment, motif_positions_to_extract)

    if motif_df.empty:
        print("Motif table is empty. This might happen if the alignment was empty or positions were out of bounds for all sequences.")
        return

    # Define site groups based on residues at 'pos28' and 'pos127' (default labels in define_site)
    # These labels must match keys in `motif_positions_to_extract`
    motif_df['site_group'] = motif_df.apply(
        lambda row: define_site(row, pos1_label='pos28', pos2_label='pos127'),
        axis=1
    )

    site_group_stats: pd.Series = motif_df['site_group'].value_counts()

    if args.stats_outfile:
        try:
            site_group_stats.to_csv(args.stats_outfile)
            print(f"\nSite group statistics saved to: {args.stats_outfile.name}")
        except IOError as e:
            print(f"\nError writing statistics to file {args.stats_outfile.name}: {e}")
            print("\nSite group statistics (fallback to console):")
            print(site_group_stats)
        finally:
            args.stats_outfile.close()
    else:
        print("\nSite group statistics:")
        print(site_group_stats)

    # Example of how res_dict could be used with seq_with_res (from modules.py)
    # This is currently commented out as it's not part of the original script's direct flow.
    # print("\n--- Example: Checking for sequences matching res_dict criteria ---")
    # from modules import seq_with_res # Keep import local to example if not broadly used
    # seq_with_res(alignment, res_dict)
    # print("--- End of Example ---")

    # NOTE: For more advanced motif analysis often done in notebooks:
    # 1. Motif Discovery: If motifs are not predefined, tools like MEME (callable externally)
    #    or Python libraries for motif discovery could be integrated.
    # 2. Motif Visualization: For the identified `site_group` or other motif patterns,
    #    generating sequence logos can be very informative.
    #    Libraries like `logomaker` (Python) or calling `weblogo` (external tool) could be used.
    #    Example placeholder for sequences belonging to a specific site_group:
    #    # sequences_for_group_X = alignment_records_for_group_X
    #    # if PLOTTING_LIBRARY_AVAILABLE and sequences_for_group_X:
    #    #   generate_sequence_logo(sequences_for_group_X, f"logo_group_X.png")

if __name__ == "__main__":
    main()