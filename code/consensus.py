
# DESCRIPTION
# This script calculates a consensus sequence from a multiple sequence alignment,
# identifies conserved regions based on specified thresholds, and optionally saves
# the consensus sequence to a FASTA file.

# IMPORT LIBRARIES
import warnings
warnings.filterwarnings("ignore") # Consider more selective warning management
import argparse
import numpy as np
from typing import List, Tuple, Dict, Any, IO
import os

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment # Corrected import
from collections import Counter

def get_consensus_scores(alignment: MultipleSeqAlignment, num_sequences: int, alignment_length: int, n_most_common: int = 0) -> List[Tuple[str, float]]:
    """
    Calculates conservation scores for each position in an alignment.

    Args:
        alignment (MultipleSeqAlignment): The input alignment.
        num_sequences (int): Number of sequences in the alignment.
        alignment_length (int): Length of the alignment.
        n_most_common (int): Index for selecting the n-th most common residue (0 for the top one).

    Returns:
        List[Tuple[str, float]]: A list of tuples, where each tuple contains the
                                 most common residue and its conservation score (percentage)
                                 for a position.
    """
    scores: List[Tuple[str, float]] = []
    if num_sequences == 0: # Avoid division by zero
        return [('-', 0.0)] * alignment_length

    for i in range(alignment_length):
        column = [record.seq[i] for record in alignment]
        residue_counts = Counter(column)

        if not residue_counts: # Handle empty column case, though unlikely in valid MSA
            top_residue = '-'
            identical_count = 0
        elif len(residue_counts) <= n_most_common: # If fewer unique residues than n_most_common
             top_residue = residue_counts.most_common(1)[0][0] # Fallback to the absolute most common
             identical_count = column.count(top_residue)
        else:
            top_residue = residue_counts.most_common(n_most_common + 1)[n_most_common][0]
            identical_count = column.count(top_residue)

        score = (identical_count / num_sequences) * 100
        scores.append((top_residue, score))
    return scores

def get_conserved_alignment(
    initial_alignment: MultipleSeqAlignment,
    consensus_sequence_str: str, # Added: pass consensus string
    conservation_scores: List[Tuple[str, float]],
    threshold: float = 80.0
) -> Tuple[MultipleSeqAlignment, List[int]]:
    """
    Extracts a new alignment containing only the conserved regions above a given threshold.

    Args:
        initial_alignment (MultipleSeqAlignment): The original alignment.
        consensus_sequence_str (str): The calculated consensus sequence string.
        conservation_scores (List[Tuple[str, float]]): List of (residue, score) tuples.
        threshold (float): Conservation threshold (percentage).

    Returns:
        Tuple[MultipleSeqAlignment, List[int]]: A tuple containing the new
                                                MultipleSeqAlignment of conserved regions
                                                and a list of 0-indexed coordinates
                                                of these conserved regions in the
                                                original alignment.
    """
    conserved_coordinates: List[int] = []
    # conserved_residues_in_consensus: List[str] = [] # Renamed from conserved_regions

    for i, score_pair in enumerate(conservation_scores):
        if score_pair[1] >= threshold:
            conserved_coordinates.append(i)
            # conserved_residues_in_consensus.append(consensus_sequence_str[i]) # Use passed consensus

    conserved_alignment_records: List[SeqRecord] = []
    for record in initial_alignment:
        conserved_sequence_str = "".join([record.seq[i] for i in conserved_coordinates])
        conserved_record = SeqRecord(Seq(conserved_sequence_str), id=record.id, description=record.description)
        conserved_alignment_records.append(conserved_record)

    final_conserved_alignment = MultipleSeqAlignment(conserved_alignment_records)
    return final_conserved_alignment, conserved_coordinates


def main():
    parser = argparse.ArgumentParser(description="Calculate consensus sequence and identify conserved regions.")
    parser.add_argument('--align_file', type=argparse.FileType('r'), required=True,
                        help='Path to the multiple sequence alignment file (e.g., Clustal format)')
    parser.add_argument('--save_consensus', action='store_true', # Changed from type=bool
                        help='Flag to save the consensus sequence to a FASTA file.')
    parser.add_argument('--consensus_outfile', type=str, default="../data/consensus.fa",
                        help='Path to save the consensus FASTA file (used if --save_consensus is active). Default: ../data/consensus.fa')
    parser.add_argument('--stats_outfile', type=argparse.FileType('w'), default=None,
                        help='Optional: Path to save the conservative statistics dictionary as a JSON file.')
    parser.add_argument('--conserved_residue_threshold', type=float, default=90.0,
                        help='Threshold for printing specific conserved residues. Default: 90.0')

    args = parser.parse_args()

    try:
        alignment = AlignIO.read(args.align_file, "clustal")
    except ValueError as e:
        print(f"Error reading alignment file '{args.align_file.name}': {e}. Ensure it's a valid Clustal file.")
        return
    finally:
        args.align_file.close()

    num_sequences = len(alignment)
    alignment_length = alignment.get_alignment_length()

    if num_sequences == 0 or alignment_length == 0:
        print("Alignment is empty. Cannot proceed.")
        return

    print(f"Number of sequences: {num_sequences}")
    print(f"Alignment length: {alignment_length}")

    # Calculate consensus scores and the consensus sequence string
    # The original script used n_res=0 for get_cons_scores, which means taking the most common.
    consensus_scores: List[Tuple[str, float]] = get_consensus_scores(alignment, num_sequences, alignment_length, n_most_common=0)

    consensus_sequence_str = "".join([score_pair[0] for score_pair in consensus_scores])
    print(f"Consensus sequence: {consensus_sequence_str}")

    if args.save_consensus:
        output_dir = os.path.dirname(args.consensus_outfile)
        if output_dir and not os.path.exists(output_dir):
            try:
                os.makedirs(output_dir)
                print(f"Created output directory: {output_dir}")
            except OSError as e:
                print(f"Error creating output directory {output_dir}: {e}")
                # Decide if to proceed or exit. For now, proceed and let file open fail.

        try:
            with open(args.consensus_outfile, "w") as f: # Use "w" to overwrite or create new
                f.write(f">consensus_{os.path.basename(args.align_file.name)}\n{consensus_sequence_str}\n")
            print(f"Consensus sequence saved to: {args.consensus_outfile}")
        except IOError as e:
            print(f"Error writing consensus sequence to {args.consensus_outfile}: {e}")


    conservative_statistics: Dict[str, Any] = {}
    thresholds_to_test = np.array((80, 90, 95, 97, 99, 100))

    print("\nCalculating conservative statistics for various thresholds...")
    for t_val in thresholds_to_test:
        # Call get_conserved_alignment once and unpack
        conserved_aln_obj, cons_coords = get_conserved_alignment(
            alignment, consensus_sequence_str, consensus_scores, threshold=t_val
        )
        conservative_statistics[f"Threshold_{t_val}%"] = {
            "conserved_length": conserved_aln_obj.get_alignment_length(),
            "conserved_coordinates": list(cons_coords), # Store as list for easier serialization if needed
            # "conserved_alignment_object": conserved_aln_obj # Storing the object might be memory intensive
        }
        print(f"  Threshold {t_val}%: Conserved Length = {conservative_statistics[f'Threshold_{t_val}%']['conserved_length']}")

    # Print or process conservative_statistics
    if args.stats_outfile:
        import json
        try:
            json.dump(conservative_statistics, args.stats_outfile, indent=2, default=str) # default=str for np.array
            print(f"\nConservative statistics saved to: {args.stats_outfile.name}")
        except IOError as e:
            print(f"\nError writing statistics to {args.stats_outfile.name}: {e}")
        except TypeError as e:
            print(f"\nError serializing statistics to JSON: {e}. Ensure all data is JSON serializable.")
        finally:
            args.stats_outfile.close()
    else:
        # Basic print if not saving to file, already printed during calculation.
        print("\nConservative Statistics Summary (already printed above during calculation).")

    # NOTE: Consider adding a plotting function here if matplotlib is available.
    # For example, plot `consensus_scores` (residue score vs position)
    # or plot lengths from `conservative_statistics` vs threshold.
    # Example placeholder:
    # if args.plot_stats and PLOTTING_LIBRARY_AVAILABLE:
    #   plot_conservation_scores(consensus_scores, f"conservation_plot_{os.path.basename(args.align_file.name)}.png")
    #   plot_conservative_lengths(conservative_statistics, f"length_vs_threshold_{os.path.basename(args.align_file.name)}.png")


    # Print residues with conservation between a specific range (e.g., >=90% and <100%)
    # Using the threshold from CLI arguments
    print(f"\nResidues with {args.conserved_residue_threshold}% <= conservation < 100% (0-indexed):")
    count_specific_residues = 0
    for i, (residue, score) in enumerate(consensus_scores):
        if args.conserved_residue_threshold <= score < 100:
            print(f"  Position {i}: Residue '{residue}', Score {score:.2f}%")
            count_specific_residues +=1
    if count_specific_residues == 0:
        print(f"  No residues found with conservation >= {args.conserved_residue_threshold}% and < 100%.")

if __name__ == "__main__":
    main()