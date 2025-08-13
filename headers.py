"""
Rewrites FASTA sequence headers based on a provided mapping file.
"""
import argparse
import csv
from typing import Dict, IO
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main():
    parser = argparse.ArgumentParser(description="Rewrite FASTA sequence headers using a mapping file.")
    parser.add_argument('--input_fasta', type=argparse.FileType('r'), required=True,
                        help='Path to the input FASTA file.')
    parser.add_argument('--mapping_file', type=argparse.FileType('r'), required=True,
                        help='Path to the mapping file (2-column CSV/TSV: old_header_id,new_header_id).')
    parser.add_argument('--output_fasta', type=argparse.FileType('w'), required=True,
                        help='Path for the output FASTA file with rewritten headers.')
    parser.add_argument('--map_delimiter', type=str, default=',',
                        help='Delimiter used in the mapping file (e.g., "," for CSV, "\\t" for TSV). Default: ","')
    # Potentially add --header_match_part later if more complex matching is needed.
    # For now, it matches based on the record.id from Biopython, which is usually the part before the first space.

    args = parser.parse_args()

    # 1. Read mapping file into a dictionary
    header_map: Dict[str, str] = {}
    try:
        reader = csv.reader(args.mapping_file, delimiter=args.map_delimiter)
        for i, row in enumerate(reader):
            if len(row) == 2:
                old_id, new_id = row[0].strip(), row[1].strip()
                if old_id in header_map:
                    print(f"Warning: Duplicate old_id '{old_id}' found in mapping file at line {i+1}. Using the last encountered mapping.")
                header_map[old_id] = new_id
            elif row: # Non-empty row that doesn't have 2 columns
                print(f"Warning: Skipping malformed row {i+1} in mapping file (expected 2 columns): {row}")
    except csv.Error as e:
        print(f"Error reading mapping file '{args.mapping_file.name}': {e}")
        args.mapping_file.close()
        args.input_fasta.close()
        args.output_fasta.close()
        return
    finally:
        args.mapping_file.close()

    if not header_map:
        print("Mapping dictionary is empty. No headers will be changed.")
    else:
        print(f"Loaded {len(header_map)} mappings.")

    # 2. Iterate through input FASTA, rewrite headers, and write to output
    records_written = 0
    records_modified = 0
    try:
        for record in SeqIO.parse(args.input_fasta, "fasta"):
            original_id = record.id
            # record.id is typically the string before the first whitespace.
            # If full description matching is needed, this part would change.

            if original_id in header_map:
                new_id = header_map[original_id]
                record.id = new_id
                record.description = "" # Keep description clean, only new ID
                records_modified += 1

            SeqIO.write(record, args.output_fasta, "fasta")
            records_written += 1

        print(f"\nProcessed {records_written} sequences.")
        print(f"Modified {records_modified} headers.")
        print(f"Output FASTA file saved to: {args.output_fasta.name}")

    except ValueError as e:
        print(f"Error parsing input FASTA file '{args.input_fasta.name}': {e}")
    except IOError as e:
        print(f"IOError during FASTA processing: {e}")
    finally:
        args.input_fasta.close()
        args.output_fasta.close()

if __name__ == "__main__":
    main()