# Bioinformatics Analysis Toolkit

This project provides a suite of Python scripts and Jupyter notebooks for performing common bioinformatics tasks, including sequence alignment analysis, taxonomic data processing, and motif discovery.

## Scripts and Notebooks

This project includes the following Python scripts and Jupyter notebooks:

*   **`code/consensus.py`**: This script takes a multiple sequence alignment file (e.g., in Clustal format) as input. It calculates the consensus sequence and identifies conserved regions based on a user-defined threshold. It can output the consensus sequence to a FASTA file.
*   **`code/taxonomy.py`**: This script processes taxonomic information. It takes a CSV file with protein IDs and NCBI taxids, retrieves their full taxonomic lineages, and can filter them (e.g., by kingdom). It can also update node names in a phylogenetic tree file (e.g., Newick format) with corresponding taxids or other taxonomic ranks.
*   **`code/motives.py`**: This script (likely intended to be `motifs.py`) analyzes aligned sequences to identify predefined amino acid motifs at specific positions. It uses helper functions from `modules.py` to define and report these motifs.
*   **`code/modules.py`**: This Python file contains various helper functions used by the other scripts. These include functions for filtering taxonomic data, retrieving NCBI lineage information, manipulating phylogenetic trees using `ete3`, and generating motif tables from alignments.
*   **Jupyter Notebooks**:
    *   `code/Conservative_analysis.ipynb`: Likely used for interactive exploration and visualization of sequence conservation analysis, complementing `consensus.py`.
    *   `code/Copy_of_ete3.ipynb`: Appears to be an exploratory or developmental notebook, possibly for testing `ete3` functionalities related to tree manipulation and visualization.
    *   `code/Motif_analysis.ipynb`: Likely provides an interactive environment for motif discovery and analysis, complementing `motives.py`.

## Data Directory

The `data/` directory serves as the default location for input files and generated output.

*   **Input Files**:
    *   Multiple sequence alignment files (e.g., `with_thg_small.clustal`, `with_thg_auto_small.trim`): Used by `consensus.py` and `motives.py`.
    *   CSV files with protein/sequence identifiers and corresponding NCBI taxids (e.g., `taxids.csv`): Used by `taxonomy.py`.
    *   Phylogenetic tree files (e.g., `tax_with_thg_small.iq.contree`, `with_thg_small.iq.contree`): Used by `taxonomy.py` for updating node labels.
*   **Output Files**:
    *   `consensus.fa`: Generated by `consensus.py`, containing the calculated consensus sequence.
    *   `tax_to_lineage.csv`: Generated by `taxonomy.py`, a CSV file mapping taxids to a specified taxonomic rank.
    *   `tax_<original_tree_filename>` (e.g., `tax_tax_with_thg_small.iq.contree`): New tree files generated by `taxonomy.py` with updated node names based on taxonomic information.

## Dependencies

The scripts rely on the following Python libraries:

*   **BioPython**: Used for sequence alignment processing, reading/writing sequence files (e.g., FASTA, Clustal).
*   **pandas**: Used for data manipulation, particularly for handling CSV files and tabular data.
*   **ete3**: Used for phylogenetic tree manipulation and visualization, including interacting with the NCBI taxonomy database.
*   **NumPy**: Used for numerical operations, especially in `consensus.py` for handling conservation scores.
*   **argparse**: Used for parsing command-line arguments in scripts like `consensus.py` and `taxonomy.py`.

You can typically install these using pip:
```bash
pip install biopython pandas ete3 numpy
```
(`argparse` is part of the Python standard library).

## Usage

### `code/consensus.py`

This script calculates a consensus sequence from a multiple sequence alignment.

**Arguments:**

*   `--align_file` (required): Path to the alignment file (e.g., Clustal format).
*   `--get_cons` (optional, boolean, default: `False`): If `True`, saves the consensus sequence to `data/consensus.fa`.

**Example:**

```bash
python code/consensus.py --align_file data/with_thg_small.clustal --get_cons True
```

This will print the number of sequences, alignment length, the consensus sequence, and conservation statistics. If `--get_cons True` is specified, it will also create/append to `data/consensus.fa`.

### `code/taxonomy.py`

This script processes taxonomic information and can update phylogenetic trees.

**Arguments:**

*   `--taxid_file` (required): Path to the CSV file containing protein IDs and taxids (tab-separated, no header).
*   `--tree_file` (optional): Path to a phylogenetic tree file (e.g., Newick format) whose node names will be updated.
*   `--get_tax` (optional, string): If specified, creates a CSV file `data/tax_to_lineage.csv` mapping taxids to the specified taxonomic rank (e.g., "phylum", "genus").

**Examples:**

1.  **Process taxids and update a tree:**
    ```bash
    python code/taxonomy.py --taxid_file data/taxids.csv --tree_file data/with_thg_small.iq.contree
    ```
    This will generate a new tree file named `data/tax_with_thg_small.iq.contree` with node names replaced by taxids where possible.

2.  **Process taxids and extract lineage information for a specific rank:**
    ```bash
    python code/taxonomy.py --taxid_file data/taxids.csv --get_tax phylum
    ```
    This will create `data/tax_to_lineage.csv` containing two columns: `taxid` and `phylum`.

### `code/motives.py`

This script identifies predefined motifs in an alignment file. It currently reads an alignment file named `tcbfprs1.aln` (hardcoded).

**To run:**

```bash
python code/motives.py
```
Ensure an alignment file named `tcbfprs1.aln` (Clustal format) exists in the same directory as `motives.py` or modify the script to point to the correct file. It will print motif statistics.

### Jupyter Notebooks

The Jupyter notebooks (`Conservative_analysis.ipynb`, `Copy_of_ete3.ipynb`, `Motif_analysis.ipynb`) can be run using a Jupyter Notebook or JupyterLab environment. They provide interactive ways to perform the analyses found in the Python scripts.
