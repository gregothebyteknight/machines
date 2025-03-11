
# IMPORT LIBRARIES
import re
import os
import ete3
import pandas as pd

def tax_filt(tax_tab):
    """
    This function filters the tax_tab DataFrame by kingdom and 
    returns a DataFrame with only the prot_id and taxid columns

    @tax_tab(pd.DataFrame): DataFrame with the taxid data
    """
    # List of patterns to remove from the taxonomy column
    patterns = [
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

    # Create a regex pattern to remove all the patterns at once
    pattern_regex = '|'.join(map(re.escape, patterns))
    # Remove the patterns from the taxonomy column
    tax_tab['taxonomy'] = tax_tab['taxonomy'].str.replace(pattern_regex, '', regex = True)
    return tax_tab

def get_lineage(tax):
  """
  This function returns the lineage of a taxon given its taxid

  @tax(int): taxid of the taxon
  """
  ncbi = ete3.NCBITaxa()
  lineage = ncbi.get_lineage(tax)
  lineage_serie = pd.Series(ncbi.get_taxid_translator(lineage), name = tax)
  return ",".join(lineage_serie[lineage].values)

def tree_tax(tree_path, tax_dict):
    tree = ete3.Tree(tree_path.name, format = 3)
    for leaf in tree:
        new_id = tax_dict.get(leaf.name, leaf.name)  # If no matching taxid found, keep the original saccver
        leaf.name = new_id

    # Write the new tree
    name, ext = os.path.splitext(os.path.basename(tree_path.name))
    tree.write(outfile = f"../data/tax_{name}{ext}", format = 3)
    print(f"tax_{name}{ext} created")

def seq_with_res(align, res_dict):

    # Iterate over the sequences in the alignment
    counter = 0
    for record in align:
        seq = str(record.seq)
        if (
            seq[pos1] == residue_1 and seq[position_2] == residue_2
            and seq[position_3] == residue_3 and seq[position_4] == residue_4
            and seq[position_5] == residue_5 and seq[position_6] == residue_6
            and seq[position_7] == residue_7 and seq[position_8] == residue_8
            and seq[position_9] == residue_9 and seq[position_10] == residue_10
            and seq[position_11] == residue_11):
            counter += 1
            print(f"Sequence: {seq}")
            print(f"Name: {record.id}")
            print(f"Description: {record.description}")
            print()  # empty line for readability
    print(f"Number of sequence fit: {counter}")


def motive_table(align):
    result = []

    for record in align:
        new_motif = {'ProteinID':record.id,

            '17': str(record.seq[17:18]),
            '60': str(record.seq[60:61]),
            '61': str(record.seq[61:62]),

            '28': str(record.seq[28:29]),
            '59': str(record.seq[59:60]),
            '127' :str(record.seq[127:128]),

            '15': str(record.seq[15:16]),
            '79': str(record.seq[79:80]),
            '99': str(record.seq[99:100]),

            '5': str(record.seq[5:6]),
            '116': str(record.seq[116:117])
                            }
    result.append(new_motif)
    motif = pd.DataFrame(result)
    return motif

def define_site(df):
    seq = "".join([df['28'], df['127']])

    if seq == 'KM' or seq == 'RM':
        return '(K|R)M_group'
    elif seq == 'KS' or seq == 'RS':
        return '(K|R)S_group'
    elif seq == 'K-' or 'R-':
        return '(K|R)-_group'
    elif seq == 'KG' or '':
        return 'KG_group'
    elif seq == 'RM':
        return 'RM_group'
    elif seq == 'RA':
        return 'RA_group'
    elif seq == 'RS':
        return 'RS_group'
    elif seq == 'RG':
        return 'RG_group'
    elif seq == 'AM':
        return 'AM_group'
    elif seq == 'AA':
        return 'AA_group'
    elif seq == 'AS':
        return 'AS_group'
    elif seq == 'A-':
        return 'A-_group'
    elif seq == 'AG':
        return 'AG_group'
    else:
        return 'other'