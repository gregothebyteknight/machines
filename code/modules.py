import ete3
import re
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