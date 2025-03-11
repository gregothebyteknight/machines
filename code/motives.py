
# IMPORT MODULES AND VARIABLES
from Bio import AlignIO
import pandas as pd
from modules import motive_table, define_site

res_dict = {
    17: 'D', # pair index: AA
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

## загружаем элаймент
align = AlignIO.read("tcbfprs1.aln", "clustal")
print("Alignment length %i" % align.get_alignment_length())

motif = motive_table(align)

motif['motif'] = motif.apply(define_site, axis = 1)

# статистика по мотивам
motif['motif'].value_counts()


