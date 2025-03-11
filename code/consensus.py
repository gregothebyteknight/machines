
# DESCRIPTION

# IMPORT LIBRARIES
import warnings
warnings.filterwarnings("ignore")
import argparse
import numpy as np

from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('--align_file', type = argparse.FileType('r'), required = True,
                    help = 'Path to the alignment file')
parser.add_argument('--get_cons', type = bool, default = False,
                    help = 'Whether save consensus.fa or not')
args = parser.parse_args()

consensus = "" # variable to store the consensus sequence
conservative_statistics = {}
count = 0

align = AlignIO.read(args.align_file, "clustal")
n_seqs = len(align)
align_len = align.get_alignment_length()
print(f"Number of sequences: {n_seqs}")
print(f"Alignment length: {align_len}")

def get_cons_scores(align, align_len, n_res = 0):
  scores = []
  for i in range(align_len):
      column = [seq[i] for seq in align]
      res_counts = Counter(column)
      top_res = res_counts.most_common(5)[n_res][0]
      identical_count = column.count(top_res)
      score = identical_count / n_seqs * 100
      scores.append((top_res, score))
  return scores

cons_scores = get_cons_scores(align, align_len)

for residue_pair in cons_scores:
  consensus += residue_pair[0]
print(f"Consensus sequence: {consensus}")

if args.get_cons:
    with open("../data/consensus.fa", "a") as f:
        f.write(f">consensus\n{consensus}")

def get_cons_aln(initial_alignment, align_len, conservation_scores, threshold = 80):
  conserved_coordinates = []
  conserved_regions = []
  conserved_alignment = []

  for i in range(align_len):
      if conservation_scores[i][1] >= threshold:
          conserved_coordinates.append(i)
          conserved_regions.append(consensus[i])

  for record in initial_alignment:
      conserved_sequence = "".join([record.seq[i] for i in conserved_coordinates])
      conserved_record = SeqRecord(Seq(conserved_sequence), id = record.id)
      conserved_alignment.append(conserved_record)

  conserved_alignment = MultipleSeqAlignment(conserved_alignment)
  return (conserved_alignment, conserved_coordinates)

for threshold in np.array((80, 90, 95, 97, 99, 100)):
    conservative_statistics[f"Threshold{threshold}"] = (get_cons_aln(align, align_len, cons_scores, 
                                                                     threshold)[0].get_alignment_length(),
                                                        np.array(get_cons_aln(align, align_len, 
                                                                              cons_scores, threshold)[1]),
                                                        get_cons_aln(align, align_len, cons_scores, threshold)[0])
threshold = 90
print(f"Residues with {threshold}% conservation:")

for i in range(len(cons_scores)):
  if 100 > cons_scores[i][1] >= threshold:
    print(cons_scores[i], count)
  count += 1