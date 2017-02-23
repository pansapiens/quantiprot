import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.sequence import SequenceSet
from quantiprot.utils.sequence import merge

# Load protein sequences from 'data/Amyload_positive.fasta':
amyload_pos_seq = load_fasta_file("data/Amyload_positive.fasta")

# Display first three sequences:
print amyload_pos_seq
for seq in amyload_pos_seq[:3]:
    print seq

# Find a sequence 'AMY438|7-13|Sup35' in 'amyload_pos_seq':
my_seq_index = amyload_pos_seq.ids().index("AMY438|7-13|Sup35")
my_seq = amyload_pos_seq[my_seq_index]
print my_seq

# And copy the sequence to a new sequence set:
my_seq_set = SequenceSet("my seq set")
my_seq_set.add(my_seq)
print my_seq_set

# Try again to add the same sequence to 'my_seq_set' with 'unique' = True:
my_seq_set.add(my_seq)
print my_seq_set

# Load another sequence sets:
amyload_neg_seq = load_fasta_file("data/Amyload_negative.fasta")

print amyload_neg_seq

# And merge with 'amyload_pos_seq'
amyload_merged_seq = merge(amyload_pos_seq, amyload_neg_seq)
print amyload_merged_seq

# Sometimes it is more convenient to merge in-place
# Here we perform two consequent mergings:
# 1) 'my_seq_set' with 'amyload_pos_seq'
# 2) the resulting set with 'amyload_neg_seq'
my_seq_set.merge_with(amyload_pos_seq).merge_with(amyload_neg_seq)
print my_seq_set
