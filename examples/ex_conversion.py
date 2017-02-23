import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.sequence import SequenceSet
from quantiprot.utils.sequence import subset, columns
from quantiprot.utils.feature import Feature, FeatureSet

# Conversions-related imports:
from quantiprot.utils.mapping import simplify
from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.metrics.basic import identity

# Load the 'data/Alphasyn.fasta' sequence set, which contains several
# peptides from alpha-synuclein deposed in the Amyload database:
alphasyn_seq = load_fasta_file("data/Alphasyn.fasta")

# Retrieve predefined mapping from aminoacids to formal charge,
# and AAindex mapping to relative frequency of occurence (entry: JOND920101)
aa2charge_map = get_aa2charge()
aa2freq_map = get_aaindex_file("JOND920101")
print aa2charge_map
print aa2freq_map

# Make Feature objects based on Mappings:
charge_feat = Feature(aa2charge_map)
freq_feat = Feature(aa2freq_map)
print charge_feat
print freq_feat

# And use them to covert 1st sequence in 'alphasyn_seq':
print charge_feat(alphasyn_seq[0])
print freq_feat(alphasyn_seq[0])

# Make a FeatureSet from a Feature and Mappings:
fs = FeatureSet("basic features")
fs.add(charge_feat)
fs.add(aa2freq_map, name="frequency")
fs.add(get_aa2hydropathy())
print fs

# And use it to convert a set of sequences:
conv_seq = fs(alphasyn_seq)
print conv_seq
for seq in conv_seq:
    print seq

# Retrieve only sequences with selected features:
sub_seq = subset(conv_seq, features=["formal_charge", "hydropathy"])
print alphasyn_seq
print conv_seq
print sub_seq

# numpy used for prettier output only.
import numpy as np

# Export data to a list of lists:
print np.matrix(columns(conv_seq, feature="hydropathy"))

#By default sequences are in columns: set 'transpose' to True to change it:
print np.matrix(columns(conv_seq, feature="hydropathy", transpose=True))

# Load a mapping from amino acids to van der Waals volume:
aa2vol_map = get_aaindex_file("FAUJ880103", default=4.04)
print aa2vol_map

# Discretize the mapping to 3 levels using the k-means algorithm:
aa2vol_km_map = simplify(aa2vol_map, "aa2vol_kmeans", k=3, method='auto', iters=10)
print aa2vol_km_map

# Discretize the mapping to 3 levels using the user-defined thresholds:
aa2vol_thresh_map = simplify(aa2vol_map, "aa2vol_thresh", thresholds=[2.5, 5.5])
print aa2vol_thresh_map

# Create other mappings by k-means clustering with customized labels 
# (text and cluster means)
aa2vol_km_text_map = simplify(aa2vol_map, "aa2vol_kmeans_text", k=3, method='auto',
                              iters=10, labels=['small','medium','large'],
                              default='unknown')
print aa2vol_km_text_map
aa2vol_km_mean_map = simplify(aa2vol_map, "aa2vol_kmeans_mean", k=3, method='auto',
                              iters=10, mean_labels=True)
print aa2vol_km_mean_map

# Create a new FeatureSet:
vol_fs = FeatureSet("aa volume levels")
vol_fs.add(identity)
vol_fs.add(aa2vol_thresh_map)
vol_fs.add(aa2vol_km_map)
vol_fs.add(aa2vol_km_text_map)
vol_fs.add(aa2vol_km_mean_map)

# And convert sequences in 'alphasyn_seq':
result_seq = vol_fs(alphasyn_seq)

for seq in result_seq:
    print seq
