import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.metrics.aaindex import get_aaindex_file
from quantiprot.metrics.basic import average

# Load data:
seq = load_fasta_file("data/Alphasyn.fasta")

# Build a feature: average polarity (Graham, 1974), AAindex entry: GRAR740102:
feat = Feature(get_aaindex_file("GRAR740102")).then(average)

# Add the feature to new feature set:
fs = FeatureSet("my set")
fs.add(feat)

# Process sequences:
res_seq = fs(seq)

# Export average polarities
res = res_seq.columns()
print res
