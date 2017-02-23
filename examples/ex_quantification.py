import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.metrics.aaindex import get_aa2charge, get_aa2hydropathy, get_aa2volume
from quantiprot.utils.mapping import simplify

# Quantification-related imports:
from quantiprot.metrics.basic import identity, average, sum_absolute, uniq_count
from quantiprot.utils.sequence import compact

# Load some data:
alphasyn_seq = load_fasta_file("data/Alphasyn.fasta")

# Prepare Features:
charge_sum_abs_feat = Feature(get_aa2charge()).then(sum_absolute)
hydropathy_average_feat = Feature(get_aa2hydropathy()).then(average)
volume_levels_feat = Feature(simplify(get_aa2volume(), name="volume levels",
                                      k=3)).then(uniq_count)

# Prepare a FeatureSet
fs = FeatureSet("simple quantification")
fs.add(hydropathy_average_feat)
fs.add(charge_sum_abs_feat)
fs.add(volume_levels_feat)

# And use it to quantify protein sequence(s):
result_seq = fs(alphasyn_seq)
print result_seq
for seq in result_seq:
    print seq

# Plain functions operating on list also work:
print Feature(len)(alphasyn_seq[0])

# Calculate the hydropathy profile smoothed over the window of length 3
hydropathy_win3_feat = Feature(get_aa2hydropathy()).then(average, window=3)
print hydropathy_win3_feat(alphasyn_seq[0])

# Buggy version of the code above:
hydropathy_feat = Feature(get_aa2hydropathy())
hydropathy_win3_feat = hydropathy_feat.then(average, window=3)
print hydropathy_feat(alphasyn_seq[0])
print hydropathy_win3_feat(alphasyn_seq[0])

# Compact multiple single-value features
compact_seq = compact(result_seq)
for seq in compact_seq:
    print seq
