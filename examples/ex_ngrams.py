import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.metrics.aaindex import get_aa2hydropathy
from quantiprot.metrics.basic import identity

# Ngram-related imports
from quantiprot.metrics.ngram import pattern_match, pattern_count
from quantiprot.analysis.ngram import ngram_count
from quantiprot.analysis.ngram import zipf_law_fit

from matplotlib import pyplot as plt

# Load some data
alphasyn_seq = load_fasta_file("data/Alphasyn.fasta")
amyload_pos_seq = load_fasta_file("data/Amyload_positive.fasta")

# Find and count matches to a pattern 'VT'
fs_aa = FeatureSet("aa patterns")
fs_aa.add(identity)
fs_aa.add(pattern_match, pattern='VT', padded=True)
fs_aa.add(pattern_count, pattern='VT')

result_seq = fs_aa(alphasyn_seq)

for seq in result_seq[:3]:
    print seq

# ...and something much more subtle:
# Map a sequence to the hydrophaty scale, and search for the pattern 0.0 - 2.0
# with the similarity radius 1.0 in the L1 norm (the 'taxi' metric).
fs_hp = FeatureSet("hydropathy patterns")
fs_hp.add(Feature(get_aa2hydropathy()))
fs_hp.add(Feature(get_aa2hydropathy()).then(pattern_match, pattern=[0.0, 2.0],
                                            metric='taxi', radius=1.0))
result_seq2 = fs_hp(alphasyn_seq)

for seq in result_seq2[:2]:
    print seq

# Calculate bigram frequencies in 'alphasyn_seq':
result_freq = ngram_count(alphasyn_seq, n=2)
print result_freq

# Fit Zipf's law for a trigram distribution in 'amyload_pos_seq':
result_fit = zipf_law_fit(amyload_pos_seq, n=3, verbose=True)

# Calculate the empirical rank-frequency plot:
counts = sorted(result_fit["ngram_counts"], reverse=True)
ranks = range(1, len(counts)+1)

# Calculate the Zipf's law-based approximation:
slope = result_fit["slope"]
harmonic_num = sum([rank**-slope for rank in ranks])
fitted_counts = [(rank**-slope) / harmonic_num * sum(counts) for rank in ranks]

# Generate the rank-frequency plot:
plt.plot(ranks, counts, 'k', label="empirical")
plt.plot(ranks, fitted_counts, 'k--',
         label="Zipf's law\nslope: {:.2f}".format((slope)))
plt.xlabel('rank')
plt.ylabel('count')
plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()
