import os
import sys
sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.metrics.aaindex import get_aa2volume, get_aa2hydropathy
from quantiprot.metrics.basic import average

# Local Fisher-test related imports:
from quantiprot.analysis.fisher import local_fisher_2d, _plot_local_fisher_2d

from matplotlib import pyplot as plt

# Load sets of amyloidogenic and non-amyloidogenic peptides:
amyload_pos_seq = load_fasta_file("data/Amyload_positive.fasta")
amyload_neg_seq = load_fasta_file("data/Amyload_negative.fasta")

# Calculate quantitive features: volume and hydropathy
mean_volume = Feature(get_aa2volume()).then(average)
mean_hydropathy = Feature(get_aa2hydropathy()).then(average)

fs = FeatureSet("volume'n'hydropathy")
fs.add(mean_volume)
fs.add(mean_hydropathy)

amyload_pos_conv_seq = fs(amyload_pos_seq)
amyload_neg_conv_seq = fs(amyload_neg_seq)

# Do local Fisher:
result = local_fisher_2d(amyload_pos_conv_seq, amyload_neg_conv_seq,
                         windows_per_frame=5, overlap_factor=5)

# Plot local Fisher:
_plot_local_fisher_2d(result, xlabel="mean volume",
                              ylabel="mean hydropathy",
                              pop1_label="amyloids",
                              pop2_label="non amyloids")
