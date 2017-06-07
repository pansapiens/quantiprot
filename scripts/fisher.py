# Copyright (c) 2016-2017 Witold Dyrka & Bogumil M. Konopka
# This software was developed in Kotulska Lab at Politechnika Wroclawska.
# This script is a command-line interface to sequence quantification
# using Quantiprot, released under the MIT license:
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
"""
This script is a command-line interface to Quantiprot for comparing two sequence
sets by performing local Fisher tests in a 2d feature space defined by a pair of
quantitative properties.
"""


from __future__ import print_function
import os
import sys
import argparse

import numpy as np

sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.sequence import compact
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.metrics.aaindex import get_aaindex_file, get_aaindex_www
from quantiprot.metrics.basic import identity
from quantiprot.metrics.basic import sum_absolute, average, average_absolute
from quantiprot.metrics.entropy import entropy
from quantiprot.metrics.rqa import recurrence, determinism, palindromism
from quantiprot.metrics.rqa import ratio_determinism, ratio_palindromism
from quantiprot.analysis.fisher import local_fisher_2d, _plot_local_fisher_2d

dsc = 'Quantiprot: feature space analysis using the Fisher test\n\n' \
      'The processing pipeline:\n' \
      ' 1. Load sequences from two fasta files\n' \
      ' 2. Convert amino acid identities to values of two quantitative properties\n' \
      ' 3. Quantify full length sequences\n' \
      '    using basic arithmetics, entropy or recurrence quantification analysis\n' \
      ' 4. Calculate local ratia of number of sequences from each set in parts of\n' \
      '    the feature space and compare with the global ratio in the whole feature' \
      '    space using the Fisher\'s exact test.\n' \
      ' 5. Output to console or to a file:\n' \
      '    - parts of the feature space with significant value of the test\n' \
      '    - list of sequences in these parts\n' \
      ' 6. Visualize results graphically.\n' \
      'Steps 2 and 6 are optional.\n\n' \
      'AAindex codes for amino acid properties are available from\n' \
      ' http://www.genome.jp/aaindex/AAindex/list_of_indices.'

parser = argparse.ArgumentParser(prog='fisher', description=dsc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i1', '--input1', required=True, help='fasta file1')
parser.add_argument('-i2', '--input2', required=True, help='fasta file2')
parser.add_argument('-o', '--output', help='tabular output file')

group_convert = parser.add_argument_group('Conversion')
group_convert.add_argument('-p1', '--property1', default=None,
                           help='AAindex code for property1')
group_convert.add_argument('-p2', '--property2', default=None,
                           help='AAindex code for property2')
group_convert.add_argument('-d1', '--default1', default=None,
                           help='default value of property1 for non-canonical amino acids')
group_convert.add_argument('-d2', '--default2', default=None,
                           help='default value of property2 for non-canonical amino acids')

group_quantify = parser.add_argument_group('Quantification')
group_quantify.add_argument('-q1', '--quantify1', required=True,
                            choices=['sum', 'sum_abs', 'avg', 'avg_abs', 'rec', 'det',
                                     'pal', 'ratio_det', 'ratio_pal', 'entropy'],
                            help='quantification function for property1')
group_quantify.add_argument('-q2', '--quantify2', required=True,
                            choices=['sum', 'sum_abs', 'avg', 'avg_abs', 'rec', 'det',
                                     'pal', 'ratio_det', 'ratio_pal', 'entropy'],
                            help='quantification function for property2')

group_rqa = parser.add_argument_group('Recurrence Quantification Analysis')
group_rqa.add_argument('-m1', '--metric1', default='identity',
                       choices=['identity', 'taxi', 'euclid', 'sup', 'inf'],
                       help='rqa: metric for property1 (default: identity)')
group_rqa.add_argument('-m2', '--metric2', default='identity',
                       choices=['identity', 'taxi', 'euclid', 'sup', 'inf'],
                       help='rqa: metric for property2 (default: identity)')
group_rqa.add_argument('-r1', '--radius1', default=0.0,
                       help='rqa: similarity radius for property1 (default: 0.0)')
group_rqa.add_argument('-r2', '--radius2', default=0.0,
                       help='rqa: similarity radius for property2 (default: 0.0)')
group_rqa.add_argument('-e1', '--dim1', default=1,
                       help='rqa: embedding dimension for property1 (default: 1)')
group_rqa.add_argument('-e2', '--dim2', default=1,
                       help='rqa: embedding dimension for property2 (default: 1)')
group_rqa.add_argument('-u1', '--tau1', default=0,
                       help='rqa: embedding delay tau for property1 (default: 0)')
group_rqa.add_argument('-u2', '--tau2', default=0,
                       help='rqa: embedding delay tau for property2 (default: 0)')
group_rqa.add_argument('-l1', '--diaglen1', default=2,
                       help='rqa: minimal diagonal length for det/pal'
                            ' for property1 (default: 2)')
group_rqa.add_argument('-l2', '--diaglen2', default=2,
                       help='rqa: minimal diagonal length for det/pal'
                            ' for property2 (default: 2)')

group_fisher = parser.add_argument_group('Local Fisher test')
group_fisher.add_argument('-b1', '--bins1', default=10,
                          help='total number of bins for property1 (default: 10)')
group_fisher.add_argument('-b2', '--bins2', default=10,
                          help='total number of bins for property2 (default: 10)')
group_fisher.add_argument('-v1', '--overlap1', default=1,
                          help='number of bins of property 1 to aggregate'
                               ' for single test (default: 1)')
group_fisher.add_argument('-v2', '--overlap2', default=1,
                          help='number of bins of property 2 to aggregate'
                               ' for single test (default: 1)')
group_fisher.add_argument('-th', '--threshold', default=0.05,
                          help='p-value threshold (default: 0.05)')
group_fisher.add_argument('-V', '--visualize', default=False, action='store_true',
                          help='visualize results (default: False)')

quantify_method = {'sum': sum,
                   'sum_abs': sum_absolute,
                   'avg': average,
                   'avg_abs': average_absolute,
                   'rec': recurrence,
                   'det': determinism,
                   'pal': palindromism,
                   'ratio_det': ratio_determinism,
                   'ratio_pal': ratio_palindromism,
                   'entropy': entropy
                  }

args = parser.parse_args()

# Load the input sequence sets
input_seq1 = load_fasta_file(args.input1)
input_seq2 = load_fasta_file(args.input2)

# Retrieve AAindex mappings for the properties if and as requested
# property1
if args.property1 is not None:
    try:
        aa_mapping1 = get_aaindex_file(args.property1)
    except ValueError:
        aa_mapping1 = get_aaindex_www(args.property1)
    try:
        aa_mapping1.default = float(args.default1)
    except (TypeError, ValueError):
        aa_mapping1.default = args.default1

    feat1 = Feature(aa_mapping1)
else:
    feat1 = Feature(identity)
# property2
if args.property2 is not None:
    try:
        aa_mapping2 = get_aaindex_file(args.property2)
    except ValueError:
        aa_mapping2 = get_aaindex_www(args.property2)
    try:
        aa_mapping2.default = float(args.default2)
    except (TypeError, ValueError):
        aa_mapping2.default = args.default2

    feat2 = Feature(aa_mapping2)
else:
    feat2 = Feature(identity)

# Make features based on mappings:
# feature1
if args.quantify1 in ['rec', 'det', 'pal', 'ratio_det', 'ratio_pal']:
    feat1 = feat1.then(quantify_method[args.quantify1],
                       metric=args.metric1, radius=float(args.radius1),
                       dim=int(args.dim1), tau=int(args.tau1),
                       det_len=int(args.diaglen1), pal_len=int(args.diaglen1))
else:
    feat1 = feat1.then(quantify_method[args.quantify1])
# feature2
if args.quantify2 in ['rec', 'det', 'pal', 'ratio_det', 'ratio_pal']:
    feat2 = feat2.then(quantify_method[args.quantify2],
                       metric=args.metric2, radius=float(args.radius2),
                       dim=int(args.dim2), tau=int(args.tau2),
                       det_len=int(args.diaglen2), pal_len=int(args.diaglen2))
else:
    feat2 = feat2.then(quantify_method[args.quantify2])

# Add the Features to a FeatureSet
fs = FeatureSet("fs")
fs.add(feat1)
fs.add(feat2)

# And use it to convert the input sets
conv_seq1 = fs(input_seq1)
conv_seq2 = fs(input_seq2)

# Do local Fisher:
result = local_fisher_2d(conv_seq1, conv_seq2,
                         windows_per_frame=(int(args.bins1)/int(args.overlap1),
                                            int(args.bins2)/int(args.overlap2)),
                         overlap_factor=(int(args.overlap1), int(args.overlap2))
                        )

# Compact converted sets for convenience
compact_conv_seq1 = compact(conv_seq1)
compact_conv_seq2 = compact(conv_seq2)

# Finally write to the output file
if args.output is not None:
    outfile = open(args.output, 'w')
else:
    outfile = sys.stdout

# Write some headers
outfile.write("#Quantiprot: feature space analysis with the Fisher test\n")
outfile.write("\n#Parameters:\n")

# Write script parameters
for key in sorted(vars(args)):
    outfile.write("#-{:}={:}\n".format(key, vars(args)[key]))

# List ranges with statistically significant discrepancies between two inputs
outfile.write("\n#Ranges with statistically significant differences:\n")

for pos in np.argwhere(result["p_value"] < float(args.threshold)):
    # Give ranges coordinates
    outfile.write("\n#plot_coordinates=[{:}, {:}]\n".format(pos[0], pos[1]))
    outfile.write("#property1=({:}, {:})\n"
                  "#property2=({:}, {:})\n".format(result["_bin_ranges_x"][pos[0]],
                                                   result["_bin_ranges_x"][pos[0]+1],
                                                   result["_bin_ranges_y"][pos[1]],
                                                   result["_bin_ranges_y"][pos[1]+1]))

    # List sequences from input1 in the range
    outfile.write("#input1_count={:}\n".format(int(result["w_counts1"][pos[0], pos[1]])))
    for s in range(len(input_seq1)): # compacted converted set preserves sequence order of the input
        if compact_conv_seq1[s].data[0] >= result["_bin_ranges_x"][pos[0]] and \
           compact_conv_seq1[s].data[0] < result["_bin_ranges_x"][pos[0]+int(args.overlap1)] and\
           compact_conv_seq1[s].data[1] >= result["_bin_ranges_y"][pos[1]] and \
           compact_conv_seq1[s].data[1] < result["_bin_ranges_y"][pos[1]+int(args.overlap2)]:
            outfile.write(input_seq1[s].fasta())

    # List sequences from input2 in the range
    outfile.write("#input2_count={:}\n".format(int(result["w_counts2"][pos[0], pos[1]])))
    for s in range(len(input_seq2)): # compacted converted set preserves sequence order of the input
        if compact_conv_seq2[s].data[0] >= result["_bin_ranges_x"][pos[0]] and \
           compact_conv_seq2[s].data[0] < result["_bin_ranges_x"][pos[0]+int(args.overlap1)] and\
           compact_conv_seq2[s].data[1] >= result["_bin_ranges_y"][pos[1]] and \
           compact_conv_seq2[s].data[1] < result["_bin_ranges_y"][pos[1]+int(args.overlap2)]:
            outfile.write(input_seq2[s].fasta())

    # List OR and p-value
    outfile.write("#OR={:}\n"
                  "#p-value={:}\n".format(result["odds_ratio"][pos[0], pos[1]],
                                          result["p_value"][pos[0], pos[1]]))

if args.output is not None:
    outfile.close()

# Plot local Fisher if requested:
if bool(args.visualize) is True:
    _plot_local_fisher_2d(result, xlabel=feat1.name,
                          ylabel=feat2.name,
                          pop1_label="input1",
                          pop2_label="input2")
