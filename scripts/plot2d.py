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
This script is a command-line interface to Quantiprot for plotting sequence sets
as points in a 2d feature space defined by a pair of quantitative properties.
"""


from __future__ import print_function
import os
import sys
import argparse

sys.path.insert(0, os.path.abspath('..'))

from matplotlib import pyplot as plt

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.metrics.aaindex import get_aaindex_file, get_aaindex_www
from quantiprot.metrics.basic import sum_absolute, average, average_absolute
from quantiprot.metrics.basic import identity
from quantiprot.metrics.entropy import entropy
from quantiprot.metrics.rqa import recurrence, determinism, palindromism
from quantiprot.metrics.rqa import ratio_determinism, ratio_palindromism

dsc = 'Quantiprot: plot sequence sets as points in a 2d feature space\n\n' \
      'The processing pipeline:\n' \
      ' 1. Load sequences from (a) fasta file(s)\n' \
      ' 2. Convert amino acid identities to values of two quantitative properties\n' \
      ' 3. Quantify full length sequences\n' \
      '    using basic arithmetics, entropy or recurrence quantification analysis\n' \
      ' 4. Plot as points in 2d feature space defined by the properties.\n' \
      'Step 2 is optional.\n\n' \
      'AAindex codes for amino acid properties are available from\n' \
      ' http://www.genome.jp/aaindex/AAindex/list_of_indices.'

parser = argparse.ArgumentParser(prog='plot2d', description=dsc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i1', '--input1', required=True,
                    help='fasta file1')
parser.add_argument('-i2', '--input2', default=None,
                    help='fasta file2')

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

# Add the features to a FeatureSet
fs = FeatureSet("fs")
fs.add(feat1)
fs.add(feat2)

# Convert and plot input1 sequences in the 2d space
input_seq1 = load_fasta_file(args.input1)
conv_seq1 = fs(input_seq1)
conv_data1_x = conv_seq1.columns(feature=feat1.name)[0]
conv_data1_y = conv_seq1.columns(feature=feat2.name)[0]
plt.plot(conv_data1_x, conv_data1_y, '.', label="input1")

# Convert and plot input1 sequences in the 2d space
if args.input2 is not None:
    input_seq2 = load_fasta_file(args.input2)
    conv_seq2 = fs(input_seq2)
    conv_data2_x = conv_seq2.columns(feature=feat1.name)[0]
    conv_data2_y = conv_seq2.columns(feature=feat2.name)[0]
    plt.plot(conv_data2_x, conv_data2_y, '.', label="input2")

# Show legend and labels
plt.xlabel(feat1.name)
plt.ylabel(feat2.name)
plt.legend()
plt.show()
