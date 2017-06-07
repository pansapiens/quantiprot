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
This script is a command-line interface to Quantiprot for quantifying sequences.
"""


from __future__ import print_function

import os
import sys
import argparse

sys.path.insert(0, os.path.abspath('..'))

from quantiprot.utils.io import load_fasta_file
from quantiprot.utils.feature import Feature, FeatureSet
from quantiprot.utils.mapping import simplify
from quantiprot.metrics.aaindex import get_aaindex_file, get_aaindex_www
from quantiprot.metrics.basic import identity
from quantiprot.metrics.basic import sum_absolute, average, average_absolute
from quantiprot.metrics.entropy import entropy
from quantiprot.metrics.rqa import recurrence, determinism, palindromism
from quantiprot.metrics.rqa import ratio_determinism, ratio_palindromism

dsc = 'Quantiprot: convert, discretize and quantify sequences\n\n' \
      'The processing pipeline:\n' \
      ' 1. Load sequences from a fasta file\n' \
      ' 2. Convert amino acid identities to quantitative properties\n' \
      ' 3. Discretize from property values to several classes\n' \
      ' 4. Quantify full length sequences or in sliding window\n' \
      '    using basic arithmetics, entropy or recurrence quantification analysis\n' \
      ' 5. Output in a tabular format to console or to a file.\n' \
      'Steps 2-4 are optional.\n\n' \
      'AAindex codes for amino acid properties are available from\n' \
      ' http://www.genome.jp/aaindex/AAindex/list_of_indices.'

parser = argparse.ArgumentParser(prog='quantify', description=dsc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', required=True, help='fasta file')
parser.add_argument('-o', '--output', default=None, help='tabular output file')

group_convert = parser.add_argument_group('Conversion')
group_convert.add_argument('-p', '--property', default=None,
                           help='AAindex property code')
group_convert.add_argument('-d', '--default', default=None,
                           help='default propert value for non-canonical amino acids')

group_simplify = parser.add_argument_group('Simplification')
group_simplify.add_argument('-s', '--simplify', default=None, choices=['linear', 'kmeans'],
                            help='discretization method')
group_simplify.add_argument('-k', '--classes', default='3',
                            help='num. of classes (default: 3)')
group_simplify.add_argument('-t', '--iterations', default='0',
                            help='num. of iterations for kmeans (default: 0)')

group_quantify = parser.add_argument_group('Quantification')
group_quantify.add_argument('-q', '--quantify', default=None,
                            choices=['sum', 'sum_abs', 'avg', 'avg_abs', 'rec', 'det',
                                     'pal', 'ratio_det', 'ratio_pal', 'entropy'],
                            help='quantification function')
group_quantify.add_argument('-w', '--window', default=0,
                            help='quantification window size'
                                 ' (default: 0 = full sequence)')

group_rqa = parser.add_argument_group('Recurrence Quantification Analysis')
group_rqa.add_argument('-m', '--metric', default='identity',
                       choices=['identity', 'taxi', 'euclid', 'sup', 'inf'],
                       help='rqa: metric (default: identity)')
group_rqa.add_argument('-r', '--radius', default=0.0,
                       help='rqa: similarity radius (default: 0.0)')
group_rqa.add_argument('-e', '--dim', default=1,
                       help='rqa: embedding dimension (default: 1)')
group_rqa.add_argument('-u', '--tau', default=0,
                       help='rqa: embedding delay tau (default: 0)')
group_rqa.add_argument('-l', '--diaglen', default=2,
                       help='rqa: minimal diagonal length for det/pal (default: 2)')

args = parser.parse_args()

# Load the 'input' sequence set
input_seq = load_fasta_file(args.input)

# Retrieve AAindex mapping for the 'property'
if args.property is not None:
    try:
        aa_mapping = get_aaindex_file(args.property)
    except ValueError:
        aa_mapping = get_aaindex_www(args.property)

    # Assign 'default' value for the Mapping
    try:
        aa_mapping.default = float(args.default)
    except (TypeError, ValueError):
        aa_mapping.default = args.default

    # Simplify if and as requested
    if args.simplify is not None:
        aa_mapping = simplify(aa_mapping, aa_mapping.__name__+"/"+args.classes,
                              method=args.simplify, k=int(args.classes),
                              iters=int(args.iterations))

    # Make a Feature from the Mapping
    feat = Feature(aa_mapping)

else:
    feat = Feature(identity)

# Order quantification if and as requested
if args.quantify is not None:
    quantify_method = {'sum': sum,
                       'sum_abs': sum_absolute,
                       'avg': average,
                       'avg_abs': average_absolute,
                       'rec': recurrence,
                       'det': determinism,
                       'pal': palindromism,
                       'ratio_det': ratio_determinism,
                       'ratio_pal': ratio_palindromism,
                       'entropy': entropy,
                      }
    if args.quantify in ['rec', 'det', 'pal', 'ratio_det', 'ratio_pal']:
        feat = feat.then(quantify_method[args.quantify], window=int(args.window),
                         metric=args.metric, radius=float(args.radius),
                         dim=int(args.dim), tau=int(args.tau),
                         det_len=int(args.diaglen), pal_len=int(args.diaglen))
    else:
        feat = feat.then(quantify_method[args.quantify], window=int(args.window))

# Add the Feature to a FeatureSet
fs = FeatureSet("fs")
fs.add(feat)

# And use it to convert the input set
conv_seq = fs(input_seq)

# Finally write to the 'output' file
if args.output is not None:
    outfile = open(args.output, 'w')
else:
    outfile = sys.stdout

for seq in conv_seq:
    outfile.write("{:s}\t{:s}".format(seq.identifier, seq.feature))
    for d in seq.data:
        outfile.write("\t{:}".format(d))
    outfile.write("\n")

if args.output is not None:
    outfile.close()
