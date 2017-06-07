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
This script is a command-line interface to Quantiprot for representing protein
sequences by n-gram vectors.
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
from quantiprot.metrics.ngram import NgramFeatureSet

dsc = 'Quantiprot: n-gram representation of sequences\n\n' \
      'The processing pipeline:\n' \
      ' 1. Load sequences from a fasta file\n' \
      ' 2. Convert amino acid identities to quantitative properties\n' \
      ' 3. Discretize from property values to several classes\n' \
      ' 4. Count occurrences of all n-grams for given\n' \
      ' 5. Output in a tabular format to console or to a file.\n' \
      'Steps 2-3 are optional.\n\n' \
      'Note that going beyond trigrams can be extremaly resource consuming.\n\n' \
      'AAindex codes for amino acid properties are available from\n' \
      ' http://www.genome.jp/aaindex/AAindex/list_of_indices.'

parser = argparse.ArgumentParser(prog='ngrams', description=dsc,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('-i', '--input', required=True, help='fasta file')
parser.add_argument('-o', '--output', default=None, help='tabular output file')

group_convert = parser.add_argument_group('Conversion')
group_convert.add_argument('-p', '--property', default=None,
                           help='AAindex property code')
group_convert.add_argument('-d', '--default', default=None,
                           help='defaul property value for non-canonical amino acids')

group_simplify = parser.add_argument_group('Simplification')
group_simplify.add_argument('-s', '--simplify', default=None,
                            choices=['linear', 'kmeans'],
                            help='discretization method (default: not specified)')
group_simplify.add_argument('-k', '--classes', default=3,
                            help='num. of classes (default: 3)')
group_simplify.add_argument('-t', '--iterations', default=0,
                            help='num. of iterations for kmeans (default: 0)')

group_ngrams = parser.add_argument_group('N-grams')
group_ngrams.add_argument('-n', '--n', default='1', help='n-gram size (default: 1)')
group_ngrams.add_argument('-m', '--metric', default='identity',
                          choices=['identity', 'taxi', 'euclid', 'sup', 'inf'],
                          help='metric for matching n-grams (default: identity)')
group_ngrams.add_argument('-r', '--radius', default=0.0,
                          help='similarity radius (default: 0.0)')

args = parser.parse_args()

# Load the 'input' sequence set
input_seq = load_fasta_file(args.input, unique=False)

# Retrieve AAindex mapping for the 'property'
if args.property is not None:
    try:
        aa_mapping = get_aaindex_file(args.property)
    except ValueError:
        aa_mapping = get_aaindex_www(args.property)

    # Simplify if and as requested
    if args.simplify is not None:
        aa_mapping = simplify(aa_mapping, aa_mapping.__name__+"/"+args.classes,
                              method=args.simplify, k=int(args.classes),
                              iters=int(args.iterations))

    # Assign 'default' value for the Mapping
    try:
        aa_mapping.default = float(args.default)
    except (TypeError, ValueError):
        aa_mapping.default = args.default

    # Make a Feature from the Mapping
    feat = Feature(aa_mapping)
else:
    feat = Feature(identity)

# Add the Feature to a FeatureSet
fs = FeatureSet("fs")
fs.add(feat)

# And use it to convert the input set
conv_seq = fs(input_seq)

# Get the alphabet of the converted set
alphabet = list(set([element for seq in conv_seq for element in seq.data]))

# Prepare the n-gram counts extractor
nfs = NgramFeatureSet('ngram_'+args.n, n=int(args.n), alphabet=alphabet,
                      mode='count', window=0, simple_names=True,
                      metric=args.metric, radius=float(args.radius))

# Extract n-grams for each sequence separately
ngram_seq = nfs(conv_seq)

# Finally write to the 'output' file
if args.output is not None:
    outfile = open(args.output, 'w')
else:
    outfile = sys.stdout

outfile.write("identifier")
for f in nfs.feature_list:
    outfile.write("\t{:}".format(f))
outfile.write("\n")

# The following code is faster than cleaner approach with compacting ngram_seq
# It takes the advantage from the fact that the inter-sequence order in ngram_seq
# preserves the order of the input set and the intra-sequence order in ngram_seq
# preserves the order of nfs feature list.
prev = None
for seq in ngram_seq:
    if seq.identifier != prev:
        if prev is not None:
            outfile.write("\n")
        outfile.write("{:s}".format(seq.identifier))
        prev = seq.identifier
    outfile.write("\t{:d}".format(int(seq.data[0])))
outfile.write("\n")

if args.output is not None:
    outfile.close()
