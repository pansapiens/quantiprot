# Copyright (c) 2016-2017 Witold Dyrka & Bogumil M. Konopka
# This software was developed in Kotulska Lab at Politechnika Wroclawska.
# This module is a part of Quantiprot, released under the MIT license:
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
This module provides functions and a class for finding patterns and ngrams
in sequences and sequence sets.

Classes:
    NgramFeatureSet: manage a set of ngram-based sequence features.

Functions:
    _list_distance: calculate distance between two lists of the same length.
    _pattern_common: do common preprocessing steps for pattern_match
                     and pattern_count.
    pattern_count: count occurrences of a given pattern.
    pattern_match: find matches of a given pattern.

"""

import itertools
import copy

import numpy as np

from quantiprot.utils.sequence import Sequence, SequenceSet
from quantiprot.utils.feature import Feature


def _list_distance(list1, list2, metric):
    """
    Calculate distance between two lists of the same length
    according to a given metric.

    Assumes that lists are of the same length.

    Args:
        list1 (list): first input list.
        list2 (list): second input list.
        metric (str): 'identity' counts identical positions,
                       'euclid' calculates the Euclidean distance (L2 norm),
                       'taxi' calculates the taxicab (Manhattan) distance
                           (L1 norm),
                       'sup' returns maximum distance between positions,
                       'inf' returns minimum distance between positions.

    Returns distance between the two lists.

    Raises:
        ValueError if unsupported distance metric used.
    """
    # Equal length of both lists is assumed
    distance = 0
    if metric == 'taxi':
        distance = sum([abs(list1[pos] - list2[pos])
                        for pos in range(len(list1))])
        return distance
    elif metric == 'euclid':
        distance = sum([(list1[pos] - list2[pos]) ** 2
                        for pos in range(len(list1))])
        return distance ** 0.5
    elif metric == 'identity':
        distance = sum([list1[pos] != list2[pos]
                        for pos in range(len(list1))])
        return distance
    elif metric == 'sup':
        distance = max([abs(list1[pos] - list2[pos])
                        for pos in range(len(list1))])
        return distance
    elif metric == 'inf':
        distance = min([abs(list1[pos] - list2[pos])
                        for pos in range(len(list1))])
        return distance
    else:
        raise ValueError("Proper distance calculation metrics are \'taxi\',\
                          \'euclid\', \'identity\', \'sup\', or \'inf\'.")


def _pattern_common(**params):
    """
    Do common preprocessing steps for pattern_match and pattern_count.

    Not really useful on its own.

    Args:
        data (list): values.
        params (kwargs):
            pattern (str or list): the pattern to be sought in data (obligatory)
            metric (str): 'identity' counts identical positions,
                           'euclid' calculates the Euclidean distance (L2 norm),
                           'taxi' calculates the taxicab (Manhattan) distance
                               (L1 norm).
                           'sup' returns maximum distance between positions,
                           'inf' returns minimum distance between positions.
                           Only 'identity' can be used with non-numerical data.
            radius (number): the similarity cutoff (non-negative)

    Returns:
        pattern (list): the pattern to be sought in data as a list.
        patlen (int): length of the pattern.
        metric (str): name of the metric for calculating match similarity.

   Raises:
        NameError when 'pattern' is not given,
        TypeError if 'pattern' is neither string nor list,
        ValueError if 'radius' is negative or unsupported distance method used.
    """

    if 'pattern' in params:
        if isinstance(params['pattern'], list):
            pattern = params['pattern']
        elif isinstance(params['pattern'], str):
            pattern = list(params['pattern'])
        else:
            raise TypeError("The pattern should be either list or a string.")
    else:
        raise NameError("No pattern provided.")

    patlen = len(pattern)

    radius = params['radius'] if 'radius' in params else 0

    if radius < 0:
        raise ValueError("Similarity radius cannot be negative.")

    metric = params['metric'] if 'metric' in params else 'identity'

    if metric not in ['identity', 'taxi', 'euclid', 'sup', 'inf']:
        raise ValueError("Unsupported distance metric.")

    return pattern, patlen, radius, metric


def pattern_count(data, **params):
    """
    Count occurrences of a given pattern.

    Args:
        data (list): values.
        params (kwargs):
            pattern (str or list): the pattern to be sought in data (obligatory)
            metric (str): 'identity' counts identical positions,
                           'euclid' calculates the Euclidean distance (L2 norm),
                           'taxi' calculates the taxicab (Manhattan) distance
                               (L1 norm).
                           'sup' returns maximum distance between positions,
                           'inf' returns minimum distance between positions.
                           Only 'identity' can be used with non-numerical data.
            radius (number): the similarity cutoff (non-negative)
            normalized (bool): whether the number of occurrences is to be
                               divided by the maximum number of occurrences.
                               (default:False)

    Returns the number of occurrences of the pattern in the data.

    Invokes internal function '_pattern_common', which raises:
        NameError when 'pattern' is not given,
        TypeError if 'pattern' is neither string nor list,
        ValueError if 'radius' is negative or unsupported distance method used.
    """

    pattern, patlen, radius, metric = _pattern_common(**params)

    normalized = params['normalized'] if 'normalized' in params else False

    counts = 0

    for pos in range(len(data) - patlen + 1):
        if _list_distance(data[pos:pos + patlen], pattern, metric) <= radius:
            counts += 1

    return counts if not normalized \
                  else 1.0 * counts / (len(data) - patlen + 1)


def pattern_match(data, **params):
    """
    Find matches of a given pattern.

    Args:
        data (list): values.
        params (kwargs):
            pattern (str or list): the pattern to be sought in data (obligatory)
            metric (str): 'identity' counts identical positions,
                           'euclid' calculates the Euclidean distance (L2 norm),
                           'taxi' calculates the taxicab (Manhattan) distance
                               (L1 norm).
                           'sup' returns maximum distance between positions,
                           'inf' returns minimum distance between positions.
                           Only 'identity' can be used with non-numerical data.
            radius (number): the similarity cutoff (non-negative).
            padding (bool): whether the rightmost positions of the data are to
                            be padded with zeros (default:False).

    Returns the list of matches (zeros and ones) to 'pattern' in the data.

    Invokes internal function '_pattern_common', which raises
        NameError when 'pattern' is not given,
        TypeError if 'pattern' is neither string nor list,
        ValueError if 'radius' is negative or unsupported distance method used.
    """
    pattern, patlen, radius, metric = _pattern_common(**params)

    padding = params['padding'] if 'padding' in params else False

    matches = []
    for pos in range(len(data) - patlen + 1):
        if _list_distance(data[pos:pos + patlen], pattern, metric) <= radius:
            matches.append(1)
        else:
            matches.append(0)
    if padding:
        for pos in range(len(data) - patlen + 1, len(data)):
            matches.append(0)

    return matches


class NgramFeatureSet(object):
    """
    Manage a set of ngram-based sequence features in efficient way.

    Attributes:
        name (str): name for the feature set.
        n (int): ngram size (positive).
        alphabet (list): alphabet to construct ngrams from.
        mode (str): 'match' finds matches to ngrams
                    'count' counts occurrences of ngrams
        window (int): length of the window over which the feature is calculated
                      (non-negative, default=0 means the whole sequence).
        simple_names (bool): whether to use 'mode' for naming n-gram features
                             (default: False) or only n-gram tuples (True).
        params (**kwargs):
            metric (str): 'identity' counts identical positions,
                           'euclid' calculates the Euclidean distance (L2 norm),
                           'taxi' calculates the taxicab (Manhattan) distance
                               (L1 norm).
                           'sup' returns maximum distance between positions,
                           'inf' returns minimum distance between positions.
                           Only 'identity' can be used with non-numerical data.
            radius (number): the similarity cutoff (non-negative, default:0)
            padding (bool): whether the rightmost positions are to be padded
                            with zeros (applicable to 'match', default:False)
            normalized (bool): whether the number of occurrences is to be
                               divided by the maximum number of occurrences.
                               (applicable to 'count', default:False)

        Raises ValueError if improper values of parameters are used.
    """

    def __init__(self, name, n, alphabet, mode, window=0, simple_names=False,
                 **params):

        if n < 1:
            raise ValueError("There are no ngrams for n < 1.")

        if mode != 'match' and mode != 'count':
            raise ValueError("Unknown mode. Available are \
                             \'match\' or \'count\'.")

        if window < 0:
            raise ValueError("Window length for a feature has to be either "
                             "positive or zero for full sequence length!")

        if 'metric' in params and \
           params['metric'] not in ['identity', 'taxi', 'euclid', 'sup', 'inf']:
            raise ValueError("Unsupported distance metric.")

        if 'radius' in params and params['radius'] < 0:
            raise ValueError("Similarity radius cannot be negative.")

        self.name = name
        self.n = n
        self.alphabet = alphabet
        self.mode = mode
        self.window = window
        self.simple_names = simple_names
        self.metric = params['metric'] if 'metric' in params else 'identity'
        self.radius = params['radius'] if 'radius' in params else 0
        self.normalized = params['normalized'] \
                          if 'normalized' in params else False
        self.padding = params['padding'] if 'padding' in params else False
        self.feature_list = list(itertools.product(alphabet, repeat=n))
        self._post_feat = None
        self._post_feat_name = ""

    def __call__(self, seq_set, window=-1):
        """
        Calculate matches or counts of all ngrams for all
        sequences in a given sequence set using a given sliding window.

        Args:
            seq_set (SequenceSet): a SequenceSet instance
            window (int): length of the sliding window. Defaults to self.window
                          (window=-1). To calculate a single value for
                          the whole sequence use window=0.

        Returns a set of sequences of quantitative values (ngram matches or
            occurences).

        Raises ValueError if unknown unsupported mode is set.
        """

        if window == -1:
            window = self.window

        result_set = SequenceSet(seq_set.name + ":" + self.name,
                                 unique=seq_set.unique)

        for seq in seq_set:

            matches = np.zeros((len(seq.data), len(self.feature_list)))

            if self.radius == 0:
                for pos in range(len(seq.data) - self.n + 1):
                    matches[pos, self.feature_list.index(
                        tuple(seq.data[pos:pos + self.n]))] = 1
            else:
                for pos in range(len(seq.data) - self.n + 1):
                    for feat_pos in range(len(self.feature_list)):
                        if _list_distance(seq.data[pos:pos + self.n],
                                          self.feature_list[feat_pos],
                                          self.metric) <= self.radius:
                            matches[pos, feat_pos] = 1

            for feat_pos in range(len(self.feature_list)):
                if window == 0:
                    if self.mode == 'match':
                        outdata = (matches[:, feat_pos]).tolist()
                        if not self.padding:
                            outdata = outdata[:-(self.n - 1)]
                    elif self.mode == 'count':
                        outdata = matches[:, feat_pos].sum()
                        if self.normalized:
                            outdata /= (len(seq.data) - self.n + 1)
                        outdata = [outdata.tolist()]
                    else:
                        raise ValueError("Unknown mode. Available are \
                                         \'match\' or \'count\'.")
                else:
                    outdata = []
                    if self.mode == 'match':
                        for win_pos in range(len(seq.data) - window + 1):
                            single_outdata = (matches[win_pos:win_pos + window,
                                                      feat_pos]).tolist()
                            if not self.padding:
                                single_outdata = single_outdata[:-(self.n - 1)]
                            outdata.append(single_outdata)
                    elif self.mode == 'count':
                        for win_pos in range(len(seq.data) - window + 1):
                            single_outdata = matches[win_pos:win_pos + window,
                                                     feat_pos].sum()
                            if self.normalized:
                                single_outdata /= (len(seq.data) - self.n + 1)
                            outdata.append(single_outdata.tolist())
                    else:
                        raise ValueError("Unknown mode. Available are \
                                         \'match\' or \'count\'.")

                my_seq = Sequence(seq.identifier, self.name, outdata)

                if self._post_feat is not None:
                    my_seq = self._post_feat(my_seq)

                my_feat_name = str(self.feature_list[feat_pos])
                if self.simple_names:
                    my_seq.feature = my_feat_name
                else:
                    my_seq.feature = my_feat_name + "_" + self.mode

                my_seq.feature += self._post_feat_name

                result_set.add(my_seq)


        return result_set

    def then_all(self, function, name=None, window=0, **params):
        """
        Define a post-processor feature.

        The method can either take a Feature or a function. In the former case
        the post-processor is a deep copy of the input Feature.

        The method modifies the self object.

        Args:
            function (function): Feature or function to serve a post-processor.
                Only when 'function' is not a Feature, the following arguments
                are taken into account:
            name (str): name for the feature. If None (default),
                        the 'function' name is used instead as 'name'.
            window (int): length of the window over which the feature is
                          calculated. Defaults to the whole sequence (window=0).
            params (**kwargs): arbitrary params to be passed to the function.

         Returns the self to allow feature chaining.
        """
        if isinstance(function, Feature):
            post_feat = copy.deepcopy(function)
        else:
            my_name = name if name is not None else function.__name__
            post_feat = Feature(function, name=my_name,
                                window=window, **params)

        if self._post_feat is None:
            self._post_feat = post_feat
        else:
            self._post_feat.then(post_feat)

        self._post_feat_name += ">" + post_feat.name

        return self

    def __len__(self):
        return len(self.feature_list)

    def __repr__(self):
        return 'Quantiprot Ngram Feature Set object {:s}'.format(self.name)

    def __str__(self):
        return "Quantiprot Ngram Feature Set object\n" + \
               "  'name': '{:s}'\n".format(self.name) + \
               "  'n'gram size: {:d}\n".format(self.n) + \
               "  'alphabet': {:}\n".format(self.alphabet) + \
               "  'mode': '{:s}'\n".format(self.mode) + \
               "  scanning 'window': {:d}\n".format(self.window) + \
               "  distance 'metric': '{:s}'\n".format(self.metric) + \
               "  similarity 'radius': {:f}\n".format(self.radius) + \
               "  'normalized' counts: {:}\n".format(self.normalized) + \
               "  sequence 'padding': {:}\n".format(self.padding) + \
               "  'simple names': '{:}'".format(self.simple_names)
