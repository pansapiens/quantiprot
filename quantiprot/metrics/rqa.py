# Copyright (c) 2016-2017 Witold Dyrka & Bogumil M. Konopka & Marta Marciniak
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
This module provides a class and functions for quantifying sequences using
measures related to Recurrence Quantification Analysis (RQA).

Classes:
    RQAFeatureSet: manage a set of RQA-based sequence features.

Functions:
    _uniq: return a list of unique elements of a list of elements.
    _recurrence_common: calculate recurrence plot for the data.
    recurrence: calculate recurrence rate of the data.
    determinism: calculate determinism of the data.
    determinism_ratio: calculate ratio of determinism and recurrence of the
        data.
    palindromism: calculate palindromism of the data.
    palindromism_ratio: calculate ratio of palindromism and recurrence of the
        data.
"""

import copy

import numpy as np

from quantiprot.utils.sequence import Sequence, SequenceSet
from quantiprot.utils.feature import Feature
from quantiprot.metrics.ngram import _list_distance


def _uniq(input_list):
    """
    Return a list of unique elements of a list of elements
    while preserving the order.

    The function adopted from Lukas comment of 13/01/2016
    at https://www.peterbe.com/plog/uniqifiers-benchmark.

    Args:
        input_list(list): the processed list of elements.

    Returns the list of unique elements of the input list.
    """
    seen = set()
    seen_add = seen.add
    return [x for x in input_list if not (x in seen or seen_add(x))]


def _recurrence_common(data, **params):
    """
    Calculate recurrence plot for the data

    Args:
        data (list): values.
        params (kwargs):

            Settings for the embedding
            dim (float):Embedding dimension (default 1)
            tau (float):Embedding delay (default 0)

            Settings for the recurrence plot.
            radius (float): Fixed threshold (default 0)
            metric (string): Distance metric in the phase space.
                Possible choices: 'identity' (default), 'euclid' and 'taxi'.

            full (bool): whether to include the main diagonal or not
                (default: False)

    Returns the lower left of the recurrence table for the data.
    """

    dim = params['dim'] if 'dim' in params else 1
    tau = params['tau'] if 'tau' in params else 1
    radius = params['radius'] if 'radius' in params else 0
    metric = params['metric'] if 'metric' in params else 'identity'
    full = params['full'] if 'full' in params else False

    len_data = len(data)
    num_points = len_data - (dim-1) * tau
    data_embed = np.zeros((num_points, dim), dtype=np.array(data).dtype)
    rec_table = np.zeros((num_points, num_points))

    if dim == 1:
        for pos in range(num_points):
            data_embed[pos, 0] = data[pos]
    elif tau == 1:
        for pos in range(num_points):
            data_embed[pos, :] = data[pos:pos + dim]
    else:
        for pos in range(num_points):
            for pos_embed in range(dim):
                data_embed[pos, pos_embed] = data[pos + pos_embed * tau]

    for pos_x in range(num_points):
        for pos_y in range(pos_x+(1*full)):
            rec_table[pos_x, pos_y] = _list_distance(data_embed[pos_x],
                                                     data_embed[pos_y],
                                                     metric=metric) <= radius
            #rec_table[pos_y, pos_x] = rec_table[pos_x, pos_y]

    return rec_table


def recurrence(data, _rec_table=None, **params):
    """
    Calculate recurrence of the data.

    Args:
        data (list): values.
        params (kwargs):
            full (bool): whether to include the main diagonal or not
                (default: False)

            Settings for the embedding
            dim (float):Embedding dimension (default 1)
            tau (float):Embedding delay (default 0)

            Settings for the recurrence plot.
            radius (float): Fixed threshold (default 0)
            metric (string): Distance metric in the phase space.
                Possible choices: 'identity' (default), 'euclid', 'taxi',
                                  'sup', and 'inf'.

    Returns the recurrence of the data.
    """

    full = params['full'] if 'full' in params else False

    rec_table = _recurrence_common(data, **params) \
                if _rec_table is None else _rec_table

    data_len = len(rec_table)
    rec_sum = np.sum(rec_table)

    if full:
        rec_sum *= 2
        rec_sum -= np.sum(rec_table.diagonal())
        rec_rate = 1.0 * rec_sum / (data_len**2)
    else:
        rec_rate = 2.0 * rec_sum / (data_len**2-data_len)

    return rec_rate.tolist()


def determinism(data, _rec_table=None, **params):
    """
    Calculate determinism of the data.

    Args:
        data (list): values.
        params (kwargs):
            det_len (int): minimum diagonal length (default: 2)
            full (bool): whether to include the main diagonal or not
                (default: False)

            Settings for the embedding
            dim (float):Embedding dimension (default 1)
            tau (float):Embedding delay (default 0)

            Settings for the recurrence plot.
            radius (float): Fixed threshold (default 0)
            metric (string): Distance metric in the phase space.
                Possible choices: 'identity' (default), 'euclid', 'taxi',
                                  'sup', and 'inf'.

    Returns the determinism of the data.
    """

    det_len = params['det_len'] if 'det_len' in params else 2
    full = params['full'] if 'full' in params else False

    rec_table = _recurrence_common(data, **params) \
                if _rec_table is None else _rec_table

    rec_sum = np.sum(rec_table)

    data_len = len(rec_table)

    diagonals = np.zeros((data_len, data_len))
    diagonals[0, :] = rec_table[0, :].copy()
    diagonals[:, 0] = rec_table[:, 0].copy()

    for pos_x in range(1, data_len):
        for pos_y in range(1, pos_x+(1*full)):
            if rec_table[pos_x, pos_y]:
                diagonals[pos_x, pos_y] = diagonals[pos_x-1, pos_y-1] \
                                        + rec_table[pos_x, pos_y]
                diagonals[pos_x-1, pos_y-1] = 0

    det_sum = np.sum(diagonals[np.where(diagonals >= det_len)])

    rec_sum = np.sum(rec_table)

    if full:
        diagonal = diagonals.diagonal()
        det_sum *= 2
        det_sum -= np.sum(diagonal[np.where(diagonal >= det_len)])
        rec_sum *= 2
        rec_sum -= np.sum(rec_table.diagonal())

    det = det_sum / rec_sum

    return det.tolist()


def ratio_determinism(data, _rec_table=None, **params):
    """
    Calculate ratio of determinism and recurrence of the data.

    Args:
        data (list): values.
        params (kwargs):
            det_len (int): minimum diagonal length (default: 2)
            full (bool): whether to include the main diagonal or not
                (default: False)

            Settings for the embedding
            dim (float):Embedding dimension (default 1)
            tau (float):Embedding delay (default 0)

            Settings for the recurrence plot.
            radius (float): Fixed threshold (default 0)
            metric (string): Distance metric in the phase space.
                Possible choices: 'identity' (default), 'euclid', 'taxi',
                                  'sup', and 'inf'.

    Returns the ratio of determinism and recurrence of the data.
    """

    rec_table = _recurrence_common(data, **params) \
                if _rec_table is None else _rec_table

    ratio_det = np.array(determinism(data, _rec_table=rec_table, **params)) / \
           np.array(recurrence(data, _rec_table=rec_table, **params))

    return ratio_det.tolist()


def palindromism(data, _rec_table=None, **params):
    """
    Calculate palindromism of the data.

    Palindromism of the data can be defined in terms of the recurrence plot
    as length of lines at antidiagonals longer than the given length divided by
    the length of all lines at antidiagonals (or all recurrence points). This
    definition counts AA, ACA, and ACCA as palindroms of lengths 2, 3 and 4,
    respectively. The function calculates this version if 'full' is True.
    Alternatively, palindromism can be defined as length of lines at anti-
    diagonals longer than the given length, s.t. do not overlap nor cross the
    main diagonal, divided by the number of all recurrence points except the
    main diagonal. This definition counts AA, ACA, and ACCA as palindroms of
    lengths 1, 1 and 2, respectively. In other words, only the reoccuring part
    of the pattern is counted. The function calculates this version if 'full'
    is False (default).

    Args:
        data (list): values.
        params (kwargs):
            pal_len (int): minimum palindromic length (see above, default: 2)
            full (bool): whether to count unpaired middle elements to
                palindromism length (default: False)

            Settings for the embedding
            dim (float): Embedding dimension (default 1)
            tau (float): Embedding delay (default 0)

            Settings for the recurrence plot
            radius (float): Fixed threshold (default 0)
            metric (string): Distance metric in the phase space.
                Possible choices: 'identity' (default), 'euclid', 'taxi',
                                  'sup', and 'inf'.

    Returns the palindromism of the data.
    """

    pal_len = params['pal_len'] if 'pal_len' in params else 2
    full = params['full'] if 'full' in params else False

    rec_table = _recurrence_common(data, **params) \
                if _rec_table is None else _rec_table

    data_len = len(rec_table)

    if full:
        rec_table += np.triu(rec_table.T, k=1)

    diagonals = np.zeros((data_len, data_len))
    diagonals[-1, :] = rec_table[-1, :].copy()
    diagonals[:, 0] = rec_table[:, 0].copy()

    for pos_x in range(data_len-2, -1, -1):
        for pos_y in range(1, pos_x+((data_len-pos_x)*full)):
            if rec_table[pos_x, pos_y]:
                diagonals[pos_x, pos_y] = diagonals[pos_x+1, pos_y-1] \
                                        + rec_table[pos_x, pos_y]
                diagonals[pos_x+1, pos_y-1] = 0

    pal_sum = np.sum(diagonals[np.where(diagonals >= pal_len)])
    rec_sum = np.sum(rec_table)

    pal = pal_sum / rec_sum

    return pal.tolist()


def ratio_palindromism(data, _rec_table=None, **params):
    """
    Calculate ratio of palindromism and recurrence of the data.

    Palindromism of the data can be defined in terms of the recurrence plot
    as length of lines at antidiagonals longer than the given length divided by
    the length of all lines at antidiagonals (or all recurrence points). This
    definition counts AA, ACA, and ACCA as palindroms of lengths 2, 3 and 4,
    respectively. The function calculates this version if 'full' is True.
    Alternatively, palindromism can be defined as length of lines at anti-
    diagonals longer than the given length, s.t. do not overlap nor cross the
    main diagonal, divided by the number of all recurrence points except the
    main diagonal. This definition counts AA, ACA, and ACCA as palindroms of
    lengths 1, 1 and 2, respectively. In other words, only the reoccuring part
    of the pattern is counted. The function calculates this version if 'full'
    is False (default).

    Args:
        data (list): values.
        params (kwargs):
            pal_len (int): minimum palindromic length (see above, default: 2)
            full (bool): whether to include the main diagonal or not
                (default: False)

            Settings for the embedding
            dim (float):Embedding dimension (default 1)
            tau (float):Embedding delay (default 0)

            Settings for the recurrence plot.
            radius (float): Fixed threshold (default 0)
            metric (string): Distance metric in the phase space.
                Possible choices: 'identity' (default), 'euclid', 'taxi',
                                  'sup', and 'inf'.

    Returns the ratio of determinism and recurrence of the data.
    """

    rec_table = _recurrence_common(data, **params) \
                if _rec_table is None else _rec_table

    ratio_pal = np.array(palindromism(data, _rec_table=rec_table, **params)) / \
           np.array(recurrence(data, _rec_table=rec_table, **params))

    return ratio_pal.tolist()


class RQAFeatureSet(object):
    """
    Manage a set of RQA-based sequence features in efficient way.

    Attributes:
        name (str): name for the feature set.
        features (list of str): selected features (default: None means all).
        window (int): length of the window over which the feature is calculated
                      (non-negative, default=0 means the whole sequence).

        params (kwargs):
            full (bool): whether to include the main diagonal or not
                (default: False)

            det_len (int): minimum diagonal length (default: 2)
            pal_len (int): minimum palindromic length (see above, default: 2)

            Settings for the embedding
            dim (float): Embedding dimension (default 1)
            tau (float): Embedding delay (default 0)

            Settings for the recurrence plot
            radius (float): Fixed threshold (default 0)
            metric (string): Distance metric in the phase space.
                Possible choices: 'identity' (default), 'euclid' and 'taxi'.

        Raises ValueError if improper window length is used.
    """

    def __init__(self, name, features=None, window=0, **params):

        if window < 0:
            raise ValueError("Window length for a feature has to be either "
                             "positive or zero for full sequence length!")

        self.name = name
        self.window = window

        self.params = {}
        self.params['dim'] = params['dim'] if 'dim' in params else 1
        self.params['tau'] = params['tau'] if 'tau' in params else 1
        self.params['radius'] = params['radius'] if 'radius' in params else 0
        self.params['metric'] = params['metric'] if 'metric' in params \
                                                 else 'identity'
        self.params['full'] = params['full'] if 'full' in params else False

        self.features = features if features is not None \
                                 else ['recurrence',
                                       'determinism', 'ratio_determinism',
                                       'palindromism', 'ratio_palindromism']

        if 'determinism' in self.features or \
           'ratio_determinism' in self.features:
            self.params['det_len'] = params['det_len'] if 'det_len' in params \
                                                       else 2
        if 'palindromism' in self.features or \
           'ratio_palindromism' in self.features:
            self.params['pal_len'] = params['pal_len'] if 'pal_len' in params \
                                                       else 2
        self._post_feat = None
        self._post_feat_name = ""

    def __call__(self, seq_set, window=-1):
        """
        Calculate RQA parameters for all sequences in a given sequence set using
        a given sliding window.

        Args:
            seq_set (SequenceSet): a SequenceSet instance
            window (int): length of the sliding window. Defaults to self.window
                          (window=-1). To calculate a single value for
                          the whole sequence use window=0.

        Returns a set of sequences of quantitative values (RQA parameters).
        """

        if window == -1:
            window = self.window

        result_set = SequenceSet(seq_set.name + ":" + self.name,
                                 unique=seq_set.unique)

        for seq in seq_set:

            rqa_output = {}

            rec_table = _recurrence_common(seq.data, **self.params)

            if window == 0:
                my_window = len(seq.data)

            if 'recurrence' in self.features or \
               'ratio_determinism' in self.features or \
               'ratio_palindromism' in self.features:
                rqa_output['recurrence'] = []
                for win_pos in range(len(seq.data) - my_window + 1):
                    rec_table_part = rec_table[win_pos:win_pos + my_window,
                                               win_pos:win_pos + my_window]
                    rqa_output['recurrence'].append(recurrence([],
                         _rec_table=rec_table_part, **self.params))

            if 'determinism' in self.features or \
               'ratio_determinism' in self.features:
                rqa_output['determinism'] = []
                for win_pos in range(len(seq.data) - my_window + 1):
                    rec_table_part = rec_table[win_pos:win_pos + my_window,
                                               win_pos:win_pos + my_window]
                    rqa_output['determinism'].append(determinism([],
                        _rec_table=rec_table_part, **self.params))

            if 'palindromism' in self.features or \
               'ratio_palindromism' in self.features:
                rqa_output['palindromism'] = []
                for win_pos in range(len(seq.data) - my_window + 1):
                    rec_table_part = rec_table[win_pos:win_pos + my_window,
                                               win_pos:win_pos + my_window]
                    rqa_output['palindromism'].append(palindromism([],
                        _rec_table=rec_table_part, **self.params))

            if 'ratio_determinism' in self.features:
                rqa_output['ratio_determinism'] = \
                    (np.array(rqa_output['determinism']) / \
                     np.array(rqa_output['recurrence'])).tolist()

            if 'ratio_palindromism' in self.features:
                rqa_output['ratio_palindromism'] = \
                    (np.array(rqa_output['palindromism']) / \
                     np.array(rqa_output['recurrence'])).tolist()

            for feat in self.features:

                my_seq = Sequence(seq.identifier, self.name, rqa_output[feat])

                if self._post_feat is not None:
                    my_seq = self._post_feat(my_seq)

                my_seq.feature = feat
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
        return len(self.features)

    def __repr__(self):
        return 'Quantiprot RQA Feature Set object {:s}'.format(self.name)

    def __str__(self):
        output = "Quantiprot RQA Feature Set object\n" + \
                 "  'name': '{:s}'\n".format(self.name) + \
                 "  features: {:}\n".format(self.features) + \
                 "  scanning 'window': {:d}".format(self.window)
        for param in self.params:
            output += "\n  '{:s}': {:}".format(param, self.params[param])

        return output
