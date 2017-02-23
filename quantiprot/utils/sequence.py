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
This module provides classes, methods and functions to store and access
sequences and sequence sets.

Classes:
    Sequence: Store and access a single sequence.
    SequenceSet: Store and access a sequence set.

Functions:
    _uniq: return a list of unique elements of a list of elements.
    merge: Merge two sequence sets.
    subset: Return a subset of the sequence set.
    compact: Compact a sequence set.
    columns: Return a list of lists of elements at selected positions from all
        sequences of selected feature(s).
"""

from warnings import warn

import itertools


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


class Sequence(object):
    """
    Store and access a single sequence.

    Attributes:
        identifier (str): identifier of the sequence.
        feature (str): type of the sequence, e.g."aminoacid".
        data (list): e.g. sequence of letters coding aminoacids
                     or numbers coding their properties.

    """

    def __init__(self, identifier, feature, data):
        self.identifier = identifier
        self.feature = feature
        self.data = data

    def __len__(self):
        return len(self.data)

    def __setitem__(self, key, val):
        self.data[key] = val

    def __getitem__(self, key):
        return self.data[key]

    def fasta(self):
        """
        Return the fasta string of the sequence.
        """
        return ">" + self.identifier + "|" + self.feature + "\n" + \
               ''.join(self.data) + "\n"

    def __repr__(self):
        return "Quantiprot Sequence object ({:s}, ".format(self.identifier) + \
               "feature={:s})".format(self.feature)

    def __str__(self):
        return "Quantiprot Sequence object\n" + \
               "  'identifier': {:s}\n".format(self.identifier) + \
               "  'feature': {:s}\n".format(self.feature) + \
               "  'data': {:}".format(self.data)


class SequenceSet(object):
    """
    Store and access a sequence set.

    Attributes:
        name (str): name of the sequence set.
        unique (bool): whether pairs (identifier, feature) are be unique
                       (default: True).
        sequence_list (list): Sequence instances.
        _sequence_unique_dict (dict): dictionary of unique (identifier, feature)
                                      pairs.
    """

    def __init__(self, name, unique=True):
        self.name = name
        self.unique = unique
        self.sequence_list = []
        self._sequence_unique_dict = {}

    def __setattr__(self, attr, value):
        """
        Set class attributes except 'unique' which is immutable.

        Args:
            attr (str): a class attribute.
            value: an atribute value.

        Raises ValueError if attempt to modify the uniqueness attribute.
        """
        if attr == "unique" and hasattr(self, 'unique'):
            raise ValueError("The uniqueness property cannot be changed.")
        self.__dict__[attr] = value

    def add(self, sequence, verbose=False):
        """
        Add a sequence to the sequence set.

        The method modifies the object.

        This method does not allow for duplication of sequences if 'self.unique'
        is set to True. In practice it verifies if a pair (sequence.identifier,
        sequence.feature) is already in the set.

        Args:
            sequence (Sequence): a sequence.
            verbose (bool): if True and the set is unique, attempts to add
                            duplicate entries are reported.
        """
        if self.unique and (sequence.identifier, sequence.feature) \
            in self._sequence_unique_dict:
            if verbose:
                warn("Quantiprot SequenceSet warning: " + \
                     "sequence ({:s}, ".format(sequence.identifier) + \
                     "feature={:s}) ".format(sequence.feature) + \
                     "already in '{:s}'".format(self.name))
            else:
                pass
        else:
            self.sequence_list.append(sequence)
            self._sequence_unique_dict[(sequence.identifier,
                                        sequence.feature)] = True


    def __getitem__(self, pos):
        return self.sequence_list[pos]

    def fasta(self):
        """
        Return the fasta string for the sequence set.

        Returns the fasta string for the sequence set.
        """
        res = ""
        for seq in self.sequence_list:
            res += seq.fasta()
        return res

    def __len__(self):
        return len(self.sequence_list)

    def ids(self, unique=False):
        """
        Return sequence identifiers in the sequence set.

        Args:
            unique: whether to report the identifier for each sequence
                    (default: False; this may result in repeated identifiers
                    e.g. if several features for the same sequence), or
                    only first time it occurs in the set (True).

        Returns the list of sequence identifiers.
        """

        if unique:
            return _uniq([seq.identifier for seq in self.sequence_list])
        else:
            return [seq.identifier for seq in self.sequence_list]

    def keys(self):
        """
        Return pairs (identifier, feature) for all sequences in the set.

        Returns the list of pairs (identifier, feature) for all sequences.
        """

        return [(seq.identifier, seq.feature) for seq in self.sequence_list]

    def merge_with(self, set2, verbose=False):
        """
        Merge the sequence set with another one.

        The method modifies the object.

        Args:
            set2 (SequenceSet): the other sequence set.
            verbose (bool): if True and the set is unique, attempts to add
                            duplicate entries are reported.

        Returns the sequence set itself.
        """
        for seq in set2.sequence_list:
            self.add(seq, verbose=verbose)
        return self

    def columns(self, positions=None, feature=None, default=None,
                transpose=False):
        """
        Return a list of lists of elements at selected positions from all
        sequences of selected feature(s).

        Args:
            positions (list or int): list of positions in a sequence
               (default: None means all positions up to the length of the
               longest sequence are used)
            feature (str): selected feature identifier (default=None
                is valid only if all sequence share the same feature).
            default: the element to be used if sequence ends
                before the positions (the default default is None)
            transpose: whether to return the list of columns from sequences
                (default: False), or the list of sequences with selected
                columns only (True).

        Returns a list of list of elements at selected positions.

        Raises ValueError if sequences in the set represent multiple features
            and none of them is selected.

        The size of the returned list of list is |positions| x |features|.
        """

        if feature is None:
            num_feat = len(set([seq.feature for seq in self.sequence_list]))
            if num_feat == 1:
                feature = self.sequence_list[0].feature
            elif num_feat > 1:
                raise ValueError("Please select one feature!")

        if positions is None:
            max_len = max([len(seq) for seq in self.sequence_list])
            positions = range(max_len)

        res_columns = []

        if not transpose:

            for pos in positions:
                column = []
                for seq in self.sequence_list:
                    if seq.feature == feature:
                        if pos < len(seq):
                            column.append(seq[pos])
                        else:
                            column.append(default)
                res_columns.append(column)

        else:

            for seq in self.sequence_list:
                if seq.feature == feature:
                    row = []
                    for pos in positions:
                        if pos < len(seq):
                            row.append(seq[pos])
                        else:
                            row.append(default)
                    res_columns.append(row)

        return res_columns

    def __repr__(self):
        return "Quantiprot SequenceSet object: {:s}".format(self.name)

    def __str__(self):
        return "Quantiprot SequenceSet object\n" + \
               "  'name': {:s}\n".format(self.name) + \
               "  'unique': {:}\n".format(self.unique) + \
               "  number of sequences: {:d}".format(self.__len__())

def merge(set1, set2, name=None, unique=True, verbose=False):
    """
    Merge two sequence sets.

    Args:
        set1 (SequenceSet): first sequence set.
        set2 (SequenceSet): second sequence set.
        name (str): name for the merged sequence set (automatic if None)
        unique (bool): whether duplicate entries (identifier, feature)
            should be omitted (default: True).
        verbose (bool): whether attempts to duplicate entry should be reported.

    Returns:
        new_set (SequenceSet): the merged sequence set.
    """

    name = set1.name + '+' + set2.name if name == None else name
    new_set = SequenceSet(name, unique=unique)
    new_set.merge_with(set1, verbose=verbose)
    new_set.merge_with(set2, verbose=verbose)
    return new_set


def subset(seq_set, features, name=None):
    """
    Return a subset of the sequence set given the (list of) features.

    Args:
        seq_set (SequenceSet): sequence set
        features (str or list): type(s) of the sequence, e.g. 'fasta'.
        name (str): name for the subset of the sequence set (automatic if None)

    Returns:
        sub_set (SequenceSet): the subset of the sequence set.
    """

    name = "subset_of:"+seq_set.name if name == None else name
    sub_set = SequenceSet(name, unique=seq_set.unique)

    if features is None:
        features = _uniq([seq.feature for seq in seq_set])

    if not hasattr(features, '__getitem__') or isinstance(features, str):
        features = [features]

    for feat in features:
        for seq in seq_set.sequence_list:
            if seq.feature == feat:
                sub_set.add(seq)

    return sub_set


def compact(seq_set, features=None, name=None):
    """
    Compact the sequence set given the (list of) features.

    The method joins values of all features for given sequence identifier into
    one sequence of values. It is intended to be used to compact sequence sets
    where each sequence is represented by single values of multiple features.

    Args:
        seq_set (SequenceSet): sequence set
        feature (str or list): type(s) of the sequence, e.g."aminoacid".
        name (str): name for the subset of the sequence set (automatic if None)

    Returns the compacted SequenceSet.
    """

    if features is None:
        features = _uniq([seq.feature for seq in seq_set])

    if not hasattr(features, '__getitem__') or isinstance(features, str):
        features = [features]

    unique_seq_ids = _uniq([seq.identifier for seq in seq_set])

    seq_data = {}
    for seq_id in unique_seq_ids:
        seq_data[seq_id] = list(itertools.repeat(None, len(features)))

    for feat_pos in range(len(features)):
        for seq in seq_set:
            if seq.feature == features[feat_pos]:
                seq_data[seq.identifier][feat_pos] = seq.data \
                         if len(seq) > 1 else seq.data[0]

    name = "compact_of:"+seq_set.name if name == None else name
    compact_set = SequenceSet(name, unique=seq_set.unique)
    for seq_id in unique_seq_ids:
        compact_set.add(Sequence(seq_id, str(features), seq_data[seq_id]))

    return compact_set


def columns(seq_set, positions=None, feature=None, default=None,
            transpose=False):
    """
    Return a list of lists of elements at selected positions from all
    sequences of selected feature(s).

    Args:
        seq_set (SequenceSet): a sequence set.
        positions (list or int): list of positions in a sequence
           (default: None means all positions up to the length of the
           longest sequence are used)
        features (str or list): selected feature identifiers (default=None).
        default: the element to be used if sequence ends
            before the positions (the default default is None)
        transpose: whether to return the list of columns from sequences
            (default: False), or the list of sequences with selected
            columns only (True).

    Returns a list of list of elements at selected positions.

    The size of the returned list of list is |positions| x |features|.
    """

    return seq_set.columns(positions=positions, feature=feature,
                           default=default, transpose=transpose)
