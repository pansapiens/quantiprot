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
This module provides classes and methods for assigning quantitative features
to sequences.

Classes:
    Feature: define sequence to feature transformation.
    FeatureSet: manage a set of features.
"""

import copy

from quantiprot.utils.sequence import Sequence, SequenceSet


class Feature(object):
    """
    Define and perform sequence to feature transformation.

    Attributes:
        function (function): function that calculates the feature.
        name (str): name for the feature. If None (default),
                        the 'function' name is used instead as 'name'.
        window (int): length of the window over which the feature is calculated.
                      Defaults to the whole sequence (window=0).
        _post_feat: 'Feature' object to be called as a post-processor.
        params (**kwargs): arbitrary parameters to be passed to the function.
    """

    def __init__(self, function, name=None, window=0, **params):
        """
        Create a quantitative feature.

        Args:
            function (function): function that calculates the feature.
            name (str): name for the feature. If None (default),
                        the 'function' name is used instead as 'name'.
            window (int): length of the window over which the feature is
                          calculated. Defaults to the whole sequence (window=0).
            params (**kwargs): arbitrary params to be passed to the function.

        Raises ValueError if the window length is negative.
        """

        if window < 0:
            raise ValueError("Window length for a feature has to be either "
                             "positive or zero for full sequence length!")

        self.function = function
        self.name = name if name is not None else self.function.__name__
        self.window = window
        self.params = params
        self._post_feat = None

    def then(self, function, name=None, window=0, **params):
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

        self.name += ">"+post_feat.name

        return self

    def __call__(self, sequence, window=-1):
        """
        Calculate values of the feature for a sequence using a given
        sliding window.

        Args:
            sequence (Sequence): a Sequence instance
            window (int): length of the sliding window. Defaults to self.window
                          (window=-1). To calculate a single value for
                          the whole sequence use window=0.

        Returns a sequence of quantitative values.
        """

        if window == -1:
            window = self.window

        if window == 0 or window == len(sequence.data):
            outdata = self.function(sequence.data, **self.params)
            if not hasattr(outdata, '__getitem__'):
                outdata = [outdata]
        else:
            outdata = []
            for i in range(len(sequence.data) - window + 1):
                outdata.append(self.function(sequence.data[i:i + window],
                                             **self.params))

        my_seq = Sequence(sequence.identifier, self.name, outdata)

        if self._post_feat is not None:
            my_seq = self._post_feat(my_seq)

        my_seq.feature = self.name

        return my_seq

    def __repr__(self):
        return "Quantiprot Feature object '{:s}'".format(self.name)

    def __str__(self):
        output = "Quantiprot Feature object\n" + \
                 "  'name': {:s}\n".format(self.name) + \
                 "  'function': {:}\n".format(self.function.__name__) + \
                 "  scanning 'window': {:d}".format(self.window)
        if len(self.params) > 0:
            output += "\n  Keyword arguments:"
            for key in self.params:
                output += "\n    '{:s}': {:}".format(key, self.params[key])

        return output


class FeatureSet(object):
    """
   Manage a set of features.

    Attributes:
        name (str): name for the feature.
        feature_list (list): list of Features.
    """

    def __init__(self, name):
        self.name = name
        self.feature_list = []

    def then_all(self, function, name=None, window=0, **params):
        """
        Define a common post-processor for all features in the set.

        The method can either take a Feature or a function. In the former case
        the post-processor is a deep copy of the input Feature.

        The method modifies the features in the set.

        Args:
            function (function): Feature or function to serve a post-processor.
                Only when 'function' is not a Feature the following arguments
                are taken into account:
            name (str): name for the feature. If None (default),
                        the 'function' name is used instead as 'name'.
            window (int): length of the window over which the feature is
                          calculated. Defaults to the whole sequence (window=0).
            params (**kwargs): arbitrary params to be passed to the function.

        Returns the self to allow feature chaining.
        """

        for feature in self.feature_list:
            feature = feature.then(function, name=name, window=window, **params)

        return self

    def add(self, function, name=None, window=0, **params):
        """
        Add a feature to the feature set.

        The method can either take a Feature or a function. In the former case
        the post-processor is a deep copy of the input Feature.

        Args:
            function (function): Feature or function to serve a post-processor.
            Only when 'function' is not a Feature then
                the following arguments are taken into account:
            name (str): name for the feature. If None (default),
                        the 'function' name is used instead as 'name'.
            window (int): length of the window over which the feature is
                          calculated. Defaults to the whole sequence (window=0).
            params (**kwargs): arbitrary params to be passed to the function.
        """

        if isinstance(function, Feature):
            my_feat = copy.deepcopy(function)
            if name is not None:
                my_feat.name = name
        else:
            my_name = name if name is not None else function.__name__
            my_feat = Feature(function, name=my_name, window=window, **params)

        self.feature_list.append(my_feat)

    def merge_with(self, set2):
        """
        Merge the feature set with another one.
        The method does not verify if there are duplicated entries.

        The method modifies the self object.

        Args:
            set2 (FeatureSet): the other feature set.

        Returns the feature set itself.
        """

        for feature in set2.feature_list:
            self.feature_list.append(feature)

        return self

    def __call__(self, seq_set, window=-1):
        """
        Calculate values of the features in the feature set for all
        sequences in a given sequence set using a given sliding window.

        Args:
            seq_set (SequenceSet): a SequenceSet instance
            window (int): length of the sliding window. Defaults to self.window
                          (window=-1). To calculate a single value for
                          the whole sequence use window=0.

        Returns a set of sequences of quantitative values.
        """

        result_set = SequenceSet(seq_set.name + ":" + self.name,
                                 unique=seq_set.unique)

        for seq in seq_set:
            for feat in self.feature_list:
                if window == -1:
                    result_set.add(feat(seq))
                else:
                    result_set.add(feat(seq, window))

        return result_set

    def __setitem__(self, key, feature):
        self.feature_list[key] = feature

    def __getitem__(self, key):
        return self.feature_list[key]

    def __len__(self):
        return len(self.feature_list)

    def __repr__(self):
        return "Quantiprot FeatureSet object '{:s}'".format(self.name)

    def __str__(self):
        return "Quantiprot FeatureSet object\n" + \
               "  'name': {:s}\n".format(self.name) + \
               "  number of features: {:d}".format(self.__len__())
