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
This module provides classes and methods for effective handling of mapping
from aminoacids (or anything else) to other features.

Classes:
    Mapping: store and handles the mapping

Functions:
    _kmeans1d: 1-d implementation of the k-means algorithm.
    simplify: reduce the mapping.
"""

import copy
import numpy as np

class Mapping(object):
    """
    Define and perform a mapping for a sequence of elements

    The Mapping object acts as a dict with operator '[]'
    or as a function operating on a list-like with operator '()'.

    Attributes:
        __name__ (str): name for the mapping.
        mapping (dict): the actual mapping.
        default: default value if an input value not in the mapping keys
                 (the default 'default' is None)
        misc: additional information (default: None)
    """

    def __init__(self, name, mapping, default=None, misc=None):
        self.__name__ = name
        self.mapping = mapping
        self.default = default
        self.misc = misc

    def __call__(self, data):
        """
        Perform mapping for a sequence of elements.

        Args:
            data (list): a sequence of elements to be mapped.

        Returns the mapped sequence.
        """
        if self.default is None:
            res = [self.mapping[d] for d in data]
        else:
            res = [self.mapping.get(d, self.default) for d in data]
        return res

    def __getitem__(self, key):
        return self.mapping[key]

    def get(self, key, default=None):
        """
        Get the mapping value for a given key.

        Args:
            key: a mapping key.
            default: default value if 'key' not in the mapping keys
                 (the default default is 'self.default')

        Returns the mapping value.
        """
        return self.mapping[key] if key in self.mapping \
            else (default if default is not None else self.default)

    def simplify(self, **params):
        """
        Reduce the mapping.

        The method modifies the self object.

        Args:
            params (**kwargs):
               'k' (int): desired number of classes (default: None).
               'thresholds' (list): list of thresholds (default: None).
               'centroids' (list): list of centroids (default: None).
               'labels' (list): list of labels (default: None).
               'method' (str): method for reduction of the mapping:
                    'linear', 'kmeans' (default: 'auto'=='kmeans')
                    Applicable if neither thresholds nor centroids given.
               'iters' (int): number of iterations if method is 'kmeans'
                    (default: 0)
               'label_from' (int): numerical label of the first class.
                    Applicable if 'mean_labels' is False.
               'mean_labels' (bool): whether class labels are to be numerical
                    (False, default) or be means of elements
                    in the class (True).
               'default': default value if an input value not in the mapping
                   keys (the default default is 'self.default' if 'mean_labels'
                   is True; otherwise None).

        Returns the reduced self.

        Raises:
            ValueError if:
                if more than one from 'thresholds', 'centroids' and 'k'
                    are provided at the same time.
                if neither 'thresholds' nor 'centroids' nor 'k'
                    are provided and also 'labels' are not provided.
                'k' < 1,
                'method' is neither 'linear' nor 'kmeans' nor 'auto',
                number of 'labels' is less than number of classes.

        Method simplifies the mapping according to provided 'thresholds'
        or 'centroids' or cluster the data to 'k' classes. If none of these
        is specified, the desired number of classes is deduced from the number
        of labels.

        Clustering using the 'linear' method simply divide the 1d space evenly.
        Clustering using the 'kmeans' method divide the 1d space using even
        percentiles of the original mapping values distribution. Then k-means
        iterations are performed if 'iters' > 0.

        If 'labels' are given, they are assigned to classes. Otherwise, if
        'mean_labels' are set to True, mean values in each class is assigned
        as the class label; otherwise, class number + 'label_from" is assigned
        as the class label. So created class labels become new mapping values
        of the modified self.

        Note that it is often required to adjust the 'default' value for the
        input values not in the mapping keys (except if new mapping labels are
        'mean_labels').
        """

        k = params['k'] if 'k' in params else None
        thresholds = params['thresholds'] if 'thresholds' in params else None
        centroids = params['centroids'] if 'centroids' in params else None
        labels = params['labels'] if 'labels' in params else None
        method = params['method'] if 'method' in params else 'auto'
        iters = params['iters'] if 'iters' in params else 0
        label_from = params['label_from'] if 'label_from' in params else 0
        mean_labels = params['mean_labels'] if 'mean_labels' in params \
                                            else False
        default = params['default'] if 'default' in params else None

        if np.sum([k is not None,
                   thresholds is not None,
                   centroids is not None]) > 1:
            raise ValueError("Specify either 'thresholds' or 'centroids'" \
                             " or number of classes 'k'!")

        mapping_values = np.transpose(np.matrix(self.mapping.values()))

        if thresholds is not None:
            classes = np.sum(np.matrix(mapping_values > np.sort(thresholds)),
                             axis=1)

        elif centroids is not None:
            classes = np.argmin(np.abs(mapping_values - centroids), axis=1)

        else:

            if k is None:
                if labels is None:
                    raise ValueError("Specify either 'thresholds'" \
                                     " or 'centroids' or number of classes " \
                                     " 'k'. Or at least 'labels'!")
                else:
                    k = len(labels)

            if k < 1:
                raise ValueError("Cannot reduce to less than one class!")

            if method == 'linear':
                thresholds = np.linspace(np.min(mapping_values),
                                         np.max(mapping_values), k+1)
                classes = np.sum(np.matrix(mapping_values > thresholds[1:-1]),
                                 axis=1)

            elif method in ['auto', 'kmeans']:
                classes = _kmeans1d(k, iters, mapping_values)

            else:
                raise ValueError("Valid methods are 'auto'," \
                                 " 'kmeans' and 'linear'!")

        if labels is not None:
            if mean_labels:
                raise ValueError("Only one method for assigning labels" \
                                 " is permitted at time!")
            if len(labels) < len(np.unique(np.array(classes))):
                raise ValueError("Insufficient number of labels!")

            mapping_labels = [labels[c] \
                              for c in np.array(classes).flatten().tolist()]
            mapping_default = None

        elif mean_labels:
            centroids = [np.mean(np.array(mapping_values[classes == c])) \
                         for c in np.unique(np.array(classes))]
            mapping_labels = [centroids[c] \
                        for c in np.array(classes).flatten().tolist()]
            mapping_default = self.default

        else:
            mapping_labels = np.array(classes+label_from).flatten().tolist()
            mapping_default = None

        self.mapping = dict(zip(self.mapping.keys(), mapping_labels))
        self.default = default if default is not None else mapping_default

        return self

    def __repr__(self):
        return "Quantiprot Mapping object '{:s}'".format(self.__name__)

    def __str__(self):
        return "Quantiprot Mapping object\n" + \
               "  '__name__': {:s}\n".format(self.__name__) + \
               "  'mapping': {:}\n".format(self.mapping) + \
               "  'default': {:}\n".format(self.default) + \
               "  'misc': {:}".format(self.misc)


def _kmeans1d(k, iters, data):
    """
    Perform k-means clustering of 1d data.

    Args:
        k (int): desired number of classes.
        iters (int): number of iterations.
        data (list): data to be clustered.

    Returns:
        classes (list): list of classes of the data.

    First, the function divide the 1d space using even percentiles of data
    distribution. The percentile interval is increased until there are 'k'
    populated clusters or number of unique data values reached.

    Then, if 'iters'>0, the standard k-means procedure is performed 'iters'
    times or until convergence.
    """

    for opt_k in range(1, len(np.unique(np.array(data)))+1):
        percs = np.linspace(0, 100, opt_k+1)+100.0/(2.0*opt_k)
        centroids = np.unique(np.array(np.percentile(data, percs[:-1])))
        classes = np.argmin(np.abs(centroids-data), axis=1)
        if len(np.unique(np.array(classes))) == k:
            break

    centroids_old = centroids

    for _ in range(1, iters):
        centroids = np.array([np.mean(np.array(data[classes == c])) \
                              for c in np.unique(np.array(classes))])

        if np.sum(np.abs(centroids-centroids_old)) == 0:
            break

        centroids_old = centroids
        classes = np.argmin(np.abs(centroids-data), axis=1)

    return classes


def simplify(mapping, name, **params):
    """
    Reduce the mapping.

    The method modifies the self object.

    Args:
        mapping (Mapping): the mapping to be reduced
        name (str): name of the reduced mapping
        params (**kwargs):
           'k' (int): desired number of classes (default: None).
           'thresholds' (list): list of thresholds (default: None).
           'centroids' (list): list of centroids (default: None).
           'labels' (list): list of labels (default: None).
           'method' (str): method for reduction of the mapping:
                'linear', 'kmeans' (default: 'auto'=='kmeans')
                Applicable if neither thresholds nor centroids given.
           'iters' (int): number of iterations if method is 'kmeans'
                (default: 0)
           'label_from' (int): numerical label of the first class.
                Applicable if 'mean_labels' is False.
           'mean_labels' (bool): whether class labels are to be numerical
                (False, default) or be means of elements
                in the class (True).
           'default': default value if an input value not in the mapping
               keys (the default default is 'self.default' if 'mean_labels'
               is True; otherwise None).

    Returns the reduced mapping.

    Raises:
        ValueError if:
            if more than one from 'thresholds', 'centroids' and 'k'
                are provided at the same time.
            if neither 'thresholds' nor 'centroids' nor 'k'
                are provided and also 'labels' are not provided.
            'k' < 1,
            'method' is neither 'linear' nor 'kmeans' nor 'auto',
            number of 'labels' is less than number of classes.

    Method simplifies the mapping according to provided 'thresholds'
    or 'centroids' or cluster the data to 'k' classes. If none of these
    is specified, the desired number of classes is deduced from the number
    of labels.

    Clustering using the 'linear' method simply divide the 1d space evenly.
    Clustering using the 'kmeans' method divide the 1d space using even
    percentiles of the original mapping values distribution. Then k-means
    iterations are performed if 'iters' > 0.

    If 'labels' are given, they are assigned to classes. Otherwise, if
    'mean_labels' are set to True, mean values in each class is assigned
    as the class label; otherwise, class number + 'label_from" is assigned
    as the class label. So created class labels become the mapping values
    of the new reduced mapping.

    Note that it is often required to adjust the 'default' value for the
    input values not in the mapping keys (except if new mapping labels are
    'mean_labels').
    """

    new_mapping = copy.deepcopy(mapping)
    new_mapping.__name__ = name
    return new_mapping.simplify(**params)
