# Copyright (c) 2016-2017 Bogumil M. Konopka & Witold Dyrka
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
This module povides functions for fitting the Zipf's law to ngrams distribution.

The module relies on the powerlaw package:
J.Alstott, E.Bullmore, D.Plenz. powerlaw: a Python package for analysis of
heavy-tailed distributions. PLoS ONE 9(1):e85777. 2014

Functions:
    ngram_count: count n-grams in a sequence set.
    zipf_law_fit: fits the Zipf's law to n-grams distribution.
"""
import numpy as np

from powerlaw import Fit as powerlaw_fit

from quantiprot.utils.sequence import compact, subset
from quantiprot.metrics.ngram import NgramFeatureSet

def ngram_count(seq_set, n=2, alphabet=None, return_all=False, **params):
    """
    Count n-grams in a sequence set.

    Args:
        seq_set (SequenceSet): input SequenceSet.
        n (int): ngram length (default 2).
        alphabet: alphabet (list): alphabet to construct ngrams from
            (default: None means alphabet is the set of sequence elements).
        return_all (bool): whether ngrams with zero occurences are to be
            returned (True) or not (default: False).
        params (**kwargs): parameters to be passed to NgramFeatureSet object.

    Returns a dictionary with ngram counts.
    """

    if alphabet is None:
        alphabet = list(set([element for seq in seq_set
                             for element in seq.data]))

    ngine = NgramFeatureSet('ngram_'+str(n), n, alphabet=alphabet,
                            mode='count', window=0, simple_names=True, **params)

    seq_ngrams = ngine(seq_set)

    ngram_counts = {}

    for seq in seq_ngrams:
        if return_all or seq.data[0]>0:
            ngram_counts[seq.feature] = seq.data[0] + ngram_counts.get(seq.feature,0)

    return ngram_counts


def zipf_law_fit(seq_set, n=2, alphabet=None,
                  discrete=True, verbose=False, **params):
    """
    Fit Zipf's law to n-gram distribution.

    The rank-frequency plot obeys Zipf's law if frequencies obey the power law.
    Thus, the function fits the power law distribution to the distribution of
    n-grams. It also performs a comparison of goodness-of-fit with respect to
    the exponential distribution. It the returns goodness-of-fit ratio and
    p-value for its significance. The null hypothesis for the test is:
    h0 - fitting with the exponential distribution is no worse than with
    the power law distribution.

    Args:
        seq_set (SequenceSet): input SequenceSet.
        n (int): ngram length (default 2).
        alphabet:  alphabet (list): alphabet to construct ngrams from
            (default: None means alphabet is the set of sequence elements).
        discrete (bool): whether the fitted values are discrete (default: True.
        verbose (bool): whether results have to be printed to std output
            (default False).
        params (**kwargs): parameters to be passed to NgramFeatureSet object.

    Returns a dictionary including:
        Parameters of the power law distribution: 'alpha', 'xmin'
            and the standard error of the power-law fit 'sigma'.
        'slope' (float): the slope of the Zipf's law:
            'slope' = 1 / ('alpha' - 1).
        'ratio' (float): the goodness-of-fit log-likelihood ratio
            power law vs exponential distribution
        'p_value' (float): the goodness-of-fit p-value.
    """

    ngram_counts = ngram_count(seq_set, n=n, alphabet=alphabet,
                               return_all=False, **params).values()

    numpy_err_params = np.geterr()
    np.seterr(divide='ignore', invalid='ignore') # otherwise a RuntimeWarning
    ngram_fit = powerlaw_fit(ngram_counts, discrete=discrete)
    ratio, p_value = ngram_fit.distribution_compare('power_law', 'exponential',
                                                    normalized_ratio=True)
    np.seterr(divide=numpy_err_params['divide'],
              invalid=numpy_err_params['invalid'])

    if verbose:
        print("Power-law fit:")
        print("  'alpha': {:f}".format(ngram_fit.alpha))
        print("  'xmin': {:f}".format(ngram_fit.xmin))
        print("  'sigma': {:f}".format(ngram_fit.sigma))
        print("Comparison to exponential:")
        print("  log-likelihood 'ratio': {:f} ".format(ratio))
        print("  'p-value': {:f}".format(p_value))
        print("Zipf's law slope: {:f}".format(1.0/(ngram_fit.alpha-1.0)))

    return {'alpha': ngram_fit.alpha, 'xmin': ngram_fit.xmin,
            'sigma': ngram_fit.sigma, 'slope': 1.0/(ngram_fit.alpha-1.0),
            'ratio':ratio, 'p_value':p_value,
            'ngram_counts': ngram_counts}
