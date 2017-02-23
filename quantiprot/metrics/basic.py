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
This module calculates basic metrics for a sequence.

Functions:
    identity: return the data itself.
    absolute: calculate the absolute values of the data.
    sum_absolute: calculate the sum of absolute values of the data.
    average: calculate the arithmetic average of the data.
    average_absolute: calculate the average of absolute values of the data.
    uniq_count: count number of unique elements in the data.
    uniq_average: calculate number of unique elements per length in the data.
    atom_count: count occurrencies of a given atomic element in the data.
    atom_freq: calculate frequency of occurrencies of a given atomic element
               in the data.
"""

def identity(data):
    """
    Return the data itself.
    Args:
        data (list): values.

    Returns the data itself.
    """
    return list(data)


def absolute(data):
    """
    Calculate the abolute values of the data.

    Args:
        data (list): values.
    Returns the absolute values of the data.
    """
    return [abs(d) for d in data]


def sum_absolute(data):
    """
    Calculate the sum of abolute values of the data.

    Args:
        data (list): values.
    Returns the sum of absolute values of the data.
    """
    return sum(absolute(data))


def average(data):
    """
    Calculate the average of values of the data.
    Args:
        data (list): values.
    Returns the average of values of the data.
    """
    return 1.0*sum(data)/len(data)


def average_absolute(data):
    """
    Calculate the average of absolute values of the data.
    Args:
        data (list): values.
    Returns the average of absolute values of the data.
    """
    return average(absolute(data))


def uniq_count(data):
    """
    Count number of unique elements in the data.
    Args:
        data (list): values.
    Returns the number of unique elements in the data.
    """
    uniq_atom_list = list(set(data))
    return len(uniq_atom_list)


def uniq_average(data):
    """
    Calculate number of unique elements per length in the data.
    Args:
        data (list): values.
    Returns the number of unique elements per length in the data.
    """
    return 1.0*uniq_count(data)/len(data)


def atom_count(data, **params):
    """
    Calculate number of occurrencies of a given atomic element in the data.
    Args:
        data (list): values.
        params (kwargs):
            atom: element for which occurencies are counted.

    Returns the number of occurencies of the atom in the data.
    """
    atom = params['atom']
    counter = sum([elem == atom for elem in data])
    return counter


def atom_freq(data, **params):
    """
    Calculate frequency of occurrencies of a given atomic element in the data.
    Args:
        data (list): values.
        params (kwargs):
            atom: element for which frequency is calculated.

    Returns the frequency of the atom in the data.
    """
    return 1.0*atom_count(data, **params)/len(data)


