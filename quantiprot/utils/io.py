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
This module provides functions that load a fasta file
and return a SequenceSet instance.

Functions:
    load_fasta_file: Load a fasta file given the filename.
    load_fasta_handle: Load a fasta file given the file handle.
"""

from quantiprot.utils.sequence import Sequence, SequenceSet

def load_fasta_file(filename, name=None, feature='fasta', unique=True):
    """
    Load a fasta file given the filename.

    Args:
        filename (str): name of the fasta file.
        name (str): name for the sequence set
            (default: None means use 'filename' as 'name')
        feature (str): name for the feature attribute of the Sequence object
            (default: 'fasta')
        unique (bool): whether sequence identifiers are to be unique
            (for the same feature, default: True)

    Returns:
        fasta_set: a SequenceSet instance with the name set to "filename".
    """

    set_name = name if name is not None else filename

    handle = open(filename, "rU")
    fasta_set = load_fasta_handle(handle, set_name, feature=feature,
                                  unique=unique)
    handle.close()
    return fasta_set


def load_fasta_handle(handle, name, feature='fasta', unique=True, ):
    """
    Load a fasta file given the file handle.

    Args:
        handle (file): handle to the fasta file.
        name (str): name for the sequence set.
        feature (str): name for the feature attribute of the Sequence object
            (default: 'fasta')
        unique (bool): whether sequence identifiers are to be unique
            (for the same feature, default: True)

    Returns:
        fasta_set: a SequenceSet instance with the name set to "name".
    """

    lines = handle.read().splitlines()
    fasta_set = SequenceSet(name, unique=unique)

    identifier = None
    data = None
    for line in lines:
        if len(line) > 0 and line[0] == '>':
            if identifier is not None:
                fasta_set.add(Sequence(identifier, feature, data))
            identifier = line[1:]
            data = []
        elif identifier is not None:
            data = data + list(line)
        else:
            pass

    if identifier is not None:
        fasta_set.add(Sequence(identifier, feature, data))

    return fasta_set

