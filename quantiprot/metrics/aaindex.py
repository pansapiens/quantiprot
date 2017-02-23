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
This module provides functions to generate the Mapping objects based on
various physico-chemical properties of aminoacids defined in the AAindex
database. In addition, generators of several common mappings are included.

Functions:
    _parse_aaindex: generate a Mapping for an AAindex entry.
    get_aaindex_file: given an index, generate a Mapping from the AAindex
        flat file.
    get_aaindex_www: given an index, generate a Mapping from the AAindex
        website.
    get_aa2numeric: generate a Mapping from aminoacids to integers.
    get_aa2charge: generate a Mapping from aminoacids to formal charges.
    get_aa2volume: generate a Mapping from aminoacids to van der Waals volumes.
    get_aa2mj: generate a Mapping from aminoacids to the Miyazawa-Jernigan
        hydrophobicity scale.
"""

import os
import re

from requests import get as requests_get

from quantiprot.utils.mapping import Mapping

def _parse_aaindex(index_id, raw_text_lines, default=None):
    """
    Generate a Mapping, for a given index_id, from the raw text lines
    in the AAindex format.

    Args:
        index_id (str): AAindex identifier.
        raw_text_lines (list of str): raw text lines in the AAindex format.
        default: default value for elements beyond the standard 20 aminoacids.

    Returns:
        A Mapping object based on the index_id entry in the AAindex database
        with 'name' set to 'index_id' and 'misc' containing dictionary
        with additional information about the index including:
        'index_id', 'description', 'authors', 'title', 'journal', 'correlations'

    Raises ValueError if 'index_id' not found in the 'raw_text_lines'.
    """

    for i in range(0, len(raw_text_lines)):
        line = raw_text_lines[i]
        if len(line) > 0 and line[0] == 'H':
            header = line[2:].strip('\n')
            if header == index_id:
                i += 1
                line = raw_text_lines[i]
                while line[0] != '/' and line[1] != '/':
                    if line[0] == 'D':
                        description = line[2:].strip('\n')
                    elif line[0] == 'A':
                        authors = line[2:].strip('\n')
                    elif line[0] == 'T':
                        title = line[2:].strip('\n')
                    elif line[0] == 'J':
                        journal = line[2:].strip('\n')
                    elif line[0] == 'C':
                        while line[0] != 'I':
                            i += 1
                            line = raw_text_lines[i]
                            #correlations = None
                        if line[0] == 'I':
                            names = ['A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', \
                                     'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', \
                                     'T', 'W', 'Y', 'V']
                            i += 1
                            line = raw_text_lines[i:i+2]
                            num = [l.split() for l in line]
                            index = [float(j) for l in num for j in l]
                            mapping = dict(zip(names, index))
                    i += 1
                    line = raw_text_lines[i]

                return Mapping(index_id, mapping,
                               misc={"index_id":index_id,
                                     "description":description,
                                     "authors":authors,
                                     "title":title,
                                     "journal":journal,
                                     #"correlations":correlations
                                    },
                               default=default)

    raise ValueError("{:s} not found in the database".format(index_id))


def get_aaindex_file(index_id, default=None, _filename='data/aaindex1'):
    """
    Generate a Mapping from the AAindex flat file for a given index.

    Args:
        index_id (str): AAindex identifier.
        default: default value for elements beyond the standard 20 aminoacids.
        _filename (str): path to the AAindex flat file (default: data/aaindex1).
    Returns:
        A Mapping object based on the index_id entry in the AAindex database
        with 'name' set to 'index_id' and 'misc' containing dictionary
        with additional information about the index including:
        'index_id', 'description', 'authors', 'title', 'journal', 'correlations'

    Invoked function _parse_aaindex raises ValueError if 'index_id' not found
    in the AAindex flat file.
    """

    file_path = os.path.dirname(os.path.abspath(__file__))+"/"+_filename
    file_id = open(file_path, 'r')
    lines = file_id.readlines()
    file_id.close()

    return _parse_aaindex(index_id, lines, default=default)


def get_aaindex_www(index_id, default=None,
                    _url='http://www.genome.jp/dbget-bin/www_bget'):
    """
    Generate a Mapping from the AAindex website.

    Args:
        index_id (str): AAindex identifier.
        default: default value for elements beyond the standard 20 aminoacids.
        _url (str): url to the online retrieval system for AAindex
                    (default: http://www.genome.jp/dbget-bin/www_bget)
    Returns:
        A Mapping object based on the index_id entry in the AAindex database
        with 'name' set to 'index_id' and 'misc' containing dictionary
        with additional information about the index including:
        'index_id', 'description', 'authors', 'title', 'journal', 'correlations'

    Invoked function _parse_aaindex raises ValueError if 'index_id' not found
    in the AAindex flat file.
    """
    index_site = _url+'?aaindex:'+index_id
    req = requests_get(index_site, auth=('user', 'pass'))
    lines = re.sub(r'<[^>]*>', '', req.content).splitlines()

    return _parse_aaindex(index_id, lines, default=default)


def get_aa2numeric(default=None):
    """
    Generate a Mapping from aminoacids to integers.

    Args:
        default: default value for elements beyond the standard 20 aminoacids.

    Returns a Mapping object from aminoacids to integers.
    """

    return Mapping("numeric", {'A': 1, 'C': 2, 'D': 3, 'E': 4,
                               'F': 5, 'G': 6, 'H': 7, 'I': 8,
                               'K': 9, 'L': 10, 'M': 11, 'N': 12,
                               'P': 13, 'Q': 14, 'R': 15, 'S': 16,
                               'T': 17, 'V': 18, 'W': 19, 'Y': 20,
                              }, default=default)


def get_aa2charge(default=None):
    """
    Generate a Mapping from aminoacids to formal charges.

    Args:
        default: default value for elements beyond the standard 20 aminoacids.

    Returns a Mapping object from aminoacids to formal charges.
    """
    return Mapping("formal_charge", {'A': 0, 'C': 0, 'D': -1, 'E': -1,
                                     'F': 0, 'G': 0, 'H': 0, 'I': 0,
                                     'K': 1, 'L': 0, 'M': 0, 'N': 0,
                                     'P': 0, 'Q': 0, 'R': 1, 'S': 0,
                                     'T': 0, 'V': 0, 'W': 0, 'Y': 0,
                                    }, default=default)


def get_aa2volume(default=None):
    """
    Generate a Mapping from aminoacids to normalized van der Waals volumes.

    The mapping is equivalent to AAindex: FAUJ880103 (Fauchere et al., 1988)

    Args:
        default: default value for elements beyond the standard 20 aminoacids.

    Returns a Mapping object from aminoacids to normalized v. der Waals volumes.
    """
    return Mapping("volume",
                   {'A': 1.00, 'C': 2.43, 'D': 2.78, 'E': 3.78,
                    'F': 5.89, 'G': 0.00, 'H': 4.66, 'I': 4.00,
                    'K': 4.77, 'L': 4.00, 'M': 4.43, 'N': 2.95,
                    'P': 2.72, 'Q': 3.95, 'R': 6.13, 'S': 1.60,
                    'T': 2.60, 'V': 3.00, 'W': 8.08, 'Y': 6.47,
                   }, default=default)


def get_aa2hydropathy(default=None):
    """
    Generate a Mapping from aminoacids to Kyte-Doolittle hydrophathy scale.

    The mapping is equivalent to AAindex: KYTJ820101 (Kyte-Doolittle, 1982).

    Args:
        default: default value for elements beyond the standard 20 aminoacids.

    Returns a Mapping object from aminoacids to Kyte-Doolitle hydropathy
        values.
    """
    return Mapping("hydropathy",
                   {'A': 1.8, 'C': 2.5, 'D': -3.5, 'E': -3.5,
                    'F': 2.8, 'G': -0.4, 'H': -3.2, 'I': 4.5,
                    'K': -3.9, 'L': 3.8, 'M': 1.9, 'N': -3.5,
                    'P': -1.6, 'Q': -3.5, 'R': -4.5, 'S': -0.8,
                    'T': -0.7, 'V': 4.2, 'W': -0.9, 'Y': -1.3,
                   }, default=default)


def get_aa2mj(default=None):
    """
    Generate a Mapping from aminoacids to Miyazwa-Jernigan hydrophobicity scale.

    The mapping is equivalent to AAindex: MIYS850101 (Miyazawa-Jernigan, 1985).

    Args:
        default: default value for elements beyond the standard 20 aminoacids.

    Returns a Mapping object from aminoacids to Miyazwa-Jernigan hydrophobicity
        values.
    """
    return Mapping("miyazawa-jernigan",
                   {'A': 2.36, 'C': 3.36, 'D': 1.67, 'E': 1.74,
                    'F': 4.37, 'G': 2.06, 'H': 2.41, 'I': 4.17,
                    'K': 1.23, 'L': 3.93, 'M': 4.22, 'N': 1.70,
                    'P': 1.89, 'Q': 1.75, 'R': 1.92, 'S': 1.81,
                    'T': 2.04, 'V': 3.49, 'W': 3.82, 'Y': 2.91,
                   }, default=default)
