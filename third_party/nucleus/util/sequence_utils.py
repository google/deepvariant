# Copyright 2018 Google LLC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
"""Utility functions for manipulating DNA sequences."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


class Error(Exception):
  """Base error class."""


def _add_lowercase(d):
  """Returns a dictionary with the lowercase keys and values entered."""
  retval = d.copy()
  retval.update({k.lower(): v.lower() for k, v in d.items()})
  return retval


STRICT_DNA_COMPLEMENT_UPPER = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
DNA_COMPLEMENT_UPPER = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
IUPAC_DNA_COMPLEMENT_UPPER = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'R': 'Y',  # R is A/G
    'Y': 'R',  # Y is C/T
    'S': 'S',  # S is C/G
    'W': 'W',  # W is A/T
    'K': 'M',  # K is G/T
    'M': 'K',  # M is A/C
    'B': 'V',  # B is C/G/T
    'V': 'B',  # V is A/C/G
    'D': 'H',  # D is A/G/T
    'H': 'D',  # H is A/C/T
    'N': 'N',  # N is any base
}

IUPAC_TO_CANONICAL_BASES_UPPER = {
    'A': ['A'],
    'T': ['T'],
    'C': ['C'],
    'G': ['G'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['C', 'G'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'V': ['A', 'C', 'G'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'N': ['A', 'C', 'G', 'T'],
}

STRICT_DNA_COMPLEMENT = _add_lowercase(STRICT_DNA_COMPLEMENT_UPPER)
DNA_COMPLEMENT = _add_lowercase(DNA_COMPLEMENT_UPPER)
IUPAC_DNA_COMPLEMENT = _add_lowercase(IUPAC_DNA_COMPLEMENT_UPPER)


STRICT_DNA_BASES_UPPER = frozenset(['A', 'C', 'G', 'T'])
STRICT_DNA_BASES = frozenset(['a', 'c', 'g', 't', 'A', 'C', 'G', 'T'])
DNA_BASES_UPPER = frozenset(['A', 'C', 'G', 'T', 'N'])
DNA_BASES = frozenset(['a', 'c', 'g', 't', 'n', 'A', 'C', 'G', 'T', 'N'])


def reverse_complement(sequence, complement_dict=None):
  """Returns the reverse complement of a DNA sequence.

  By default this will successfully reverse complement sequences comprised
  solely of A, C, G, and T letters. Other complement dictionaries can be
  passed in for more permissive matching.

  Args:
    sequence: str. The input sequence to reverse complement.
    complement_dict: dict[str, str]. The lookup dictionary holding the
      complement base pairs.

  Returns:
    The reverse complement DNA sequence.

  Raises:
    Error: The sequence contains letters not present in complement_dict.
  """
  if complement_dict is None:
    complement_dict = STRICT_DNA_COMPLEMENT_UPPER

  try:
    return ''.join(complement_dict[nt] for nt in reversed(sequence))
  except KeyError:
    raise Error('Unknown base in {}, cannot reverse complement using {}'.format(
        sequence, str(complement_dict)))
