# Copyright 2018 Google Inc.
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

"""Utility functions for working with reads."""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import numpy as np

from third_party.nucleus.util import cigar
from third_party.nucleus.util import ranges


def read_range(read):
  """Creates a Range proto from the alignment of Read.

  Args:
    read: the read to calculate range

  Returns:
    A third_party.nucleus.protos.Range for read.
  """
  start = read.alignment.position.position
  end = start + cigar.alignment_length(read.alignment.cigar)
  return ranges.make_range(read.alignment.position.reference_name, start, end)


def reservoir_sample(iterable, k, random=None):
  """Samples k elements with uniform probability from an iterable.

  Selects a subset of k elements from n input elements with uniform probability
  without needing to hold all n elements in memory at the same time. This
  implementation has max space complexity O(min(k, n)), i.e., we allocate up to
  min(k, n) elements to store the samples. This means that we only use ~n
  elements when n is smaller than k, which can be important when k is large. If
  n elements are added to this sampler, and n <= k, all n elements will be
  retained. If n > k, each added element will be retained with a uniform
  probability of k / n.

  The order of the k retained samples from our n elements is undefined. In
  particular that means that the elements in the returned list can occur in a
  different order than they appeared in the iterable.

  More details about reservoir sampling (and the specific algorithm used here
  called Algorithm R) can be found on wikipedia:

  https://en.wikipedia.org/wiki/Reservoir_sampling#Algorithm_R

  Args:
    iterable: Python iterable. The iterable to sample from.
    k: int. The number of elements to sample.
    random: A random number generator or None.
  Returns:
    A list containing the k sampled elements.
  Raises:
    ValueError: If k is negative.
  """
  if k < 0:
    raise ValueError('k must be nonnegative, but got {}'.format(k))
  if random is None:
    random = np.random
  sample = []
  for i, item in enumerate(iterable):
    if len(sample) < k:
      sample.append(item)
    else:
      j = random.randint(0, i + 1)
      if j < k:
        sample[j] = item
  return sample
