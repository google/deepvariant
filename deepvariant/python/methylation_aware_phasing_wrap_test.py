# Copyright 2026 Google LLC.
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
"""Tests for the methylation-aware phasing Python binding."""

from absl.testing import absltest
import numpy as np

from deepvariant.protos import deepvariant_pb2
from deepvariant.python import methylation_aware_phasing
from third_party.nucleus.protos import reads_pb2


class MethylationAwarePhasingWrapTest(absltest.TestCase):

  def test_empty_inputs(self):
    phases, p_values = methylation_aware_phasing.phase([], [], [])
    self.assertEmpty(phases)
    self.assertEmpty(p_values)

  def test_reads_without_methylated_sites_preserve_initial_phases(self):
    reads = [
        reads_pb2.Read(fragment_name=f'read_{i}', read_number=1)
        for i in range(3)
    ]
    phases, p_values = methylation_aware_phasing.phase(
        reads, [1, 0, 2], []
    )
    self.assertEqual(phases, [1, 0, 2])
    self.assertEmpty(p_values)

  def test_accepts_numpy_candidate_array(self):
    # RegionProcessor passes a NumPy array after selecting reference sites.
    site = deepvariant_pb2.DeepVariantCall(methylation_p_value=0.25)
    reads = [reads_pb2.Read(fragment_name='read', read_number=1)]
    phases, p_values = methylation_aware_phasing.phase(
        reads, [1], np.asarray([site], dtype=object)
    )
    self.assertEqual(phases, [1])
    self.assertEqual(p_values, [0.25])

  def test_accepts_empty_numpy_candidate_array(self):
    phases, p_values = methylation_aware_phasing.phase(
        [], [], np.asarray([], dtype=object)
    )
    self.assertEmpty(phases)
    self.assertEmpty(p_values)


if __name__ == '__main__':
  absltest.main()
