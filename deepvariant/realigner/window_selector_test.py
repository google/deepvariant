# Copyright 2017 Google Inc.
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
"""Tests for deepvariant.realigner.window_selector."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest

from deepvariant import test_utils
from deepvariant.core.genomics import range_pb2
from deepvariant.protos import realigner_pb2
from deepvariant.realigner.window_selector import WindowSelector


class WindowSelectorTest(absltest.TestCase):

  def test_ws_config(self):
    return realigner_pb2.RealignerOptions.WindowSelectorOptions(
        min_num_supporting_reads=2,
        max_num_supporting_reads=10,
        min_mapq=20,
        min_base_quality=20,
        min_windows_distance=4)

  def test_process_read(self):
    """Test WindowSelector.process_read()."""
    window = WindowSelector(self.test_ws_config())

    ref = 'A' * 100

    read_1 = test_utils.make_read(
        'AAGA', start=10, cigar='4M', quals=[64] * 4, name='read_1')
    read_2 = test_utils.make_read(
        'AAGTA', start=10, cigar='2M2I1M', quals=[64] * 5, name='read_2')
    read_3 = test_utils.make_read(
        'AAA', start=10, cigar='2M2D1M', quals=[64] * 3, name='read_3')
    read_4 = test_utils.make_read(
        'TGATAC', start=10, cigar='2S3M1S', quals=[64] * 6, name='read_4')
    read_5 = test_utils.make_read(
        'AAGA', start=10, cigar='2M1X1M', quals=[64] * 4, name='read_5')

    self.assertEqual(list(window.process_read(ref, read_1)), [12])
    self.assertEqual(list(window.process_read(ref, read_2)), [10, 11, 12, 13])
    self.assertEqual(list(window.process_read(ref, read_3)), [12, 13])
    self.assertEqual(list(window.process_read(ref, read_4)), [8, 9, 11, 13])
    self.assertEqual(list(window.process_read(ref, read_5)), [12])

  def test_candidate_pos_low_qual(self):
    """Test WindowSelector.process_read() with reads of low quality."""
    window = WindowSelector(self.test_ws_config())

    ref = 'A' * 100

    read_1 = test_utils.make_read(
        'AAGA', start=10, cigar='4M', quals=[64, 64, 10, 30], name='read_1')
    read_2 = test_utils.make_read(
        'AAGTA',
        start=10,
        cigar='2M2I1M',
        quals=[64, 64, 10, 30, 64],
        name='read_2')
    read_3 = test_utils.make_read(
        'TGATAC',
        start=10,
        cigar='2S3M1S',
        quals=[64, 10, 64, 64, 64, 64],
        name='read_3')
    read_4 = test_utils.make_read(
        'AAGA', start=10, cigar='2M1X1M', quals=[64, 64, 30, 10], name='read_4')

    self.assertEqual(list(window.process_read(ref, read_1)), [])
    self.assertEqual(list(window.process_read(ref, read_2)), [11, 13])
    self.assertEqual(list(window.process_read(ref, read_3)), [8, 11, 13])
    self.assertEqual(list(window.process_read(ref, read_4)), [12])

  def test_windows(self):
    """Test WindowSelector.windows()."""
    window = WindowSelector(self.test_ws_config())
    candidates = {0: 2, 2: 4, 3: 11, 8: 3}

    self.assertEqual(
        list(window.windows(candidates, 'ref', 0)), [
            range_pb2.Range(reference_name='ref', start=-4, end=6),
            range_pb2.Range(reference_name='ref', start=4, end=12)
        ])


if __name__ == '__main__':
  absltest.main()
