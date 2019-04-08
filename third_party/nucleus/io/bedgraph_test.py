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
"""Tests for third_party.nucleus.io.bedgraph."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import bedgraph
from third_party.nucleus.protos import bedgraph_pb2
from third_party.nucleus.testing import test_utils


class BedGraphTests(parameterized.TestCase):

  @parameterized.parameters('test_regions.bedgraph', 'test_regions.bedgraph.gz')
  def test_iterate_bedgraph_reader(self, bedgraph_path):
    bedgraph_path = test_utils.genomics_core_testdata(bedgraph_path)
    expected = [('chr1', 10, 20, 100), ('chr1', 100, 200, 250),
                ('chr1', 300, 400, 150.1), ('chr1', 500, 501, 20.13)]
    with bedgraph.BedGraphReader(bedgraph_path) as reader:
      records = list(reader.iterate())
    self.assertLen(records, 4)
    self.assertEqual(
        [(r.reference_name, r.start, r.end, r.data_value) for r in records],
        expected)

  @parameterized.parameters('test_regions.bedgraph', 'test_regions.bedgraph.gz')
  def test_roundtrip_writer(self, bedgraph_path):
    output_path = test_utils.test_tmpfile(bedgraph_path)
    input_path = test_utils.genomics_core_testdata(bedgraph_path)
    records = []
    with bedgraph.BedGraphReader(input_path) as reader:
      records = list(reader.iterate())

    with bedgraph.BedGraphWriter(output_path) as writer:
      for record in records:
        writer.write(record)

    with bedgraph.BedGraphReader(output_path) as reader:
      v2_records = list(reader.iterate())

    self.assertLen(records, 4)
    self.assertEqual(records, v2_records)


if __name__ == '__main__':
  absltest.main()
