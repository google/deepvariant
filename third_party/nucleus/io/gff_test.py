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
"""Tests for third_party.nucleus.io.gff."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import gff
from third_party.nucleus.io import tfrecord
from third_party.nucleus.protos import gff_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges

# Names of testdata GFF files; we also reuse these basenames for output files
# in the tmp directory.
TEXT_GFF_FILES = ('test_features.gff', 'test_features.gff.gz')
TFRECORD_GFF_FILES = ('test_features.gff.tfrecord',
                      'test_features.gff.tfrecord.gz')
ALL_GFF_FILES = TEXT_GFF_FILES + TFRECORD_GFF_FILES

EXPECTED_GFF_VERSION = 'gff-version 3.2.1'


class GffReaderTests(parameterized.TestCase):

  @parameterized.parameters(*ALL_GFF_FILES)
  def test_iterate_gff_reader(self, gff_filename):
    gff_path = test_utils.genomics_core_testdata(gff_filename)
    expected = [('ctg123', 999, 9000), ('ctg123', 999, 1012)]

    with gff.GffReader(gff_path) as reader:
      records = list(reader.iterate())
    self.assertLen(records, 2)
    self.assertEqual(
        [(r.range.reference_name, r.range.start, r.range.end) for r in records],
        expected)

  @parameterized.parameters(*TEXT_GFF_FILES)
  def test_native_gff_header(self, gff_filename):
    gff_path = test_utils.genomics_core_testdata(gff_filename)
    with gff.GffReader(gff_path) as reader:
      self.assertEqual(EXPECTED_GFF_VERSION, reader.header.gff_version)
    with gff.NativeGffReader(gff_path) as native_reader:
      self.assertEqual(EXPECTED_GFF_VERSION, native_reader.header.gff_version)


class GffWriterTests(parameterized.TestCase):
  """Tests for GffWriter."""

  def setUp(self):
    tfrecord_file = test_utils.genomics_core_testdata(
        'test_features.gff.tfrecord')
    self.records = list(
        tfrecord.read_tfrecords(tfrecord_file, proto=gff_pb2.GffRecord))
    self.header = gff_pb2.GffHeader(
        sequence_regions=[ranges.make_range('ctg123', 0, 1497228)])

  @parameterized.parameters(*ALL_GFF_FILES)
  def test_roundtrip_writer(self, filename):
    output_path = test_utils.test_tmpfile(filename)
    with gff.GffWriter(output_path, header=self.header) as writer:
      for record in self.records:
        writer.write(record)

    with gff.GffReader(output_path) as reader:
      v2_records = list(reader.iterate())

    self.assertEqual(self.records, v2_records)


if __name__ == '__main__':
  absltest.main()
