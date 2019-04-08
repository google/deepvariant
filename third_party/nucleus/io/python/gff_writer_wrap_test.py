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
"""Tests for GffWriter CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import tfrecord
from third_party.nucleus.io.python import gff_writer
from third_party.nucleus.protos import gff_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges

_DOUBLE_CLOSE_ERROR = 'Cannot close an already closed GffWriter'
_WRITE_TO_CLOSED_ERROR = 'Cannot write to closed GFF stream'


class WrapGffWriterTest(parameterized.TestCase):

  def setUp(self):
    out_fname = test_utils.test_tmpfile('output.gff')
    self.writer = gff_writer.GffWriter.to_file(out_fname, gff_pb2.GffHeader(),
                                               gff_pb2.GffWriterOptions())
    self.expected_gff_content = open(
        test_utils.genomics_core_testdata('test_features.gff')).readlines()
    self.header = gff_pb2.GffHeader(
        sequence_regions=[ranges.make_range('ctg123', 0, 1497228)])
    self.record = gff_pb2.GffRecord(
        range=ranges.make_range('ctg123', 1000, 1100))

  def test_writing_canned_records(self):
    """Tests writing all the records that are 'canned' in our tfrecord file."""
    # This file is in TFRecord format.
    tfrecord_file = test_utils.genomics_core_testdata(
        'test_features.gff.tfrecord')
    writer_options = gff_pb2.GffWriterOptions()
    gff_records = list(
        tfrecord.read_tfrecords(tfrecord_file, proto=gff_pb2.GffRecord))
    out_fname = test_utils.test_tmpfile('output.gff')
    with gff_writer.GffWriter.to_file(out_fname, self.header,
                                      writer_options) as writer:
      for record in gff_records:
        writer.write(record)

    with open(out_fname) as f:
      self.assertEqual(f.readlines(), self.expected_gff_content)

  def test_context_manager(self):
    with self.writer:
      # Writing within the context manager succeeds.
      self.assertIsNone(self.writer.write(self.record))

    # self.writer should be closed, so writing again will fail.
    with self.assertRaisesRegexp(ValueError, _WRITE_TO_CLOSED_ERROR):
      self.writer.write(self.record)

  def test_double_close(self):
    with self.writer:
      # Writing within the context manager succeeds.
      self.assertIsNone(self.writer.write(self.record))

    with self.assertRaisesRegexp(ValueError, _DOUBLE_CLOSE_ERROR):
      # Entering the closed writer should be fine.
      with self.writer:
        pass  # We want to raise an error on exit, so nothing to do in context.


if __name__ == '__main__':
  absltest.main()
