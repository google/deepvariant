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
"""Tests for BedWriter CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from third_party.nucleus.io.python import bed_writer
from third_party.nucleus.protos import bed_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import io_utils

_DOUBLE_CLOSE_ERROR = 'Cannot close an already closed BedWriter'
_WRITE_TO_CLOSED_ERROR = 'Cannot write to closed BED stream'


class WrapBedWriterTest(parameterized.TestCase):

  def setUp(self):
    out_fname = test_utils.test_tmpfile('output.bed')
    self.writer = bed_writer.BedWriter.to_file(
        out_fname, bed_pb2.BedHeader(num_fields=12), bed_pb2.BedWriterOptions())
    self.expected_bed_content = [
        'chr1\t10\t20\tfirst\t100\t+\t12\t18\t255,124,1\t3\t2,6,2\t10,12,18\n',
        'chr1\t100\t200\tsecond\t250\t.\t120\t180\t252,122,12\t2\t35,40\t'
        '100,160\n'
    ]
    self.record = bed_pb2.BedRecord(
        reference_name='chr1', start=20, end=30, name='r')

  def test_writing_canned_records(self):
    """Tests writing all the records that are 'canned' in our tfrecord file."""
    # This file is in TFRecord format.
    tfrecord_file = test_utils.genomics_core_testdata(
        'test_regions.bed.tfrecord')

    header = bed_pb2.BedHeader(num_fields=12)
    writer_options = bed_pb2.BedWriterOptions()
    bed_records = list(
        io_utils.read_tfrecords(tfrecord_file, proto=bed_pb2.BedRecord))
    out_fname = test_utils.test_tmpfile('output.bed')
    with bed_writer.BedWriter.to_file(out_fname, header,
                                      writer_options) as writer:
      for record in bed_records:
        writer.write(record)

    with tf.gfile.GFile(out_fname, 'r') as f:
      self.assertEqual(f.readlines(), self.expected_bed_content)

  def test_context_manager(self):
    with self.writer:
      # Writing within the context manager succeeds.
      self.assertIsNone(self.writer.write(self.record))

    # self.writer should be closed, so writing again will fail.
    with self.assertRaisesRegexp(ValueError, _WRITE_TO_CLOSED_ERROR):
      self.writer.write(self.record)

  def test_double_context_manager(self):
    with self.writer:
      # Writing within the context manager succeeds.
      self.assertIsNone(self.writer.write(self.record))

    with self.assertRaisesRegexp(ValueError, _DOUBLE_CLOSE_ERROR):
      # Entering the closed writer should be fine.
      with self.writer:
        pass  # We want to raise an error on exit, so nothing to do in context.


if __name__ == '__main__':
  absltest.main()
