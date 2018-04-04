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
"""Tests for FastqWriter CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from third_party.nucleus.io import fastq
from third_party.nucleus.io.python import fastq_writer
from third_party.nucleus.protos import fastq_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import io_utils

_DOUBLE_CLOSE_ERROR = 'Cannot close an already closed FastqWriter'
_WRITE_TO_CLOSED_ERROR = 'Cannot write to closed FASTQ stream'


class WrapFastqWriterTest(parameterized.TestCase):

  def setUp(self):
    writer_options = fastq_pb2.FastqWriterOptions()
    out_fname = test_utils.test_tmpfile('output.fastq')
    self.writer = fastq_writer.FastqWriter.to_file(out_fname, writer_options)
    self.expected_fastq_content = [
        '@NODESC:header\n', 'GATTACA\n', '+\n', 'BB>B@FA\n',
        '@M01321:49:000000000-A6HWP:1:1101:17009:2216 1:N:0:1\n',
        'CGTTAGCGCAGGGGGCATCTTCACACTGGTGACAGGTAACCGCCGTAGTAAAGGTTCCGCCTTTCACT\n',
        '+\n',
        'AAAAABF@BBBDGGGG?FFGFGHBFBFBFABBBHGGGFHHCEFGGGGG?FGFFHEDG3EFGGGHEGHG\n',
        '@FASTQ contains multiple spaces in description\n',
        'CGGCTGGTCAGGCTGACATCGCCGCCGGCCTGCAGCGAGCCGCTGC\n', '+\n',
        'FAFAF;F/9;.:/;999B/9A.DFFF;-->.AAB/FC;9-@-=;=.\n'
    ]
    self.record = fastq_pb2.FastqRecord(
        id='ID', description='desc', sequence='ACGTAC', quality='ABCDEF')

  def test_writing_canned_records(self):
    """Tests writing all the variants that are 'canned' in our tfrecord file."""
    # This file is in TFRecord format.
    tfrecord_file = test_utils.genomics_core_testdata('test_reads.tfrecord')

    writer_options = fastq_pb2.FastqWriterOptions()
    fastq_records = list(
        io_utils.read_tfrecords(tfrecord_file, proto=fastq_pb2.FastqRecord))
    out_fname = test_utils.test_tmpfile('output.fastq')
    with fastq_writer.FastqWriter.to_file(out_fname, writer_options) as writer:
      for record in fastq_records:
        writer.write(record)

    with tf.gfile.GFile(out_fname, 'r') as f:
      self.assertEqual(f.readlines(), self.expected_fastq_content)

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


class WrapFastqWriterRoundTripTests(parameterized.TestCase):

  @parameterized.parameters('test_reads.fastq', 'test_reads.fastq.gz')
  def test_round_trip_fastq(self, test_datum_name):
    # Round-trip FASTQ records through writing and reading:
    # 1. Read records v1 from FastqReader;
    # 2. Write v1 to fastq using our FastqWriter;
    # 3. Read back in using FastqReader -- v2;
    # 4. compare v1 and v2.
    in_file = test_utils.genomics_core_testdata(test_datum_name)
    out_file = test_utils.test_tmpfile('output_' + test_datum_name)

    v1_reader = fastq.FastqReader(in_file)
    v1_records = list(v1_reader.iterate())
    self.assertTrue(v1_records, 'Reader failed to find records')

    writer_options = fastq_pb2.FastqWriterOptions()

    with fastq_writer.FastqWriter.to_file(out_file, writer_options) as writer:
      for record in v1_records:
        writer.write(record)

    v2_reader = fastq.FastqReader(out_file)
    v2_records = list(v2_reader.iterate())
    self.assertEqual(v1_records, v2_records,
                     'Round-tripped FASTQ files not as expected')


if __name__ == '__main__':
  absltest.main()
