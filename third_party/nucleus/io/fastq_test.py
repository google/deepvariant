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
"""Tests for third_party.nucleus.io.fastq."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import fastq
from third_party.nucleus.protos import fastq_pb2
from third_party.nucleus.testing import test_utils


class FastqReaderTests(parameterized.TestCase):

  @parameterized.parameters('test_reads.fastq', 'test_reads.fastq.gz',
                            'test_reads.tfrecord', 'test_reads.tfrecord.gz')
  def test_iterate_fastq_reader(self, fastq_filename):
    fastq_path = test_utils.genomics_core_testdata(fastq_filename)
    expected_ids = [
        'NODESC:header', 'M01321:49:000000000-A6HWP:1:1101:17009:2216', 'FASTQ'
    ]
    with fastq.FastqReader(fastq_path) as reader:
      records = list(reader.iterate())
    self.assertLen(records, 3)
    self.assertEqual([r.id for r in records], expected_ids)


class FastqWriterTests(parameterized.TestCase):
  """Tests for FastqWriter."""

  def setUp(self):
    self.records = [
        fastq_pb2.FastqRecord(id='id1', sequence='ACGTG', quality='ABCDE'),
        fastq_pb2.FastqRecord(id='id2', sequence='ATTT', quality='ABC@'),
        fastq_pb2.FastqRecord(
            id='ID3',
            description='multi desc',
            sequence='GATAC',
            quality='ABC@!'),
    ]

  @parameterized.parameters('test_raw.fastq', 'test_zipped.fastq.gz',
                            'test_raw.tfrecord', 'test_zipped.tfrecord.gz')
  def test_roundtrip_writer(self, filename):
    output_path = test_utils.test_tmpfile(filename)
    with fastq.FastqWriter(output_path) as writer:
      for record in self.records:
        writer.write(record)

    with fastq.FastqReader(output_path) as reader:
      v2_records = list(reader.iterate())

    self.assertEqual(self.records, v2_records)


if __name__ == '__main__':
  absltest.main()
