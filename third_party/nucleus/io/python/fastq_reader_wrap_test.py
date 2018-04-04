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
"""Tests for fastq_reader CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import clif_postproc
from third_party.nucleus.io.python import fastq_reader
from third_party.nucleus.protos import fastq_pb2
from third_party.nucleus.testing import test_utils


class FastqReaderTest(parameterized.TestCase):

  def setUp(self):
    self.fastq = test_utils.genomics_core_testdata('test_reads.fastq')
    self.zipped_fastq = test_utils.genomics_core_testdata('test_reads.fastq.gz')
    self.options = fastq_pb2.FastqReaderOptions()

  def test_fastq_iterate(self):
    with fastq_reader.FastqReader.from_file(self.fastq, self.options) as reader:
      iterable = reader.iterate()
      self.assertIsInstance(iterable, clif_postproc.WrappedCppIterable)
      self.assertEqual(test_utils.iterable_len(iterable), 3)

    zreader = fastq_reader.FastqReader.from_file(
        self.zipped_fastq,
        fastq_pb2.FastqReaderOptions(
            compression_type=fastq_pb2.FastqReaderOptions.GZIP))
    with zreader:
      ziterable = zreader.iterate()
      self.assertIsInstance(ziterable, clif_postproc.WrappedCppIterable)
      self.assertEqual(test_utils.iterable_len(ziterable), 3)

  def test_from_file_raises_with_missing_fastq(self):
    with self.assertRaisesRegexp(ValueError,
                                 'Not found: Could not open missing.fastq'):
      fastq_reader.FastqReader.from_file('missing.fastq', self.options)

  def test_ops_on_closed_reader_raise(self):
    reader = fastq_reader.FastqReader.from_file(self.fastq, self.options)
    with reader:
      pass
    # At this point the reader is closed.
    with self.assertRaisesRegexp(ValueError, 'Cannot Iterate a closed'):
      reader.iterate()

  @parameterized.parameters('malformed.fastq', 'malformed2.fastq')
  def test_fastq_iterate_raises_on_malformed_record(self, filename):
    malformed = test_utils.genomics_core_testdata(filename)
    reader = fastq_reader.FastqReader.from_file(malformed, self.options)
    iterable = iter(reader.iterate())
    self.assertIsNotNone(next(iterable))
    with self.assertRaises(ValueError):
      list(iterable)


if __name__ == '__main__':
  absltest.main()
