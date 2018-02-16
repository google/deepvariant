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
"""Tests for deepvariant.util.io."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.util.io import sam
from deepvariant.util.genomics import reads_pb2
from deepvariant.util.genomics import reference_pb2
from deepvariant.util import io_utils
from deepvariant.util import test_utils


class ReadWriterTests(parameterized.TestCase):
  """Tests for sam.SamWriter."""

  def setUp(self):
    self.read1 = test_utils.make_read(
        bases='ACCGT',
        chrom='chr1',
        start=10,
        cigar='5M',
        mapq=50,
        quals=range(30, 35),
        name='read1')
    self.read2 = test_utils.make_read(
        bases='AACCTT',
        chrom='chr2',
        start=15,
        cigar='7M',
        mapq=40,
        quals=range(20, 26),
        name='read2')
    self.contigs = [
        reference_pb2.ContigInfo(name='chr1'),
        reference_pb2.ContigInfo(name='chr2'),
    ]

  def test_make_read_writer_tfrecords(self):
    outfile = test_utils.test_tmpfile('test.tfrecord')
    writer = sam.SamWriter(outfile, contigs=[])

    # Test that the writer is a context manager and that we can write a read to
    # it.
    with writer:
      writer.write(self.read1)
      writer.write(self.read2)

    # Our output should have exactly one read in it.
    self.assertEqual([self.read1, self.read2],
                     list(
                         io_utils.read_tfrecords(outfile,
                                                 proto=reads_pb2.Read)))

  def test_make_read_writer_bam_fails_with_not_implemented_error(self):
    with self.assertRaises(NotImplementedError):
      sam.SamWriter('test.bam', contigs=self.contigs)


if __name__ == '__main__':
  absltest.main()
