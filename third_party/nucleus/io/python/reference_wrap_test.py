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
"""Tests for GenomeReference CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
import sys
if 'google' in sys.modules and 'google.protobuf' not in sys.modules:
  del sys.modules['google']


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io.python import reference
from third_party.nucleus.protos import fasta_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges


class WrapReferenceTest(parameterized.TestCase):

  @parameterized.parameters(
      ('test.fasta', False, 'TAACC'), ('test.fasta', True, 'TaaCC'),
      ('test.fasta.gz', False, 'TAACC'), ('test.fasta.gz', True, 'TaaCC'))
  def test_wrap(self, fasta_filename, keep_true_case, expected_bases):
    chr_names = ['chrM', 'chr1', 'chr2']
    chr_lengths = [100, 76, 121]
    fasta = test_utils.genomics_core_testdata(fasta_filename)
    fai = test_utils.genomics_core_testdata(fasta_filename + '.fai')
    options = fasta_pb2.FastaReaderOptions(keep_true_case=keep_true_case)
    with reference.IndexedFastaReader.from_file(fasta, fai, options) as ref:
      self.assertEqual(ref.contig_names, chr_names)
      self.assertEqual(
          ref.bases(ranges.make_range('chrM', 22, 27)), expected_bases)

      self.assertTrue(ref.is_valid_interval(ranges.make_range('chrM', 1, 10)))
      self.assertFalse(
          ref.is_valid_interval(ranges.make_range('chrM', 1, 100000)))

      self.assertLen(ref.contigs, 3)
      self.assertEqual([c.name for c in ref.contigs], chr_names)
      self.assertEqual([c.n_bases for c in ref.contigs], chr_lengths)
      for contig in ref.contigs:
        self.assertEqual(ref.contig(contig.name), contig)
        self.assertTrue(ref.has_contig(contig.name))
        self.assertFalse(ref.has_contig(contig.name + '.unknown'))

  @parameterized.parameters(
      # The fasta and the FAI are both missing.
      ('missing.fasta', 'missing.fasta.fai'),
      # The fasta is present but the FAI is missing.
      ('test.fasta', 'missing.fasta.fai'),
      # The fasta is missing but the FAI is present.
      ('missing.fasta', 'test.fasta.fai'),
  )
  def test_from_file_raises_with_missing_inputs(self, fasta_filename,
                                                fai_filename):
    fasta = test_utils.genomics_core_testdata(fasta_filename)
    fai = test_utils.genomics_core_testdata(fai_filename)
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegexp(
        ValueError, 'could not load fasta and/or fai for fasta ' + fasta):
      reference.IndexedFastaReader.from_file(fasta, fai,
                                             fasta_pb2.FastaReaderOptions())


if __name__ == '__main__':
  absltest.main()
