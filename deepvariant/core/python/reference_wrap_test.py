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
"""Tests for GenomeReference CLIF python wrappers."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.core import ranges
from deepvariant.core import test_utils
from deepvariant.core.python import reference_fai


class WrapReferenceTest(parameterized.TestCase):

  @parameterized.parameters('test.fasta', 'test.fasta.gz')
  def test_wrap(self, fasta_filename):
    chr_names = ['chrM', 'chr1', 'chr2']
    chr_lengths = [100, 76, 121]
    fasta = test_utils.genomics_core_testdata(fasta_filename)
    fai = test_utils.genomics_core_testdata(fasta_filename + '.fai')
    with reference_fai.GenomeReferenceFai.from_file(fasta, fai) as ref:
      self.assertEqual(ref.n_contigs, 3)
      self.assertIn(fasta, ref.fasta_path)
      self.assertIn('GenomeReference backed by htslib FAI index', str(ref))
      self.assertEqual(ref.contig_names, chr_names)
      self.assertEqual(ref.n_bp, sum(chr_lengths))
      self.assertEqual(ref.bases(ranges.make_range('chrM', 1, 10)), 'ATCACAGGT')

      self.assertTrue(ref.is_valid_interval(ranges.make_range('chrM', 1, 10)))
      self.assertFalse(
          ref.is_valid_interval(ranges.make_range('chrM', 1, 100000)))

      self.assertEqual(len(ref.contigs), 3)
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
    with self.assertRaisesRegexp(
        ValueError,
        'Not found: could not load fasta and/or fai for fasta ' + fasta):
      reference_fai.GenomeReferenceFai.from_file(fasta, fai)


if __name__ == '__main__':
  absltest.main()
