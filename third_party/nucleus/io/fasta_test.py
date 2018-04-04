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

"""Tests for third_party.nucleus.io.fasta."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import fasta
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges


class RefFastaReaderTests(parameterized.TestCase):

  @parameterized.parameters('test.fasta', 'test.fasta.gz')
  def test_make_ref_reader_default(self, fasta_filename):
    fasta_path = test_utils.genomics_core_testdata(fasta_filename)
    with fasta.RefFastaReader(fasta_path) as reader:
      self.assertEqual(reader.query(ranges.make_range('chrM', 1, 6)), 'ATCAC')

  @parameterized.parameters('test.fasta', 'test.fasta.gz')
  def test_make_ref_reader_cache_specified(self, fasta_filename):
    fasta_path = test_utils.genomics_core_testdata(fasta_filename)
    with fasta.RefFastaReader(fasta_path, cache_size=10) as reader:
      self.assertEqual(reader.query(ranges.make_range('chrM', 1, 5)), 'ATCA')


class InMemoryRefReaderTests(parameterized.TestCase):

  @classmethod
  def setUpClass(cls):
    cls.fasta_reader = fasta.RefFastaReader(
        test_utils.genomics_core_testdata('test.fasta'))

    cls.in_mem = fasta.InMemoryRefReader(
        [(contig.name, 0,
          cls.fasta_reader.query(
              ranges.make_range(contig.name, 0, contig.n_bases)))
         for contig in cls.fasta_reader.header.contigs])

  def test_non_zero_start_query(self):
    """Checks all of the ways we can construct an InMemoryRefReader."""
    bases = 'ACGTAACCGGTT'
    for start in range(len(bases)):
      reader = fasta.InMemoryRefReader([('1', start, bases[start:])])
      self.assertEqual(reader.header.contigs[0].name, '1')
      self.assertEqual(reader.header.contigs[0].n_bases, len(bases))

      # Check that our query operation works as expected with a start position.
      for end in range(start, len(bases)):
        self.assertEqual(reader.query(ranges.make_range('1', start, end)),
                         bases[start:end])

  @parameterized.parameters(
      # Start is 10, so this raises because it's before the bases starts.
      dict(start=0, end=1),
      # Spans into the start of the bases; make sure it detects it's bad.
      dict(start=8, end=12),
      # Spans off the end of the bases.
      dict(start=12, end=15),
  )
  def test_bad_query_with_start(self, start, end):
    reader = fasta.InMemoryRefReader([('1', 10, 'ACGT')])
    with self.assertRaises(ValueError):
      reader.query(ranges.make_range('1', start, end))

  def test_contigs(self):
    # Our contigs can have a different order, descriptions are dropped, etc so
    # we need to check specific fields by hand.
    fasta_contigs = {
        contig.name: contig for contig in self.fasta_reader.header.contigs
    }
    mem_contigs = {contig.name: contig for contig in self.in_mem.header.contigs}

    self.assertEqual(fasta_contigs.keys(), mem_contigs.keys())
    for name, fasta_contig in fasta_contigs.iteritems():
      self.assertContigsAreEqual(mem_contigs[name], fasta_contig)

  def assertContigsAreEqual(self, actual, expected):
    self.assertEqual(actual.name, expected.name)
    self.assertEqual(actual.n_bases, expected.n_bases)
    self.assertEqual(actual.pos_in_fasta, expected.pos_in_fasta)

  def test_iterate(self):
    with self.assertRaises(NotImplementedError):
      self.in_mem.iterate()

  @parameterized.parameters(
      dict(region=ranges.make_range('chr1', 0, 10), expected=True),
      dict(region=ranges.make_range('chr1', 10, 50), expected=True),
      dict(region=ranges.make_range('chr1', 10, 500), expected=False),
      dict(region=ranges.make_range('chr3', 10, 20), expected=False),
  )
  def test_is_valid(self, region, expected):
    self.assertEqual(self.in_mem.is_valid(region), expected)
    self.assertEqual(self.fasta_reader.is_valid(region), expected)

  def test_known_contig(self):
    for contig in self.fasta_reader.header.contigs:
      self.assertContigsAreEqual(self.in_mem.contig(contig.name), contig)

  def test_str_and_repr(self):
    self.assertIsInstance(str(self.in_mem), basestring)
    self.assertIsInstance(repr(self.in_mem), basestring)

  def test_unknown_contig(self):
    for reader in [self.fasta_reader, self.in_mem]:
      with self.assertRaises(ValueError):
        reader.contig('unknown')

  def test_good_query(self):
    for contig in self.fasta_reader.header.contigs:
      for start in range(contig.n_bases):
        for end in range(start, contig.n_bases):
          region = ranges.make_range(contig.name, start, end)
          self.assertEqual(
              self.in_mem.query(region), self.fasta_reader.query(region))

  @parameterized.parameters(
      ranges.make_range('chr1', -1, 10),     # bad start.
      ranges.make_range('chr1', 10, 1),      # end < start.
      ranges.make_range('chr1', 0, 1000),    # off end of chromosome.
      ranges.make_range('unknown', 0, 10),   # unknown chromosome.
  )
  def test_bad_query(self, region):
    for reader in [self.fasta_reader, self.in_mem]:
      with self.assertRaises(ValueError):
        reader.query(region)


if __name__ == '__main__':
  absltest.main()
