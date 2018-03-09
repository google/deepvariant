"""Tests for deepvariant.util.in_memory_vcf_reader."""


from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.util.io import vcf
from deepvariant.util.genomics import reference_pb2
from deepvariant.util import in_memory_vcf_reader
from deepvariant.util import ranges
from deepvariant.util import test_utils


class InMemoryVcfReaderTests(parameterized.TestCase):
  """Test the functionality provided by vcf.InMemoryVcfReader."""

  def setUp(self):
    self.variants = [
        test_utils.make_variant(chrom='1', start=10),
        test_utils.make_variant(chrom='1', start=20),
        test_utils.make_variant(chrom='1', start=30),
        test_utils.make_variant(chrom='2', start=25),
        test_utils.make_variant(chrom='2', start=55),
        test_utils.make_variant(chrom='3', start=10),
    ]
    self.header = vcf.VcfHeader(
        contigs=[
            reference_pb2.ContigInfo(name='1', n_bases=100),
            reference_pb2.ContigInfo(name='2', n_bases=100),
            reference_pb2.ContigInfo(name='3', n_bases=100),
            reference_pb2.ContigInfo(name='4', n_bases=100),
        ],
        filters=[],
        samples=['NA12878'])
    self.reader = in_memory_vcf_reader.InMemoryVcfReader(
        self.variants, self.header)

  def test_iterate(self):
    """Tests that iterate returns an iterable containing our variants."""
    self.assertEqual(list(self.reader.iterate()), self.variants)

  def test_header(self):
    """Tests that the reader provides us back the header we gave it."""
    self.assertEqual(self.reader.header, self.header)

  @parameterized.parameters(
      dict(query='1', expected_variant_indices=[0, 1, 2]),
      dict(query='2', expected_variant_indices=[3, 4]),
      dict(query='3', expected_variant_indices=[5]),
      dict(query='4', expected_variant_indices=[]),
      dict(query='1:1-15', expected_variant_indices=[0]),
      dict(query='1:1-25', expected_variant_indices=[0, 1]),
      dict(query='1:1-35', expected_variant_indices=[0, 1, 2]),
      dict(query='1:15-35', expected_variant_indices=[1, 2]),
      dict(query='1:25-35', expected_variant_indices=[2]),
  )
  def test_query(self, query, expected_variant_indices):
    range1 = ranges.parse_literal(query, ranges.contigs_dict(
        self.header.contigs))
    self.assertEqual(
        list(self.reader.query(range1)),
        [self.variants[i] for i in expected_variant_indices])


if __name__ == '__main__':
  absltest.main()
