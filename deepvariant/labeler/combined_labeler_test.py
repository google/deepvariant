# Copyright 2025 Google LLC.
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

from absl.testing import absltest
from absl.testing import parameterized

from deepvariant.labeler import combined_labeler
from third_party.nucleus.io import fasta
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import ranges


def _test_variant(start=10, alleles=('A', 'C'), gt=None):
  variant = variants_pb2.Variant(
      reference_name='chr1',
      start=start,
      end=start + len(alleles[0]),
      reference_bases=alleles[0],
      alternate_bases=alleles[1:],
  )

  if gt:
    variant.calls.add(genotype=gt)

  return variant


class CombinedLabelerTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self.confident_regions = ranges.RangeSet(
        [ranges.make_range('chr1', 0, 10000000)]
    )

  def _make_labeler(self, truths, ref_reader):
    return combined_labeler.CombinedLabeler(
        truth_vcf_reader=vcf.InMemoryVcfReader(truths),
        ref_reader=ref_reader,
        confident_regions=self.confident_regions,
    )

  def assertGetsCorrectLabels(
      self,
      candidates,
      true_variants,
      ref,
      expected_genotypes,
  ):
    labeler = self._make_labeler(true_variants, fasta.InMemoryFastaReader(ref))
    labeled_variants = list(
        labeler.label_variants(
            candidates, ranges.make_range('chr1', 279767, 279782)
        )
    )
    self.assertIsNotNone(labeled_variants)

    # Check that the genotypes of our labeled variants are the ones we expect.
    self.assertEqual(
        [variant_label.genotype for variant_label in labeled_variants],
        [tuple(x) for x in expected_genotypes],
    )

  # This is a real case from HG002 sample where the haplotype labeler cannot
  # correctly label 1 variant while positional labeler resolves it correctly.
  def test_haplotype_labeler_over_positional(self):
    """Tests that haplotype labeler is used when positional is hom-ref."""
    self.assertGetsCorrectLabels(
        candidates=[
            _test_variant(start=3512103, alleles=['A', 'AT']),
            _test_variant(start=3512104, alleles=['TA', 'T']),
            _test_variant(start=3512105, alleles=['A', 'T']),
            _test_variant(start=3512107, alleles=['A', 'T']),
        ],
        true_variants=[
            _test_variant(start=3512105, alleles=['A', 'T'], gt=[1, 0]),
            _test_variant(start=3512107, alleles=['A', 'T'], gt=[1, 1]),
        ],
        ref=[('chr1', 3512100, 'TATATATATTTTTTTTTT')],
        expected_genotypes=[
            [0, 0],
            [0, 0],
            [0, 1],
            [1, 1],
        ],
    )


if __name__ == '__main__':
  absltest.main()
