# Copyright 2017 Google LLC.
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
"""Tests for deepvariant .variant_labeler."""

import collections


from absl.testing import absltest
from absl.testing import parameterized

from third_party.nucleus.io import vcf
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from deepvariant import testdata
from deepvariant.labeler import customized_classes_labeler
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import vcf_constants


def setUpModule():
  testdata.init()


FakeVCFObject = collections.namedtuple('FakeVCFObject', ['field_access_cache'])
CUSTOMIZED_INFO_FIELD_NAME = 'type'
CUSTOMIZED_CLASSES_LIST = 'ref,class1,class2'


def _add_class_to_variant(variant, class_status):
  if class_status is None:
    return variant

  header = variants_pb2.VcfHeader(infos=[
      variants_pb2.VcfInfo(
          id=CUSTOMIZED_INFO_FIELD_NAME,
          number='A',
          type=vcf_constants.STRING_TYPE,
          description='Customized class label for the variant.')
  ])
  my_cache = vcf.VcfHeaderCache(header)
  vcf_object = FakeVCFObject(field_access_cache=my_cache)
  variant_utils.set_info(
      variant, CUSTOMIZED_INFO_FIELD_NAME, class_status, vcf_object=vcf_object)
  return variant


class CustomizedClassesVariantLabelerTest(parameterized.TestCase):

  # Confident variants: SNP, deletion, and multi-allelic.
  snp_class1 = _add_class_to_variant(
      test_utils.make_variant(start=10, alleles=['A', 'C'], gt=[0, 1]),
      class_status='class1')

  snp_class2 = _add_class_to_variant(
      test_utils.make_variant(start=20, alleles=['ACG', 'A'], gt=[1, 1]),
      class_status='class2')

  multiallelic = _add_class_to_variant(
      test_utils.make_variant(
          start=30, alleles=['ACT', 'ACTGT', 'A'], gt=[1, 2]),
      class_status='class2')

  # Outside our confident regions.
  non_confident = _add_class_to_variant(
      test_utils.make_variant(start=200, alleles=['A', 'C'], gt=[0, 1]),
      class_status='class1')

  filtered = _add_class_to_variant(
      test_utils.make_variant(start=40, filters='FAILED', gt=[0, 1]),
      class_status='class1')

  # TODO: check the following cases:
  # no_class_status
  # invalid_class_status
  # (Value error should be produced in both cases)

  variants = [snp_class1, snp_class2, multiallelic, non_confident, filtered]

  def _make_labeler(self, variants, confident_regions):
    return customized_classes_labeler.CustomizedClassesVariantLabeler(
        truth_vcf_reader=vcf.InMemoryVcfReader(variants),
        confident_regions=confident_regions,
        classes_list=CUSTOMIZED_CLASSES_LIST,
        info_field_name=CUSTOMIZED_INFO_FIELD_NAME)

  @parameterized.parameters(
      # Simple tests: we get back our matching variants in the confident regions
      dict(
          candidate=snp_class1,
          expected_confident=True,
          expected_truth=snp_class1,
          expected_label=1),
      dict(
          candidate=snp_class2,
          expected_confident=True,
          expected_truth=snp_class2,
          expected_label=2),
      # For multiallelic variants, we default to class 0.
      dict(
          candidate=multiallelic,
          expected_confident=True,
          expected_truth=multiallelic,
          expected_label=2),

      # Test the behavior outside of our confident regions.
      # If we provide a variant outside the confident regions (non_confident) we
      # don't get back any expected_truth variants.
      dict(
          candidate=non_confident,
          expected_confident=False,
          expected_truth=None,
          expected_label=0),
      # No matching variant, so we get a None as well as False.
      dict(
          candidate=test_utils.make_variant(start=300, alleles=['A', 'C']),
          expected_confident=False,
          expected_truth=None,
          expected_label=0),

      # This variant doesn't have any match but we're confident in it.
      dict(
          candidate=test_utils.make_variant(start=15, alleles=['C', 'A']),
          expected_confident=True,
          expected_label=0,
          expected_truth=test_utils.make_variant(
              start=15, alleles=['C', 'A'], gt=[0, 0])),

      # These variant start at our SNP but has a different allele. We are
      # confident and we get back the true snp variant.
      # However, we are on a different allele, so its status is unknown.
      # TODO: Confirm this case.
      dict(
          candidate=test_utils.make_variant(
              start=snp_class1.start, alleles=['A', 'G']),
          expected_confident=True,
          expected_label=0,
          expected_truth=snp_class1),
      # TODO: checking this assumption is correct:
      # If the alleles don't match, return class 0?
      dict(
          candidate=test_utils.make_variant(
              start=snp_class1.start, alleles=['AC', 'C']),
          expected_confident=True,
          expected_label=0,
          expected_truth=snp_class1),
      dict(
          candidate=test_utils.make_variant(
              start=snp_class1.start, alleles=['A', 'CA']),
          expected_confident=True,
          expected_label=0,
          expected_truth=snp_class1),
      # Checks that we don't match against the filtered truth variant in our
      # database. This means that we return not the filtered variant but one
      # with a (0, 0) genotype.
      dict(
          candidate=test_utils.make_variant(start=filtered.start),
          expected_confident=True,
          expected_label=0,
          expected_truth=test_utils.make_variant(
              start=filtered.start, gt=(0, 0))),
      # These variant start at our SNP but has a different first alt allele 'G'.
      # The second alt ('C') matches snp_class1, so we still got back the
      # expected_label.
      dict(
          candidate=test_utils.make_variant(
              start=snp_class1.start, alleles=['A', 'G', 'C']),
          expected_confident=True,
          expected_label=1,
          expected_truth=snp_class1,
          variant_alt_alleles_indices=[1]),
      # And, even if the variant_alt_alleles_indices is a composite one ([0,1]),
      # We still label it as long as one of them matches the truth_alt.
      dict(
          candidate=test_utils.make_variant(
              start=snp_class1.start, alleles=['A', 'G', 'C']),
          expected_confident=True,
          expected_label=1,
          expected_truth=snp_class1,
          variant_alt_alleles_indices=[0, 1]),
      # ... But we won't label it if the alt_allele_indices does not cover the
      # truth alt.
      dict(
          candidate=test_utils.make_variant(
              start=snp_class1.start, alleles=['A', 'G', 'C']),
          expected_confident=True,
          expected_label=0,
          expected_truth=snp_class1,
          variant_alt_alleles_indices=[0]),
  )
  def test_label_variants(self,
                          candidate,
                          expected_confident,
                          expected_truth,
                          expected_label=None,
                          variant_alt_alleles_indices=None):
    if variant_alt_alleles_indices is None:
      variant_alt_alleles_indices = [0]
    labeler = self._make_labeler(
        self.variants,
        ranges.RangeSet(
            [ranges.make_range(self.snp_class1.reference_name, 10, 100)]))

    # Call _match so we can compare our expected truth with the actual one.
    is_confident, truth_variant = labeler._match(candidate)
    self.assertEqual(expected_truth, truth_variant)
    self.assertEqual(is_confident, expected_confident)

    # Now call label_variants to exercise the higher-level API.
    classes_dict = (
        customized_classes_labeler.CustomizedClassesVariantLabel.classes_dict)
    if expected_label is None and expected_truth is not None:
      expected_class_str = expected_truth.info[
          customized_classes_labeler.CustomizedClassesVariantLabel
          .info_field_name].values[0].string_value
      expected_label = classes_dict[expected_class_str]

    labels = list(labeler.label_variants([candidate]))
    self.assertLen(labels, 1)
    self.assertEqual(candidate, labels[0].variant)
    self.assertEqual(expected_confident, labels[0].is_confident)
    self.assertEqual(
        expected_label,
        labels[0].label_for_alt_alleles(variant_alt_alleles_indices))

  def test_match_selects_variant_by_start(self):
    # Tests that match() selects the variant at the same start even if that
    # variant doesn't have the same alleles at candidate and there's an
    # overlapping with the same alleles.
    overlapping = [
        test_utils.make_variant(start=20, alleles=['CC', 'A'], gt=[1, 1]),
        test_utils.make_variant(start=21, alleles=['AAA', 'A'], gt=[0, 1]),
        test_utils.make_variant(start=22, alleles=['AA', 'A'], gt=[1, 1]),
    ]
    candidate = test_utils.make_variant(start=21, alleles=['CC', 'A'])

    labeler = self._make_labeler(
        overlapping,
        ranges.RangeSet(
            [ranges.make_range(overlapping[0].reference_name, 0, 100)]))
    is_confident, truth_variant = labeler._match(candidate)
    self.assertEqual(is_confident, True)
    self.assertEqual(truth_variant, overlapping[1])


if __name__ == '__main__':
  absltest.main()
