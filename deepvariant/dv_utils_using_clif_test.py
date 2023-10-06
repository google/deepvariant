# Copyright 2023 Google LLC.
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
"""Tests for deepvariant.dv_utils_using_clif."""

from unittest import mock



from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from deepvariant import dv_utils
from deepvariant import dv_utils_using_clif
from third_party.nucleus.protos import variants_pb2


class DVUtilsUsingClifTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self.alts = ['A']
    self.variant = variants_pb2.Variant(
        reference_name='1',
        start=10,
        end=11,
        reference_bases='C',
        alternate_bases=self.alts,
    )
    self.encoded_image = b'encoded_image_data'
    self.default_shape = [5, 5, 7]

  def testMakeExample(self):
    example = dv_utils_using_clif.make_example(
        self.variant, self.alts, self.encoded_image, self.default_shape
    )

    self.assertEqual(
        self.encoded_image, dv_utils.example_encoded_image(example)
    )
    self.assertEqual(self.variant, dv_utils.example_variant(example))
    self.assertEqual(b'1:11-11', dv_utils.example_locus(example))
    self.assertEqual([0], dv_utils.example_alt_alleles_indices(example))
    self.assertEqual('1:11:C->A', dv_utils.example_key(example))
    self.assertEqual(
        dv_utils_using_clif.EncodedVariantType.SNP.value,
        dv_utils.example_variant_type(example),
    )
    pass

  def testMakeExampleMultiAllelic(self):
    alts = ['AA', 'CC', 'GG']
    self.variant.alternate_bases[:] = alts
    # Providing GG, AA checks that we're sorting the indices.
    example = dv_utils_using_clif.make_example(
        self.variant, ['GG', 'AA'], b'foo', self.default_shape
    )
    self.assertEqual([0, 2], dv_utils.example_alt_alleles_indices(example))
    self.assertEqual(['AA', 'GG'], dv_utils.example_alt_alleles(example))
    self.assertEqual('1:11:C->AA/GG', dv_utils.example_key(example))
    self.assertEqual(
        dv_utils_using_clif.EncodedVariantType.INDEL.value,
        dv_utils.example_variant_type(example),
    )

  def testAltAllelesWithVariant(self):
    alts = list(self.variant.alternate_bases)
    example = dv_utils_using_clif.make_example(
        self.variant, alts, b'foo', self.default_shape
    )
    self.assertEqual([0], dv_utils.example_alt_alleles_indices(example))
    with mock.patch(
        'deepvariant.dv_utils.example_variant'
    ) as mock_ex_variant:
      # Providing variant directly avoids the call to example_variant().
      self.assertEqual(
          alts, dv_utils.example_alt_alleles(example, variant=self.variant)
      )
      mock_ex_variant.assert_not_called()

      # Checks that we load the variant if needed and that our mock is working.
      mock_ex_variant.return_value = self.variant
      self.assertEqual(alts, dv_utils.example_alt_alleles(example))
      mock_ex_variant.assert_called_once_with(example)

  def assertIsNotAFeature(self, label, example):
    self.assertNotIn(label, example.features.feature)

  def testExampleSetLabel(self):
    example = dv_utils_using_clif.make_example(
        self.variant, self.alts, self.encoded_image, self.default_shape
    )

    self.assertIsNotAFeature('label', example)
    for label in [0, 1, 2]:
      dv_utils.example_set_label(example, label)
      self.assertEqual(label, dv_utils.example_label(example))

  def testExampleImageShape(self):
    example = dv_utils_using_clif.make_example(
        self.variant, self.alts, self.encoded_image, self.default_shape
    )
    self.assertEqual(self.default_shape, dv_utils.example_image_shape(example))


if __name__ == '__main__':
  tf.compat.v1.disable_eager_execution()
  absltest.main()
