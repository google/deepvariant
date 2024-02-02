# Copyright 2020 Google LLC.
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
"""Tests for deepvariant .show_examples."""

import glob
import os

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from deepvariant import show_examples
from deepvariant import testdata
from third_party.nucleus.testing import test_utils

FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


class ShowExamplesTest(parameterized.TestCase):

  def test_create_region_filter(self):
    region_flag = 'chr1:10-20'
    variants_inside = [test_utils.make_variant(start=15, alleles=['A', 'G'])]
    variants_outside = [
        test_utils.make_variant(start=4, alleles=['A', 'G']),
        test_utils.make_variant(chrom='chr2', start=4, alleles=['A', 'G']),
    ]

    region_filter_function = show_examples.create_region_filter(
        region_flag_string=region_flag
    )
    for v in variants_inside:
      self.assertTrue(
          region_filter_function(v),
          msg='Variant at {} should pass filter {}'.format(
              show_examples.get_short_id(v, [0]), region_flag
          ),
      )
    for v in variants_outside:
      self.assertFalse(
          region_filter_function(v),
          msg='Variant at {} should NOT pass filter {}'.format(
              show_examples.get_short_id(v, [0]), region_flag
          ),
      )

  @parameterized.parameters([
      # The first allele is the ref, the rest are alts.
      dict(alleles=['A', 'G'], indices=[0], expected='chr1:10_A->G'),
      dict(alleles=['A', 'AA', 'AC'], indices=[0], expected='chr1:10_A->AA'),
      dict(
          alleles=['A', 'AA', 'AC'], indices=[0, 1], expected='chr1:10_A->AA|AC'
      ),
  ])
  def test_get_full_id(self, alleles, indices, expected):
    variant = test_utils.make_variant(start=10, alleles=alleles)
    output = show_examples.get_full_id(variant, indices)
    self.assertEqual(expected, output)

  @parameterized.parameters([
      # The first allele is the ref, the rest are alts.
      dict(alleles=['A', 'G'], indices=[0], expected='chr1:10_A->G'),
      dict(alleles=['A', 'AA', 'AC'], indices=[0], expected='chr1:10_A->AA'),
      dict(
          alleles=['A', 'AA', 'AC'], indices=[0, 1], expected='chr1:10_A->AA|AC'
      ),
  ])
  def test_get_short_id_for_small_variants(self, alleles, indices, expected):
    variant = test_utils.make_variant(start=10, alleles=alleles)
    output = show_examples.get_short_id(variant, indices)
    self.assertEqual(expected, output)
    full_id = show_examples.get_full_id(variant, indices)
    self.assertEqual(output, full_id)

  @parameterized.parameters([
      # The first allele is the ref, the rest are alts.
      dict(alleles=['ACGTACGT', 'A'], indices=[0], expected='chr1:10_DEL7bp'),
      dict(alleles=['A', 'ACGTACGT'], indices=[0], expected='chr1:10_INS7bp'),
      dict(
          alleles=['A', 'ACGTACGT', 'ACGTACGTACGTACGT'],
          indices=[0, 1],
          expected='chr1:10_INS7bp|INS15bp',
      ),
      dict(
          alleles=['A', 'ACGTACGT', 'AAAAAAAA'],
          indices=[0],
          expected='chr1:10_alt0INS7bp',
      ),
      dict(
          alleles=['A', 'ACGTACGT', 'AAAAAAAA'],
          indices=[1],
          expected='chr1:10_alt1INS7bp',
      ),
      dict(
          alleles=['A', 'ACGTACGT', 'AAAAAAAA'],
          indices=[0, 1],
          expected='chr1:10_alt0INS7bp|alt1INS7bp',
      ),
  ])
  def test_get_short_id_for_longer_variants(self, alleles, indices, expected):
    variant = test_utils.make_variant(start=10, alleles=alleles)
    output = show_examples.get_short_id(variant, indices)
    self.assertEqual(expected, output)


class ShowExamplesEnd2EndTest(absltest.TestCase):

  @flagsaver.flagsaver
  def test_show_examples_end2end_calling_examples(self):
    output_prefix = test_utils.test_tmpfile('calling')
    FLAGS.examples = testdata.GOLDEN_CALLING_EXAMPLES
    FLAGS.output = output_prefix
    show_examples.run()
    ls = glob.glob('{}*'.format(output_prefix))
    filenames = [os.path.basename(path) for path in ls]
    self.assertTrue(all(['calling' in filename for filename in filenames]))
    self.assertTrue(all([filename.endswith('.png') for filename in filenames]))
    self.assertFalse(
        any(['label' in filename for filename in filenames]),
        msg='Calling examples should NOT produce labeled images.',
    )

  @flagsaver.flagsaver
  def test_show_examples_end2end_training_examples(self):
    output_prefix = test_utils.test_tmpfile('training')
    FLAGS.examples = testdata.GOLDEN_TRAINING_EXAMPLES
    FLAGS.output = output_prefix
    show_examples.run()
    ls = glob.glob('{}*'.format(output_prefix))
    filenames = [os.path.basename(path) for path in ls]
    self.assertTrue(all(['training' in filename for filename in filenames]))
    self.assertTrue(all([filename.endswith('.png') for filename in filenames]))
    self.assertTrue(
        all(['label' in filename for filename in filenames]),
        msg='Training examples should produce labeled images.',
    )

  @flagsaver.flagsaver
  def test_show_examples_end2end_all_optional_parameters(self):
    # Set all the optional parameters to check that they all work together.
    output_prefix = test_utils.test_tmpfile('kitchen_sink')
    FLAGS.examples = testdata.GOLDEN_TRAINING_EXAMPLES
    FLAGS.example_info_json = FLAGS.examples + '.example_info.json'
    FLAGS.output = output_prefix
    FLAGS.annotation = False
    FLAGS.regions = 'chr20:10,003,650-10,005,000'
    FLAGS.vcf = testdata.TRUTH_VARIANTS_VCF
    FLAGS.image_type = 'both'
    FLAGS.num_records = 5
    FLAGS.verbose = False
    FLAGS.truth_labels = False  # On by default for training examples.

    show_examples.run()
    ls = glob.glob('{}*'.format(output_prefix))
    filenames = [os.path.basename(path) for path in ls]

    self.assertTrue(
        any(['kitchen_sink_chr20' in filename for filename in filenames]),
        msg='image_type=both, so there should be images without "rgb"',
    )
    self.assertTrue(
        any(['rgb' in filename for filename in filenames]),
        msg='image_type=both, so there should be images with "rgb"',
    )
    self.assertTrue(all([filename.endswith('.png') for filename in filenames]))
    self.assertFalse(
        any(['label' in filename for filename in filenames]),
        msg='Should be no "label" when truth_labels=False',
    )
    self.assertLen(
        filenames,
        10,
        msg='Should be 10 filenames, i.e. 5 records with channels+rgb for each',
    )
    # Despite "Count", this checks that the items are the same, unordered.
    self.assertCountEqual(
        filenames,
        [
            'kitchen_sink_chr20:10004146_A->G.png',
            'kitchen_sink_chr20:10004146_A->G.rgb.png',
            'kitchen_sink_chr20:10004093_A->C.png',
            'kitchen_sink_chr20:10004093_A->C.rgb.png',
            'kitchen_sink_chr20:10003831_G->A.png',
            'kitchen_sink_chr20:10003831_G->A.rgb.png',
            'kitchen_sink_chr20:10003691_A->G.png',
            'kitchen_sink_chr20:10003691_A->G.rgb.png',
            'kitchen_sink_chr20:10003650_T->C.png',
            'kitchen_sink_chr20:10003650_T->C.rgb.png',
        ],
        msg=(
            'Specific examples and their output filenames should be the same '
            'if the inputs are the same.'
        ),
    )

  @flagsaver.flagsaver
  def test_show_examples_raises_on_wrong_column_labels(self):
    output_prefix = test_utils.test_tmpfile('column_labels')
    FLAGS.examples = testdata.GOLDEN_TRAINING_EXAMPLES
    FLAGS.output = output_prefix
    FLAGS.column_labels = 'read base,base quality,mapping quality,strand'
    with self.assertRaisesRegex(ValueError, '--column_labels'):
      show_examples.run()

    # With 6 channel names, it should run without error:
    FLAGS.column_labels = '1,2,3,4,5,6,19'
    show_examples.run()


if __name__ == '__main__':
  absltest.main()
