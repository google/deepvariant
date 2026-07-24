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
import gzip

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from deepvariant import testdata
from deepvariant.labeler import labeled_examples_to_vcf
from third_party.nucleus.io import vcf
from third_party.nucleus.testing import test_utils

FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


class ExamplesToVCFUnitTest(parameterized.TestCase):

  @flagsaver.flagsaver
  def test_end2end(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.examples = testdata.GOLDEN_TRAINING_EXAMPLES + '@3'  # Sharded.
    FLAGS.output_vcf = test_utils.test_tmpfile('examples_to_vcf.vcf.gz')

    labeled_examples_to_vcf.main(0)

    with gzip.open(FLAGS.output_vcf, 'rt') as f:
      vcf_lines = f.readlines()
    with open(
        testdata.deepvariant_testdata('golden.training_examples.vcf')
    ) as f:
      golden_lines = f.readlines()
    self.assertEqual(vcf_lines, golden_lines)

  @flagsaver.flagsaver
  def test_sample_name_flag(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.examples = testdata.GOLDEN_TRAINING_EXAMPLES
    FLAGS.sample_name = 'sample_name'
    FLAGS.output_vcf = test_utils.test_tmpfile('no_sample_name.vcf.gz')

    labeled_examples_to_vcf.main(0)

    with vcf.VcfReader(FLAGS.output_vcf) as vcf_reader:
      self.assertEqual(
          list(vcf_reader.header.sample_names), [FLAGS.sample_name]
      )

  @flagsaver.flagsaver
  def test_raises_for_unlabeled_examples(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.examples = testdata.GOLDEN_CALLING_EXAMPLES
    FLAGS.output_vcf = test_utils.test_tmpfile('unlabeled.vcf')

    with self.assertRaisesRegex(
        ValueError,
        (
            'Variant .* does not have any genotypes. This tool only works with '
            'variants that have been labeled'
        ),
    ):
      labeled_examples_to_vcf.main(0)

  @parameterized.parameters(
      # Case 1: Simple bi-allelic, no pruning, no simplification
      dict(
          alleles=['A', 'C'],
          gt=[0, 1],
          ad=[10, 15],
          expected_alleles=['A', 'C'],
          expected_gt=[0, 1],
          expected_ad=[10, 15],
      ),
      # Case 2: Multi-allelic, prune Alt1, simplify
      # Ref: CAA, Alt1: C, Alt2: CA
      # GT: (0, 2) i.e. Ref and Alt2 (CA)
      # Alt1 (C) is pruned. Remaining: CAA, CA
      # Simplifies to CA, C.
      # GT becomes (0, 1).
      dict(
          alleles=['CAA', 'C', 'CA'],
          gt=[0, 2],
          ad=[10, 5, 8],
          expected_alleles=['CA', 'C'],
          expected_gt=[0, 1],
          expected_ad=[10, 8],
      ),
      # Case 3: Multi-allelic, hom-ref (all pruned except first)
      # Ref: CAA, Alt1: C, Alt2: CA
      # GT: (0, 0)
      # Keep first Alt1 (C). Alt2 (CA) pruned.
      # Remaining: CAA, C.
      # Cannot simplify CAA, C.
      dict(
          alleles=['CAA', 'C', 'CA'],
          gt=[0, 0],
          ad=[10, 5, 8],
          expected_alleles=['CAA', 'C'],
          expected_gt=[0, 0],
          expected_ad=[10, 5],
      ),
      # Case 4: Multi-allelic, no pruning (all active), simplify
      # Ref: ATT, Alt1: TTT, Alt2: CTT
      # GT: (1, 2)
      # Keep both. Simplifies to A, T, C.
      dict(
          alleles=['ATT', 'TTT', 'CTT'],
          gt=[1, 2],
          ad=[10, 5, 8],
          expected_alleles=['A', 'T', 'C'],
          expected_gt=[1, 2],
          expected_ad=[10, 5, 8],
      ),
  )
  def test_prune_and_simplify_variant(
      self, alleles, gt, ad, expected_alleles, expected_gt, expected_ad
  ):
    variant = test_utils.make_variant(
        chrom='chr1', start=10, alleles=alleles, gt=gt, ad=ad
    )
    result = labeled_examples_to_vcf.prune_and_simplify_variant(variant)

    self.assertEqual(result.reference_bases, expected_alleles[0])
    self.assertEqual(list(result.alternate_bases), expected_alleles[1:])

    call = result.calls[0]
    self.assertEqual(list(call.genotype), expected_gt)

    # Check AD remapping
    self.assertEqual([v.int_value for v in call.info['AD'].values], expected_ad)
    # Check DP remapping (should be sum of expected_ad)
    self.assertEqual(
        [v.int_value for v in call.info['DP'].values], [sum(expected_ad)]
    )


if __name__ == '__main__':
  absltest.main()
