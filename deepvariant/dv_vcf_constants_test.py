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
import itertools
from absl.testing import absltest
from absl.testing import parameterized
from deepvariant import dv_vcf_constants
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import variantcall_utils


class DvVcfConstantsTest(parameterized.TestCase):

  @parameterized.parameters(
      dict(contigs=[], sample_names=[]),
      dict(
          contigs=[reference_pb2.ContigInfo(name='chr1')],
          sample_names=['single_sample'],
      ),
      dict(
          contigs=[
              reference_pb2.ContigInfo(name='1'),
              reference_pb2.ContigInfo(name='2'),
          ],
          sample_names=['multiple', 'samples'],
      ),
  )
  def test_deepvariant_header(self, contigs, sample_names):
    header = dv_vcf_constants.deepvariant_header(
        contigs=contigs, sample_names=sample_names
    )
    self.assertCountEqual(header.contigs, contigs)
    self.assertCountEqual(header.sample_names, sample_names)
    self.assertNotEmpty(header.filters)
    self.assertNotEmpty(header.infos)
    self.assertNotEmpty(header.formats)

  def test_compute_filter_fields(self):
    # This generates too many tests as a parameterized test.
    for qual, min_qual in itertools.product(range(100), range(100)):
      # First test with no alleleic depth.
      variant = variants_pb2.Variant()
      variant.quality = qual
      expected = []
      expected.append(dv_vcf_constants.DEEP_VARIANT_NO_CALL)
      self.assertEqual(
          dv_vcf_constants.compute_filter_fields(variant, min_qual),
          expected,
      )
      # Now add hom ref genotype and AD --> qual shouldn't affect filter field
      del variant.filter[:]
      variant.calls.add(genotype=[0, 0])
      variantcall_utils.set_ad(variant.calls[0], [1, 1])
      expected = []
      expected.append(dv_vcf_constants.DEEP_VARIANT_REF_FILTER)
      self.assertEqual(
          dv_vcf_constants.compute_filter_fields(variant, min_qual),
          expected,
      )
      # Now add variant genotype --> qual filter should matter again
      del variant.filter[:]
      del variant.calls[:]
      variant.calls.add(genotype=[0, 1])
      variantcall_utils.set_ad(variant.calls[0], [1, 1])
      expected = []
      expected.append(
          dv_vcf_constants.DEEP_VARIANT_PASS
          if qual >= min_qual
          else dv_vcf_constants.DEEP_VARIANT_QUAL_FILTER
      )
      self.assertEqual(
          dv_vcf_constants.compute_filter_fields(variant, min_qual),
          expected,
      )


if __name__ == '__main__':
  absltest.main()
