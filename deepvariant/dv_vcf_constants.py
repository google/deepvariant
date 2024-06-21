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
"""Library for generating VCF information created by DeepVariant."""

from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import vcf_constants

# Current release version of DeepVariant.
DEEP_VARIANT_VERSION = '1.6.1'

# FILTER field IDs.
DEEP_VARIANT_PASS = 'PASS'
DEEP_VARIANT_REF_FILTER = 'RefCall'
DEEP_VARIANT_QUAL_FILTER = 'LowQual'
DEEP_VARIANT_NO_CALL = 'NoCall'
DEEP_VARIANT_GERMLINE = 'GERMLINE'
DEEP_VARIANT_PON = 'PON'

# FORMAT field IDs.
DEEP_VARIANT_MIN_DP_FORMAT = 'MIN_DP'
DEEP_VARIANT_MED_DP_FORMAT = 'MED_DP'
DEEP_VARIANT_VAF_FORMAT = 'VAF'
DEEP_VARIANT_MODEL_ID_FORMAT = 'MID'

# Genotype codes:
UNCALLED_GENOTYPE = -1


def deepvariant_header(
    contigs,
    sample_names,
    add_info_candidates=False,
    include_med_dp=True,
    include_model_id=False,
):
  """Returns a VcfHeader used for writing VCF output.

  This function fills out the FILTER, INFO, FORMAT, and extra header information
  created by the DeepVariant pipeline using consistent fields that DeepVariant
  creates. The `contigs` and `sample_names` fields are unique depending on the
  input data used, so are required inputs.

  Args:
    contigs: list(ContigInfo). The list of contigs on which variants were
      called.
    sample_names: list(str). The list of samples present in the run.
    add_info_candidates: Adds the 'CANDIDATES' info field for debugging
      purposes.
    include_med_dp: boolean. If True, we will include MED_DP.
    include_model_id: boolean, If True, tag variants by model that called them.

  Returns:
    A nucleus.genomics.v1.VcfHeader proto with known fixed headers and the given
    samples and contigs populated.
  """
  version = variants_pb2.VcfExtra(
      key='DeepVariant_version', value=DEEP_VARIANT_VERSION
  )

  info_fields = [
      vcf_constants.reserved_info_field('END'),
  ]
  formats = [
      vcf_constants.reserved_format_field('GT'),
      vcf_constants.reserved_format_field('GQ'),
      vcf_constants.reserved_format_field('DP'),
      variants_pb2.VcfFormatInfo(
          id=DEEP_VARIANT_MIN_DP_FORMAT,
          number='1',
          type='Integer',
          description='Minimum DP observed within the GVCF block.',
      ),
      vcf_constants.reserved_format_field('AD'),
      variants_pb2.VcfFormatInfo(
          id=DEEP_VARIANT_VAF_FORMAT,
          number='A',
          type='Float',
          description='Variant allele fractions.',
      ),
      vcf_constants.reserved_format_field('PL'),
  ]
  if add_info_candidates:
    info_fields.append(
        variants_pb2.VcfInfo(
            id='CANDIDATES',
            number='1',
            type=vcf_constants.STRING_TYPE,
            description='pipe-delimited candidate alleles.',
        )
    )

  if include_med_dp:
    formats.append(
        variants_pb2.VcfFormatInfo(
            id=DEEP_VARIANT_MED_DP_FORMAT,
            number='1',
            type='Integer',
            description=(
                'Median DP observed within the GVCF block '
                'rounded to the nearest integer.'
            ),
        )
    )
  if include_model_id:
    formats.append(
        variants_pb2.VcfFormatInfo(
            id=DEEP_VARIANT_MODEL_ID_FORMAT,
            number='1',
            type=vcf_constants.STRING_TYPE,
            description='Identifies which model called this variant.',
        )
    )

  return variants_pb2.VcfHeader(
      fileformat='VCFv4.2',
      filters=[
          vcf_constants.reserved_filter_field(DEEP_VARIANT_PASS),
          variants_pb2.VcfFilterInfo(
              id=DEEP_VARIANT_REF_FILTER,
              description='Genotyping model thinks this site is reference.',
          ),
          variants_pb2.VcfFilterInfo(
              id=DEEP_VARIANT_QUAL_FILTER,
              description=(
                  'Confidence in this variant being real is below '
                  'calling threshold.'
              ),
          ),
          variants_pb2.VcfFilterInfo(
              id=DEEP_VARIANT_NO_CALL,
              description='Site has depth=0 resulting in no call.',
          ),
      ],
      infos=info_fields,
      formats=formats,
      contigs=contigs,
      sample_names=sample_names,
      extras=[version],
  )
