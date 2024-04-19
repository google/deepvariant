# Copyright 2018 Google LLC.
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
"""VariantCall utilities."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from third_party.nucleus.protos import struct_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import struct_utils
from third_party.nucleus.util import vcf_constants

# Special-cased FORMAT fields that are first-class fields in the VariantCall.
_GL = 'GL'
_GT = 'GT'


def set_format(variant_call, field_name, value, vcf_object=None):
  """Sets a field of the info map of the `VariantCall` to the given value(s).

  `variant_call.info` is analogous to the FORMAT field of a VCF call.

  Example usage:
  with vcf.VcfReader('/path/to/my.vcf') as vcf_reader:
    for variant in vcf_reader:
      first_call = variant.calls[0]
      # Type can be inferred for reserved VCF fields.
      set_format(first_call, 'AD', 25)
      # Specify the reader explicitly for unknown fields.
      set_format(first_call, 'MYFIELD', 30, vcf_reader)

  Args:
    variant_call: VariantCall proto. The VariantCall to modify.
    field_name: str. The name of the field to set.
    value: A single value or list of values to update the VariantCall with.
      The type of the value is determined by the `vcf_object` if one is given,
      otherwise is looked up based on the reserved FORMAT fields in the VCF
      specification.
    vcf_object: (Optional) A VcfReader or VcfWriter object. If not None, the
      type of the field is inferred from the associated VcfReader or VcfWriter
      based on its name. Otherwise, the type is inferred if it is a reserved
      field.
  """
  if field_name == _GL:
    set_gl(variant_call, value)
    return
  if field_name == _GT:
    set_gt(variant_call, value)
    return

  if vcf_object is None:
    set_field_fn = vcf_constants.reserved_format_field_set_fn(field_name)
  else:
    set_field_fn = vcf_object.field_access_cache.format_field_set_fn(field_name)
  set_field_fn(variant_call.info, field_name, value)


def get_format(variant_call, field_name, vcf_object=None):
  """Returns the value of the `field_name` FORMAT field.

  The `vcf_object` is used to determine the type of the resulting value. If it
  is a single value or a Flag, that single value will be returned. Otherwise,
  the list of values is returned.

  Args:
    variant_call: VariantCall proto. The VariantCall of interest.
    field_name: str. The name of the field to retrieve values from.
    vcf_object: (Optional) A VcfReader or VcfWriter object. If not None, the
      type of the field is inferred from the associated VcfReader or VcfWriter
      based on its name. Otherwise, the type is inferred if it is a reserved
      field.
  """
  if field_name == _GL:
    return get_gl(variant_call)
  if field_name == _GT:
    return get_gt(variant_call)

  if vcf_object is None:
    get_field_fn = vcf_constants.reserved_format_field_get_fn(field_name)
  else:
    get_field_fn = vcf_object.field_access_cache.format_field_get_fn(field_name)
  return get_field_fn(variant_call.info, field_name)


# The following functions are convenience methods for getting/setting some
# reserved FORMAT fields of a VariantCall as well as some non-reserved FORMAT
# fields used by DeepVariant. Note that these functions will use the types of
# each field as defined by the VCF 4.3 specification, mirrored in
# vcf_constants.py, so if you have redefined any of these fields to have
# different types these functions will not do what you want.
def set_ad(variant_call, ad):
  """Sets the allele depth of the VariantCall."""
  set_format(variant_call, 'AD', ad)


def get_ad(variant_call):
  """Gets the allele depth of the VariantCall."""
  return get_format(variant_call, 'AD')


def set_gl(variant_call, gl):
  """Sets the genotype likelihoods of the VariantCall.

  Args:
    variant_call: VariantCall proto. The VariantCall to modify.
    gl: list(float). The list of genotype likelihoods for the VariantCall.
  """
  # Note: genotype_likelihood is extracted to a first-class field within
  # VariantCall. Consequently, we just set its value directly here.
  variant_call.genotype_likelihood[:] = gl


def get_gl(variant_call):
  """Returns the genotype likelihoods of the VariantCall.

  Args:
    variant_call: VariantCall proto. The VariantCall for which to return GLs.

  Returns:
    A list of floats representing the genotype likelihoods of this call.
  """
  return variant_call.genotype_likelihood


def set_gt(variant_call, gt):
  """Sets the genotypes of the VariantCall.

  Args:
    variant_call: VariantCall proto. The VariantCall to modify.
    gt: list(int). The list of genotypes for the VariantCall.
  """
  # Note: genotype is extracted to a first-class field within
  # VariantCall. Consequently, we just set its value directly here.
  variant_call.genotype[:] = gt


def get_gt(variant_call):
  """Returns the genotypes of the VariantCall.

  Args:
    variant_call: VariantCall proto. The VariantCall for which to return GTs.

  Returns:
    A list of ints representing the genotype indices of this call.
  """
  return variant_call.genotype


def set_gq(variant_call, gq):
  """Sets the genotype quality of the VariantCall."""
  set_format(variant_call, 'GQ', gq)


def get_gq(variant_call):
  """Gets the genotype quality of the VariantCall."""
  return get_format(variant_call, 'GQ')


def set_med_dp(variant_call, med_dp):
  """Sets the 'MED_DP' field of the VariantCall."""
  struct_utils.set_int_field(variant_call.info, 'MED_DP', med_dp)


def get_med_dp(variant_call):
  """Gets the 'MED_DP' field of the VariantCall."""
  return struct_utils.get_int_field(
      variant_call.info, 'MED_DP', is_single_field=True)


def set_min_dp(variant_call, min_dp):
  """Sets the 'MIN_DP' field of the VariantCall."""
  struct_utils.set_int_field(variant_call.info, 'MIN_DP', min_dp)


def set_model_id(variant_call, model_id):
  """Sets the 'MID' field of the VariantCall."""
  set_format(variant_call, 'MID', model_id)


def get_min_dp(variant_call):
  """Gets the 'MIN_DP' field of the VariantCall."""
  return struct_utils.get_int_field(
      variant_call.info, 'MIN_DP', is_single_field=True)


def set_bam_fname(variant_call, bam_fname):
  """Sets 'BAM_FNAME' field of the VariantCall."""
  return struct_utils.set_string_field(variant_call.info, 'BAM_FNAME',
                                       bam_fname)


def has_genotypes(variant_call):
  """Returns True iff the VariantCall has one or more called genotypes.

  Args:
    variant_call: VariantCall proto. The VariantCall to evaluate.

  Returns:
    True if the VariantCall has one or more called genotypes, False otherwise.
  """
  return any(gt >= 0 for gt in variant_call.genotype)


def has_full_genotypes(variant_call):
  """Returns True iff the VariantCall has only known genotypes.

  Args:
    variant_call: VariantCall proto. The VariantCall to evaluate.

  Returns:
    True if all `genotype` fields are known genotypes.
  """
  return all(gt >= 0 for gt in variant_call.genotype)


def ploidy(variant_call):
  """Returns the ploidy of the VariantCall.

  Args:
    variant_call: VariantCall proto. The VariantCall to evaluate.

  Returns:
    The ploidy of the call (a non-negative integer).
  """
  # Unknown genotypes are represented as -1 in VariantCall protos. When
  # a VCF is parsed that contains multiple ploidies in different samples,
  # a separate padding value of -2**30 - 1 is inserted into the calls.
  return sum(gt >= -1 for gt in variant_call.genotype)


def has_variation(variant_call):
  """Returns True if and only if the call has a non-reference genotype.

  Args:
    variant_call: VariantCall proto. The VariantCall to evaluate.

  Returns:
    True if and only if the call has a non-reference genotype.
  """
  return any(gt > 0 for gt in variant_call.genotype)


def is_heterozygous(variant_call):
  """Returns True if and only if the call is heterozygous.

  Args:
    variant_call: VariantCall proto. The VariantCall to evaluate.

  Returns:
    True if and only if the call is heterozygous.
  """
  return len({gt for gt in variant_call.genotype if gt >= 0}) >= 2
