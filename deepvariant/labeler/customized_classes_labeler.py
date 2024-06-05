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
"""variant_labeler for DeepVariant."""



from deepvariant.labeler import positional_labeler
from deepvariant.labeler import variant_labeler
from third_party.nucleus.util import struct_utils


# ---------------------------------------------------------------------------
# CustomizedClassesVariantLabel
#
class CustomizedClassesVariantLabel(variant_labeler.VariantLabel):
  """Dataclass containing information about a label assigned to a variant.

  Attributes:
    is_confident: bool. True if we could confidently assign a label to this
      variant, False otherwise.
    variant: nucleus.protos.Variant proto that we assigned a label for.
    class_status: string. One of the keys in classes_dict
  """

  classes_dict = None
  info_field_name = None

  def __init__(
      self, is_confident, variant, truth_variant, classes_list, info_field_name
  ):
    self.info_field_name = info_field_name
    self.classes_dict = {k: v for v, k in enumerate(classes_list.split(','))}
    self.is_confident = is_confident
    self.variant = variant
    self.truth_variant = truth_variant

  def get_class(self):
    """Returns the class of the label."""
    try:
      return self.classes_dict[self.get_class_status(self.truth_variant.info)]
    except ValueError:
      return 0

  def label_for_alt_alleles(self, alt_alleles_indices):
    """Computes the label value for an example.

    This function computes the TensorFlow label value (0, 1, 2, .. N-1) we train
    DeepVariant to predict.
    The `alt_alleles_indices` being passed in is from the candidates (not
    truth), so they could still have multiple alts. If any of the alt alleles
    matches the truth, we'll return the label of the truth.
    TODO: Fix multi-allelic cases. Add corresponding unit test cases.
    Note that this function currently doesn't handle multi-allelic cases
    correctly. For example it assumes `truth_alt` is the first one.

    Args:
      alt_alleles_indices: list[int]. A list of the alt_allele_indices.

    Returns:
      int >= 0. Label for the classes in `classes_dict`.
    """
    if not self.truth_variant:
      return 0

    if self.truth_variant.calls:
      if self.truth_variant.calls[0].genotype == [0, 0]:
        return 0

    # If the ref of the candidate and the truth doesn't match, return 0 (ref).
    if self.truth_variant.reference_bases != self.variant.reference_bases:
      return 0
    true_class_status = self.get_class_status(self.truth_variant.info)
    truth_alt = self.truth_variant.alternate_bases[0]
    # Default is label 0. Usually reference.
    label = 0
    # Note that this logic below might not be the best when
    # `alt_alleles_indices` is a composite one, like [0, 1]. For now we'll
    # return the corresponding label if any of them matches truth_alt.
    for ind in alt_alleles_indices:
      if self.variant.alternate_bases[ind] == truth_alt:
        # allele in called variant is the same as truth_alt
        label = self.classes_dict[true_class_status]

    return label

  def get_class_status(self, info_field):
    """Extract class status from nucleus.protos.Variant.info.

    Args:
      info_field: INFO field of nucleus.protos.Variant proto to extract the
        classes status from. Must contain `info_field_name` field which is set
        to one of self.classes_dict.keys().

    Returns:
      string. Class status. Has to be one of the keys of `classes_dict`.

    Raises:
      ValueError: if type is missing in info_field
      ValueError: if type is not in self.classes_dict.keys()
    """

    if self.info_field_name not in info_field.keys():
      raise ValueError(
          'Cannot create class labels: '
          + 'VCF file does not contain INFO/{} field'.format(
              self.info_field_name
          )
      )

    class_status = struct_utils.get_string_field(
        info_field, self.info_field_name, True
    )

    if class_status not in self.classes_dict.keys():
      raise ValueError(
          'class_status status unknown: {}. Known status: {}'.format(
              class_status, self.classes_dict.keys()
          )
      )
    return class_status


# ---------------------------------------------------------------------------
# CustomizedClassesVariantLabeler
#
class CustomizedClassesVariantLabeler(
    positional_labeler.PositionalVariantLabeler
):
  """Extracts the class of the variant (possible values are keys in

  `classes_dict`) from INFO/`info_field_name` field in VCF file.
  """

  def __init__(
      self, truth_vcf_reader, confident_regions, classes_list, info_field_name
  ):
    """Creates a new CustomizedClassesVariantLabeler.

    Args:
      truth_vcf_reader: a VcfReader object that points to our truth variant set.
      confident_regions: A RangeSet containing all of the confidently called
        regions. A variant that falls outside of one of these regions will be
        receive a special not-confident marker.
      classes_list: A common-separated string of classes.
      info_field_name: the name in INFO field where we should get the customized
        field from.

    Raises:
      ValueError: if vcf_reader is None.
    """
    super(CustomizedClassesVariantLabeler, self).__init__(
        truth_vcf_reader=truth_vcf_reader, confident_regions=confident_regions
    )
    self.classes_list = classes_list
    self.info_field_name = info_field_name

  def label_variants(self, variants, region=None):
    """Gets label information for each variant in variants.

    This is the primary API for assigning labels to variants. This function
    takes and iterable of variants and yield a VariantLabel object for each
    variant. The VariantLabel can be used to determine the variant type label
    for each variant suitable for training a DeepVariant model. The API accepts
    an iterable of Variants because, in the general case, the labeling of
    variants aren't independent, in that the label assigned to one variant may
    impact the label we assign to a nearby variant.

    Args:
      variants: iterable[nucleus.protos.Variant]: An iterable of variants to
        label. The variants should be in coordinate-sorted order and all on the
        same chromosome.
      region: A nucleus.genomics.v1.Range object specifying the region over
        which we are labeling variants. This should span at least the span of
        variants, but may be larger. Statistics about the labeling will be
        computed over region.

    Yields:
      A VariantLabel object for each variant in variants, in order.
    """
    for variant in variants:
      is_confident, truth_variant = self._match(variant)

      yield CustomizedClassesVariantLabel(
          is_confident=is_confident,
          variant=variant,
          truth_variant=truth_variant,
          classes_list=self.classes_list,
          info_field_name=self.info_field_name,
      )
