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
"""Postprocess output from call_variants to produce a VCF file."""
# TODO: Add type annotations to this module
import collections
import functools
import itertools
import os
import tempfile
import time
from typing import Iterable, Iterator, Sequence

from absl import flags
from absl import logging
import numpy as np
import pysam
import tensorflow as tf

from deepvariant import calling_regions_utils
from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import dv_vcf_constants
from deepvariant import haplotypes
from deepvariant import logging_level
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import postprocess_variants as postprocess_variants_lib
from absl import app
import multiprocessing
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tabix
from third_party.nucleus.io import tfrecord
from third_party.nucleus.io import vcf
from third_party.nucleus.io.python import merge_variants
from third_party.nucleus.io.python import vcf_concat
from third_party.nucleus.protos import range_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.util import errors
from third_party.nucleus.util import genomics_math
from third_party.nucleus.util import proto_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils
from third_party.nucleus.util.struct_utils import add_string_field


_INFILE = flags.DEFINE_string(
    'infile',
    None,
    (
        'Required. Path(s) to CallVariantOutput protos in TFRecord format to '
        'postprocess. These should be the complete set of outputs for '
        'call_variants.py.'
    ),
)
_OUTFILE = flags.DEFINE_string(
    'outfile',
    None,
    (
        'Required. Destination path where we will write output variant calls in'
        ' VCF format.'
    ),
)
_REF = flags.DEFINE_string(
    'ref',
    None,
    (
        'Required. Genome reference in FAI-indexed FASTA format. Used to'
        ' determine the sort order for the emitted variants and the VCF header.'
    ),
)
_SMALL_MODEL_CVO_RECORDS = flags.DEFINE_string(
    'small_model_cvo_records',
    None,
    (
        'Optional. Path(s) to CallVariantOutput protos in TFRecord format that'
        ' were called by the small model to include in postprocess .'
    ),
)
_QUAL_FILTER = flags.DEFINE_float(
    'qual_filter',
    1.0,
    'Any variant with QUAL < qual_filter will be filtered in the VCF file.',
)
_CNN_HOMREF_CALL_MIN_GQ = flags.DEFINE_float(
    'cnn_homref_call_min_gq',
    20.0,
    (
        'All CNN RefCalls whose GQ is less than this value will have ./.'
        ' genotype instead of 0/0.'
    ),
)
_MULT_ALLELIC_QUAL_FILTER = flags.DEFINE_float(
    'multi_allelic_qual_filter',
    1.0,
    'The qual value below which to filter multi-allelic variants.',
)
_NONVARIANT_SITE_TFRECORD_PATH = flags.DEFINE_string(
    'nonvariant_site_tfrecord_path',
    None,
    (
        'Optional. Path(s) to the non-variant sites protos in TFRecord format'
        ' to convert to gVCF file. This should be the complete set of outputs'
        ' from the --gvcf flag of make_examples.py.'
    ),
)
_GVCF_OUTFILE = flags.DEFINE_string(
    'gvcf_outfile',
    None,
    'Optional. Destination path where we will write the Genomic VCF output.',
)
_GROUP_VARIANTS = flags.DEFINE_boolean(
    'group_variants',
    True,
    (
        'If using vcf_candidate_importer and multi-allelic '
        'sites are split across multiple lines in VCF, set to False so that '
        'variants are not grouped when transforming CallVariantsOutput to '
        'Variants.'
    ),
)
_VCF_STATS_REPORT = flags.DEFINE_boolean(
    'vcf_stats_report',
    False,
    'Deprecated. Use vcf_stats_report.py instead.',
)
_SAMPLE_NAME = flags.DEFINE_string(
    'sample_name',
    None,
    (
        'Optional. If set, this will only be used if the sample name cannot be '
        'determined from the CallVariantsOutput or non-variant sites protos.'
    ),
)
_USE_MULTIALLELIC_MODEL = flags.DEFINE_boolean(
    'use_multiallelic_model',
    False,
    (
        'If True, use a specialized model for genotype resolution of'
        ' multiallelic cases with two alts.'
    ),
)
_DEBUG_OUTPUT_ALL_CANDIDATES = flags.DEFINE_enum(
    'debug_output_all_candidates',
    None,
    ['ALT', 'INFO'],
    (
        'Outputs all candidates considered by DeepVariant as additional ALT'
        ' alleles  or as an INFO field. For ALT, filtered candidates are'
        ' assigned a GL=0 and added as ALTs alleles, but do not appear in any'
        ' sample genotypes. This flag is useful for debugging purposes.'
        ' ALT-mode is incompatible with the multiallelic caller.'
    ),
)
_ONLY_KEEP_PASS = flags.DEFINE_boolean(
    'only_keep_pass', False, 'If True, only keep PASS calls.'
)
_HAPLOID_CONTIGS = flags.DEFINE_string(
    'haploid_contigs',
    None,
    (
        'Optional list of non autosomal chromosomes. For all listed chromosomes'
        'HET probabilities are not considered. The list can be either comma '
        'or space-separated.'
    ),
)
_CPUS = flags.DEFINE_integer(
    'cpus',
    multiprocessing.cpu_count(),
    'Number of worker processes to use. Set --cpus < 2 to disable parallel'
    ' processing.',
    short_name='j',
    required=False,
)
_NUM_PARTITIONS = flags.DEFINE_integer(
    'num_partitions',
    0,
    'Number of partitions to use for parallel or sequential processing. --Set'
    ' --num_partitions > --cpus to trade runtime for lower memory usage. Set'
    ' --num_partitions < 2 and --cpus < 2 to disable partitioning.',
    required=False,
)
_PAR_REGIONS = flags.DEFINE_string(
    'par_regions_bed',
    None,
    (
        'Optional BED file containing Human Pseudoautosomal Region (PAR) '
        'regions.'
        'Variants within this region are unaffected by genotype reallocation '
        'applied on regions supplied by --haploid_contigs flag.'
    ),
)
_REGIONS = flags.DEFINE_string(
    'regions',
    '',
    (
        'Optional. Space-separated list of regions we want to process. Elements'
        ' can be region literals (e.g., chr20:10-20) or paths to BED/BEDPE'
        ' files. This should match the flag passed to make_examples.py.'
    ),
)
_PROCESS_SOMATIC = flags.DEFINE_boolean(
    'process_somatic',
    False,
    'Optional. If specified the input is treated as somatic.',
)
_PON_FILTERING = flags.DEFINE_string(
    'pon_filtering',
    None,
    (
        'Optional. Only used if --process_somatic is true. '
        'A VCF file with Panel of Normals (PON) data.'
        'If set, the output VCF will be filtered: any variants that appear in '
        'PON will be marked with a PON filter, and PASS filter value will be '
        'removed.'
    ),
)

# Some format fields are indexed by alt allele, such as AD (depth by allele).
# These need to be cleaned up if we remove any alt alleles. Any info field
# listed here will be have its values cleaned up if we've removed any alt
# alleles.
# Each tuple contains: field name, ref_is_zero.
_ALT_ALLELE_INDEXED_FORMAT_FIELDS = frozenset(
    [('AD', True), ('VAF', False), ('MF', True), ('MD', True)]
)

# The number of places past the decimal point to round QUAL estimates to.
_QUAL_PRECISION = 7

# When this was set, it's about 20 seconds per log.
_LOG_EVERY_N = 100000

# When outputting all alt alleles, use placeholder value to indicate genotype
# will be soft-filtered.
_FILTERED_ALT_PROB = -9.0


def _extract_single_sample_name(
    record: deepvariant_pb2.CallVariantsOutput,
) -> str:
  """Returns the name of the single sample within the CallVariantsOutput file.

  Args:
    record: A deepvariant_pb2.CallVariantsOutput record.

  Returns:
    The name of the single individual in the first proto in the file.

  Raises:
    ValueError: There is not exactly one VariantCall in the proto or the
        call_set_name of the VariantCall is not populated.
  """
  variant = record.variant
  call = variant_utils.only_call(variant)
  name = call.call_set_name
  if not name:
    raise ValueError(
        'Error extracting name: no call_set_name set: {}'.format(record)
    )

  return name


def compute_filter_fields(
    variant: variants_pb2.Variant, min_quality: float
) -> list[str]:
  """Computes the filter fields for this variant.

  Variant filters are generated based on its quality score value and particular
  genotype call.

  Args:
    variant: Variant to filter.
    min_quality: Minimum acceptable phred scaled variant detection probability.

  Returns:
    Filter field strings to be added to the variant.
  """
  if variant_utils.genotype_type(variant) == variant_utils.GenotypeType.no_call:
    return [dv_vcf_constants.DEEP_VARIANT_NO_CALL]
  if variant_utils.genotype_type(variant) == variant_utils.GenotypeType.hom_ref:
    return [dv_vcf_constants.DEEP_VARIANT_REF_FILTER]
  elif variant.quality < min_quality:
    return [dv_vcf_constants.DEEP_VARIANT_QUAL_FILTER]
  else:
    return [dv_vcf_constants.DEEP_VARIANT_PASS]


def _pysam_resolve_file_path(file_path: str) -> str:
  """Prepends a prefix to the file_path when accessing Google files.

  Args:
    file_path: str. Full path pointing a specific file to access with pysam.

  Returns:
    str. The full configured file path for pysam to open.
  """
  # BEGN_INTERNAL
  if (
      file_path.startswith('/cns/')
      or file_path.startswith('/placer/')
      or file_path.startswith('/readahead/')
      or file_path.startswith('/bigstore/')
  ):
    return f'google:{file_path}'
  # END_INTERNAL
  return file_path


def most_likely_genotype(
    predictions: Sequence[float], ploidy: int = 2, n_alleles: int = 2
) -> tuple[int, list[int]]:
  """Gets the most likely genotype from predictions.

  From https://samtools.github.io/hts-specs/VCFv4.3.pdf:

  Genotype Ordering. In general case of ploidy P and N alternate alleles (0 is
  the REF and 1..N the alternate alleles), the ordering of genotypes for the
  likelihoods can be expressed by the following pseudocode with as many nested
  loops as ploidy:

  * Note that we use inclusive for loop boundaries.
  for a_P = 0 . . . N
    for a_P-1 = 0 . . . aP
      . . .
      for a_1 = 0 . . . a2
        println a1 a2 . . . aP

  Alternatively, the same can be achieved recursively with the following
  pseudocode:

    Ordering (P , N , suffix =""):
      for a in 0 . . . N
        if (P == 1) println str (a) + suffix
        if (P > 1) Ordering (P -1 , a, str (a) + suffix)

  Examples:
  * for P=2 and N=1, the ordering is 00,01,11
  * for P=2 and N=2, the ordering is 00,01,11,02,12,22
  * for P=3 and N=2, the ordering is 000,001,011,111,002,012,112,022,122,222
  * for P=1, the index of the genotype a is a
  * for P=2, the index of the genotype "a/b", where a <= b, is b(b + 1)/2 + a
  * for P=2 and arbitrary N, the ordering can be easily derived from a
    triangular matrix:
      b / a 0 1 2 3
      0     0
      1     1 2
      2     3 4 5
      3     6 7 8 9

  Args:
    predictions: N element array-like. The real-space probabilities of each
      genotype state for this variant. The number of elements in predictions is
      related to ploidy and n_alleles is given by N = choose(ploidy + n_alleles
      - 1, n_alleles -1) for more information see:
      http://genome.sph.umich.edu/wiki/Relationship_between_Ploidy,_Alleles_and_Genotypes
    ploidy: int >= 1. The ploidy (e.g., number of chromosomes) of this sample.
    n_alleles: int >= 2. The number of alleles (ref + n_alts).

  Returns:
    Two values. The first is the index of the most likely prediction in
    predictions. The second is a list of P elements with the VCF-style genotype
    indices corresponding to this index. For example, with P = 2 and an index of
    1, this returns the value (1, [0, 1]).

  Raises:
    NotImplementedError: if ploidy != 2 as this not yet implemented.
    ValueError: If n_alleles < 2.
    ValueError: If we cannot determine the genotype given prediction, n_alts,
      and ploidy.
  """
  # TODO: This can be memoized for efficiency.
  if ploidy != 2:
    raise NotImplementedError('Ploidy != 2 not yet implemented.')
  if n_alleles < 2:
    raise ValueError('n_alleles must be >= 2 but got', n_alleles)
  # TODO: would be nice to add test that predictions has the right
  # number of elements. But that would involve calculating the binomial
  # coefficient of n_alleles and ploidy, which would be expensive. Probably
  # need to memoize the whole function if we are going to add this.
  index_of_max = np.argmax(predictions)
  # This is the general case solution for fixed ploidy of 2 and arbitrary
  # n_alleles. We should generalize this code to the arbitrary ploidy case when
  # needed and memoize the mapping here.
  index = 0
  for h1 in range(0, n_alleles + 1):
    for h2 in range(0, h1 + 1):
      if index == index_of_max:
        return index, [h2, h1]
      index += 1
  raise ValueError('No corresponding GenotypeType for predictions', predictions)


def uncall_gt_if_no_ad(variant: variants_pb2.Variant) -> None:
  """Converts genotype to "./." if sum(AD)=0."""
  vcall = variant_utils.only_call(variant)
  if sum(variantcall_utils.get_ad(vcall)) == 0:
    # Set GT to ./.; GLs set to 0; GQ=0
    vcall.genotype[:] = [-1, -1]
    vcall.genotype_likelihood[:] = [0, 0]
    variantcall_utils.set_gq(vcall, 0)


def uncall_homref_gt_if_lowqual(
    variant: variants_pb2.Variant, min_homref_gq: float
) -> None:
  """Converts genotype to "./." if variant is CNN RefCall and has low GQ.

  If the variant has "RefCall" filter (which means an example was created for
  this site but CNN didn't call this as variant) and if the GQ is less than
  the given min_homref_gq threshold, set the genotype of the variant proto
  to "./.". See http://internal for more info.

  Args:
    variant: third_party.nucleus.protos.Variant proto.
    min_homref_gq: float.
  """
  vcall = variant_utils.only_call(variant)
  if (
      variant.filter == [dv_vcf_constants.DEEP_VARIANT_REF_FILTER]
      and variantcall_utils.get_gq(vcall) < min_homref_gq
  ):
    vcall.genotype[:] = [-1, -1]


def maybe_phase_genotype(
    variant: variants_pb2.Variant,
    genotype: list[int],
) -> tuple[bool, list[int]]:
  """Phases the genotype if phase information is available.

  The `ALT_PS` field contains phases assigned to each allele in range [0..2],
  for HP tags 0,1,2. The length of this array is number of alleles + 1, since
  REF is added implicitly.

  For example, [1,2] means that REF allele is assigned phase 1 and ALT_1 allele
  is assigned a phase 2. [2,1] means that REF allele is assigned phase 2 and
  ALT_1 allele is assigned phase 1. [2,2,1,1] means REF is phase 2, ALT_1 is
  phase 2, ALT_2 is phase 1, and ALT_3 is phase 1, etc.

  Args:
    variant: third_party.nucleus.protos.Variant proto.
    genotype: list of ints. The genotype indices to be written to the VCF.

  Returns:
    is_phased: bool. Whether it was possible to phase the genotype.
    genotype: genotype in the correct phasing order, if phased.
  """
  if not (
      variant_utils.get_info(variant, 'PS_CONTIG')
      and variant_utils.get_info(variant, 'ALT_PS')
  ):
    return False, genotype
  phase_info = [p.int_value for p in variant.info['ALT_PS'].values]
  if max(genotype) >= len(phase_info):
    logging.warning(
        (
            'Genotype %s is out of range for phase info %s for variant %s. '
            'Phasing was not applied.'
        ),
        genotype,
        phase_info,
        variant,
    )
    return False, genotype

  allele_1_haplotype = phase_info[genotype[0]]
  allele_2_haplotype = phase_info[genotype[1]]
  is_phased = (
      0 not in (allele_1_haplotype, allele_2_haplotype)
      and allele_1_haplotype != allele_2_haplotype
  )
  if is_phased:
    genotype = [
        genotype[allele_1_haplotype - 1],
        genotype[allele_2_haplotype - 1],
    ]
  return is_phased, genotype


def add_call_to_variant(
    variant: variants_pb2.Variant,
    predictions: Sequence[float],
    qual_filter: float = 0,
    sample_name: str | None = None,
) -> variants_pb2.Variant:
  """Fills in Variant record using the prediction probabilities.

  This functions sets the call[0].genotype, call[0].info['GQ'],
  call[0].genotype_probabilities, variant.filter, and variant.quality fields of
  variant based on the genotype likelihoods in predictions.

  Args:
    variant: third_party.nucleus.protos.Variant protobuf to be filled in with
      info derived from predictions.
    predictions: N element array-like. The real-space probabilities of each
      genotype state for this variant.
    qual_filter: float. If predictions implies that this isn't a reference call
      and the QUAL of the prediction isn't larger than qual_filter variant will
      be marked as FILTERed.
    sample_name: str. The name of the sample to assign to the Variant proto
      call_set_name field.

  Returns:
    A Variant record.

  Raises:
    ValueError: If variant doesn't have exactly one variant.call record.
  """
  call = variant_utils.only_call(variant)
  n_alleles = len(variant.alternate_bases) + 1
  index, genotype = most_likely_genotype(predictions, n_alleles=n_alleles)
  gq, variant.quality = compute_quals(predictions, index)
  call.call_set_name = sample_name
  call.is_phased, genotype = maybe_phase_genotype(variant, genotype)
  if call.is_phased:
    variantcall_utils.set_ps(call, variant_utils.get_info(variant, 'PS_CONTIG'))
  if is_methylated(call):
    mf = variantcall_utils.get_mf(call)
    mt = variantcall_utils.determine_methylation_type(mf)
    variantcall_utils.set_mt(call, mt)
  variantcall_utils.set_gt(call, genotype)
  variantcall_utils.set_gq(call, gq)
  gls = [genomics_math.perror_to_bounded_log10_perror(gp) for gp in predictions]
  variantcall_utils.set_gl(call, gls)
  uncall_gt_if_no_ad(variant)
  variant.filter[:] = compute_filter_fields(variant, qual_filter)
  uncall_homref_gt_if_lowqual(variant, _CNN_HOMREF_CALL_MIN_GQ.value)
  return variant


def compute_quals(
    predictions: Sequence[float], prediction_index: int
) -> tuple[int, int]:
  """Computes GQ and QUAL values from a set of prediction probabilities.

  Prediction probabilities are represented as a probability distribution over
  the N genotype states (e.g., for 3 genotype states {HOM_REF, HET, HOM_VAR}).
  Genotype Quality (or GQ) represents the PHRED scaled confidence in the
  particular genotype assignment. Likewise the QUAL representes the PHRED scaled
  confidence in variant as compared to reference, that is, P(NON_REF) / P(ALL)
  which in the diploid genotype case is P(HET) + P(HOM_VAR) / P(ALL). These
  quality scores are capped by _MAX_CONFIDENCE.

  Args:
    predictions: N element array-like. The real-space probabilities of each
      genotype state for this variant.
    prediction_index: int. The actual called genotype from the distribution.

  Returns:
    GQ and QUAL values for output in a Variant record.
  """
  # GQ is prob(genotype) / prob(all genotypes)
  # GQ is rounded to the nearest integer to comply with the VCF spec.
  gq = int(
      np.around(
          genomics_math.ptrue_to_bounded_phred(predictions[prediction_index])
      )
  )
  # QUAL is prob(variant genotype) / prob(all genotypes)
  # Taking the min to avoid minor numerical issues than can push sum > 1.0.
  # TODO: this is equivalent to the likely better implementation:
  #   genomics_math.perror_to_phred(max(predictions[0], min_ref_confidence))
  # where min_ref_confidence is roughly 1.25e-10 (producing a qual of 99).
  qual = genomics_math.ptrue_to_bounded_phred(min(sum(predictions[1:]), 1.0))
  rounded_qual = round(qual, _QUAL_PRECISION)
  return gq, rounded_qual


def expected_alt_allele_indices(num_alternate_bases: int) -> list[list[int]]:
  """Returns (sorted) expected list of alt_allele_indices, given #alt bases."""
  num_alleles = num_alternate_bases + 1
  alt_allele_indices_list = [
      sorted(list(set(x) - {0}))
      for x in itertools.combinations(range(num_alleles), 2)
  ]
  # alt_allele_indices starts from 0, where 0 refers to the first alt allele.
  # pylint: disable=g-complex-comprehension
  return sorted([
      [i - 1 for i in alt_allele_indices]
      for alt_allele_indices in alt_allele_indices_list
  ])
  # pylint: enable=g-complex-comprehension


def _check_alt_allele_indices(
    call_variants_outputs: Sequence[deepvariant_pb2.CallVariantsOutput],
) -> bool:
  """Returns True if and only if the alt allele indices are valid."""
  all_alt_allele_indices = sorted([
      list(call_variants_output.alt_allele_indices.indices)
      for call_variants_output in call_variants_outputs
  ])
  if all_alt_allele_indices != expected_alt_allele_indices(
      len(call_variants_outputs[0].variant.alternate_bases)
  ):
    logging.warning(
        (
            'Alt allele indices found from call_variants_outputs for '
            'variant %s is %s, which is invalid.'
        ),
        call_variants_outputs[0].variant,
        all_alt_allele_indices,
    )
    return False
  return True


def is_valid_call_variants_outputs(
    call_variants_outputs: Sequence[deepvariant_pb2.CallVariantsOutput],
) -> bool:
  """Returns True if the call_variants_outputs follows our assumptions.

  Args:
    call_variants_outputs: list of CallVariantsOutput to check.

  Returns:
    True if the sanity check passes.
  """

  if not call_variants_outputs:
    return True  # An empty list is a degenerate case.

  if not _check_alt_allele_indices(call_variants_outputs):
    return False

  first_call, other_calls = call_variants_outputs[0], call_variants_outputs[1:]
  # Sanity check that all call_variants_outputs have the same `variant`.
  for call_to_check in other_calls:
    if first_call.variant != call_to_check.variant:
      logging.warning(
          (
              'Expected all inputs to merge_predictions to have the '
              'same `variant`, but getting %s and %s.'
          ),
          first_call.variant,
          call_to_check.variant,
      )
      return False
  return True


def convert_call_variants_outputs_to_probs_dict(
    canonical_variant: variants_pb2.Variant,
    call_variants_outputs: Sequence[deepvariant_pb2.CallVariantsOutput],
    alt_alleles_to_remove: set[str],
    debug_output_all_candidates: str | None = None,
) -> dict[tuple[str, str], list[float]]:
  """Converts a list of CallVariantsOutput to an internal allele probs dict.

  Args:
    canonical_variant: variants_pb2.Variant.
    call_variants_outputs: list of CallVariantsOutput.
    alt_alleles_to_remove: set of strings. Alleles to remove.
    debug_output_all_candidates: If 'ALT', set low qual alleles to be
      soft-filtered.

  Returns:
    Dictionary of {(allele1, allele2): list of probabilities},
    where allele1 and allele2 are strings.
  """
  flattened_dict = collections.defaultdict(list)
  if not call_variants_outputs:
    return flattened_dict

  for call_variants_output in call_variants_outputs:
    allele_set1 = frozenset([canonical_variant.reference_bases])
    allele_set2 = frozenset(
        canonical_variant.alternate_bases[index]
        for index in call_variants_output.alt_allele_indices.indices
    )
    has_alleles_to_rm = bool(alt_alleles_to_remove.intersection(allele_set2))
    if has_alleles_to_rm and debug_output_all_candidates != 'ALT':
      continue
    if has_alleles_to_rm:
      # This block is run when debug_output_all_candidates=ALT
      # It sets genotype likelihood to a placeholder value,
      # which is later used to set GL=1.0 (prob=0).
      p11, p12, p22 = (
          _FILTERED_ALT_PROB,
          _FILTERED_ALT_PROB,
          _FILTERED_ALT_PROB,
      )
    else:
      p11, p12, p22 = call_variants_output.genotype_probabilities
    for set1, set2, p in [
        (allele_set1, allele_set1, p11),
        (allele_set1, allele_set2, p12),
        (allele_set2, allele_set2, p22),
    ]:
      for indices in itertools.product(set1, set2):
        flattened_dict[indices].append(p)
  return flattened_dict


def get_alt_alleles_to_remove(
    call_variants_outputs: Sequence[deepvariant_pb2.CallVariantsOutput],
    qual_filter: float,
) -> set[str]:
  """Returns all the alt alleles with quality below qual_filter.

  Quality is defined as (1-p(ref/ref)). This removes all alt alleles whose
  quality is below the filter value, with the exception that if the set of
  removed alt alleles covers everything in the alternate_bases, the single alt
  allele where the 1-p(ref/ref) is the highest is retained.

  Args:
    call_variants_outputs: list of CallVariantsOutput.
    qual_filter: double. The qual value below which to filter variants.

  Returns:
    Set of strings: alt alleles to remove.
  """
  alt_alleles_to_remove = set()  # first alt is represented as 0.
  if not qual_filter or not call_variants_outputs:
    return alt_alleles_to_remove

  max_qual, max_qual_allele = None, None
  canonical_variant = call_variants_outputs[0].variant
  for call_variants_output in call_variants_outputs:
    # Go through the ones where alt_allele_indices has
    # exactly one element. There are the pileup images that contains information
    # like:
    #    p00, p01, p11
    # or p00, p02, p22
    # ...p00, p0N, pNN
    if len(call_variants_output.alt_allele_indices.indices) == 1:
      # From here, we want to see which ones of these alt alleles (1-N) that we
      # can skip. We can use the concept of QUAL in VCF, and filter out ones
      # where QUAL < FLAGS.qual_filter. This is because if QUAL is too low,
      # it means it is unlikely this has a variant genotype.
      _, qual = compute_quals(
          call_variants_output.genotype_probabilities, prediction_index=0
      )
      alt_allele_index = call_variants_output.alt_allele_indices.indices[0]
      # Keep track of one alt allele with the highest qual score.
      if max_qual is None or max_qual < qual:
        max_qual, max_qual_allele = (
            qual,
            canonical_variant.alternate_bases[alt_allele_index],
        )
      if qual < qual_filter:
        alt_alleles_to_remove.add(
            canonical_variant.alternate_bases[alt_allele_index]
        )

  # If all alt alleles are below `qual_filter`, keep at least one.
  if len(alt_alleles_to_remove) == len(canonical_variant.alternate_bases):
    alt_alleles_to_remove -= set([max_qual_allele])
  return alt_alleles_to_remove


def is_methylated(call: variants_pb2.VariantCall) -> bool:
  """Determines if a VariantCall is methylated.

  A variant is considered methylated if any of its methylation fractions
  (MF) for the reference or alternate alleles is greater than 0.

  Args:
    call: A `variants_pb2.VariantCall` object.

  Returns:
    bool: True if any `methylation_fraction` value is greater than 0,
          otherwise False.
  """
  if 'MF' not in call.info:
    return False

  mf_values = variantcall_utils.get_mf(call)
  if any(mf > 0 for mf in mf_values):
    return True

  return False


class AlleleRemapper:
  """Facilitates removing alt alleles from a Variant.

  This class provides a one-to-shop for managing the information needed to
  remove alternative alleles from Variant. It provides functions and properties
  to get the original alts, the new alts, and asking if alleles (strings) or
  indices (integers) should be retained or eliminated.
  """

  def __init__(
      self, original_alt_alleles: Sequence[str], alleles_to_remove: set[str]
  ):
    self.original_alts = list(original_alt_alleles)
    self.alleles_to_remove = set(alleles_to_remove)

  def keep_index(
      self, allele_index: int, ref_is_zero: bool | str = False
  ) -> bool:
    if ref_is_zero:
      return True if allele_index == 0 else self.keep_index(allele_index - 1)
    else:
      return self.original_alts[allele_index] not in self.alleles_to_remove

  def retained_alt_alleles(self) -> Sequence[str]:
    return [
        alt for alt in self.original_alts if alt not in self.alleles_to_remove
    ]

  def reindex_allele_indexed_fields(
      self, variant: variants_pb2.Variant, fields: frozenset[tuple[str, bool]]
  ) -> None:
    """Updates variant.call fields indexed by ref + alt_alleles.

    Args:
      variant: Variant proto. We will update the info fields of the Variant.call
        protos.
      fields: Iterable of string. Each string should provide a key to an
        alternative allele indexed field in VariantCall.info fields. Each field
        specified here will be updated to remove values associated with alleles
        no longer wanted according to this remapper object.
    """
    for field_info in fields:
      field = field_info[0]
      ref_is_zero = field_info[1]
      for call in variant.calls:
        if field in call.info:
          entry = call.info[field]
          updated = [
              v
              for i, v in enumerate(entry.values)
              if self.keep_index(i, ref_is_zero=ref_is_zero)
          ]
          # We cannot do entry.values[:] = updated as the ListValue type "does
          # not support assignment" so we have to do this grossness.
          del entry.values[:]
          entry.values.extend(updated)


def prune_alleles(
    variant: variants_pb2.Variant, alt_alleles_to_remove: set[str]
) -> variants_pb2.Variant:
  """Remove the alt alleles in alt_alleles_to_remove from canonical_variant.

  Args:
    variant: variants_pb2.Variant.
    alt_alleles_to_remove: iterable of str. Alt alleles to remove from variant.

  Returns:
    variants_pb2.Variant with the alt alleles removed from alternate_bases.
  """
  # If we aren't removing any alt alleles, just return the unmodified variant.
  if not alt_alleles_to_remove:
    return variant

  new_variant = variants_pb2.Variant()
  new_variant.CopyFrom(variant)

  # Cleanup any VariantCall.info fields indexed by alt allele.
  remapper = AlleleRemapper(variant.alternate_bases, alt_alleles_to_remove)
  remapper.reindex_allele_indexed_fields(
      new_variant, _ALT_ALLELE_INDEXED_FORMAT_FIELDS
  )
  new_variant.alternate_bases[:] = remapper.retained_alt_alleles()

  return new_variant


def get_multiallelic_distributions(
    call_variants_outputs: Sequence[deepvariant_pb2.CallVariantsOutput],
    pruned_alleles: set[str],
) -> np.ndarray:
  """Return 9 values for 3 distributions from given multiallelic CVOs.

  This function is only called for sites with two alt alleles remaining after
  pruning. However, call_variants_outputs contains CVOs from pruned and unpruned
  alleles, so we ignore the CVOs containing alleles that were pruned.

  Args:
    call_variants_outputs: list of CVOs for a multiallelic site with exactly two
      alts after pruning. For such a site, we would expect 3 CVOs (alt1, alt2,
      alt1/2). However, there may be more than 3 CVOs if some alleles were
      pruned at this site.
    pruned_alleles: set of strings corresponding to pruned alleles. Used to
      filter CVOs for pruned alleles.

  Returns:
    final_probs: array of shape (1, 9). The 9 values correspond to three model
      output distributions. The first is from the image containing alt1, the
      second is from the image for alt2, the third is from the image with both
      alt1 and alt2.
  """

  alt_allele_indices_to_probs = {}
  first_alt_index = None
  second_alt_index = None

  # Find the CVOs with two alts, corresponding to the image with alt1 and alt2.
  for cvo in call_variants_outputs:
    indices = cvo.alt_allele_indices.indices[:]
    curr_alleles = [cvo.variant.alternate_bases[i] for i in indices]
    curr_alleles_pruned = any([a in pruned_alleles for a in curr_alleles])

    # Ignore CVOs containing pruned alleles.
    if len(indices) == 2 and not curr_alleles_pruned:
      first_alt_index = min(indices)
      second_alt_index = max(indices)
      probs = cvo.genotype_probabilities[:]
      alt_allele_indices_to_probs[(first_alt_index, second_alt_index)] = probs

  # Find the single alt CVOs.
  for cvo in call_variants_outputs:
    if len(cvo.alt_allele_indices.indices[:]) == 1:
      index = cvo.alt_allele_indices.indices[0]
      if index == first_alt_index or index == second_alt_index:
        probs = cvo.genotype_probabilities[:]
        alt_allele_indices_to_probs[index] = probs

  assert len(alt_allele_indices_to_probs) == 3
  # Concatenate all probabilities into one array.
  final_probs = np.array([
      alt_allele_indices_to_probs[first_alt_index]
      + alt_allele_indices_to_probs[second_alt_index]
      + alt_allele_indices_to_probs[(first_alt_index, second_alt_index)]
  ])
  return final_probs


@functools.lru_cache
def get_multiallelic_model(
    use_multiallelic_model: bool,
) -> tf.keras.Model | None:
  """Loads and returns the model, which must be in saved model format.

  Args:
    use_multiallelic_model: if True, use a specialized model for genotype
      resolution of multiallelic cases with two alts.

  Returns:
    A keras model instance if use_multiallelic_model, else None.
  """

  if not use_multiallelic_model:
    return None

  curr_dir = os.path.dirname(__file__)
  multiallelic_model_path = os.path.join(curr_dir, 'multiallelic_model')


  return tf.keras.models.load_model(multiallelic_model_path, compile=False)


def normalize_predictions(predictions: Sequence[float]) -> Sequence[float]:
  """Normalize predictions and handle soft-filtered alt alleles."""
  if sum(predictions) == 0:
    predictions = [1.0] * len(predictions)
  denominator = (
      sum([i if i != _FILTERED_ALT_PROB else 0.0 for i in predictions]) or 1.0
  )
  normalized_predictions = [
      i / denominator if i != _FILTERED_ALT_PROB else 0.0 for i in predictions
  ]
  return normalized_predictions


def correct_nonautosome_probabilities(
    probabilities: list[float],
    variant: variants_pb2.Variant,
) -> Sequence[float]:
  """Recalculate probabilities for non-autosome heterozygous calls."""
  n_alleles = len(variant.alternate_bases) + 1

  # It is assumed that probabilities are stored in the specific order. See
  # most_likely_genotype for details.
  # Each heterozyhous probability is zeroed. For example, for biallelic case
  # the probability of 0/1 genotype becomes zero.
  index = 0
  for h1 in range(0, n_alleles):
    for h2 in range(0, h1 + 1):
      if h2 != h1:
        if len(probabilities) <= index:
          raise ValueError("Probabilties array doesn't match alt alleles.")
        probabilities[index] = 0
      index += 1

  new_sum = sum(probabilities) or 1.0
  return list(map(lambda p: p / new_sum, probabilities))


def is_non_autosome(variant: variants_pb2.Variant) -> bool:
  """Returns True if variant is non_autosome."""
  haploid_contigs_str = _HAPLOID_CONTIGS.value or ''
  parts = haploid_contigs_str.split(',')
  # pylint: disable=g-complex-comprehension
  haploid_contigs = [item for part in parts for item in part.split()]
  # pylint: enable=g-complex-comprehension
  return haploid_contigs and variant.reference_name in haploid_contigs  # pytype: disable=bad-return-type


def is_in_regions(
    variant: variants_pb2.Variant, regions: ranges.RangeSet
) -> bool:
  """Returns True of variant overlaps one of the regions."""
  if regions:
    return regions.variant_overlaps(variant)
  else:
    return False


@functools.lru_cache
def get_par_regions() -> ranges.RangeSet | None:
  """Returns the cached par regions if specified, else None."""
  if _PAR_REGIONS.value:
    return ranges.RangeSet.from_bed(_PAR_REGIONS.value, enable_logging=False)
  else:
    return None


def merge_predictions(
    call_variants_outputs: Sequence[deepvariant_pb2.CallVariantsOutput],
    qual_filter: float | None = None,
    multiallelic_model: tf.keras.Model | None = None,
    debug_output_all_candidates: str | None = None,
) -> tuple[variants_pb2.Variant, Sequence[float]]:
  """Merges the predictions from the multi-allelic calls."""
  # See the logic described in the class PileupImageCreator pileup_image.py
  #
  # Because of the logic above, this function expects all cases above to have
  # genotype_predictions that we can combine from.

  par_regions = get_par_regions()

  if not call_variants_outputs:
    raise ValueError('Expected 1 or more call_variants_outputs.')

  if not is_valid_call_variants_outputs(call_variants_outputs):
    raise ValueError('`call_variants_outputs` did not pass sanity check.')

  first_call, other_calls = call_variants_outputs[0], call_variants_outputs[1:]
  canonical_variant = first_call.variant
  if not other_calls:
    canonical_variant = variant_utils.simplify_variant_alleles(
        canonical_variant
    )
    if is_non_autosome(canonical_variant) and not is_in_regions(
        canonical_variant, par_regions
    ):
      return canonical_variant, correct_nonautosome_probabilities(
          list(first_call.genotype_probabilities), canonical_variant
      )
    return canonical_variant, first_call.genotype_probabilities

  # Special handling of multiallelic variants
  alt_alleles_to_remove = get_alt_alleles_to_remove(
      call_variants_outputs, qual_filter
  )

  # flattened_probs_dict is only used with the multiallelic model
  flattened_probs_dict = convert_call_variants_outputs_to_probs_dict(
      canonical_variant,
      call_variants_outputs,
      alt_alleles_to_remove,
      debug_output_all_candidates,
  )

  if debug_output_all_candidates == 'INFO':
    add_string_field(
        canonical_variant.info,
        'CANDIDATES',
        '|'.join(canonical_variant.alternate_bases),
    )
  if debug_output_all_candidates != 'ALT':
    canonical_variant = prune_alleles(canonical_variant, alt_alleles_to_remove)
  # Run alternate model for multiallelic cases.
  num_alts = len(canonical_variant.alternate_bases)
  if num_alts == 2 and multiallelic_model is not None:
    # We have 3 CVOs for 2 alts. In this case, there are 6 possible genotypes.
    cvo_probs = get_multiallelic_distributions(
        call_variants_outputs, alt_alleles_to_remove
    )
    normalized_predictions = multiallelic_model(cvo_probs).numpy().tolist()[0]
  else:

    def min_alt_filter(probs):
      return min([x for x in probs if x != _FILTERED_ALT_PROB] or [0])

    predictions = [
        min_alt_filter(flattened_probs_dict[(m, n)])
        for _, _, m, n in variant_utils.genotype_ordering_in_likelihoods(
            canonical_variant
        )
    ]
    if sum(predictions) == 0:
      predictions = [1.0] * len(predictions)
    normalized_predictions = normalize_predictions(predictions)

  # Note the simplify_variant_alleles call *must* happen after the predictions
  # calculation above. flattened_probs_dict is indexed by alt allele, and
  # simplify can change those alleles so we cannot simplify until afterwards.
  canonical_variant = variant_utils.simplify_variant_alleles(canonical_variant)
  if is_non_autosome(canonical_variant) and not is_in_regions(
      canonical_variant, par_regions
  ):
    return canonical_variant, correct_nonautosome_probabilities(
        normalized_predictions, canonical_variant
    )
  else:
    return canonical_variant, normalized_predictions


def should_filter(
    variant: variants_pb2.Variant,
    pon_vcf_reader: vcf.VcfReader,
    padding_bases: int = 0,
) -> bool:
  """Returns True if the variant should be filtered based on PON."""
  if pon_vcf_reader is None:
    return False
  query_region = ranges.make_range(
      chrom=variant.reference_name,
      start=variant.start - padding_bases,
      end=variant.end + padding_bases,
  )
  pon_variants = list(pon_vcf_reader.query(query_region))
  if not pon_variants:
    return False
  # TODO: Consider improving this logic to directly match the
  # contig name, the position, and REF and ALT directly.
  variant_key = variant_utils.variant_key(variant)
  for pon_variant in pon_variants:
    if variant_key == variant_utils.variant_key(pon_variant):
      return True
  return False


def add_pon_filter(
    variant_generator: Iterator[variants_pb2.Variant],
    pon_vcf_reader: vcf.VcfReader,
) -> Iterator[variants_pb2.Variant]:
  for variant in variant_generator:
    if dv_vcf_constants.DEEP_VARIANT_PASS in variant.filter and should_filter(
        variant, pon_vcf_reader
    ):
      variant.filter.remove(dv_vcf_constants.DEEP_VARIANT_PASS)
      variant.filter.append(dv_vcf_constants.DEEP_VARIANT_PON)
    yield variant


def write_variants_to_vcf(
    variant_iterable: Iterator[variants_pb2.Variant],
    output_vcf_path: str,
    header: variants_pb2.VcfHeader,
):
  """Writes Variant protos to a VCF file.

  Args:
    variant_iterable: iterable. An iterable of sorted Variant protos.
    output_vcf_path: str. Output file in VCF format.
    header: VcfHeader proto. The VCF header to use for writing the variants.
  """
  logging.info('Writing output to VCF file: %s', output_vcf_path)
  with vcf.VcfWriter(
      output_vcf_path, header=header, round_qualities=True
  ) as writer:
    count = 0
    for variant in variant_iterable:
      if not _ONLY_KEEP_PASS.value or variant.filter == [
          dv_vcf_constants.DEEP_VARIANT_PASS
      ]:
        count += 1
        if _PROCESS_SOMATIC.value:
          writer.write_somatic(variant)
        else:
          writer.write(variant)
        logging.log_every_n(
            logging.INFO, '%s variants written.', _LOG_EVERY_N, count
        )


def _sort_grouped_variants(group: Sequence[deepvariant_pb2.CallVariantsOutput]):
  return sorted(group, key=lambda x: sorted(x.alt_allele_indices.indices))


def _transform_call_variant_group_to_output_variant(
    call_variant_group: Sequence[deepvariant_pb2.CallVariantsOutput],
    qual_filter: float,
    multi_allelic_qual_filter: float,
    sample_name: str,
    use_multiallelic_model: bool,
    debug_output_all_candidates: str | None,
) -> variants_pb2.Variant:
  """Transforms a group of CalVariantOutput to VariantOutput.

  The group of CVOs present in the call_variants_group are converted to the
  Variant proto, with the following filters applied: 1) variants are omitted
  if their quality is lower than the `qual_filter` threshold. 2) multi-allelic
  variants omit individual alleles whose qualities are lower than the
  `multi_allelic_qual_filter` threshold.

  Args:
    call_variant_group: list[CVO]. A group of CallVariantsOutput protos.
    qual_filter: double. The qual value below which to filter variants.
    multi_allelic_qual_filter: double. The qual value below which to filter
      multi-allelic variants.
    sample_name: str. Sample name to write to VCF file.
    use_multiallelic_model: if True, use a specialized model for genotype
      resolution of multiallelic cases with two alts.
    debug_output_all_candidates: if 'ALT', output all alleles considered by
      DeepVariant as ALT alleles.

  Returns:
    Variant proto representing the group of CallVariantsOutput protos.
  """
  multiallelic_model = get_multiallelic_model(
      use_multiallelic_model=use_multiallelic_model
  )
  outputs = _sort_grouped_variants(call_variant_group)
  canonical_variant, predictions = merge_predictions(
      outputs,
      multi_allelic_qual_filter,
      multiallelic_model=multiallelic_model,
      debug_output_all_candidates=debug_output_all_candidates,
  )
  return add_call_to_variant(
      canonical_variant,
      predictions,
      qual_filter=qual_filter,
      sample_name=sample_name,
  )


def _transform_call_variants_output_to_variants(
    input_sorted_tfrecord_path: str,
    sample_name: str,
) -> Iterator[variants_pb2.Variant]:
  """Yields Variant protos in sorted order from CallVariantsOutput protos.

  Args:
    input_sorted_tfrecord_path: str. TFRecord format file containing sorted
      CallVariantsOutput protos.
    sample_name: str. Sample name use in the output VCF and gVCF.

  Yields:
    Variant protos in sorted order representing the CallVariantsOutput calls.
  """
  for call_variant_group in group_call_variants_outputs(
      input_sorted_tfrecord_path, _GROUP_VARIANTS.value
  ):
    yield _transform_call_variant_group_to_output_variant(
        call_variant_group,
        _QUAL_FILTER.value,
        _MULT_ALLELIC_QUAL_FILTER.value,
        sample_name,
        _USE_MULTIALLELIC_MODEL.value,
        _DEBUG_OUTPUT_ALL_CANDIDATES.value,
    )


def dump_variants_to_temp_file(
    variant_protos: Iterator[variants_pb2.Variant],
) -> tempfile._TemporaryFileWrapper:
  temp = tempfile.NamedTemporaryFile()
  tfrecord.write_tfrecords(variant_protos, temp.name)
  return temp


def group_call_variants_outputs(
    input_sorted_tfrecord_path: str, group_variants: bool
) -> Iterator[Sequence[deepvariant_pb2.CallVariantsOutput]]:
  """Yields CallVariantOutputs grouped by their variant range.

  Args:
    input_sorted_tfrecord_path: str. TFRecord format file containing sorted
      CallVariantsOutput protos.
    group_variants: bool. If true, group variants that have same start and end
      position.
  """
  group_fn = None
  if group_variants:
    group_fn = lambda x: variant_utils.variant_range(x.variant)
  for _, group in itertools.groupby(
      tfrecord.read_tfrecords(
          input_sorted_tfrecord_path, proto=deepvariant_pb2.CallVariantsOutput
      ),
      group_fn,
  ):
    yield list(group)


def _concat_vcf(
    output_file: str, temp_vcf_files: Sequence[tempfile._TemporaryFileWrapper]
) -> None:
  """Concatenates a set of temp (g)VCF files."""
  vcf_files_to_concat = [f.name for f in temp_vcf_files]
  vcf_concat.concat(output_file, vcf_files_to_concat)


def process_contiguous_partition(
    contiguous_range_set: Sequence[range_pb2.Range],
    contigs: Sequence[reference_pb2.ContigInfo],
    cvo_paths: Sequence[str],
    temp_file_name: str,
    sample_name: str,
) -> Iterator[variants_pb2.Variant]:
  """Postprocess all CVOs in the given partition and returns an iterator.

  Args:
    contiguous_range_set: set of contiguous ranges to load and transform.
    contigs: all contigs from ref
    cvo_paths: paths to all CVO files
    temp_file_name: path to temp file to variants to.
    sample_name: the sample name to use for the output VCF and gVCF.

  Returns:
    An iterator of processed variants.
  """
  start_time = time.time()
  num_cvo_records = postprocess_variants_lib.process_single_sites_tfrecords(
      contigs,
      cvo_paths,
      temp_file_name,
      contiguous_range_set,
  )
  if contiguous_range_set:
    logging.info(
        'Processing region %s:%s-%s:%s',
        contiguous_range_set[0].reference_name,
        contiguous_range_set[0].start,
        contiguous_range_set[-1].reference_name,
        contiguous_range_set[-1].end,
    )
  logging.info('CVO sorting took %s minutes', (time.time() - start_time) / 60)
  if num_cvo_records == 0:
    return iter([])

  logging.info('Transforming call_variants_output to variants.')
  independent_variants = _transform_call_variants_output_to_variants(
      input_sorted_tfrecord_path=temp_file_name,
      sample_name=sample_name,
  )
  variant_generator = haplotypes.maybe_resolve_conflicting_variants(
      independent_variants
  )
  logging.info('Processed %s variants.', num_cvo_records)
  return variant_generator


def _yield_variants_from_temp_files(
    temp_files: Sequence[tempfile._TemporaryFileWrapper],
) -> Iterable[variants_pb2.Variant]:
  """Yields variants from all the temp files in order.

  Args:
    temp_files: a list of NamedTemporaryFiles objects

  Yields:
    variants read in order from the given temp files.
  """
  for temp_file in temp_files:
    for variant in tfrecord.read_tfrecords(
        temp_file.name, proto=variants_pb2.Variant
    ):
      yield variant


def _decide_to_use_csi(contigs: Sequence[reference_pb2.ContigInfo]) -> bool:
  """Return True if CSI index is to be used over tabix index format.

  If the length of any reference chromosomes exceeds 512M
  (here we use 5e8 to keep a safety margin), we will choose csi
  as the index format. Otherwise we use tbi as default.

  Args:
    contigs: list of contigs.

  Returns:
    A boolean variable indicating if the csi format is to be used or not.
  """
  max_chrom_length = max([c.n_bases for c in contigs])
  return max_chrom_length > 5e8


def build_index(vcf_file: str, csi: bool = False) -> None:
  """A helper function for indexing VCF files.

  Args:
    vcf_file: string. Path to the VCF file to be indexed.
    csi: bool. If true, index using the CSI format.
  """

  if csi:
    tabix.build_csi_index(vcf_file, min_shift=14)
  else:
    tabix.build_index(vcf_file)


def get_cvo_paths(cvo_file_spec: str) -> list[str]:
  """Returns sharded filenames for the `cvo_file_spec` parameter."""
  if sharded_file_utils.is_sharded_file_spec(cvo_file_spec):
    # Input is already sharded, so dynamic sharding check is disabled.
    paths = sharded_file_utils.maybe_generate_sharded_filenames(cvo_file_spec)
  else:
    # Input is expected to be dynamically sharded.
    filename_resolver = cvo_file_spec.replace('.tfrecord.gz', '*')
    all_files = sharded_file_utils.glob_list_sharded_file_patterns(
        filename_resolver
    )
    filename_pattern = cvo_file_spec.replace(
        '.tfrecord.gz', '@' + str(len(all_files)) + '.tfrecord.gz'
    )
    paths = sharded_file_utils.maybe_generate_sharded_filenames(
        filename_pattern
    )
    # This check is to make sure all files we glob is exactly the same as the
    # paths we create, otherwise we have multiple file patterns.
    if sorted(all_files) != sorted(paths):
      raise ValueError(
          'Found multiple file patterns in input filename space: ',
          cvo_file_spec,
      )
  return paths


def get_first_cvo_record(
    paths: Sequence[str],
) -> deepvariant_pb2.CallVariantsOutput | None:
  """Returns the first record from the given paths."""
  return dv_utils.get_one_example_from_examples_path(
      ','.join(paths), proto=deepvariant_pb2.CallVariantsOutput
  )


def get_sample_name(cvo_paths: Sequence[str]) -> str:
  """Determines the sample name to be used for the output VCF and gVCF.

  We check the following sources to determine the sample name and use the first
  name available:
    1) CallVariantsOutput
    2) nonvariant site TFRecords
    3) --sample_name flag
    4) default sample name

  Args:
    cvo_paths: file paths to all CVO files.

  Returns:
    sample_name used when writing the output VCF and gVCF.
  """

  record = get_first_cvo_record(cvo_paths)
  gvcf_record = None
  if _NONVARIANT_SITE_TFRECORD_PATH.value:
    gvcf_record = dv_utils.get_one_example_from_examples_path(
        _NONVARIANT_SITE_TFRECORD_PATH.value, proto=variants_pb2.Variant
    )

  if record is not None:
    sample_name = _extract_single_sample_name(record)
    logging.info(
        'Using sample name from call_variants output. Sample name: %s',
        sample_name,
    )
    if _SAMPLE_NAME.value:
      logging.info('--sample_name is set but was not used.')

  elif (
      _NONVARIANT_SITE_TFRECORD_PATH.value and gvcf_record and gvcf_record.calls
  ):
    sample_name = gvcf_record.calls[0].call_set_name
    logging.info(
        (
            'call_variants output is empty, so using sample name from TFRecords'
            ' at --nonvariant_site_tfrecord_path. Sample name: %s'
        ),
        sample_name,
    )
    if _SAMPLE_NAME.value:
      logging.info('--sample_name is set but was not used.')

  elif _SAMPLE_NAME.value:
    sample_name = _SAMPLE_NAME.value
    logging.info(
        (
            'call_variants output and nonvariant TFRecords are empty. Using'
            ' sample name set with --sample_name. Sample name: %s'
        ),
        sample_name,
    )

  else:
    sample_name = dv_constants.DEFAULT_SAMPLE_NAME
    logging.info(
        (
            'Could not determine sample name and --sample_name is unset. Using'
            ' the default sample name. Sample name: %s'
        ),
        sample_name,
    )
  return sample_name


def run_postprocess_variants_on_region(
    output_vcf: str,
    output_gvcf: str,
    partition: Sequence[range_pb2.Range],
    contigs: Sequence[reference_pb2.ContigInfo],
    all_cvo_paths: Sequence[str],
    header: variants_pb2.VcfHeader,
    is_empty: bool,
    sample_name: str,
) -> None:
  """Runs postprocess_variants on the given partition.

  If the partition is empty, we process all CVO records. If the partition is not
  empty, we process only the CVO records in the partition.

  Args:
    output_vcf: path to the output VCF file.
    output_gvcf: path to the output gVCF file.
    partition: a list of nucleus.genomics.v1.Range protos.
    contigs: all contigs from ref
    all_cvo_paths: paths to all CVO files
    header: the VCF header
    is_empty: if the partition is empty.
    sample_name: the sample name to use for the output VCF and gVCF.

  Returns:
    None (the output is written to the output_vcf and output_gvcf files).
  """
  temp = tempfile.NamedTemporaryFile()
  start_time = time.time()
  if not is_empty:
    variant_generator = process_contiguous_partition(
        partition,
        contigs,
        all_cvo_paths,
        temp.name,
        sample_name,
    )
    pon_reader = (
        vcf.VcfReader(_PON_FILTERING.value) if _PON_FILTERING.value else None
    )
    variant_generator = add_pon_filter(variant_generator, pon_reader)
  else:
    logging.info('call_variants_output is empty. Writing out empty VCF.')
    variant_generator = iter([])
  logging.info(
      'Processing variants (and writing to temporary files) took %s minutes',
      (time.time() - start_time) / 60,
  )

  start_time = time.time()
  if not _NONVARIANT_SITE_TFRECORD_PATH.value:
    if _PROCESS_SOMATIC.value:
      logging.info('Writing variants to somatic VCF.')
    else:
      logging.info('Writing variants to VCF.')
    write_variants_to_vcf(
        variant_iterable=variant_generator,
        output_vcf_path=output_vcf,
        header=header,
    )
    logging.info(
        'VCF creation took %s minutes', (time.time() - start_time) / 60
    )
  else:
    tmp_variant_file = dump_variants_to_temp_file(variant_generator)
    merge_variants.merge_and_write_variants_and_nonvariants(
        _ONLY_KEEP_PASS.value,
        tmp_variant_file.name,
        tfrecord.expanded_paths_if_sharded(
            _NONVARIANT_SITE_TFRECORD_PATH.value
        ),
        _REF.value,
        output_vcf,
        output_gvcf,
        header,
        partition,
        _PROCESS_SOMATIC.value,
    )
    tmp_variant_file.close()
    logging.info(
        'VCF and gVCF creation took %s minutes.',
        (time.time() - start_time) / 60,
    )
  temp.close()


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: postprocess_variants does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)),
          errors.CommandLineError,
      )
    del argv  # Unused.

    if (not _NONVARIANT_SITE_TFRECORD_PATH.value) != (not _GVCF_OUTFILE.value):
      errors.log_and_raise(
          (
              'gVCF creation requires both nonvariant_site_tfrecord_path and '
              'gvcf_outfile flags to be set.'
          ),
          errors.CommandLineError,
      )

    if (
        _USE_MULTIALLELIC_MODEL.value
        and _DEBUG_OUTPUT_ALL_CANDIDATES.value == 'ALT'
    ):
      errors.log_and_raise(
          (
              'debug_output_all_candidates=ALT is incompatible with the '
              'multiallelic model. Use INFO instead.'
          ),
          errors.CommandLineError,
      )

    if _NUM_PARTITIONS.value > 0 and _NUM_PARTITIONS.value < _CPUS.value:
      logging.warning(
          '--num_partitions is less than --cpus. Setting --num_partitions to'
          ' --cpus=%s',
          _CPUS.value,
      )

    proto_utils.uses_fast_cpp_protos_or_die()
    logging_level.set_from_flag()

    fasta_reader = pysam.FastaFile(
        filename=_pysam_resolve_file_path(_REF.value)
    )
    contigs = []
    for reference_index in range(fasta_reader.nreferences):
      contigs.append(
          reference_pb2.ContigInfo(
              name=fasta_reader.references[reference_index],
              n_bases=fasta_reader.lengths[reference_index],
              pos_in_fasta=reference_index,
          )
      )

    cvo_paths = get_cvo_paths(_INFILE.value)
    small_model_cvo_paths = []
    if _SMALL_MODEL_CVO_RECORDS.value:
      small_model_cvo_paths = get_cvo_paths(_SMALL_MODEL_CVO_RECORDS.value)
    all_cvo_paths = cvo_paths + small_model_cvo_paths

    sample_name = get_sample_name(all_cvo_paths)
    header = dv_vcf_constants.deepvariant_header(
        contigs=contigs,
        sample_names=[sample_name],
        add_info_candidates=_DEBUG_OUTPUT_ALL_CANDIDATES.value == 'INFO',
        include_model_id=_SMALL_MODEL_CVO_RECORDS.value is not None,
    )
    if _PROCESS_SOMATIC.value:
      header.filters.append(
          variants_pb2.VcfFilterInfo(
              id=dv_vcf_constants.DEEP_VARIANT_GERMLINE,
              description='Non somatic variants',
          )
      )
    if _PON_FILTERING.value:
      if not _PROCESS_SOMATIC.value:
        raise ValueError(
            'PON filtering is only supported for somatic variant calling.'
        )
      header.filters.append(
          variants_pb2.VcfFilterInfo(
              id=dv_vcf_constants.DEEP_VARIANT_PON,
              description='Filtered by Panel of Normals (PON)',
          )
      )

    is_empty = get_first_cvo_record(all_cvo_paths) is None
    # Run sequentially in the absence of multiple CPUs or partitions.
    if _CPUS.value < 1 and _NUM_PARTITIONS.value < 1:
      logging.info(
          'Running postprocess_variants without parallelism or partitions.'
      )
      run_postprocess_variants_on_region(
          _OUTFILE.value,
          _GVCF_OUTFILE.value,
          [],
          contigs,
          all_cvo_paths,
          header,
          is_empty,
          sample_name,
      )
    else:
      # Otherwise run over multiple partitions, in parallel if cpus > 1.
      calling_regions = calling_regions_utils.build_calling_regions(
          contigs=contigs,
          regions_to_include=calling_regions_utils.parse_regions_flag(
              _REGIONS.value
          ),
          regions_to_exclude=[],
          ref_n_regions=[],
      )
      num_partitions = max(_NUM_PARTITIONS.value, _CPUS.value)
      partitions = calling_regions_utils.partition_calling_regions(
          calling_regions, num_partitions=num_partitions
      )
      temp_vcf_files = [
          tempfile.NamedTemporaryFile(suffix='.gz') for _ in partitions
      ]
      temp_gvcf_files = [
          tempfile.NamedTemporaryFile(suffix='.gz') for _ in partitions
      ]

      if _CPUS.value > 1:
        logging.info(
            'Running postprocess_variants with parallelism using %s CPUs over'
            ' %s partitions.',
            _CPUS.value,
            num_partitions,
        )
        async_results = []
        with multiprocessing.Pool(_CPUS.value) as pool:
          for task_id in range(num_partitions):
            async_results.append(
                pool.apply_async(
                    run_postprocess_variants_on_region,
                    (
                        temp_vcf_files[task_id].name,
                        temp_gvcf_files[task_id].name,
                        partitions[task_id],
                        contigs,
                        all_cvo_paths,
                        header,
                        is_empty,
                        sample_name,
                    ),
                )
            )
          pool.close()
          pool.join()
      else:
        logging.info(
            'Running postprocess_variants sequentially over %s partitions.',
            num_partitions,
        )
        for task_id in range(num_partitions):
          run_postprocess_variants_on_region(
              temp_vcf_files[task_id].name,
              temp_gvcf_files[task_id].name,
              partitions[task_id],
              contigs,
              all_cvo_paths,
              header,
              is_empty,
              sample_name,
          )

      _concat_vcf(_OUTFILE.value, temp_vcf_files)
      if _NONVARIANT_SITE_TFRECORD_PATH.value:
        _concat_vcf(_GVCF_OUTFILE.value, temp_gvcf_files)
      for temp_vcf_file in temp_vcf_files:
        temp_vcf_file.close()
      for temp_gvcf_file in temp_gvcf_files:
        temp_gvcf_file.close()

    start_time = time.time()
    use_csi = _decide_to_use_csi(contigs)
    if str(_OUTFILE.value).endswith('.gz'):
      build_index(_OUTFILE.value, use_csi)
    if _NONVARIANT_SITE_TFRECORD_PATH.value and str(
        _GVCF_OUTFILE.value
    ).endswith('.gz'):
      build_index(_GVCF_OUTFILE.value, use_csi)
    logging.info(
        'Indexing VCF and gVCF took %s minutes.',
        (time.time() - start_time) / 60,
    )


if __name__ == '__main__':
  flags.mark_flags_as_required(['infile', 'outfile', 'ref'])
  logging.set_verbosity(logging.INFO)
  logging.get_absl_logger().setLevel(logging.INFO)
  app.run(main)
