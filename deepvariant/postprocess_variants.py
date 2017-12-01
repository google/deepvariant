# Copyright 2017 Google Inc.
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

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import collections
import itertools
import tempfile


import numpy as np
import tensorflow as tf

from absl import logging

from deepvariant import logging_level
from deepvariant.core import errors
from deepvariant.core import genomics_io
from deepvariant.core import io_utils
from deepvariant.core import math
from deepvariant.core import proto_utils
from deepvariant.core import variantutils
from deepvariant.core.genomics import variants_pb2
from deepvariant.core.protos import core_pb2
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import postprocess_variants as postprocess_variants_lib

FLAGS = tf.flags.FLAGS

tf.flags.DEFINE_string(
    'infile', None,
    'Required. Path(s) to CallVariantOutput protos in TFRecord format to '
    'postprocess. These should be the complete set of outputs for '
    'call_variants.py.')
tf.flags.DEFINE_string(
    'outfile', None,
    'Required. Destination path where we will write output variant calls in '
    'VCF format.')
tf.flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference in FAI-indexed FASTA format. Used to determine '
    'the sort order for the emitted variants and the VCF header.')
tf.flags.DEFINE_float(
    'qual_filter', 1.0,
    'Any variant with QUAL < qual_filter will be filtered in the VCF file.')
tf.flags.DEFINE_float(
    'multi_allelic_qual_filter', 1.0,
    'The qual value below which to filter multi-allelic variants.')

# The filter field strings to add to variants created by this method.
DEEP_VARIANT_REF_FILTER = 'RefCall'
DEEP_VARIANT_QUAL_FILTER = 'LowQual'
DEEP_VARIANT_PASS = 'PASS'
DEEP_VARIANT_ALL_FILTER_VALUES = frozenset(
    [DEEP_VARIANT_REF_FILTER, DEEP_VARIANT_QUAL_FILTER, DEEP_VARIANT_PASS])

FILTERS = [
    core_pb2.VcfFilterInfo(
        id=DEEP_VARIANT_REF_FILTER,
        description='Genotyping model thinks this site is reference.'),
    core_pb2.VcfFilterInfo(
        id=DEEP_VARIANT_QUAL_FILTER,
        description=
        'Confidence in this variant being real is below calling threshold.'),
]

# Some format fields are indexed by alt allele, such as AD (depth by allele).
# These need to be cleaned up if we remove any alt alleles. Any info field
# listed here will be have its values cleaned up if we've removed any alt
# alleles.
# Each tuple contains: field name, ref_is_zero.
_ALT_ALLELE_INDEXED_FORMAT_FIELDS = frozenset([('AD', True), ('VAF', False)])

# The number of places past the decimal point to round QUAL estimates to.
_QUAL_PRECISION = 7


def _extract_single_sample_name(record):
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
  if len(variant.calls) != 1:
    raise ValueError(
        'Error extracting name: expected one call in {}, not {}'.format(
            record, len(variant.calls)))
  name = variant.calls[0].call_set_name
  if not name:
    raise ValueError(
        'Error extracting name: no call_set_name set: {}'.format(record))

  return name


def compute_filter_fields(variant, min_quality):
  """Computes the filter fields for this variant.

  Variant filters are generated based on its quality score value and particular
  genotype call.

  Args:
    variant: Variant to filter.
    min_quality: Minimum acceptable phred scaled variant detection probability.

  Returns:
    Filter field strings to be added to the variant.
  """
  if variantutils.genotype_type(variant) == variantutils.GenotypeType.hom_ref:
    return [DEEP_VARIANT_REF_FILTER]
  elif variant.quality < min_quality:
    return [DEEP_VARIANT_QUAL_FILTER]
  else:
    return [DEEP_VARIANT_PASS]


def most_likely_genotype(predictions, ploidy=2, n_alleles=2):
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
      related to ploidy and n_alleles is given by
        N = choose(ploidy + n_alleles - 1, n_alleles -1)
      for more information see:
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
  # redacted
  if ploidy != 2:
    raise NotImplementedError('Ploidy != 2 not yet implemented.')
  if n_alleles < 2:
    raise ValueError('n_alleles must be >= 2 but got', n_alleles)
  # redacted
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


def add_call_to_variant(variant, predictions, qual_filter=0, sample_name=None):
  """Fills in Variant record using the prediction probabilities.

  This functions sets the call[0].genotype, call[0].info['GQ'],
  call[0].genotype_probabilities, variant.filter, and variant.quality fields of
  variant based on the genotype likelihoods in predictions.

  Args:
    variant: learning.genomics.deepvariant.core.genomics.Variant protobuf
      to be filled in with info derived from predictions.
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
  if len(variant.calls) != 1:
    raise ValueError('Variant must have exactly one VariantCall record',
                     variant)
  n_alleles = len(variant.alternate_bases) + 1
  index, genotype = most_likely_genotype(predictions, n_alleles=n_alleles)
  gq, variant.quality = compute_quals(predictions, index)
  call = variant.calls[0]
  call.call_set_name = sample_name
  call.genotype[:] = genotype
  variantutils.set_variantcall_gq(call, gq)
  call.genotype_likelihood[:] = [
      math.perror_to_bounded_log10_perror(gp) for gp in predictions
  ]
  variant.filter[:] = compute_filter_fields(variant, qual_filter)
  return variant


def compute_quals(predictions, prediction_index):
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
  gq = np.around(math.ptrue_to_bounded_phred(predictions[prediction_index]))
  # QUAL is prob(variant genotype) / prob(all genotypes)
  # Taking the min to avoid minor numerical issues than can push sum > 1.0.
  # redacted
  #   math.perror_to_phred(max(predictions[0], min_ref_confidence))
  # where min_ref_confidence is something like 1e-15 (producing a qual of 150).
  qual = math.ptrue_to_bounded_phred(min(sum(predictions[1:]), 1.0))
  rounded_qual = round(qual, _QUAL_PRECISION)
  return gq, rounded_qual


def expected_alt_allele_indices(num_alternate_bases):
  """Returns (sorted) expected list of alt_allele_indices, given #alt bases."""
  num_alleles = num_alternate_bases + 1
  alt_allele_indices_list = [
      sorted(list(set(x) - {0}))
      for x in itertools.combinations(range(num_alleles), 2)
  ]
  # alt_allele_indices starts from 0, where 0 refers to the first alt allele.
  return sorted([[i - 1
                  for i in alt_allele_indices]
                 for alt_allele_indices in alt_allele_indices_list])


def _check_alt_allele_indices(call_variants_outputs):
  all_alt_allele_indices = sorted([
      list(call_variants_output.alt_allele_indices.indices)
      for call_variants_output in call_variants_outputs
  ])
  if all_alt_allele_indices != expected_alt_allele_indices(
      len(call_variants_outputs[0].variant.alternate_bases)):
    return False
  return True


def is_valid_call_variants_outputs(call_variants_outputs):
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
      logging.warning('Expected all inputs to merge_predictions to have the '
                      'same `variant`, but getting %s and %s.',
                      first_call.variant, call_to_check.variant)
      return False
  return True


def convert_call_variants_outputs_to_probs_dict(
    canonical_variant, call_variants_outputs, alt_alleles_to_remove):
  """Converts a list of CallVariantsOutput to an internal allele probs dict.

  Args:
    canonical_variant: variants_pb2.Variant.
    call_variants_outputs: list of CallVariantsOutput.
    alt_alleles_to_remove: set of strings. Alleles to remove.

  Returns:
    Dictionary of {(allele1, allele2): list of probabilities},
    where allele1 and allele2 are strings.
  """
  flattened_dict = collections.defaultdict(lambda: [])
  if not call_variants_outputs:
    return flattened_dict

  for call_variants_output in call_variants_outputs:
    allele_set1 = frozenset([canonical_variant.reference_bases])
    allele_set2 = frozenset(
        canonical_variant.alternate_bases[index]
        for index in call_variants_output.alt_allele_indices.indices)
    if alt_alleles_to_remove.intersection(allele_set2):
      continue
    p11, p12, p22 = call_variants_output.genotype_probabilities
    for (set1, set2, p) in [(allele_set1, allele_set1, p11), (allele_set1,
                                                              allele_set2, p12),
                            (allele_set2, allele_set2, p22)]:
      for indices in itertools.product(set1, set2):
        flattened_dict[indices].append(p)
  return flattened_dict


def get_alt_alleles_to_remove(call_variants_outputs, qual_filter):
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
          call_variants_output.genotype_probabilities, prediction_index=0)
      alt_allele_index = call_variants_output.alt_allele_indices.indices[0]
      # Keep track of one alt allele with the highest qual score.
      if max_qual is None or max_qual < qual:
        max_qual, max_qual_allele = (
            qual, canonical_variant.alternate_bases[alt_allele_index])
      if qual < qual_filter:
        alt_alleles_to_remove.add(
            canonical_variant.alternate_bases[alt_allele_index])

  # If all alt alleles are below `qual_filter`, keep at least one.
  if len(alt_alleles_to_remove) == len(canonical_variant.alternate_bases):
    alt_alleles_to_remove -= set([max_qual_allele])
  # redacted
  return alt_alleles_to_remove


# redacted
class AlleleRemapper(object):
  """Faciliates removing alt alleles from a Variant.

  This class provides a one-to-shop for managing the information needed to
  remove alternative alleles from Variant. It provides functions and properties
  to get the original alts, the new alts, and asking if alleles (strings) or
  indices (integers) should be retained or eliminated.
  """

  # redacted

  def __init__(self, original_alt_alleles, alleles_to_remove):
    self.original_alts = list(original_alt_alleles)
    self.alleles_to_remove = set(alleles_to_remove)

  def keep_index(self, allele_index, ref_is_zero=False):
    if ref_is_zero:
      return True if allele_index == 0 else self.keep_index(allele_index - 1)
    else:
      return self.original_alts[allele_index] not in self.alleles_to_remove

  def retained_alt_alleles(self):
    return [
        alt for alt in self.original_alts if alt not in self.alleles_to_remove
    ]

  def reindex_allele_indexed_fields(self, variant, fields):
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
              v for i, v in enumerate(entry.values)
              if self.keep_index(i, ref_is_zero=ref_is_zero)
          ]
          # We cannot do entry.values[:] = updated as the ListValue type "does
          # not support assignment" so we have to do this grossness.
          del entry.values[:]
          entry.values.extend(updated)


def prune_alleles(variant, alt_alleles_to_remove):
  """Remove the alt alleles in alt_alleles_to_remove from canonical_variant.

  Args:
    variant: variants_pb2.Variant.
    alt_alleles_to_remove: iterable of str. Alt alleles to remove from
                           variant.
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
  remapper.reindex_allele_indexed_fields(new_variant,
                                         _ALT_ALLELE_INDEXED_FORMAT_FIELDS)
  new_variant.alternate_bases[:] = remapper.retained_alt_alleles()

  return new_variant


def simplify_alleles(variant):
  """Replaces the alleles in variants with their simplified versions.

  This function takes a variant and replaces its ref and alt alleles with those
  produced by a call to variantutils.simplify_alleles() to remove common postfix
  bases in the alleles that may be present due to pruning away alleles.

  Args:
    variant: learning.genomics.genomics.Variant proto we want to simplify.

  Returns:
    variant with its ref and alt alleles replaced with their simplified
      equivalents.
  """
  simplified_alleles = variantutils.simplify_alleles(variant.reference_bases,
                                                     *variant.alternate_bases)
  variant.reference_bases = simplified_alleles[0]
  variant.alternate_bases[:] = simplified_alleles[1:]
  return variant


def merge_predictions(call_variants_outputs, qual_filter=None):
  """Merges the predictions from the multi-allelic calls."""
  # See the logic described in the class PileupImageCreator pileup_image.py
  #
  # Because of the logic above, this function expects all cases above to have
  # genotype_predictions that we can combine from.
  if not call_variants_outputs:
    raise ValueError('Expected 1 or more call_variants_outputs.')

  if not is_valid_call_variants_outputs(call_variants_outputs):
    raise ValueError('`call_variants_outputs` did not pass sanity check.')

  first_call, other_calls = call_variants_outputs[0], call_variants_outputs[1:]
  canonical_variant = first_call.variant
  if not other_calls:
    return canonical_variant, first_call.genotype_probabilities

  alt_alleles_to_remove = get_alt_alleles_to_remove(call_variants_outputs,
                                                    qual_filter)
  flattened_probs_dict = convert_call_variants_outputs_to_probs_dict(
      canonical_variant, call_variants_outputs, alt_alleles_to_remove)

  canonical_variant = prune_alleles(canonical_variant, alt_alleles_to_remove)
  predictions = [
      min(flattened_probs_dict[(m, n)])
      for _, _, m, n in variantutils.genotype_ordering_in_likelihoods(
          canonical_variant)
  ]
  denominator = sum(predictions)
  # Note the simplify_alleles call *must* happen after the predictions
  # calculation above. flattened_probs_dict is indexed by alt allele, and
  # simplify can change those alleles so we cannot simplify until afterwards.
  canonical_variant = simplify_alleles(canonical_variant)
  return canonical_variant, [i / denominator for i in predictions]


def write_call_variants_output_to_vcf(contigs, input_sorted_tfrecord_path,
                                      output_vcf_path, qual_filter,
                                      multi_allelic_qual_filter, sample_name):
  """Reads CallVariantsOutput protos and writes to a VCF file.

  Variants present in the input TFRecord are converted to VCF format, with the
  following filters applied: 1) variants are omitted if their quality is lower
  than the `qual_filter` threshold. 2) multi-allelic variants omit individual
  alleles whose qualities are lower than the `multi_allelic_qual_filter`
  threshold.

  Args:
    contigs: list(ContigInfo). A list of the reference genome contigs for
      writers that need contig information.
    input_sorted_tfrecord_path: str. TFRecord format file containing sorted
      CallVariantsOutput protos.
    output_vcf_path: str. Output file in VCF format.
    qual_filter: double. The qual value below which to filter variants.
    multi_allelic_qual_filter: double. The qual value below which to filter
      multi-allelic variants.
    sample_name: str. Sample name to write to VCF file.
  """
  logging.info('Writing calls to VCF file: %s', output_vcf_path)
  sync_writer, writer_fn = genomics_io.make_variant_writer(
      output_vcf_path, contigs, samples=[sample_name], filters=FILTERS)
  with sync_writer, io_utils.AsyncWriter(writer_fn) as writer:
    for _, group in itertools.groupby(
        io_utils.read_tfrecords(
            input_sorted_tfrecord_path,
            proto=deepvariant_pb2.CallVariantsOutput),
        lambda x: variantutils.variant_range(x.variant)):
      outputs = list(group)
      canonical_variant, predictions = merge_predictions(
          outputs, multi_allelic_qual_filter)
      variant = add_call_to_variant(
          canonical_variant,
          predictions,
          qual_filter=qual_filter,
          sample_name=sample_name)
      writer.write(variant)


def main(argv=()):
  with errors.clean_commandline_error_exit():
    if len(argv) > 1:
      errors.log_and_raise(
          'Command line parsing failure: postprocess_variants does not accept '
          'positional arguments but some are present on the command line: '
          '"{}".'.format(str(argv)), errors.CommandLineError)
    del argv  # Unused.
    proto_utils.uses_fast_cpp_protos_or_die()

    logging_level.set_from_flag()

    with genomics_io.make_ref_reader(FLAGS.ref) as reader:
      contigs = reader.contigs
    paths = io_utils.maybe_generate_sharded_filenames(FLAGS.infile)
    with tempfile.NamedTemporaryFile() as temp:
      postprocess_variants_lib.process_single_sites_tfrecords(
          contigs, paths, temp.name)
      # Read one CallVariantsOutput record and extract the sample name from it.
      # Note that this assumes that all CallVariantsOutput protos in the infile
      # contain a single VariantCall within their constituent Variant proto, and
      # that the call_set_name is identical in each of the records.
      record = next(
          io_utils.read_tfrecords(
              paths[0], proto=deepvariant_pb2.CallVariantsOutput,
              max_records=1))
      sample_name = _extract_single_sample_name(record)
      write_call_variants_output_to_vcf(
          contigs=contigs,
          input_sorted_tfrecord_path=temp.name,
          output_vcf_path=FLAGS.outfile,
          qual_filter=FLAGS.qual_filter,
          multi_allelic_qual_filter=FLAGS.multi_allelic_qual_filter,
          sample_name=sample_name)


if __name__ == '__main__':
  tf.flags.mark_flags_as_required(['infile', 'outfile', 'ref'])
  tf.app.run()
