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
# pylint: disable=line-too-long
r"""Converts labeled DeepVariant examples protos into a VCF file.

By default, the GT for each of the VCF entries will be parsed from the variant
in the DeepVariant tf.Example. If the variant doesn't have the GT field, we'll
use the `label` in the example to fill the GT field.

There is an optional --allow_unlabeled_examples flag which will make any
unlabeled examples with ./. as GT. Default for --allow_unlabeled_examples is
false, which means the code will crash if any examples are unlabeled (no GT in
variant AND also no label in tf.Example.)
"""
# pylint: enable=line-too-long

import itertools

from absl import app
from absl import flags
from absl import logging

from deepvariant import dv_utils
from deepvariant import dv_vcf_constants
from third_party.nucleus.io import fasta
from third_party.nucleus.io import tabix
from third_party.nucleus.io import tfrecord
from third_party.nucleus.io import vcf
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils

_ALLOW_UNLABELED_EXAMPLES = flags.DEFINE_bool(
    'allow_unlabeled_examples',
    None,
    'If True, allow unlabeled examples as input and output ./. as the GT.',
)
_REF = flags.DEFINE_string(
    'ref',
    None,
    (
        'Required. Genome reference. Used to get the reference contigs for the '
        'VCF file.'
    ),
)
_EXAMPLES = flags.DEFINE_string(
    'examples',
    None,
    'Required. Path to labeled, DeepVariant tf.Example protos.',
)
_OUTPUT_VCF = flags.DEFINE_string(
    'output_vcf', None, 'Required. Path where we will write out output VCF.'
)
_SAMPLE_NAME = flags.DEFINE_string(
    'sample_name',
    '',
    (
        'The sample name to write into the VCF. By default this is None, '
        'indicating we will use the call_set_name of the sample encoded in the '
        'example variant.'
    ),
)
_MAX_RECORDS = flags.DEFINE_integer(
    'max_records',
    -1,
    (
        'If provided, we will only read in at most max_record examples for '
        'conversion to VCF.'
    ),
)
_LOG_EVERY = flags.DEFINE_integer(
    'log_every',
    10000,
    (
        'How frequently should we provide updates on the conversion process? We'
        ' will log our conversion of every `log_every` variants.'
    ),
)


def _example_sort_key(example):
  return variant_utils.variant_range_tuple(dv_utils.example_variant(example))


def prune_and_simplify_variant(variant):
  """Prunes unused alt alleles and simplifies the variant representation."""
  call = variant_utils.only_call(variant)
  if not call or not variantcall_utils.has_genotypes(call):
    return variant

  gt = call.genotype
  # If genotype is missing (./.), do nothing.
  if any(g < 0 for g in gt):
    return variant

  # Genotype indices: 0 is Ref, 1 is Alt1, 2 is Alt2, etc.
  # Active alt indices (0-based for alternate_bases)
  active_alt_indices = {g - 1 for g in gt if g > 0}

  # If hom-ref (no active alts), keep the first alt.
  if not active_alt_indices:
    active_alt_indices = {0}

  alt_alleles_to_remove = set()
  for i, alt in enumerate(variant.alternate_bases):
    if i not in active_alt_indices:
      alt_alleles_to_remove.add(alt)

  if not alt_alleles_to_remove:
    # Still simplify even if no alleles were pruned.
    return variant_utils.simplify_variant_alleles(variant)

  # Use dv_utils.prune_alleles to prune and remap FORMAT fields.
  pruned_variant = dv_utils.prune_alleles(variant, alt_alleles_to_remove)

  # Remap genotype indices.
  retained_alts = pruned_variant.alternate_bases
  index_map = {0: 0}  # Ref maps to Ref
  new_idx = 1
  for old_idx_0based, alt in enumerate(variant.alternate_bases):
    old_idx = old_idx_0based + 1
    if alt in retained_alts:
      index_map[old_idx] = new_idx
      new_idx += 1
    else:
      index_map[old_idx] = -1  # Pruned

  new_gt = []
  for g in gt:
    new_g = index_map.get(g, -1)
    if new_g == -1:
      raise ValueError(
          f'Genotype index {g} was pruned but it was active in GT {gt}'
      )
    new_gt.append(new_g)

  pruned_call = variant_utils.only_call(pruned_variant)
  variantcall_utils.set_gt(pruned_call, new_gt)

  if 'AD' in pruned_call.info and 'DP' in pruned_call.info:
    new_ad_values = [v.int_value for v in pruned_call.info['AD'].values]
    variantcall_utils.set_format(pruned_call, 'DP', [sum(new_ad_values)])

  # Finally, simplify the alleles.
  simplified_variant = variant_utils.simplify_variant_alleles(pruned_variant)
  return simplified_variant


def examples_to_variants(examples_path, max_records=None):
  """Yields Variant protos from the examples in examples_path.

  This function reads in tf.Examples produced by DeepVariant from examples_path,
  which may contain a sharded spec, sorts them, selects a representive example
  when there are multiple versions representing different alt_alleles, and
  yields the example_variant field from those examples.

  Args:
    examples_path: str. Path, or sharded spec, to labeled tf.Examples produced
      by DeepVariant in training mode.
    max_records: int or None. Maximum number of records to read, or None, to
      read all of the records.

  Yields:
    nucleus.protos.Variant protos in coordinate-sorted order.

  Raises:
    ValueError: if we find a Variant in any example that doesn't have genotypes.
  """
  examples = tfrecord.read_tfrecords(
      examples_path,
      max_records=max_records,
      compression_type='GZIP',
  )
  variants_and_labels = sorted(
      (
          (dv_utils.example_variant(example), dv_utils.example_label(example))
          for example in examples
      ),
      key=lambda x: variant_utils.variant_range_tuple(x[0]),
  )
  for _, group in itertools.groupby(
      variants_and_labels, lambda x: variant_utils.variant_range_tuple(x[0])
  ):
    sorted_group = sorted(
        group, key=lambda x: (x[1] if x[1] is not None else -1), reverse=True
    )
    variant, label = sorted_group[0]
    if not variantcall_utils.has_genotypes(variant_utils.only_call(variant)):
      if label is not None:
        logging.log_every_n(
            logging.INFO,
            'Variant in the example does not have GT. Use label to fill GT.',
            _LOG_EVERY.value,
        )
        if label == 0:
          gt = (0, 0)
        elif label == 1:
          gt = (0, 1)
        elif label == 2:
          gt = (1, 1)
        else:
          raise ValueError(
              'Variant {} has an invalid label: {}. Label must be 0, 1, or 2.'
              .format(variant_utils.variant_key(variant), label)
          )
        call = variant.calls[0] if variant.calls else variant.calls.add()
        variantcall_utils.set_gt(call, gt)
      elif _ALLOW_UNLABELED_EXAMPLES.value:
        call = variant.calls[0] if variant.calls else variant.calls.add()
        variantcall_utils.set_gt(call, (-1, -1))
      else:
        raise ValueError(
            (
                'Variant {} does not have any genotypes. This tool only works '
                'with variants that have been labeled.'
            ).format(variant_utils.variant_key(variant))
        )
    variant = prune_and_simplify_variant(variant)
    yield variant


def peek_sample_name(variants_iter):
  """Gets the call_set_name from the first Variant of variants_iter.

  Args:
    variants_iter: iterable[nucleus.protos.Variant]. Our source of variants.

  Returns:
    tuple of (str, iterable[Variant]). The first element is the call_set_name of
    the first variant of variants_iter, or 'UNKNOWN' if the iterable is empty.
    The second is a new iterable that yields the same elements of variant_iter,
    in the same order, which is necessary to return as we need to peek into
    the original iterator.
  """
  try:
    first = next(variants_iter)
    return first.calls[0].call_set_name, itertools.chain([first], variants_iter)
  except StopIteration:
    # No variants, just return a placeholder value.
    return 'UNKNOWN', iter([])


def main(argv):
  del argv

  contigs = fasta.IndexedFastaReader(_REF.value).header.contigs
  max_records = _MAX_RECORDS.value if _MAX_RECORDS.value >= 0 else None
  variants_iter = examples_to_variants(_EXAMPLES.value, max_records=max_records)

  if not _SAMPLE_NAME.value:
    sample_name, variants_iter = peek_sample_name(variants_iter)
  else:
    sample_name = _SAMPLE_NAME.value
  header = dv_vcf_constants.deepvariant_header(
      contigs=contigs, sample_names=[sample_name]
  )
  with vcf.VcfWriter(_OUTPUT_VCF.value, header=header) as writer:
    for variant in variants_iter:
      variant.calls[0].call_set_name = sample_name
      logging.log_every_n(
          logging.INFO,
          'Converted %s',
          _LOG_EVERY.value,
          variant_utils.variant_key(variant),
      )
      writer.write(variant)
  if _OUTPUT_VCF.value.endswith('.gz'):
    tabix.build_index(_OUTPUT_VCF.value)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'ref',
      'output_vcf',
  ])
  app.run(main)
