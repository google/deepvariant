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
# pylint: disable=line-too-long
r"""Converts labeled DeepVariant examples protos into a VCF file.

./blaze-bin/learning/genomics/deepvariant/labeler/labeled_examples_to_vcf \
  --ref $(pwd)/learning/genomics/deepvariant/testdata/ucsc.hg19.chr20.unittest.fasta.gz \
  --examples $(pwd)/learning/genomics/deepvariant/testdata/golden.training_examples.tfrecord \
  --output_vcf /tmp/golden_training.vcf
"""
# pylint: enable=line-too-long

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import itertools

from absl import app
from tensorflow import flags
from absl import logging

from third_party.nucleus.io import fasta
from third_party.nucleus.io import vcf
from third_party.nucleus.util import io_utils
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils

from deepvariant import dv_vcf_constants
from deepvariant import tf_utils

FLAGS = flags.FLAGS

flags.DEFINE_string(
    'ref', None,
    'Required. Genome reference. Used to get the reference contigs for the '
    'VCF file.')
flags.DEFINE_string(
    'examples', None,
    'Required. Path to labeled, DeepVariant tf.Example protos.'
)
flags.DEFINE_string('output_vcf', None,
                    'Required. Path where we will write out output VCF.')
flags.DEFINE_string(
    'sample_name', '',
    'The sample name to write into the VCF. By default this is None, '
    'indicating we will use the call_set_name of the sample encoded in the '
    'example variant.')
flags.DEFINE_integer(
    'max_records', -1,
    'If provided, we will only read in at most max_record examples for '
    'conversion to VCF.'
)
flags.DEFINE_integer(
    'log_every', 10000,
    'How frequently should we provide updates on the conversion process? We '
    'will log our conversion of every `log_every` variants.'
)


def _example_sort_key(example):
  return variant_utils.variant_range_tuple(tf_utils.example_variant(example))


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
  examples = io_utils.read_tfrecords(examples_path, max_records=max_records)
  variants = sorted(
      (tf_utils.example_variant(example) for example in examples),
      key=variant_utils.variant_range_tuple)

  for _, group in itertools.groupby(variants,
                                    variant_utils.variant_range_tuple):
    variant = next(group)
    if not variantcall_utils.has_genotypes(variant_utils.only_call(variant)):
      raise ValueError((
          'Variant {} does not have any genotypes. This tool only works with '
          'variants that have been labeled.').format(
              variant_utils.variant_key(variant)))
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
    # No variants, just return a dummy value.
    return 'UNKNOWN', iter([])


def main(argv):
  del argv

  contigs = fasta.RefFastaReader(FLAGS.ref).header.contigs
  max_records = FLAGS.max_records if FLAGS.max_records >= 0 else None
  variants_iter = examples_to_variants(FLAGS.examples, max_records=max_records)

  if not FLAGS.sample_name:
    sample_name, variants_iter = peek_sample_name(variants_iter)
  else:
    sample_name = FLAGS.sample_name
  header = dv_vcf_constants.deepvariant_header(
      contigs=contigs, sample_names=[sample_name])
  with vcf.VcfWriter(FLAGS.output_vcf, header=header) as writer:
    for variant in variants_iter:
      variant.calls[0].call_set_name = sample_name
      logging.log_every_n(logging.INFO, 'Converted %s', FLAGS.log_every,
                          variant_utils.variant_key(variant))
      writer.write(variant)


if __name__ == '__main__':
  flags.mark_flags_as_required([
      'examples',
      'ref',
      'output_vcf',
  ])
  app.run(main)
