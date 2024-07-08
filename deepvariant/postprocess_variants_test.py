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
"""Tests for deepvariant .postprocess_variants."""

import copy
import errno
import gzip
import io
import itertools
import os
import shutil
import sys
from unittest import mock



from absl import flags
from absl import logging
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
import numpy as np
import pysam
import tensorflow as tf

from deepvariant import dv_constants
from deepvariant import dv_vcf_constants
from deepvariant import postprocess_variants
from deepvariant import testdata
from deepvariant.protos import deepvariant_pb2
from third_party.nucleus.io import fasta
from third_party.nucleus.io import tfrecord
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import range_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import genomics_math
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils
from third_party.nucleus.util import vcf_constants
from third_party.nucleus.util.struct_utils import get_string_field


FLAGS = flags.FLAGS

_DEFAULT_SAMPLE_NAME = 'NA12878'

_GVCF_ALT_ALLELE_GL = -99

_FILTERED_ALT_PROB = postprocess_variants._FILTERED_ALT_PROB


def placeholder_reference_reader():
  return fasta.InMemoryFastaReader(
      chromosomes=[
          ('1', 0, 'AACCGGTTACGTTCGATTTTAAAACCCCGGGG'),
          ('2', 0, 'GCAGTGACGTAGCGATGACGTAGACGCTTACG'),
      ]
  )


def placeholder_reference_pysam_reader():
  # Start: 9995000, End: 63025519
  return pysam.FastaFile(testdata.CHR20_FASTA)


def setUpModule():
  testdata.init()


class MockVcfWriter(object):
  """A mock VcfWriter that records the variants written in a list."""

  def __init__(self):
    self.variants_written = []

  def write(self, proto):
    self.variants_written.append(copy.deepcopy(proto))


def _create_variant(
    ref_name,
    start,
    ref_base,
    alt_bases,
    qual,
    filter_field,
    genotype,
    gq,
    likelihoods,
    ad=None,
):
  """Creates a Variant record for testing.

  Args:
    ref_name: reference name for this variant
    start: start position on the contig
    ref_base: reference base(s)
    alt_bases: list(str). alternate base(s)
    qual: PHRED scaled detection probability
    filter_field: filter string for this variant
    genotype: list of integers corresponding to the called genotype
    gq: PHRED scaled genotype quality
    likelihoods: genotype likelihoods for this variant
    ad: list of integers corresponding to allelic depths.

  Returns:
    A Variant record created with the specified arguments.
  """
  return test_utils.make_variant(
      chrom=ref_name,
      start=start,
      alleles=[ref_base] + alt_bases,
      qual=qual,
      filters=filter_field,
      gt=genotype,
      gq=gq,
      gls=likelihoods,
      sample_name=_DEFAULT_SAMPLE_NAME,
      ad=ad,
  )


def _create_variant_with_alleles(ref=None, alts=None, start=0):
  """Creates a Variant record with specified alternate_bases."""
  return variants_pb2.Variant(
      reference_bases=ref,
      alternate_bases=alts,
      start=start,
      calls=[variants_pb2.VariantCall(call_set_name=_DEFAULT_SAMPLE_NAME)],
  )


def _create_call_variants_output(
    indices, probabilities=None, ref=None, alts=None, variant=None
):
  if alts is None != variant is None:
    raise ValueError('Exactly one of either `alts` or `variant` should be set.')
  if not variant:
    variant = _create_variant_with_alleles(ref=ref, alts=alts)
  return deepvariant_pb2.CallVariantsOutput(
      genotype_probabilities=probabilities,
      alt_allele_indices=deepvariant_pb2.CallVariantsOutput.AltAlleleIndices(
          indices=indices
      ),
      variant=variant,
  )


def _read_contents(path, decompress=False):
  with tf.io.gfile.GFile(path, 'rb') as fin:
    contents = fin.read()
    if decompress:
      contents = gzip.GzipFile(path, fileobj=io.BytesIO(contents)).read()
    return contents


def _create_nonvariant(ref_name, start, end, ref_base):
  """Creates a non-variant Variant record for testing.

  Args:
    ref_name: str. Reference name for this variant.
    start: int. start position on the contig [0-based, half open).
    end: int. end position on the contig [0-based, half open).
    ref_base: str. reference base at the start position.

  Returns:
    A non-variant Variant record created with the specified arguments.
  """
  return test_utils.make_variant(
      chrom=ref_name,
      start=start,
      end=end,
      alleles=[ref_base, vcf_constants.GVCF_ALT_ALLELE],
      gt=[0, 0],
      gls=[-0.001, -5, -10],
  )


def make_golden_dataset(compressed_inputs=False):
  if compressed_inputs:
    # Call variants now produce sharded outputs so the golden test has been
    # changed to have sharded input.
    source_path = test_utils.test_tmpfile(
        'golden.postprocess_single_site_input-00000-of-00001.tfrecord.gz'
    )
    tfrecord.write_tfrecords(
        tfrecord.read_tfrecords(
            testdata.GOLDEN_POSTPROCESS_INPUT_SHARDED,
            proto=deepvariant_pb2.CallVariantsOutput,
            compression_type='GZIP',
        ),
        source_path,
        compression_type='GZIP',
    )
    # However, the input filename should still not be in sharded pattern.
    # So replacing the shard pattern to input pattern.
    source_path = testdata.GOLDEN_POSTPROCESS_INPUT
  else:
    source_path = testdata.GOLDEN_POSTPROCESS_INPUT
  return source_path


def create_outfile(file_name, compressed_outputs=False, only_keep_pass=False):
  if only_keep_pass:
    file_name += '.pass_only'
  if compressed_outputs:
    file_name += '.gz'
  return test_utils.test_tmpfile(file_name)


class AlleleRemapperTest(parameterized.TestCase):

  @parameterized.parameters(
      (list('C'), [], [True]),
      (list('C'), ['C'], [False]),
      (list('CT'), [], [True, True]),
      (list('CT'), ['C'], [False, True]),
      (list('CT'), ['T'], [True, False]),
      (list('CT'), ['C', 'T'], [False, False]),
      (list('CTG'), ['C'], [False, True, True]),
      (list('CTG'), ['T'], [True, False, True]),
      (list('CTG'), ['G'], [True, True, False]),
      (list('CTG'), ['C', 'G'], [False, True, False]),
  )
  def test_basic(self, alt_alleles, remove, keep_index_expected):
    remapper = postprocess_variants.AlleleRemapper(alt_alleles, remove)
    self.assertEqual(remapper.original_alts, alt_alleles)
    self.assertEqual(remapper.alleles_to_remove, set(remove))
    self.assertEqual(
        keep_index_expected,
        [remapper.keep_index(i) for i in range(len(alt_alleles))],
    )
    # When our i is 1 for the first alt allele, we expect that we will get back
    # our keep_index_expected but also that keep_index(i==0) is True for the
    # reference allele.
    self.assertEqual(
        [True] + keep_index_expected,
        [
            remapper.keep_index(i, ref_is_zero=True)
            for i in range(len(alt_alleles) + 1)
        ],
    )

  def test_makes_copy_of_inputs(self):
    alt_alleles = ['A', 'B']
    removes = {'B'}
    remapper = postprocess_variants.AlleleRemapper(alt_alleles, removes)
    del alt_alleles[0]
    removes -= {'B'}
    self.assertEqual(remapper.original_alts, ['A', 'B'])
    self.assertEqual(remapper.alleles_to_remove, {'B'})


class PostprocessVariantsTest(parameterized.TestCase):

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters(
      (compressed_inputs_and_outputs, only_keep_pass)
      for compressed_inputs_and_outputs in [False, True]
      for only_keep_pass in [False, True]
  )
  # pylint: enable=g-complex-comprehension
  @flagsaver.flagsaver
  def test_call_end2end(self, compressed_inputs_and_outputs, only_keep_pass):
    FLAGS.infile = make_golden_dataset(compressed_inputs_and_outputs)
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.outfile = create_outfile(
        'calls.vcf', compressed_inputs_and_outputs, only_keep_pass
    )
    FLAGS.nonvariant_site_tfrecord_path = testdata.GOLDEN_POSTPROCESS_GVCF_INPUT
    FLAGS.gvcf_outfile = create_outfile(
        'gvcf_calls.vcf', compressed_inputs_and_outputs, only_keep_pass
    )
    FLAGS.only_keep_pass = only_keep_pass
    FLAGS.cpus = 0
    postprocess_variants.main(['postprocess_variants.py'])

    if only_keep_pass:
      vcf_output = testdata.GOLDEN_POSTPROCESS_OUTPUT_PASS_ONLY
    else:
      vcf_output = testdata.GOLDEN_POSTPROCESS_OUTPUT
    self.assertEqual(
        _read_contents(FLAGS.outfile, compressed_inputs_and_outputs),
        _read_contents(vcf_output),
    )
    self.assertEqual(
        _read_contents(FLAGS.gvcf_outfile, compressed_inputs_and_outputs),
        _read_contents(testdata.GOLDEN_POSTPROCESS_GVCF_OUTPUT),
    )

    if compressed_inputs_and_outputs:
      self.assertTrue(tf.io.gfile.exists(FLAGS.outfile + '.tbi'))
      self.assertTrue(tf.io.gfile.exists(FLAGS.gvcf_outfile + '.tbi'))

  @flagsaver.flagsaver
  def test_haploid_contigs(self):
    FLAGS.infile = make_golden_dataset()
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.outfile = create_outfile('haploid_contigs.calls.vcf')
    FLAGS.nonvariant_site_tfrecord_path = testdata.GOLDEN_POSTPROCESS_GVCF_INPUT
    FLAGS.gvcf_outfile = create_outfile('haploid_contigs.gvcf_calls.vcf')
    FLAGS.haploid_contigs = 'chr20'
    FLAGS.cpus = 0
    postprocess_variants.main(['postprocess_variants.py'])
    vcf_output = testdata.GOLDEN_POSTPROCESS_OUTPUT_HAPLOID
    self.assertEqual(
        _read_contents(FLAGS.outfile),
        _read_contents(vcf_output),
    )
    self.assertEqual(
        _read_contents(FLAGS.gvcf_outfile),
        _read_contents(testdata.GOLDEN_POSTPROCESS_GVCF_OUTPUT_HAPLOID),
    )

  @flagsaver.flagsaver
  def test_group_variants(self):
    FLAGS.infile = testdata.GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_INPUT
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.outfile = create_outfile('calls.vcf')
    FLAGS.cpus = 0

    FLAGS.group_variants = True
    with self.assertRaisesRegex(
        ValueError, '`call_variants_outputs` did not pass sanity check.'
    ):
      postprocess_variants.main(['postprocess_variants.py'])

    FLAGS.group_variants = False
    postprocess_variants.main(['postprocess_variants.py'])
    self.assertEqual(
        _read_contents(FLAGS.outfile),
        _read_contents(
            testdata.GOLDEN_VCF_CANDIDATE_IMPORTER_POSTPROCESS_OUTPUT
        ),
    )

  @parameterized.parameters(False, True)
  def test_build_index(self, use_csi):
    vcf_file_gz = os.path.join(
        absltest.get_default_test_tmpdir(), 'call_test_id_%s.vcf.gz' % use_csi
    )
    shutil.copy(testdata.GOLDEN_POSTPROCESS_OUTPUT_COMPRESSED, vcf_file_gz)
    postprocess_variants.build_index(vcf_file_gz, use_csi)

    if use_csi:
      self.assertFalse(tf.io.gfile.exists(vcf_file_gz + '.tbi'))
      self.assertTrue(tf.io.gfile.exists(vcf_file_gz + '.csi'))
    else:
      self.assertFalse(tf.io.gfile.exists(vcf_file_gz + '.csi'))
      self.assertTrue(tf.io.gfile.exists(vcf_file_gz + '.tbi'))

  @flagsaver.flagsaver
  def test_reading_sharded_input_with_empty_shards_does_not_crash(self):
    valid_variants = tfrecord.read_tfrecords(
        testdata.GOLDEN_POSTPROCESS_INPUT_SHARDED,
        proto=deepvariant_pb2.CallVariantsOutput,
        compression_type='GZIP',
    )
    empty_shard_one = test_utils.test_tmpfile(
        'reading_empty_shard-00000-of-00002.tfrecord.gz'
    )
    none_empty_shard_two = test_utils.test_tmpfile(
        'reading_empty_shard-00001-of-00002.tfrecord.gz'
    )
    tfrecord.write_tfrecords([], empty_shard_one, compression_type='GZIP')
    tfrecord.write_tfrecords(
        valid_variants, none_empty_shard_two, compression_type='GZIP'
    )
    FLAGS.infile = test_utils.test_tmpfile('reading_empty_shard.tfrecord.gz')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.outfile = test_utils.test_tmpfile('calls_reading_empty_shard.vcf')
    FLAGS.cpus = 0

    postprocess_variants.main(['postprocess_variants.py'])

  @parameterized.parameters(
      (
          [
              test_utils.make_variant(
                  chrom='chr20',
                  start=3,
                  end=4,
                  alleles=['A', 'T'],
                  sample_name='vcf_sample_name',
                  gt=[0, 1],
                  gls=[-0.001, -5, -10],
                  ad=[0, 1, 2],
              )
          ],
          [],
          'flag_sample_name',
          'vcf_sample_name',
      ),
      (
          [],
          [
              test_utils.make_variant(
                  chrom='chr20',
                  start=0,
                  end=1,
                  alleles=['A', vcf_constants.GVCF_ALT_ALLELE],
                  sample_name='gvcf_sample_name',
                  gt=[0, 0],
                  gls=[-0.001, -5, -10],
                  ad=[0, 1, 2],
              )
          ],
          'flag_sample_name',
          'gvcf_sample_name',
      ),
      # flag_sample_name only used when no CVOs or nonvariant TFRecords present.
      ([], [], 'flag_sample_name', 'flag_sample_name'),
      ([], [], None, dv_constants.DEFAULT_SAMPLE_NAME),
  )
  @flagsaver.flagsaver
  def test_sample_name_set_correctly(
      self, variants, nonvariants, sample_name_flag, expected_sample_name
  ):
    shard = test_utils.test_tmpfile('records.cvo.tfrecord-00000-of-00001')
    cvos = [
        _create_call_variants_output(
            indices=[0], probabilities=[0.19, 0.75, 0.06], variant=variant
        )
        for variant in variants
    ]

    tfrecord.write_tfrecords(cvos, shard, compression_type='GZIP')
    FLAGS.infile = test_utils.test_tmpfile('records.cvo.tfrecord@1')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.outfile = test_utils.test_tmpfile('records.vcf')
    FLAGS.sample_name = sample_name_flag
    FLAGS.cpus = 0

    FLAGS.nonvariant_site_tfrecord_path = test_utils.test_tmpfile(
        'records.postprocess_gvcf_input.tfrecord.gz'
    )
    tfrecord.write_tfrecords(
        nonvariants,
        FLAGS.nonvariant_site_tfrecord_path,
        compression_type='GZIP',
    )
    FLAGS.gvcf_outfile = test_utils.test_tmpfile('records.g.vcf')
    postprocess_variants.main(['postprocess_variants.py'])

    fasta_reader = fasta.IndexedFastaReader(FLAGS.ref)
    contigs = fasta_reader.header.contigs
    expected_vcf_header = dv_vcf_constants.deepvariant_header(
        contigs=contigs, sample_names=[expected_sample_name]
    )
    vcf_reader = vcf.VcfReader(FLAGS.outfile)
    self.assertEqual(vcf_reader.header, expected_vcf_header)
    self.assertLen(list(vcf_reader), len(variants))

    expected_gvcf = dv_vcf_constants.deepvariant_header(
        contigs=contigs, sample_names=[expected_sample_name]
    )
    gvcf_reader = vcf.VcfReader(FLAGS.gvcf_outfile)
    self.assertEqual(gvcf_reader.header, expected_gvcf)
    self.assertLen(list(gvcf_reader), len(nonvariants) + len(variants))

  def test_extract_single_variant_name(self):
    record = _create_call_variants_output(
        indices=[0], probabilities=[0.19, 0.75, 0.06], ref='A', alts=['C', 'T']
    )
    expected = _DEFAULT_SAMPLE_NAME
    actual = postprocess_variants._extract_single_sample_name(record)
    self.assertEqual(actual, expected)

  @parameterized.parameters(
      ([],),
      (['multiple', 'names'],),
  )
  def test_exception_extract_single_variant_name(self, names):
    variant_calls = [
        variants_pb2.VariantCall(call_set_name=name) for name in names
    ]
    variant = variants_pb2.Variant(calls=variant_calls)
    record = deepvariant_pb2.CallVariantsOutput(variant=variant)
    with self.assertRaisesRegex(ValueError, 'Expected exactly one VariantCal'):
      postprocess_variants._extract_single_sample_name(record)

  @parameterized.parameters(
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.19, 0.75, 0.06],
                  ref='A',
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.03, 0.93, 0.04],
                  ref='A',
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.92, 0.05],
                  ref='A',
                  alts=['C', 'T'],
              ),
          ],
          set(),
          False,
          {
              ('A', 'A'): [0.19, 0.03, 0.03],
              ('A', 'C'): [0.75, 0.92],
              ('C', 'C'): [0.06, 0.05],
              ('A', 'T'): [0.93, 0.92],
              ('T', 'T'): [0.04, 0.05],
              ('C', 'T'): [0.05],
              ('T', 'C'): [0.05],
          },
      ),
      # Example where all alt alleles are below qual_filter, but we keep one
      # where the qual is highest among the ones filtered out ('T')
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.19, 0.75, 0.06],
                  ref='A',
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.03, 0.93, 0.04],
                  ref='A',
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.92, 0.05],
                  ref='A',
                  alts=['C', 'T'],
              ),
          ],
          set(['C']),
          False,
          {
              ('A', 'A'): [0.03],
              ('A', 'T'): [0.93],
              ('T', 'T'): [0.04],
          },
      ),
      # Test debug_output_all_candidates=ALT
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.03, 0.05, 0.85],
                  ref='A',
                  alts=['ACC', 'ACCC'],
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.03, 0.04, 0.92],
                  ref='A',
                  alts=['ACC', 'ACCC'],
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.05, 0.92],
                  ref='A',
                  alts=['ACC', 'ACCC'],
              ),
          ],
          set(['ACC']),
          'ALT',
          {
              ('A', 'A'): [_FILTERED_ALT_PROB, 0.03, _FILTERED_ALT_PROB],
              ('A', 'ACC'): [_FILTERED_ALT_PROB, _FILTERED_ALT_PROB],
              ('ACC', 'ACC'): [_FILTERED_ALT_PROB, _FILTERED_ALT_PROB],
              ('A', 'ACCC'): [0.04, _FILTERED_ALT_PROB],
              ('ACC', 'ACCC'): [_FILTERED_ALT_PROB],
              ('ACCC', 'ACC'): [_FILTERED_ALT_PROB],
              ('ACCC', 'ACCC'): [0.92, _FILTERED_ALT_PROB],
          },
      ),
      # Test debug_output_all_candidates=None
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.03, 0.05, 0.85],
                  ref='A',
                  alts=['ACC', 'ACCC'],
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.03, 0.04, 0.92],
                  ref='A',
                  alts=['ACC', 'ACCC'],
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.05, 0.92],
                  ref='A',
                  alts=['ACC', 'ACCC'],
              ),
          ],
          set(['ACC']),
          None,
          {
              ('A', 'A'): [0.03],
              ('A', 'ACCC'): [0.04],
              ('ACCC', 'ACCC'): [0.92],
          },
      ),
  )
  def test_convert_call_variants_outputs_to_probs_dict(
      self,
      call_variants_outputs,
      alt_alleles_to_remove,
      debug_output_all_candidates,
      expected_probs_dict,
  ):
    # In the current code, all call_variants_outputs have the same variant
    # field.
    canonical_variant = call_variants_outputs[0].variant
    self.assertEqual(
        postprocess_variants.convert_call_variants_outputs_to_probs_dict(
            canonical_variant,
            call_variants_outputs,
            alt_alleles_to_remove,
            debug_output_all_candidates,
        ),
        expected_probs_dict,
    )

  @parameterized.parameters(
      # Example with 2 alternate_bases:
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[1], probabilities=[0.03, 0.93, 0.04], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.92, 0.05],
                  alts=['C', 'T'],
              ),
          ],
          [0.03, 0.75, 0.05, 0.92, 0.05, 0.04],
      ),
      # One more example with 2 alternate_bases:
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
      (
          [
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.978, 0.03, 0.002],
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.992, 0.007, 0.001],
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.99997, 0.00002, 0.00001],
                  alts=['C', 'T'],
              ),
          ],
          [0.978, 0.00002, 0.00001, 0.007, 0.001, 0.001],
      ),
      # An extreme case where our logic could result in ZeroDivisionError if
      # we don't handle this special case.
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.0, 1.0, 0.0], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[1], probabilities=[0.00, 1.0, 0.0], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[0, 1], probabilities=[1.0, 0.0, 0.0], alts=['C', 'T']
              ),
          ],
          [1.0 / 6] * 6,
      ),
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1.
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['A']
              ),
          ],
          [0.19, 0.75, 0.06],
      ),
      # expected_unnormalized_probs is min of
      # 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, 0/3, 1/3, 2/3, 3/3.
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.999, 0.001, 0],
                  alts=['C', 'G', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0, 1], probabilities=[0, 1, 0], alts=['C', 'G', 'T']
              ),
              _create_call_variants_output(
                  indices=[0, 2],
                  probabilities=[0.0001, 0.9996, 0.0003],
                  alts=['C', 'G', 'T'],
              ),
              _create_call_variants_output(
                  indices=[1], probabilities=[0, 1, 0], alts=['C', 'G', 'T']
              ),
              _create_call_variants_output(
                  indices=[1, 2],
                  probabilities=[0.0001, 0.0002, 0.9997],
                  alts=['C', 'G', 'T'],
              ),
              _create_call_variants_output(
                  indices=[2],
                  probabilities=[0.00004, 0.9999, 0.00006],
                  alts=['C', 'G', 'T'],
              ),
          ],
          [0, 0.001, 0, 0.0002, 0, 0, 0.0002, 0.0003, 0.9997, 0.00006],
      ),
  )
  def test_merge_predictions_probs(self, inputs, expected_unnormalized_probs):
    denominator = sum(expected_unnormalized_probs)
    for permuted_inputs in itertools.permutations(inputs):
      _, predictions = postprocess_variants.merge_predictions(permuted_inputs)
      np.testing.assert_almost_equal(
          predictions, [x / denominator for x in expected_unnormalized_probs]
      )

  @parameterized.parameters(
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.999, 0.001, 0],
                  alts=['C', 'G', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0, 1], probabilities=[0, 1, 0], alts=['C', 'G', 'T']
              ),
              _create_call_variants_output(
                  indices=[0, 2],
                  probabilities=[0.0001, 0.9996, 0.0003],
                  alts=['C', 'G', 'T'],
              ),
              _create_call_variants_output(
                  indices=[1], probabilities=[0, 1, 0], alts=['C', 'G', 'T']
              ),
              _create_call_variants_output(
                  indices=[1, 2],
                  probabilities=[0.0001, 0.0002, 0.9997],
                  alts=['C', 'G', 'T'],
              ),
              _create_call_variants_output(
                  indices=[2],
                  probabilities=[0.00004, 0.9999, 0.00006],
                  alts=['C', 'G', 'T'],
              ),
          ],
          ['C', 'G', 'T'],
      ),
  )
  def test_candidate_info_field(self, inputs, expected_candidates):
    merge_pred = postprocess_variants.merge_predictions
    variant, _ = merge_pred(inputs, debug_output_all_candidates='INFO')
    obs_candidates = get_string_field(variant.info, 'CANDIDATES')[0].split('|')
    self.assertSameElements(obs_candidates, expected_candidates)

  @parameterized.parameters((
      [0.001, 0.017, 0.30, _FILTERED_ALT_PROB, 0.327],
      [0.0015504, 0.0263566, 0.4651163, 0.0, 0.5069767],
  ))
  def test_normalize_predictions(self, predictions, expected_predictions):
    norm_predictions = postprocess_variants.normalize_predictions(predictions)
    np.testing.assert_almost_equal(norm_predictions, expected_predictions)

  @parameterized.parameters(
      # Example with 2 alternate_bases:
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[1], probabilities=[0.03, 0.93, 0.04], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.92, 0.05],
                  alts=['C', 'T'],
              ),
          ],
          [0.033062, 0.10498016, 0.00496365, 0.5842303, 0.2543793, 0.01838462],
      ),
      # One more example with 2 alternate_bases:
      # expected_unnormalized_probs is min of 0/0, 0/1, 1/1, 0/2, 1/2, 2/2
      (
          [
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.978, 0.03, 0.002],
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.992, 0.007, 0.001],
                  alts=['C', 'T'],
              ),
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.99997, 0.00002, 0.00001],
                  alts=['C', 'T'],
              ),
          ],
          [
              9.3330729e-01,
              1.5126608e-02,
              6.1836297e-04,
              4.9650513e-02,
              2.9180625e-05,
              1.2679433e-03,
          ],
      ),
      # An extreme case where our logic could result in ZeroDivisionError if
      # we don't handle this special case.
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.0, 1.0, 0.0], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[1], probabilities=[0.0, 1.0, 0.0], alts=['C', 'T']
              ),
              _create_call_variants_output(
                  indices=[0, 1], probabilities=[1.0, 0.0, 0.0], alts=['C', 'T']
              ),
          ],
          [
              1.3300395e-03,
              9.5756045e-03,
              1.9776919e-05,
              7.6043198e-04,
              9.3802148e-01,
              5.0292656e-02,
          ],
      ),
      # Example where all alt alleles are below qual_filter, but we keep one
      # where the qual is highest among the ones filtered out.
      (
          [
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T']),
              ),
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.99, 0.01, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T']),
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T']),
              ),
          ],
          [0.99, 0.01, 0.0],
          6,
      ),
  )
  @flagsaver.flagsaver
  def test_merge_predictions_multiallelics_probs(
      self, inputs, expected_unnormalized_probs, qual_filter=None
  ):
    FLAGS.use_multiallelic_model = True
    multiallelic_model = postprocess_variants.get_multiallelic_model(
        use_multiallelic_model=FLAGS.use_multiallelic_model
    )
    denominator = sum(expected_unnormalized_probs)
    for permuted_inputs in itertools.permutations(inputs):
      _, predictions = postprocess_variants.merge_predictions(
          permuted_inputs,
          multiallelic_model=multiallelic_model,
          qual_filter=qual_filter,
      )
      np.testing.assert_almost_equal(
          predictions,
          [x / denominator for x in expected_unnormalized_probs],
          decimal=5,
      )

  @parameterized.parameters(
      # With 1 alt allele, we expect to see 1 alt_allele_indices: [0].
      (
          [
              _create_call_variants_output(
                  indices=[1], probabilities=[0.19, 0.75, 0.06], alts=['A']
              ),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      # With 2 alt alleles, we expect to see 3 alt_allele_indices.
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['G', 'T']
              ),
              _create_call_variants_output(
                  indices=[1], probabilities=[0.03, 0.93, 0.04], alts=['G', 'T']
              ),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      # With 2 alt alleles, we expect to see 3 alt_allele_indices:
      # [0], [1], [0, 1].
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['G', 'T']
              ),
              _create_call_variants_output(
                  indices=[0], probabilities=[0.03, 0.93, 0.04], alts=['G', 'T']
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.03, 0.93, 0.04],
                  alts=['G', 'T'],
              ),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      (
          [
              _create_call_variants_output(
                  indices=[0], probabilities=[0.19, 0.75, 0.06], alts=['G', 'T']
              ),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      ([], 'Expected 1 or more call_variants_outputs.'),
      # With 3 alt alleles, we expect to see 6 alt_allele_indices.
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.999, 0.001, 0],
                  alts=['AA', 'T', 'AAA'],
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0, 1, 0],
                  alts=['AA', 'T', 'AAA'],
              ),
              _create_call_variants_output(
                  indices=[0, 2],
                  probabilities=[0.0001, 0.9996, 0.0003],
                  alts=['AA', 'T', 'AAA'],
              ),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      # reference_bases have to be exactly the same.
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.999, 0.001, 0],
                  variant=_create_variant_with_alleles(
                      ref='A', alts=['T', 'C']
                  ),
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.2, 0.8, 0],
                  variant=_create_variant_with_alleles(
                      ref='A', alts=['T', 'C']
                  ),
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.2, 0.8, 0],
                  variant=_create_variant_with_alleles(
                      ref='G', alts=['T', 'C']
                  ),
              ),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
      # alternate_bases have to be exactly the same. Different orders are
      # not acceptable either.
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.999, 0.001, 0],
                  variant=_create_variant_with_alleles(alts=['T', 'C']),
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[0.2, 0.8, 0],
                  variant=_create_variant_with_alleles(alts=['T', 'C']),
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.2, 0.8, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T']),
              ),
          ],
          '`call_variants_outputs` did not pass sanity check.',
      ),
  )
  def test_exception_merge_predictions(self, inputs, text):
    with self.assertRaisesRegex(ValueError, text):
      postprocess_variants.merge_predictions(inputs)

  @parameterized.parameters(
      (
          _create_variant(
              '1',
              1,
              'A',
              ['C'],
              20.0,
              dv_vcf_constants.DEEP_VARIANT_PASS,
              [1, 1],
              15,
              [-2.0, -9.90308995105826, -0.0043648054],
              [10, 10],
          ),
          [1, 1],
      ),
      (
          _create_variant(
              'GL000220.1',
              10000210,
              'C',
              ['T'],
              50.0,
              dv_vcf_constants.DEEP_VARIANT_NO_CALL,
              [1, 1],
              25,
              [-2.0, -9.90308995105826, -0.0043648054],
              [0, 0],
          ),
          [-1, -1],
      ),
      (
          _create_variant(
              'GL000220.1',
              10000210,
              'C',
              ['T'],
              5.0,
              dv_vcf_constants.DEEP_VARIANT_NO_CALL,
              [1, 1],
              25,
              [-2.0, -9.90308995105826, -0.0043648054],
              [0, 0],
          ),
          [-1, -1],
      ),
  )
  def test_uncall_gt_if_no_ad(self, variant, expected_gt):
    postprocess_variants.uncall_gt_if_no_ad(variant)
    self.assertEqual(variant.calls[0].genotype, expected_gt)

  @parameterized.parameters(
      (
          _create_variant(
              'X',
              10000210,
              'CACA',
              ['C'],
              0.04364805402,
              dv_vcf_constants.DEEP_VARIANT_REF_FILTER,
              [0, 0],
              19,
              [-0.0043648054, -2.30102999566, -2.30102999566],
              [1, 0],
          ),
          [-1, -1],
      ),
      (
          _create_variant(
              'chrY',
              10000210,
              'C',
              ['T'],
              0.00043431619,
              dv_vcf_constants.DEEP_VARIANT_REF_FILTER,
              [0, 0],
              20,
              [-0.00004343161, -4.0, -9.90308995105826],
              [0, 1],
          ),
          [0, 0],
      ),
      (
          _create_variant(
              'X',
              10000210,
              'CACA',
              ['C', 'A'],
              0.0217691925,
              dv_vcf_constants.DEEP_VARIANT_REF_FILTER,
              [0, 0],
              13,
              [-0.00217691925, -3, -3, -3, -3, -3],
              [1, 1, 1],
          ),
          [-1, -1],
      ),
  )
  def test_uncall_homref_gt_if_lowqual(self, variant, expected_gt):
    postprocess_variants.uncall_homref_gt_if_lowqual(variant, 20)
    self.assertEqual(variant.calls[0].genotype, expected_gt)

  @parameterized.parameters(
      (
          [0.01, 0.0, 0.99],
          _create_variant(
              'GL000220.1',
              1,
              'A',
              ['.'],
              20.0,
              dv_vcf_constants.DEEP_VARIANT_PASS,
              [1, 1],
              20,
              [-2.0, -9.90308995105826, -0.0043648054],
              [1, 0],
          ),
      ),
      (
          [0.01, 0.0, 0.99],
          _create_variant(
              'GL000220.1',
              10000210,
              'C',
              ['T'],
              20.0,
              dv_vcf_constants.DEEP_VARIANT_PASS,
              [1, 1],
              20,
              [-2.0, -9.90308995105826, -0.0043648054],
              [0, 2],
          ),
      ),
      (
          [0.001, 0.999, 0.0],
          _create_variant(
              '20',
              10000210,
              'C',
              ['CT'],
              30.0,
              dv_vcf_constants.DEEP_VARIANT_PASS,
              [0, 1],
              30,
              [-3.0, -0.00043451177, -9.90308995105826],
          ),
      ),
      (
          [0.0001, 0.0, 0.9999],
          _create_variant(
              '1',
              1,
              'C',
              ['T'],
              40.0,
              dv_vcf_constants.DEEP_VARIANT_PASS,
              [1, 1],
              40,
              [-4.0, -9.90308995105826, -0.00004343161],
          ),
      ),
      (
          [0.1, 0.90, 0.0],
          _create_variant(
              '20',
              10000210,
              'A',
              ['T'],
              10.0,
              dv_vcf_constants.DEEP_VARIANT_PASS,
              [0, 1],
              10,
              [-1.0, -0.04575749056, -9.90308995105826],
          ),
      ),
      (
          [0.99, 0.005, 0.005],
          _create_variant(
              'X',
              10000210,
              'CACA',
              ['C'],
              0.04364805402,
              dv_vcf_constants.DEEP_VARIANT_REF_FILTER,
              [0, 0],
              20,
              [-0.0043648054, -2.30102999566, -2.30102999566],
          ),
      ),
      (
          [0.9999, 0.0001, 0.0],
          _create_variant(
              'chrY',
              10000210,
              'C',
              ['T'],
              0.00043431619,
              dv_vcf_constants.DEEP_VARIANT_REF_FILTER,
              [0, 0],
              40,
              [-0.00004343161, -4.0, -9.90308995105826],
          ),
      ),
      # Multi-allelic test examples.
      (
          [0.995, 0.001, 0.001, 0.001, 0.001, 0.001],
          _create_variant(
              'X',
              10000210,
              'CACA',
              ['C', 'A'],
              0.0217691925,
              dv_vcf_constants.DEEP_VARIANT_REF_FILTER,
              [0, 0],
              23,
              [-0.00217691925, -3, -3, -3, -3, -3],
          ),
      ),
      (
          [0.001, 0.001, 0.001, 0.995, 0.001, 0.001],
          _create_variant(
              'X',
              10000210,
              'CACA',
              ['C', 'A'],
              30,
              dv_vcf_constants.DEEP_VARIANT_PASS,
              [0, 2],
              23,
              [-3, -3, -3, -0.00217691925, -3, -3],
          ),
      ),
  )
  def test_add_call_to_variant(self, probs, expected):
    raw_variant = variants_pb2.Variant(
        reference_name=expected.reference_name,
        reference_bases=expected.reference_bases,
        alternate_bases=expected.alternate_bases,
        start=expected.start,
        end=expected.end,
        calls=[variants_pb2.VariantCall(call_set_name=_DEFAULT_SAMPLE_NAME)],
    )
    variantcall_utils.set_ad(raw_variant.calls[0], [1, 1])
    variant = postprocess_variants.add_call_to_variant(
        variant=raw_variant, predictions=probs, sample_name=_DEFAULT_SAMPLE_NAME
    )
    self.assertEqual(variant.reference_bases, expected.reference_bases)
    self.assertEqual(variant.alternate_bases, expected.alternate_bases)
    self.assertEqual(variant.reference_name, expected.reference_name)
    self.assertEqual(variant.start, expected.start)
    self.assertEqual(variant.end, expected.end)
    self.assertAlmostEqual(variant.quality, expected.quality, places=6)
    self.assertEqual(variant.filter, expected.filter)
    self.assertLen(variant.calls, 1)
    self.assertLen(expected.calls, 1)
    self.assertEqual(variant.calls[0].genotype, expected.calls[0].genotype)
    self.assertEqual(variant.calls[0].info['GQ'], expected.calls[0].info['GQ'])
    for gl, expected_gl in zip(
        variant.calls[0].genotype_likelihood,
        expected.calls[0].genotype_likelihood,
    ):
      self.assertAlmostEqual(gl, expected_gl, places=6)

  @parameterized.parameters(
      (
          0,
          [0, 0],
      ),
      (
          1,
          [0, 1],
      ),
      (
          2,
          [1, 1],
      ),
      (
          3,
          [0, 2],
      ),
      (
          4,
          [1, 2],
      ),
      (
          5,
          [2, 2],
      ),
  )
  def test_triallelic_genotype_in_add_call_to_variant(
      self, highest_prob_position, expected_best_genotype
  ):
    """Ensures the order of GLs are interpreted correctly for triallelics."""
    raw_variant = _create_variant_with_alleles(ref='CACA', alts=['C', 'A'])
    # Create a probability with 6 elements, one of them 0.995 (best genotype),
    # and the rest 0.001.
    probs = [0.001] * 6
    assert 0 <= highest_prob_position <= len(probs)
    probs[highest_prob_position] = 0.995
    variantcall_utils.set_ad(raw_variant.calls[0], [1, 1, 1])
    variant = postprocess_variants.add_call_to_variant(
        variant=raw_variant, predictions=probs, sample_name=_DEFAULT_SAMPLE_NAME
    )
    self.assertEqual(variant.calls[0].genotype, expected_best_genotype)

  @parameterized.parameters(
      # Q20 tests
      ([0.01, 0.0, 0.99], 0, 0, 20.0),
      ([0.01, 0.0, 0.99], 1, 0, 20.0),
      ([0.01, 0.0, 0.99], 2, 20, 20.0),
      # Q30 tests
      ([0.001, 0.0, 0.999], 0, 0, 30.0),
      ([0.001, 0.0, 0.999], 1, 0, 30.0),
      ([0.001, 0.0, 0.999], 2, 30, 30.0),
      # Q40 tests
      ([0.0001, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.9999], 1, 0, 40.0),
      ([0.0001, 0.0, 0.9999], 2, 40, 40.0),
      # Test that this works with any sized genotype vector.
      ([0.0001, 0.0, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.0, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.0, 0.0, 0.0, 0.9999], 0, 0, 40.0),
      ([0.0001, 0.0, 0.0, 0.0, 0.0, 0.9999], 5, 40, 40.0),
      ([0.0001, 0.0, 0.0, 0.0, 0.0, 0.0, 0.9999], 0, 0, 40.0),
      # Test that probabilities more extreme than genomics_math._MAX_CONFIDENCE
      # are appropriately rounded.
      ([1e-11, 1 - 1e-11, 0.0], 0, 0, 99.03089987),
      ([1e-11, 1 - 1e-11, 0.0], 1, 99, 99.03089987),
      ([1e-11, 1 - 1e-11, 0.0], 2, 0, 99.03089987),
      ([1e-15, 1 - 1e-15, 0.0], 0, 0, 99.03089987),
      ([1e-15, 1 - 1e-15, 0.0], 1, 99, 99.03089987),
      ([1e-15, 1 - 1e-15, 0.0], 2, 0, 99.03089987),
  )
  def test_compute_quals(self, probs, call, expected_gq, expected_qual):
    gq, qual = postprocess_variants.compute_quals(probs, call)
    self.assertEqual(gq, expected_gq)
    self.assertAlmostEqual(qual, expected_qual, places=6)

  @parameterized.parameters(
      # Make sure code is robust to minor numerical issues where the sum of
      # the vector isn't exactly 1.0.
      # This vector sums to ~1.0, minus any parsing uncertainty, and will
      # return a GQ of 40 but a qual of MAX_QUAL.
      ([0.0, 0.0001, 0.9999], 2, 40),
      # This vector sums to >1.0, but we should get back a max qual value
      # instead of throwing an exception.
      ([0.0, 0.00011, 0.9999], 2, 40),
  )
  def test_compute_quals_numerical_stability(self, probs, call, expected_gq):
    max_qual = round(
        genomics_math.ptrue_to_bounded_phred(1.0),
        postprocess_variants._QUAL_PRECISION,
    )
    gq, qual = postprocess_variants.compute_quals(probs, call)
    self.assertEqual(expected_gq, gq)
    self.assertEqual(max_qual, qual)

  @parameterized.parameters(
      # Standard diploid case with 1 alt allele.
      ([1, 0, 0], 1, (0, [0, 0])),
      ([0, 1, 0], 1, (1, [0, 1])),
      ([0, 0, 1], 1, (2, [1, 1])),
      # Diploid case with 2 alt alleles.
      ([1, 0, 0, 0, 0, 0], 2, (0, [0, 0])),
      ([0, 1, 0, 0, 0, 0], 2, (1, [0, 1])),
      ([0, 0, 1, 0, 0, 0], 2, (2, [1, 1])),
      ([0, 0, 0, 1, 0, 0], 2, (3, [0, 2])),
      ([0, 0, 0, 0, 1, 0], 2, (4, [1, 2])),
      ([0, 0, 0, 0, 0, 1], 2, (5, [2, 2])),
  )
  def test_most_likely_genotype(self, probs, n_alts, expected):
    del n_alts
    self.assertEqual(expected, postprocess_variants.most_likely_genotype(probs))

  @parameterized.parameters(1, 3, 4)
  def test_most_likely_genotype_raises_with_non2_ploidy(self, bad_ploidy):
    with self.assertRaises(NotImplementedError):
      postprocess_variants.most_likely_genotype([1, 0, 0], ploidy=bad_ploidy)

  def test_compute_filter_fields(self):
    # This generates too many tests as a parameterized test.
    for qual, min_qual in itertools.product(range(100), range(100)):
      # First test with no alleleic depth.
      variant = variants_pb2.Variant()
      variant.quality = qual
      expected = []
      expected.append(dv_vcf_constants.DEEP_VARIANT_NO_CALL)
      self.assertEqual(
          postprocess_variants.compute_filter_fields(variant, min_qual),
          expected,
      )
      # Now add hom ref genotype and AD --> qual shouldn't affect filter field
      del variant.filter[:]
      variant.calls.add(genotype=[0, 0])
      variantcall_utils.set_ad(variant.calls[0], [1, 1])
      expected = []
      expected.append(dv_vcf_constants.DEEP_VARIANT_REF_FILTER)
      self.assertEqual(
          postprocess_variants.compute_filter_fields(variant, min_qual),
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
          postprocess_variants.compute_filter_fields(variant, min_qual),
          expected,
      )

  @parameterized.parameters(
      (
          [
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.01, 0.98, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT']),
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT']),
              ),
              _create_call_variants_output(
                  indices=[2],
                  probabilities=[0.01, 0.97, 0.02],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT']),
              ),
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[0.01, 0.98, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT']),
              ),
              _create_call_variants_output(
                  indices=[0, 2],
                  probabilities=[0.04, 0.95, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT']),
              ),
              _create_call_variants_output(
                  indices=[1, 2],
                  probabilities=[0.01, 0.98, 0.01],
                  variant=_create_variant_with_alleles(alts=['C', 'T', 'TT']),
              ),
          ],
          6,
          set(['T']),
      ),
      # Example where all alt alleles are below qual_filter, but we keep one
      # where the qual is highest among the ones filtered out.
      (
          [
              _create_call_variants_output(
                  indices=[0, 1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T']),
              ),
              _create_call_variants_output(
                  indices=[0],
                  probabilities=[0.99, 0.01, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T']),
              ),
              _create_call_variants_output(
                  indices=[1],
                  probabilities=[1, 0, 0],
                  variant=_create_variant_with_alleles(alts=['C', 'T']),
              ),
          ],
          6,
          set(['T']),
      ),
  )
  def test_get_alt_alleles_to_remove(
      self, call_variants_outputs, qual_filter, expected_output
  ):
    self.assertEqual(
        postprocess_variants.get_alt_alleles_to_remove(
            call_variants_outputs, qual_filter
        ),
        expected_output,
    )

  @parameterized.parameters(
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set([]),
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['C']),
          _create_variant_with_alleles(alts=['T', 'TT']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['T']),
          _create_variant_with_alleles(alts=['C', 'TT']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['TT']),
          _create_variant_with_alleles(alts=['C', 'T']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['C', 'TT']),
          _create_variant_with_alleles(alts=['T']),
      ),
      (
          _create_variant_with_alleles(alts=['C', 'T', 'TT']),
          set(['C', 'T', 'TT']),
          _create_variant_with_alleles(alts=[]),
      ),
  )
  def test_prune_alleles(
      self, canonical_variant, alt_alleles_to_remove, expected_variant
  ):
    self.assertEqual(
        postprocess_variants.prune_alleles(
            canonical_variant, alt_alleles_to_remove
        ),
        expected_variant,
    )

  @parameterized.parameters(
      # Check that we are simplifying alleles and that the simplification deps
      # on the alleles we've removed.
      dict(
          alleles=['CAA', 'CA', 'C'],
          start=5,
          alt_alleles_to_remove=[],
          expected_alleles=['CAA', 'CA', 'C'],
          expected_end=8,
      ),
      # Removing the C allele allows us to simplify CAA + CA => CA + C.
      dict(
          alleles=['CAA', 'CA', 'C'],
          start=4,
          alt_alleles_to_remove=['C'],
          expected_alleles=['CA', 'C'],
          expected_end=6,
      ),
      # Removing the CA allele doesn't allow any simplification.
      dict(
          alleles=['CAA', 'CA', 'C'],
          start=3,
          alt_alleles_to_remove=['CA'],
          expected_alleles=['CAA', 'C'],
          expected_end=6,
      ),
      # Make sure we keep at least one anchor base when pruning.
      dict(
          alleles=['CCA', 'CA', 'T'],
          start=2,
          alt_alleles_to_remove=['T'],
          expected_alleles=['CC', 'C'],
          expected_end=4,
      ),
  )
  def test_prune_and_simplify_alleles(
      self,
      alleles,
      start,
      alt_alleles_to_remove,
      expected_alleles,
      expected_end,
  ):
    """Test that prune_alleles + simplify_variant_alleles works as expected."""
    variant = _create_variant_with_alleles(
        ref=alleles[0], alts=alleles[1:], start=start
    )
    pruned = postprocess_variants.prune_alleles(variant, alt_alleles_to_remove)
    simplified = variant_utils.simplify_variant_alleles(pruned)
    self.assertEqual(simplified.reference_bases, expected_alleles[0])
    self.assertEqual(simplified.alternate_bases, expected_alleles[1:])
    self.assertEqual(simplified.start, start)
    self.assertEqual(simplified.end, expected_end)

  def test_merge_predictions_simplifies_alleles(self):
    """Checks that merge_predictions simplifies alleles."""
    ref, alts = 'CCA', ['CA', 'C']
    inputs = [
        _create_call_variants_output(
            ref=ref, indices=[0], probabilities=[0.0, 1.0, 0.0], alts=alts
        ),
        _create_call_variants_output(
            ref=ref, indices=[1], probabilities=[1.0, 0.0, 0.0], alts=alts
        ),
        _create_call_variants_output(
            ref=ref, indices=[0, 1], probabilities=[0.0, 1.0, 0.0], alts=alts
        ),
    ]

    for permuted_inputs in itertools.permutations(inputs):
      # qual_filter=2 is needed so we remove our middle 'C' allele.
      variant, probs = postprocess_variants.merge_predictions(
          permuted_inputs, qual_filter=2
      )
      np.testing.assert_almost_equal(probs, [0.0, 1.0, 0.0])
      self.assertEqual(variant.reference_bases, 'CC')
      self.assertEqual(variant.alternate_bases, ['C'])

  @parameterized.parameters(
      (['A'], {}, [1, 2], [1, 2]),
      (['A'], {'A'}, [1, 2], [1]),
      (['A', 'C'], {}, [1, 2, 3], [1, 2, 3]),
      (['A', 'C'], {'A'}, [1, 2, 3], [1, 3]),
      (['A', 'C'], {'C'}, [1, 2, 3], [1, 2]),
      (['A', 'C'], {'A', 'C'}, [1, 2, 3], [1]),
  )
  def test_prune_alleles_handles_format_fields(
      self, alts, to_remove, orig_ad, expected_ad
  ):
    variant = _create_variant_with_alleles(alts=alts)
    test_utils.set_list_values(variant.calls[0].info['AD'], orig_ad)
    actual = postprocess_variants.prune_alleles(variant, to_remove)
    self.assertEqual(
        [v.int_value for v in actual.calls[0].info['AD'].values], expected_ad
    )

  @parameterized.parameters(
      (1, [[0]]),
      (2, [[0], [0, 1], [1]]),
      (3, [[0], [0, 1], [0, 2], [1], [1, 2], [2]]),
      (4, [[0], [0, 1], [0, 2], [0, 3], [1], [1, 2], [1, 3], [2], [2, 3], [3]]),
      (
          8,
          [
              [0],
              [0, 1],
              [0, 2],
              [0, 3],
              [0, 4],
              [0, 5],
              [0, 6],
              [0, 7],
              [1],
              [1, 2],
              [1, 3],
              [1, 4],
              [1, 5],
              [1, 6],
              [1, 7],
              [2],
              [2, 3],
              [2, 4],
              [2, 5],
              [2, 6],
              [2, 7],
              [3],
              [3, 4],
              [3, 5],
              [3, 6],
              [3, 7],
              [4],
              [4, 5],
              [4, 6],
              [4, 7],
              [5],
              [5, 6],
              [5, 7],
              [6],
              [6, 7],
              [7],
          ],
      ),
  )
  def test_expected_alt_allele_indices(
      self, num_alternate_bases, expected_indices
  ):
    self.assertEqual(
        postprocess_variants.expected_alt_allele_indices(num_alternate_bases),
        expected_indices,
    )

  @flagsaver.flagsaver
  def test_catches_bad_argv(self):
    # Define valid flags to ensure raise occurs due to argv issues.
    FLAGS.infile = make_golden_dataset(False)
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.outfile = test_utils.test_tmpfile('nonempty_outfile.vcf')
    with mock.patch.object(logging, 'error') as mock_logging,\
        mock.patch.object(sys, 'exit') as mock_exit:
      postprocess_variants.main(['postprocess_variants.py', 'extra_arg'])
    mock_logging.assert_called_once_with(
        'Command line parsing failure: postprocess_variants does not accept '
        'positional arguments but some are present on the command line: '
        '"[\'postprocess_variants.py\', \'extra_arg\']".')
    mock_exit.assert_called_once_with(errno.ENOENT)

  @flagsaver.flagsaver
  def test_catches_bad_flags(self):
    FLAGS.infile = make_golden_dataset(False)
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.outfile = 'nonempty_outfile.vcf'
    FLAGS.nonvariant_site_tfrecord_path = testdata.GOLDEN_POSTPROCESS_GVCF_INPUT
    # This is the bad flag.
    FLAGS.gvcf_outfile = ''
    with mock.patch.object(logging, 'error') as mock_logging, \
        mock.patch.object(sys, 'exit') as mock_exit:
      postprocess_variants.main(['postprocess_variants.py'])
    mock_logging.assert_called_once_with(
        'gVCF creation requires both nonvariant_site_tfrecord_path and '
        'gvcf_outfile flags to be set.')
    mock_exit.assert_called_once_with(errno.ENOENT)

  @parameterized.parameters(
      dict(
          call=_create_call_variants_output(
              indices=[0],
              probabilities=[0.02, 0.98, 0],
              variant=_create_variant_with_alleles(ref='CA', alts=['C']),
          ),
          expected_probabilities=[1.0, 0, 0],
      ),
      dict(
          call=_create_call_variants_output(
              indices=[0],
              probabilities=[0.98, 0.02, 0],
              variant=_create_variant_with_alleles(ref='CA', alts=['C']),
          ),
          expected_probabilities=[1.0, 0, 0],
      ),
      dict(
          call=_create_call_variants_output(
              indices=[0],
              probabilities=[0.2, 0.5, 0.3],
              variant=_create_variant_with_alleles(ref='CA', alts=['C']),
          ),
          expected_probabilities=[0.4, 0, 0.6],
      ),
      dict(
          call=_create_call_variants_output(
              indices=[0],
              probabilities=[0.0, 1.0, 0.0],
              variant=_create_variant_with_alleles(ref='CA', alts=['C']),
          ),
          expected_probabilities=[0, 0, 0],
      ),
      dict(
          call=_create_call_variants_output(
              indices=[0, 1],
              probabilities=[0.02, 0.03, 0.45, 0.07, 0.3, 0.13],
              variant=_create_variant_with_alleles(ref='CA', alts=['C', 'CAA']),
          ),
          expected_probabilities=[0.033, 0, 0.75, 0, 0, 0.216],
      ),
  )
  def test_correct_nonautosome_probabilities(
      self, call, expected_probabilities
  ):
    self.assertSequenceAlmostEqual(
        expected_probabilities,
        postprocess_variants.correct_nonautosome_probabilities(
            call.genotype_probabilities, call.variant
        ),
        places=2,
    )


class MergeVcfAndGvcfTest(parameterized.TestCase):

  # TODO use itertools.permutations to improve the test.
  def test_sort_grouped_variants(self):
    group = [
        _create_call_variants_output(
            indices=[0, 1],
            probabilities=[0.98, 0.02, 0],
            variant=_create_variant_with_alleles(ref='CA', alts=['C', 'CAA']),
        ),
        _create_call_variants_output(
            indices=[1],
            probabilities=[0.995, 0.005, 0],
            variant=_create_variant_with_alleles(ref='CA', alts=['C', 'CAA']),
        ),
        _create_call_variants_output(
            indices=[0],
            probabilities=[0.95, 0.05, 0],
            variant=_create_variant_with_alleles(ref='CA', alts=['C', 'CAA']),
        ),
    ]
    output = postprocess_variants._sort_grouped_variants(group)
    # In sorted output, 1st has indices=[0].
    self.assertEqual(output[0], group[2])
    self.assertEqual(output[0].alt_allele_indices.indices, [0])
    # In sorted output, 2nd has indices=[0, 1].
    self.assertEqual(output[1], group[0])
    self.assertEqual(output[1].alt_allele_indices.indices, [0, 1])
    # In sorted output, 3rd has indices=[1].
    self.assertEqual(output[2], group[1])
    self.assertEqual(output[2].alt_allele_indices.indices, [1])


class PartitionContigsTest(parameterized.TestCase):
  # Total bps: 2100
  CONTIGS = [
      reference_pb2.ContigInfo(
          name='chr1',
          n_bases=1000,
          pos_in_fasta=0,
      ),
      reference_pb2.ContigInfo(
          name='chr2',
          n_bases=500,
          pos_in_fasta=1,
      ),
      reference_pb2.ContigInfo(
          name='chr3',
          n_bases=300,
          pos_in_fasta=2,
      ),
      reference_pb2.ContigInfo(
          name='chr4',
          n_bases=200,
          pos_in_fasta=3,
      ),
      reference_pb2.ContigInfo(
          name='chr5',
          n_bases=100,
          pos_in_fasta=4,
      ),
  ]

  @parameterized.parameters(
      dict(
          num_partitions=1,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=1000),
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=2,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=1000),
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=3,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=700),
                  range_pb2.Range(reference_name='chr1', start=700, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=4,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=525),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=525, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=500),
                  range_pb2.Range(reference_name='chr3', start=0, end=300),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=10,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=210),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=210, end=420),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=420, end=630),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=630, end=840),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=840, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=210),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=210, end=420),
                  range_pb2.Range(reference_name='chr2', start=420, end=500),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=0, end=210),
                  range_pb2.Range(reference_name='chr3', start=210, end=300),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=0, end=200),
              ],
              [
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
      dict(
          num_partitions=11,
          expected_partition_groups=[
              [
                  range_pb2.Range(reference_name='chr1', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=190, end=380),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=380, end=570),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=570, end=760),
              ],
              [
                  range_pb2.Range(reference_name='chr1', start=760, end=950),
                  range_pb2.Range(reference_name='chr1', start=950, end=1000),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=190, end=380),
              ],
              [
                  range_pb2.Range(reference_name='chr2', start=380, end=500),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr3', start=190, end=300),
                  range_pb2.Range(reference_name='chr4', start=0, end=190),
              ],
              [
                  range_pb2.Range(reference_name='chr4', start=190, end=200),
                  range_pb2.Range(reference_name='chr5', start=0, end=100),
              ],
          ],
      ),
  )
  def test_split_into_contiguous_partitions(
      self, num_partitions, expected_partition_groups
  ):
    partition_groups = postprocess_variants.split_into_contiguous_partitions(
        self.CONTIGS, num_partitions=num_partitions
    )
    self.assertLen(partition_groups, num_partitions)
    self.assertEqual(partition_groups, expected_partition_groups)


if __name__ == '__main__':
  absltest.main()
