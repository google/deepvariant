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
"""Tests for deepvariant.make_examples."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import copy
import errno
import sys



from tensorflow import flags
from absl.testing import absltest
from absl.testing import parameterized
import mock

from absl import logging

from third_party.nucleus.io import vcf
from third_party.nucleus.protos import reads_pb2
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import io_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils
from third_party.nucleus.util import vcf_constants
from deepvariant import make_examples
from deepvariant import testdata
from deepvariant import tf_utils
from deepvariant.labeler import variant_labeler
from deepvariant.protos import deepvariant_pb2
from deepvariant.protos import realigner_pb2
from deepvariant.testing import flagsaver

FLAGS = flags.FLAGS

# Dictionary mapping keys to decoders for decode_example function.
_EXAMPLE_DECODERS = {
    'locus': tf_utils.example_locus,
    'alt_allele_indices/encoded': tf_utils.example_alt_alleles_indices,
    'image/encoded': tf_utils.example_encoded_image,
    'variant/encoded': tf_utils.example_variant,
    'label': tf_utils.example_label,
    'image/format': tf_utils.example_image_format,
    'image/shape': tf_utils.example_image_shape,
}


def decode_example(example):
  """Decodes a tf.Example from DeepVariant into a dict of Pythonic structures.

  Args:
    example: tf.Example proto. The example to make into a dictionary.

  Returns:
    A python dictionary with key/value pairs for each of the fields of example,
    with each value decoded as needed into Python structures like protos, list,
    etc.

  Raises:
    KeyError: If example contains a feature without a known decoder.
  """
  as_dict = {}
  for key in example.features.feature:
    if key not in _EXAMPLE_DECODERS:
      raise KeyError('Unexpected example key', key)
    as_dict[key] = _EXAMPLE_DECODERS[key](example)
  return as_dict


def setUpModule():
  testdata.init()


def _make_contigs(specs):
  """Makes ContigInfo protos from specs.

  Args:
    specs: A list of 2- or 3-tuples. All tuples should be of the same length.
      If 2-element, these should be the name and length in basepairs of each
      contig, and their pos_in_fasta will be set to their index in the list. If
      the 3-element, the tuple should contain name, length, and pos_in_fasta.

  Returns:
    A list of ContigInfo protos, one for each spec in specs.
  """
  if specs and len(specs[0]) == 3:
    return [
        reference_pb2.ContigInfo(name=name, n_bases=length, pos_in_fasta=i)
        for name, length, i in specs
    ]
  else:
    return [
        reference_pb2.ContigInfo(name=name, n_bases=length, pos_in_fasta=i)
        for i, (name, length) in enumerate(specs)
    ]


def _from_literals_list(literals, contig_map=None):
  """Makes a list of Range objects from literals."""
  return ranges.parse_literals(literals, contig_map)


def _from_literals(literals, contig_map=None):
  """Makes a RangeSet of intervals from literals."""
  return ranges.RangeSet.from_regions(literals, contig_map)


def _sharded(basename, num_shards):
  if num_shards:
    return basename + '@' + str(num_shards)
  else:
    return basename


@parameterized.parameters(
    dict(mode='calling', num_shards=0),
    dict(mode='calling', num_shards=3),
    dict(mode='training', num_shards=0, labeler_algorithm='positional_labeler'),
    dict(mode='training', num_shards=0, labeler_algorithm='haplotype_labeler'),
    dict(mode='training', num_shards=3, labeler_algorithm='haplotype_labeler'),
)
class MakeExamplesEnd2EndTest(parameterized.TestCase):

  @flagsaver.FlagSaver
  def test_make_examples_end2end(self, mode, num_shards,
                                 labeler_algorithm=None):
    self.maxDiff = None
    self.assertIn(mode, {'calling', 'training'})
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.candidates = test_utils.test_tmpfile(
        _sharded('vsc.tfrecord', num_shards))
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('examples.tfrecord', num_shards))
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.partition_size = 1000
    FLAGS.mode = mode
    FLAGS.gvcf_gq_binsize = 5
    if labeler_algorithm is not None:
      FLAGS.labeler_algorithm = labeler_algorithm

    if mode == 'calling':
      FLAGS.gvcf = test_utils.test_tmpfile(
          _sharded('gvcf.tfrecord', num_shards))
    else:
      FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
      FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED

    for task_id in range(max(num_shards, 1)):
      FLAGS.task = task_id
      options = make_examples.default_options(add_flags=True)
      make_examples.make_examples_runner(options)

    # Test that our candidates are reasonable, calling specific helper functions
    # to check lots of properties of the output.
    candidates = sorted(
        io_utils.read_tfrecords(
            FLAGS.candidates, proto=deepvariant_pb2.DeepVariantCall),
        key=lambda c: variant_utils.variant_range_tuple(c.variant))
    self.verify_deepvariant_calls(candidates, options)
    self.verify_variants(
        [call.variant for call in candidates], region, options, is_gvcf=False)

    # Verify that the variants in the examples are all good.
    examples = self.verify_examples(
        FLAGS.examples, region, options, verify_labels=mode == 'training')
    example_variants = [tf_utils.example_variant(ex) for ex in examples]
    self.verify_variants(example_variants, region, options, is_gvcf=False)

    # Verify the integrity of the examples and then check that they match our
    # golden labeled examples. Note we expect the order for both training and
    # calling modes to produce deterministic order because we fix the random
    # seed.
    if mode == 'calling':
      golden_file = _sharded(testdata.GOLDEN_CALLING_EXAMPLES, num_shards)
    else:
      golden_file = _sharded(testdata.GOLDEN_TRAINING_EXAMPLES, num_shards)
    self.assertDeepVariantExamplesEqual(
        examples, list(io_utils.read_tfrecords(golden_file)))

    if mode == 'calling':
      nist_reader = vcf.VcfReader(testdata.TRUTH_VARIANTS_VCF)
      nist_variants = list(nist_reader.query(region))
      self.verify_nist_concordance(example_variants, nist_variants)

      # Check the quality of our generated gvcf file.
      gvcfs = variant_utils.sorted_variants(
          io_utils.read_tfrecords(FLAGS.gvcf, proto=variants_pb2.Variant))
      self.verify_variants(gvcfs, region, options, is_gvcf=True)
      self.verify_contiguity(gvcfs, region)
      gvcf_golden_file = _sharded(testdata.GOLDEN_POSTPROCESS_GVCF_INPUT,
                                  num_shards)
      expected_gvcfs = list(
          io_utils.read_tfrecords(gvcf_golden_file, proto=variants_pb2.Variant))
      self.assertItemsEqual(gvcfs, expected_gvcfs)

  def verify_nist_concordance(self, candidates, nist_variants):
    # Tests that we call all of the real variants (according to NIST's Genome
    # in a Bottle callset for NA12878) in our candidate callset.
    # Tests that we don't have an enormous number of FP calls. We should have
    # no more than 5x (arbitrary) more candidate calls than real calls. If we
    # have more it's likely due to some major pipeline problem.
    self.assertLess(len(candidates), 5 * len(nist_variants))
    for nist_variant in nist_variants:
      self.assertVariantIsPresent(nist_variant, candidates)

  def assertDeepVariantExamplesEqual(self, actual, expected):
    """Asserts that actual and expected tf.Examples from DeepVariant are equal.

    Args:
      actual: iterable of tf.Examples from DeepVariant. DeepVariant examples
        that we want to check.
      expected: iterable of tf.Examples. Expected results for actual.
    """
    self.assertEqual(len(actual), len(expected))
    for i in range(len(actual)):
      self.assertEqual(decode_example(actual[i]), decode_example(expected[i]))

  def assertVariantIsPresent(self, to_find, variants):

    def variant_key(v):
      return (v.reference_bases, v.start, v.end, v.reference_bases)

    # Finds a call in our actual call set for each NIST variant, asserting
    # that we found exactly one.
    matches = [
        variant for variant in variants
        if variant_key(to_find) == variant_key(variant)
    ]
    self.assertEqual(len(matches), 1, 'to_find failed {}'.format(to_find))

    # Verify that every alt allele appears in the call (but the call might)
    # have more than just those calls.
    for alt in matches[0].alternate_bases:
      self.assertIn(alt, to_find.alternate_bases)

  def verify_variants(self, variants, region, options, is_gvcf):
    # Verifies simple properties of the Variant protos in variants. For example,
    # checks that the reference_name() is our expected chromosome. The flag
    # is_gvcf determines how we check the VariantCall field of each variant,
    # enforcing expectations for gVCF records if true or variant calls if false.
    for variant in variants:
      self.assertEqual(variant.reference_name, region.reference_name)
      self.assertNotEqual(variant.reference_bases, '')
      self.assertGreater(len(variant.alternate_bases), 0)
      self.assertGreaterEqual(variant.start, region.start)
      self.assertLessEqual(variant.start, region.end)
      self.assertEqual(len(variant.calls), 1)

      call = variant_utils.only_call(variant)
      self.assertEqual(call.call_set_name,
                       options.variant_caller_options.sample_name)
      if is_gvcf:
        # GVCF records should have 0/0 genotypes as they are reference sites,
        # have genotype likelihoods and a GQ value.
        self.assertEqual(call.genotype, [0, 0])
        self.assertEqual(len(call.genotype_likelihood), 3)
        self.assertGreaterEqual(variantcall_utils.get_gq(call), 0)

  def verify_contiguity(self, contiguous_variants, region):
    """Verifies region is fully covered by gvcf records."""
    # We expect that the intervals cover every base, so the first variant should
    # be at our interval start and the last one should end at our interval end.
    self.assertGreater(len(contiguous_variants), 0)
    self.assertEqual(region.start, contiguous_variants[0].start)
    self.assertEqual(region.end, contiguous_variants[-1].end)

    # After this loop completes successfully we know that together the GVCF and
    # Variants form a fully contiguous cover of our calling interval, as
    # expected.
    for v1, v2 in zip(contiguous_variants, contiguous_variants[1:]):
      # Sequential variants should be contiguous, meaning that v2.start should
      # be v1's end, as the end is exclusive and the start is inclusive.
      if v1.start == v2.start and v1.end == v2.end:
        # Skip duplicates here as we may have multi-allelic variants turning
        # into multiple bi-allelic variants at the same site.
        continue
      # We expect to immediately follow the end of a gvcf record but to occur
      # at the base immediately after a variant, since the variant's end can
      # span over a larger interval when its a deletion and we still produce
      # gvcf records under the deletion.
      expected_start = v1.end if v1.alternate_bases == ['<*>'] else v1.start + 1
      self.assertEqual(v2.start, expected_start)

  def verify_deepvariant_calls(self, dv_calls, options):
    # Verifies simple structural properties of the DeepVariantCall objects
    # emitted by the VerySensitiveCaller, such as that the AlleleCount and
    # Variant both have the same position.
    for call in dv_calls:
      for alt_allele in call.variant.alternate_bases:
        # Skip ref calls.
        if alt_allele == vcf_constants.NO_ALT_ALLELE:
          continue
        # Make sure allele appears in our allele_support field and that at
        # least our min number of reads to call an alt allele are present in
        # the supporting reads list for that allele.
        self.assertIn(alt_allele, list(call.allele_support))
        self.assertGreaterEqual(
            len(call.allele_support[alt_allele].read_names),
            options.variant_caller_options.min_count_snps)

  def verify_examples(self, examples_filename, region, options, verify_labels):
    # Do some simple structural checks on the tf.Examples in the file.
    expected_features = [
        'variant/encoded', 'locus', 'image/format', 'image/encoded',
        'alt_allele_indices/encoded'
    ]
    if verify_labels:
      expected_features += ['label']

    examples = list(io_utils.read_tfrecords(examples_filename))
    for example in examples:
      for label_feature in expected_features:
        self.assertIn(label_feature, example.features.feature)
      # pylint: disable=g-explicit-length-test
      self.assertGreater(len(tf_utils.example_alt_alleles_indices(example)), 0)

    # Check that the variants in the examples are good.
    variants = [tf_utils.example_variant(x) for x in examples]
    self.verify_variants(variants, region, options, is_gvcf=False)

    return examples


class MakeExamplesUnitTest(parameterized.TestCase):

  @parameterized.parameters(
      dict(
          flag_value='CALLING',
          expected=deepvariant_pb2.DeepVariantOptions.CALLING,
      ),
      dict(
          flag_value='TRAINING',
          expected=deepvariant_pb2.DeepVariantOptions.TRAINING,
      ),
  )
  def test_parse_proto_enum_flag(self, flag_value, expected):
    enum_pb2 = deepvariant_pb2.DeepVariantOptions.Mode
    self.assertEqual(
        make_examples.parse_proto_enum_flag(enum_pb2, flag_value), expected)

  def test_parse_proto_enum_flag_error_handling(self):
    with self.assertRaisesRegexp(
        ValueError,
        'Unknown enum option "foo". Allowed options are CALLING,TRAINING'):
      make_examples.parse_proto_enum_flag(
          deepvariant_pb2.DeepVariantOptions.Mode, 'foo')

  @flagsaver.FlagSaver
  def test_default_options_with_training_random_emit_ref_sites(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''

    FLAGS.training_random_emit_ref_sites = 0.3
    options = make_examples.default_options(add_flags=True)
    self.assertAlmostEqual(
        options.variant_caller_options.fraction_reference_sites_to_emit, 0.3)

  @flagsaver.FlagSaver
  def test_default_options_without_training_random_emit_ref_sites(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''

    options = make_examples.default_options(add_flags=True)
    # In proto3, there is no way to check presence of scalar field:
    # redacted
    # As an approximation, we directly check that the value should be exactly 0.
    self.assertEqual(
        options.variant_caller_options.fraction_reference_sites_to_emit, 0.0)

  def test_extract_sample_name_from_reads_single_sample(self):
    mock_sample_reader = mock.Mock()
    mock_sample_reader.header = reads_pb2.SamHeader(
        read_groups=[reads_pb2.ReadGroup(sample_id='sample_name')])
    self.assertEqual(
        make_examples.extract_sample_name_from_sam_reader(mock_sample_reader),
        'sample_name')

  @parameterized.parameters(
      # No samples could be found in the reads.
      dict(
          samples=[],
          expected_error_message='No non-empty sample name found in the input '
          'reads. Please provide the name of the sample with the --sample_name '
          'argument.'),
      # Check that we detect an empty sample name and raise an exception.
      dict(
          samples=[''],
          expected_error_message='No non-empty sample name found in the input '
          'reads. Please provide the name of the sample '
          'with the --sample_name argument.'),
      # We have more than one sample in the reads.
      dict(
          samples=['sample1', 'sample2'],
          expected_error_message=
          r'Multiple samples \(sample1, sample2\) were found in the input '
          'reads. DeepVariant can only call variants from a BAM file '
          'containing a single sample.'),
  )
  def test_extract_sample_name_from_reads_detects_bad_samples(
      self, samples, expected_error_message):
    mock_sample_reader = mock.Mock()
    mock_sample_reader.header = reads_pb2.SamHeader(read_groups=[
        reads_pb2.ReadGroup(sample_id=sample) for sample in samples
    ])
    with self.assertRaisesRegexp(ValueError, expected_error_message):
      make_examples.extract_sample_name_from_sam_reader(mock_sample_reader)

  @flagsaver.FlagSaver
  def test_confident_regions(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''

    options = make_examples.default_options(add_flags=True)
    confident_regions = make_examples.read_confident_regions(options)

    # Our expected intervals, inlined from CONFIDENT_REGIONS_BED.
    expected = _from_literals_list([
        'chr20:10000847-10002407', 'chr20:10002521-10004171',
        'chr20:10004274-10004964', 'chr20:10004995-10006386',
        'chr20:10006410-10007800', 'chr20:10007825-10008018',
        'chr20:10008044-10008079', 'chr20:10008101-10008707',
        'chr20:10008809-10008897', 'chr20:10009003-10009791',
        'chr20:10009934-10010531'
    ])
    # Our confident regions should be exactly those found in the BED file.
    self.assertCountEqual(expected, list(confident_regions))

  @parameterized.parameters(
      ({
          'examples': ('foo', 'foo')
      },),
      ({
          'examples': ('foo', 'foo'),
          'gvcf': ('bar', 'bar')
      },),
      ({
          'examples': ('foo@10', 'foo-00000-of-00010')
      },),
      ({
          'task': (0, 0),
          'examples': ('foo@10', 'foo-00000-of-00010')
      },),
      ({
          'task': (1, 1),
          'examples': ('foo@10', 'foo-00001-of-00010')
      },),
      ({
          'task': (1, 1),
          'examples': ('foo@10', 'foo-00001-of-00010'),
          'gvcf': ('bar@10', 'bar-00001-of-00010')
      },),
      ({
          'task': (1, 1),
          'examples': ('foo@10', 'foo-00001-of-00010'),
          'gvcf': ('bar@10', 'bar-00001-of-00010'),
          'candidates': ('baz@10', 'baz-00001-of-00010')
      },),
  )
  @flagsaver.FlagSaver
  def test_sharded_outputs1(self, settings):
    # Set all of the requested flag values.
    for name, (flag_val, _) in settings.items():
      setattr(FLAGS, name, flag_val)

    FLAGS.mode = 'training'
    FLAGS.reads = ''
    FLAGS.ref = ''
    options = make_examples.default_options(add_flags=True)

    # Check all of the flags.
    for name, option_val in [('examples', options.examples_filename),
                             ('candidates', options.candidates_filename),
                             ('gvcf', options.gvcf_filename)]:
      expected = settings[name][1] if name in settings else ''
      self.assertEqual(expected, option_val)

  @flagsaver.FlagSaver
  def test_gvcf_output_enabled_is_false_without_gvcf_flag(self):
    FLAGS.mode = 'training'
    FLAGS.gvcf = ''
    FLAGS.reads = ''
    FLAGS.ref = ''
    FLAGS.examples = ''
    options = make_examples.default_options(add_flags=True)
    self.assertFalse(make_examples.gvcf_output_enabled(options))

  @flagsaver.FlagSaver
  def test_gvcf_output_enabled_is_true_with_gvcf_flag(self):
    FLAGS.mode = 'training'
    FLAGS.gvcf = '/tmp/foo.vcf'
    FLAGS.reads = ''
    FLAGS.ref = ''
    FLAGS.examples = ''
    options = make_examples.default_options(add_flags=True)
    self.assertTrue(make_examples.gvcf_output_enabled(options))

  def test_validate_ref_contig_coverage(self):
    ref_contigs = _make_contigs([('1', 100), ('2', 100)])

    # Fully covered reference contigs don't trigger an error.
    for threshold in [0.5, 0.9, 1.0]:
      self.assertIsNone(
          make_examples.validate_reference_contig_coverage(
              ref_contigs, ref_contigs, threshold))

    # No common contigs always blows up.
    for threshold in [0.0, 0.1, 0.5, 0.9, 1.0]:
      with self.assertRaisesRegexp(ValueError, 'span 200'):
        make_examples.validate_reference_contig_coverage(
            ref_contigs, [], threshold)

    # Dropping either contig brings up below our 0.9 threshold.
    with self.assertRaisesRegexp(ValueError, 'span 200'):
      make_examples.validate_reference_contig_coverage(ref_contigs,
                                                       _make_contigs(
                                                           [('1', 100)]), 0.9)

    with self.assertRaisesRegexp(ValueError, 'span 200'):
      make_examples.validate_reference_contig_coverage(ref_contigs,
                                                       _make_contigs(
                                                           [('2', 100)]), 0.9)

    # Our actual overlap is 50%, so check that we raise when appropriate.
    with self.assertRaisesRegexp(ValueError, 'span 200'):
      make_examples.validate_reference_contig_coverage(ref_contigs,
                                                       _make_contigs(
                                                           [('2', 100)]), 0.6)
    self.assertIsNone(
        make_examples.validate_reference_contig_coverage(
            ref_contigs, _make_contigs([('2', 100)]), 0.4))

  @parameterized.parameters(
      # all intervals are shared.
      ([[('chrM', 10)], [('chrM', 10)]], [('chrM', 10)]),
      # No common intervals.
      ([[('chrM', 10)], [('chr1', 10)]], []),
      # The names are the same but sizes are different, so not common.
      ([[('chrM', 10)], [('chrM', 20)]], []),
      # One common interval and one not.
      ([[('chrM', 10), ('chr1', 20)], [('chrM', 10),
                                       ('chr2', 30)]], [('chrM', 10)]),
      # Check that the order doesn't matter.
      ([[('chr1', 20), ('chrM', 10)], [('chrM', 10),
                                       ('chr2', 30)]], [('chrM', 10, 1)]),
      # Three-way merges.
      ([
          [('chr1', 20), ('chrM', 10)],
          [('chrM', 10), ('chr2', 30)],
          [('chr2', 30), ('chr3', 30)],
      ], []),
  )
  def test_common_contigs(self, contigs_list, expected):
    self.assertEqual(
        _make_contigs(expected),
        make_examples.common_contigs(
            [_make_contigs(contigs) for contigs in contigs_list]))

  @parameterized.parameters(
      # Note that these tests aren't so comprehensive as we are trusting that
      # the intersection code logic itself is good and well-tested elsewhere.
      # Here we are focusing on some basic tests and handling of missing
      # calling_region and confident_region data.
      (['1:1-10'], ['1:1-10']),
      (['1:1-100'], ['1:1-100']),
      (['1:50-150'], ['1:50-100']),
      (None, ['1:1-100', '2:1-200']),
      (['1:20-50'], ['1:20-50']),
      # Chr3 isn't part of our contigs; make sure we tolerate it.
      (['1:20-30', '1:40-60', '3:10-50'], ['1:20-30', '1:40-60']),
      # Check that we handle overlapping calling or confident regions.
      (['1:25-30', '1:20-40'], ['1:20-40']),
  )
  def test_regions_to_process(self, calling_regions, expected):
    contigs = _make_contigs([('1', 100), ('2', 200)])
    self.assertCountEqual(
        _from_literals_list(expected),
        make_examples.regions_to_process(
            contigs, 1000, calling_regions=_from_literals(calling_regions)))

  @parameterized.parameters(
      (50, None, [
          '1:1-50', '1:51-100', '2:1-50', '2:51-76', '3:1-50', '3:51-100',
          '3:101-121'
      ]),
      (120, None, ['1:1-100', '2:1-76', '3:1-120', '3:121']),
      (500, None, ['1:1-100', '2:1-76', '3:1-121']),
      (10, ['1:1-20', '1:30-35'], ['1:1-10', '1:11-20', '1:30-35']),
      (8, ['1:1-20', '1:30-35'], ['1:1-8', '1:9-16', '1:17-20', '1:30-35']),
  )
  def test_regions_to_process_partition(self, max_size, calling_regions,
                                        expected):
    contigs = _make_contigs([('1', 100), ('2', 76), ('3', 121)])
    self.assertCountEqual(
        _from_literals_list(expected),
        make_examples.regions_to_process(
            contigs, max_size, calling_regions=_from_literals(calling_regions)))

  @parameterized.parameters(
      dict(includes=[], excludes=[], expected=['1:1-100', '2:1-200']),
      dict(includes=['1'], excludes=[], expected=['1:1-100']),
      # Check that excludes work as expected.
      dict(includes=[], excludes=['1'], expected=['2:1-200']),
      dict(includes=[], excludes=['2'], expected=['1:1-100']),
      dict(includes=[], excludes=['1', '2'], expected=[]),
      # Check that excluding pieces works. The main checks on taking the
      # difference between two RangeSets live in ranges.py so here we are just
      # making sure some basic logic works.
      dict(includes=['1'], excludes=['1:1-10'], expected=['1:11-100']),
      # Check that includes and excludes work together.
      dict(
          includes=['1', '2'],
          excludes=['1:5-10', '1:20-50', '2:10-20'],
          expected=['1:1-4', '1:11-19', '1:51-100', '2:1-9', '2:21-200']),
      dict(
          includes=['1'],
          excludes=['1:5-10', '1:20-50', '2:10-20'],
          expected=['1:1-4', '1:11-19', '1:51-100']),
      dict(
          includes=['2'],
          excludes=['1:5-10', '1:20-50', '2:10-20'],
          expected=['2:1-9', '2:21-200']),
      # A complex example of including and excluding.
      dict(
          includes=['1:10-20', '2:50-60', '2:70-80'],
          excludes=['1:1-13', '1:19-50', '2:10-65'],
          expected=['1:14-18', '2:70-80']),
  )
  def test_build_calling_regions(self, includes, excludes, expected):
    contigs = _make_contigs([('1', 100), ('2', 200)])
    actual = make_examples.build_calling_regions(contigs, includes, excludes)
    self.assertCountEqual(actual, _from_literals_list(expected))

  def test_regions_to_process_sorted_within_contig(self):
    # These regions are out of order but within a single contig.
    contigs = _make_contigs([('z', 100)])
    in_regions = _from_literals(['z:15', 'z:20', 'z:6', 'z:25-30', 'z:3-4'])
    sorted_regions = _from_literals_list(
        ['z:3-4', 'z:6', 'z:15', 'z:20', 'z:25-30'])
    actual_regions = list(
        make_examples.regions_to_process(
            contigs, 100, calling_regions=in_regions))
    # The assertEqual here is checking the order is exactly what we expect.
    self.assertEqual(sorted_regions, actual_regions)

  def test_regions_to_process_sorted_contigs(self):
    # These contig names are out of order lexicographically.
    contigs = _make_contigs([('z', 100), ('a', 100), ('n', 100)])
    in_regions = _from_literals(['a:10', 'n:1', 'z:20', 'z:5'])
    sorted_regions = _from_literals_list(['z:5', 'z:20', 'a:10', 'n:1'])
    actual_regions = list(
        make_examples.regions_to_process(
            contigs, 100, calling_regions=in_regions))
    # The assertEqual here is checking the order is exactly what we expect.
    self.assertEqual(sorted_regions, actual_regions)

  @parameterized.parameters([2, 3, 4, 5, 50])
  def test_regions_to_process_sharding(self, num_shards):
    """Makes sure we deterministically split up regions."""

    def get_regions(task_id, num_shards):
      return make_examples.regions_to_process(
          contigs=_make_contigs([('z', 100), ('a', 100), ('n', 100)]),
          partition_size=5,
          task_id=task_id,
          num_shards=num_shards)

    # Check that the regions are the same unsharded vs. sharded.
    unsharded_regions = get_regions(0, 0)
    sharded_regions = []
    for task_id in range(num_shards):
      task_regions = get_regions(task_id, num_shards)
      sharded_regions.extend(task_regions)
    self.assertCountEqual(unsharded_regions, sharded_regions)

  @parameterized.parameters(
      # Providing one of task id and num_shards but not the other is bad.
      (None, 0),
      (None, 2),
      (2, None),
      (0, None),
      # Negative values are illegal.
      (-1, 2),
      (0, -2),
      # task_id >= num_shards is bad.
      (2, 2),
      (3, 2),
  )
  def test_regions_to_process_fails_with_bad_shard_args(self, task, num_shards):
    with self.assertRaises(ValueError):
      make_examples.regions_to_process(
          contigs=_make_contigs([('z', 100), ('a', 100), ('n', 100)]),
          partition_size=10,
          task_id=task,
          num_shards=num_shards)

  def test_catches_bad_argv(self):
    with mock.patch.object(logging, 'error') as mock_logging,\
        mock.patch.object(sys, 'exit') as mock_exit:
      make_examples.main(['make_examples.py', 'extra_arg'])
    mock_logging.assert_called_once_with(
        'Command line parsing failure: make_examples does not accept '
        'positional arguments but some are present on the command line: '
        '"[\'make_examples.py\', \'extra_arg\']".')
    mock_exit.assert_called_once_with(errno.ENOENT)

  @flagsaver.FlagSaver
  def test_catches_bad_flags(self):
    # Set all of the requested flag values.
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.candidates = test_utils.test_tmpfile('vsc.tfrecord')
    FLAGS.examples = test_utils.test_tmpfile('examples.tfrecord')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.partition_size = 1000
    FLAGS.mode = 'training'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    # This is the bad flag.
    FLAGS.confident_regions = ''

    with mock.patch.object(logging, 'error') as mock_logging,\
        mock.patch.object(sys, 'exit') as mock_exit:
      make_examples.main(['make_examples.py'])
    mock_logging.assert_called_once_with(
        'confident_regions is required when in training mode.')
    mock_exit.assert_called_once_with(errno.ENOENT)

  @parameterized.parameters(
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2', '3'],
          vcf_names=None,
          names_to_exclude=[],
          min_coverage_fraction=1.0,
          expected_names=['1', '2', '3']),
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=None,
          names_to_exclude=[],
          min_coverage_fraction=0.66,
          expected_names=['1', '2']),
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=['1', '3'],
          names_to_exclude=[],
          min_coverage_fraction=0.33,
          expected_names=['1']),
      dict(
          ref_names=['1', '2', '3', '4', '5'],
          sam_names=['1', '2', '3'],
          vcf_names=None,
          names_to_exclude=['4', '5'],
          min_coverage_fraction=1.0,
          expected_names=['1', '2', '3']),
  )
  def test_ensure_consistent_contigs(self, ref_names, sam_names, vcf_names,
                                     names_to_exclude, min_coverage_fraction,
                                     expected_names):
    ref_contigs = _make_contigs([(name, 100) for name in ref_names])
    sam_contigs = _make_contigs([(name, 100) for name in sam_names])
    if vcf_names is not None:
      vcf_contigs = _make_contigs([(name, 100) for name in vcf_names])
    else:
      vcf_contigs = None
    actual = make_examples._ensure_consistent_contigs(
        ref_contigs, sam_contigs, vcf_contigs, names_to_exclude,
        min_coverage_fraction)
    self.assertEqual([a.name for a in actual], expected_names)

  @parameterized.parameters(
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=None,
          names_to_exclude=[],
          min_coverage_fraction=0.67),
      dict(
          ref_names=['1', '2', '3'],
          sam_names=['1', '2'],
          vcf_names=['1', '3'],
          names_to_exclude=[],
          min_coverage_fraction=0.34),
  )
  def test_ensure_inconsistent_contigs(self, ref_names, sam_names, vcf_names,
                                       names_to_exclude, min_coverage_fraction):
    ref_contigs = _make_contigs([(name, 100) for name in ref_names])
    sam_contigs = _make_contigs([(name, 100) for name in sam_names])
    if vcf_names is not None:
      vcf_contigs = _make_contigs([(name, 100) for name in vcf_names])
    else:
      vcf_contigs = None
    with self.assertRaisesRegexp(ValueError, 'Reference contigs span'):
      make_examples._ensure_consistent_contigs(ref_contigs, sam_contigs,
                                               vcf_contigs, names_to_exclude,
                                               min_coverage_fraction)

  @flagsaver.FlagSaver
  def test_regions_and_exclude_regions_flags(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.regions = 'chr20:10,000,000-11,000,000'
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.exclude_regions = 'chr20:10,010,000-10,100,000'

    options = make_examples.default_options(add_flags=True)
    self.assertCountEqual(
        list(
            ranges.RangeSet(
                make_examples.processing_regions_from_options(options))),
        _from_literals_list(
            ['chr20:10,000,000-10,009,999', 'chr20:10,100,001-11,000,000']))


class RegionProcessorTest(parameterized.TestCase):

  def setUp(self):
    self.region = ranges.parse_literal('chr20:10,000,000-10,000,100')

    FLAGS.reads = ''
    self.options = make_examples.default_options(add_flags=False)
    self.options.reference_filename = testdata.CHR20_FASTA
    self.options.reads_filename = testdata.CHR20_BAM
    self.options.truth_variants_filename = testdata.TRUTH_VARIANTS_VCF
    self.options.mode = deepvariant_pb2.DeepVariantOptions.TRAINING

    self.processor = make_examples.RegionProcessor(self.options)
    self.mock_init = self.add_mock('_initialize')
    self.default_shape = [5, 5, 7]
    self.default_format = 'raw'

  def add_mock(self, name, retval='dontadd', side_effect='dontadd'):
    patcher = mock.patch.object(self.processor, name, autospec=True)
    self.addCleanup(patcher.stop)
    mocked = patcher.start()
    if retval != 'dontadd':
      mocked.return_value = retval
    if side_effect != 'dontadd':
      mocked.side_effect = side_effect
    return mocked

  def test_on_demand_initialization_called_if_not_initialized(self):
    candidates = ['Candidates']
    self.assertFalse(self.processor.initialized)
    self.processor.in_memory_sam_reader = mock.Mock()
    mock_rr = self.add_mock('region_reads', retval=[])
    mock_cir = self.add_mock('candidates_in_region', retval=(candidates, []))
    mock_lc = self.add_mock('label_candidates', retval=[])
    self.processor.process(self.region)
    test_utils.assert_called_once_workaround(self.mock_init)
    mock_rr.assert_called_once_with(self.region)
    self.processor.in_memory_sam_reader.replace_reads.assert_called_once_with(
        [])
    mock_cir.assert_called_once_with(self.region)
    mock_lc.assert_called_once_with(candidates)

  def test_on_demand_initialization_not_called_if_initialized(self):
    self.processor.initialized = True
    self.assertTrue(self.processor.initialized)
    self.processor.in_memory_sam_reader = mock.Mock()
    mock_rr = self.add_mock('region_reads', retval=[])
    mock_cir = self.add_mock('candidates_in_region', retval=([], []))
    mock_lc = self.add_mock('label_candidates', retval=[])
    self.processor.process(self.region)
    test_utils.assert_not_called_workaround(self.mock_init)
    mock_rr.assert_called_once_with(self.region)
    self.processor.in_memory_sam_reader.replace_reads.assert_called_once_with(
        [])
    mock_cir.assert_called_once_with(self.region)
    test_utils.assert_called_once_workaround(mock_lc)

  def test_process_calls_no_candidates(self):
    self.processor.in_memory_sam_reader = mock.Mock()
    mock_rr = self.add_mock('region_reads', retval=[])
    mock_cir = self.add_mock('candidates_in_region', retval=([], []))
    mock_cpe = self.add_mock('create_pileup_examples', retval=[])
    mock_lc = self.add_mock('label_candidates')
    self.assertEqual(([], [], []), self.processor.process(self.region))
    mock_rr.assert_called_once_with(self.region)
    self.processor.in_memory_sam_reader.replace_reads.assert_called_once_with(
        [])
    mock_cir.assert_called_once_with(self.region)
    test_utils.assert_not_called_workaround(mock_cpe)
    mock_lc.assert_called_once_with([])

  @parameterized.parameters([
      deepvariant_pb2.DeepVariantOptions.TRAINING,
      deepvariant_pb2.DeepVariantOptions.CALLING
  ])
  def test_process_calls_with_candidates(self, mode):
    self.processor.options.mode = mode

    self.processor.in_memory_sam_reader = mock.Mock()
    mock_read = mock.MagicMock()
    mock_candidate = mock.MagicMock()
    mock_example = mock.MagicMock()
    mock_label = mock.MagicMock()
    mock_rr = self.add_mock('region_reads', retval=[mock_read])
    mock_cir = self.add_mock(
        'candidates_in_region', retval=([mock_candidate], []))
    mock_cpe = self.add_mock('create_pileup_examples', retval=[mock_example])
    mock_lc = self.add_mock(
        'label_candidates', retval=[(mock_candidate, mock_label)])
    mock_alte = self.add_mock('add_label_to_example', retval=mock_example)
    self.assertEqual(([mock_candidate], [mock_example], []),
                     self.processor.process(self.region))
    mock_rr.assert_called_once_with(self.region)
    self.processor.in_memory_sam_reader.replace_reads.assert_called_once_with(
        [mock_read])
    mock_cir.assert_called_once_with(self.region)
    mock_cpe.assert_called_once_with(mock_candidate)

    if mode == deepvariant_pb2.DeepVariantOptions.TRAINING:
      mock_lc.assert_called_once_with([mock_candidate])
      mock_alte.assert_called_once_with(mock_example, mock_label)
    else:
      # In training mode we don't label our candidates.
      test_utils.assert_not_called_workaround(mock_lc)
      test_utils.assert_not_called_workaround(mock_alte)

  @parameterized.parameters([
      deepvariant_pb2.DeepVariantOptions.TRAINING,
      deepvariant_pb2.DeepVariantOptions.CALLING
  ])
  def test_process_keeps_ordering_of_candidates_and_examples(self, mode):
    self.processor.options.mode = mode

    r1, r2 = mock.Mock(), mock.Mock()
    c1, c2 = mock.Mock(), mock.Mock()
    l1, l2 = mock.Mock(), mock.Mock()
    e1, e2, e3 = mock.Mock(), mock.Mock(), mock.Mock()
    self.processor.in_memory_sam_reader = mock.Mock()
    self.add_mock('region_reads', retval=[r1, r2])
    self.add_mock('candidates_in_region', retval=([c1, c2], []))
    mock_cpe = self.add_mock(
        'create_pileup_examples', side_effect=[[e1], [e2, e3]])
    mock_lc = self.add_mock('label_candidates', retval=[(c1, l1), (c2, l2)])
    mock_alte = self.add_mock('add_label_to_example', side_effect=[e1, e2, e3])
    self.assertEqual(([c1, c2], [e1, e2, e3], []),
                     self.processor.process(self.region))
    self.processor.in_memory_sam_reader.replace_reads.assert_called_once_with(
        [r1, r2])
    # We don't try to label variants when in calling mode.
    self.assertEqual([mock.call(c1), mock.call(c2)], mock_cpe.call_args_list)

    if mode == deepvariant_pb2.DeepVariantOptions.CALLING:
      # In calling mode, we never try to label.
      test_utils.assert_not_called_workaround(mock_lc)
      test_utils.assert_not_called_workaround(mock_alte)
    else:
      mock_lc.assert_called_once_with([c1, c2])
      self.assertEqual([
          mock.call(e1, l1),
          mock.call(e2, l2),
          mock.call(e3, l2),
      ], mock_alte.call_args_list)

  def test_process_with_realigner(self):
    self.processor.options.mode = deepvariant_pb2.DeepVariantOptions.CALLING
    self.processor.options.realigner_options.CopyFrom(
        realigner_pb2.RealignerOptions())
    self.processor.realigner = mock.Mock()
    self.processor.realigner.realign_reads.return_value = [], []

    self.processor.sam_reader = mock.Mock()
    self.processor.sam_reader.query.return_value = []
    self.processor.in_memory_sam_reader = mock.Mock()

    c1, c2 = mock.Mock(), mock.Mock()
    e1, e2, e3 = mock.Mock(), mock.Mock(), mock.Mock()
    self.add_mock('candidates_in_region', retval=([c1, c2], []))
    mock_cpe = self.add_mock(
        'create_pileup_examples', side_effect=[[e1], [e2, e3]])
    mock_lc = self.add_mock('label_candidates')

    self.assertEqual(([c1, c2], [e1, e2, e3], []),
                     self.processor.process(self.region))
    self.processor.sam_reader.query.assert_called_once_with(self.region)
    self.processor.realigner.realign_reads.assert_called_once_with([],
                                                                   self.region)
    self.processor.in_memory_sam_reader.replace_reads.assert_called_once_with(
        [])
    self.assertEqual([mock.call(c1), mock.call(c2)], mock_cpe.call_args_list)
    test_utils.assert_not_called_workaround(mock_lc)

  def test_candidates_in_region_no_reads(self):
    self.processor.in_memory_sam_reader = mock.Mock()
    self.processor.in_memory_sam_reader.query.return_value = []
    mock_ac = self.add_mock('_make_allele_counter_for_region')

    self.assertEqual(([], []), self.processor.candidates_in_region(self.region))

    self.processor.in_memory_sam_reader.query.assert_called_once_with(
        self.region)
    # A region with no reads should return out without making an AlleleCounter.
    test_utils.assert_not_called_workaround(mock_ac)

  @parameterized.parameters(True, False)
  def test_candidates_in_region(self, include_gvcfs):
    self.options.gvcf_filename = 'foo.vcf' if include_gvcfs else ''
    self.processor.in_memory_sam_reader = mock.Mock()
    reads = ['read1', 'read2']
    self.processor.in_memory_sam_reader.query.return_value = reads
    # Setup our make_allele_counter and other mocks.
    mock_ac = mock.Mock()
    mock_make_ac = self.add_mock(
        '_make_allele_counter_for_region', retval=mock_ac)
    # Setup our make_variant_caller and downstream mocks.
    mock_vc = mock.Mock()
    expected_calls = (['variant'], ['gvcf'] if include_gvcfs else [])
    mock_vc.calls_from_allele_counter.return_value = expected_calls
    self.processor.variant_caller = mock_vc

    actual = self.processor.candidates_in_region(self.region)

    # Make sure we're getting our reads for the region.
    self.processor.in_memory_sam_reader.query.assert_called_once_with(
        self.region)

    # Make sure we're creating an AlleleCounter once and adding each of our
    # reads to it.
    mock_make_ac.assert_called_once_with(self.region)
    self.assertEqual([mock.call(r) for r in reads], mock_ac.add.call_args_list)

    # Make sure we call CallVariant for each of the counts returned by the
    # allele counter.
    mock_vc.calls_from_allele_counter.assert_called_once_with(
        mock_ac, include_gvcfs)

    # Finally, our actual result should be the single 'variant' and potentially
    # the gvcf records.
    self.assertEqual(expected_calls, actual)

  def test_create_pileup_examples_handles_none(self):
    self.processor.pic = mock.Mock()
    dv_call = mock.Mock()
    self.processor.pic.create_pileup_images.return_value = None
    self.assertEqual([], self.processor.create_pileup_examples(dv_call))
    self.processor.pic.create_pileup_images.assert_called_once_with(dv_call)

  def test_create_pileup_examples(self):
    self.processor.pic = mock.Mock()
    self.add_mock(
        '_encode_tensor',
        side_effect=[('tensor1', self.default_shape, self.default_format),
                     ('tensor2', self.default_shape, self.default_format)])
    dv_call = mock.Mock()
    dv_call.variant = test_utils.make_variant(start=10, alleles=['A', 'C', 'G'])
    ex = mock.Mock()
    alt1, alt2 = ['C'], ['G']
    self.processor.pic.create_pileup_images.return_value = [(alt1, 'tensor1'),
                                                            (alt2, 'tensor2')]

    actual = self.processor.create_pileup_examples(dv_call)

    self.processor.pic.create_pileup_images.assert_called_once_with(dv_call)

    self.assertEquals(len(actual), 2)
    for ex, (alt, img) in zip(actual, [(alt1, 'tensor1'), (alt2, 'tensor2')]):
      self.assertEqual(tf_utils.example_alt_alleles(ex), alt)
      self.assertEqual(tf_utils.example_variant(ex), dv_call.variant)
      self.assertEqual(tf_utils.example_encoded_image(ex), img)
      self.assertEqual(tf_utils.example_image_shape(ex), self.default_shape)
      self.assertEqual(tf_utils.example_image_format(ex), self.default_format)

  @parameterized.parameters(
      # Test that a het variant gets a label value of 1 assigned to the example.
      dict(
          label=variant_labeler.VariantLabel(
              is_confident=True,
              variant=test_utils.make_variant(start=10, alleles=['A', 'C']),
              genotype=(0, 1)),
          expected_label_value=1,
      ),
      # Test that a reference variant gets a label value of 0 in the example.
      dict(
          label=variant_labeler.VariantLabel(
              is_confident=True,
              variant=test_utils.make_variant(start=10, alleles=['A', '.']),
              genotype=(0, 0)),
          expected_label_value=0,
      ),
  )
  def test_add_label_to_example(self, label, expected_label_value):
    example = self._example_for_variant(label.variant)
    labeled = copy.deepcopy(example)
    actual = self.processor.add_label_to_example(labeled, label)

    # The add_label_to_example command modifies labeled and returns it.
    self.assertIs(actual, labeled)

    # Check that all keys from example are present in labeled.
    for key, value in example.features.feature.iteritems():
      if key != 'variant/encoded':  # Special case tested below.
        self.assertEqual(value, labeled.features.feature[key])

    # The genotype of our example_variant should be set to the true genotype
    # according to our label.
    self.assertEqual(expected_label_value, tf_utils.example_label(labeled))
    labeled_variant = tf_utils.example_variant(labeled)
    call = variant_utils.only_call(labeled_variant)
    self.assertEqual(tuple(call.genotype), label.genotype)

    # The original variant and labeled_variant from out tf.Example should be
    # equal except for the genotype field, since this is set by
    # add_label_to_example.
    label.variant.calls[0].genotype[:] = []
    call.genotype[:] = []
    self.assertEqual(label.variant, labeled_variant)

  def test_label_variant_raises_for_non_confident_variant(self):
    label = variant_labeler.VariantLabel(
        is_confident=False,
        variant=test_utils.make_variant(start=10, alleles=['A', 'C']),
        genotype=(0, 1))
    example = self._example_for_variant(label.variant)
    with self.assertRaisesRegexp(
        ValueError, 'Cannot add a non-confident label to an example'):
      self.processor.add_label_to_example(example, label)

  def _example_for_variant(self, variant):
    return tf_utils.make_example(variant, list(variant.alternate_bases), 'foo',
                                 self.default_shape, self.default_format)


if __name__ == '__main__':
  absltest.main()
