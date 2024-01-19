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
"""Tests for deepvariant.make_examples."""

import enum
import errno
import json
import platform
import sys
from unittest import mock



from absl import flags
from absl import logging
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
from etils import epath
import numpy as np

from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import make_examples
from deepvariant import make_examples_core
from deepvariant import testdata
from deepvariant.protos import deepvariant_pb2
from tensorflow.python.platform import gfile
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils
from third_party.nucleus.util import vcf_constants
from third_party.nucleus.util import vis

FLAGS = flags.FLAGS

# Dictionary mapping keys to decoders for decode_example function.
_EXAMPLE_DECODERS = {
    'locus': dv_utils.example_locus,
    'alt_allele_indices/encoded': dv_utils.example_alt_alleles_indices,
    'image/encoded': dv_utils.example_encoded_image,
    'variant/encoded': dv_utils.example_variant,
    'variant_type': dv_utils.example_variant_type,
    'label': dv_utils.example_label,
    'image/shape': dv_utils.example_image_shape,
    'sequencing_type': dv_utils.example_sequencing_type,
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
  logging.set_verbosity(logging.FATAL)
  testdata.init()


def _sharded(basename, num_shards=None):
  if num_shards:
    return basename + '@' + str(num_shards)
  else:
    return basename


class TestConditions(enum.Enum):
  """Enum capturing what the test condition we're using."""

  USE_BAM = 1
  USE_CRAM = 2
  USE_MULTI_BAMS = 3


class MakeExamplesEnd2EndTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self._saved_flags = flagsaver.save_flag_values()

  def tearDown(self):
    super().tearDown()
    flagsaver.restore_flag_values(self._saved_flags)

  @flagsaver.flagsaver
  def test_make_examples_check_diff_channels_ordering(self):
    """Confirms the channels with diff_channels are ordered as expected."""
    FLAGS.reads = testdata.CHR20_BAM
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.candidates = test_utils.test_tmpfile(
        'check_pb_channels.vsc.tfrecord.gz'
    )
    FLAGS.examples = test_utils.test_tmpfile('check_pb_channels.ex.tfrecord.gz')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.add_hp_channel = True
    FLAGS.alt_aligned_pileup = 'diff_channels'
    FLAGS.max_reads_per_partition = 600
    FLAGS.min_mapping_quality = 1
    FLAGS.parse_sam_aux_fields = True
    FLAGS.partition_size = 25000
    FLAGS.phase_reads = True
    FLAGS.pileup_image_width = 199
    FLAGS.realign_reads = False
    FLAGS.sort_by_haplotypes = True
    FLAGS.track_ref_reads = True
    FLAGS.vsc_min_fraction_indels = 0.12
    FLAGS.mode = 'calling'
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    example_info_json = dv_utils.get_example_info_json_filename(
        FLAGS.examples, None
    )
    example_info = json.load(gfile.GFile(example_info_json, 'r'))
    self.assertEqual(example_info['channels'], [1, 2, 3, 4, 5, 6, 7, 9, 10])

  @flagsaver.flagsaver
  def test_make_examples_compare_realignment_modes(self):
    def _run_with_realignment_mode(enable_joint_realignment):
      FLAGS.enable_joint_realignment = enable_joint_realignment
      num_shards = 1
      FLAGS.reads = testdata.CHR20_BAM
      region = ranges.parse_literal('chr20:10,000,000-10,010,000')
      FLAGS.ref = testdata.CHR20_FASTA
      FLAGS.candidates = test_utils.test_tmpfile(
          _sharded(f'jr-{enable_joint_realignment}.vsc.tfrecord', num_shards)
      )
      FLAGS.examples = test_utils.test_tmpfile(
          _sharded(f'jr-{enable_joint_realignment}.ex.tfrecord', num_shards)
      )
      FLAGS.regions = [ranges.to_literal(region)]
      FLAGS.partition_size = 1000
      FLAGS.mode = 'calling'
      FLAGS.gvcf_gq_binsize = 5
      FLAGS.gvcf = test_utils.test_tmpfile(
          _sharded('compare_realignment_modes.gvcf.tfrecord', num_shards)
      )
      FLAGS.task = 0
      options = make_examples.default_options(add_flags=True)
      # We need to overwrite bam_fname for USE_CRAM test since Golden Set
      # generated from BAM file. BAM filename is stored in candidates. If we
      # don't overwrite default_options variants won't match and test fail.
      options.bam_fname = 'NA12878_S1.chr20.10_10p1mb.bam'
      make_examples_core.make_examples_runner(options)

      # Test that our candidates are reasonable, calling specific helper
      # functions to check lots of properties of the output.
      candidates = sorted(
          tfrecord.read_tfrecords(
              FLAGS.candidates,
              proto=deepvariant_pb2.DeepVariantCall,
          ),
          key=lambda c: variant_utils.variant_range_tuple(c.variant),
      )
      self.verify_deepvariant_calls(candidates, options)
      self.verify_variants(
          [call.variant for call in candidates], region, options, is_gvcf=False
      )

      # Verify that the variants in the examples are all good.
      examples = self.verify_examples(
          FLAGS.examples, region, options, verify_labels=False
      )
      example_variants = [dv_utils.example_variant(ex) for ex in examples]
      self.verify_variants(example_variants, region, options, is_gvcf=False)
      return examples

    examples1 = _run_with_realignment_mode(False)
    examples2 = _run_with_realignment_mode(True)
    # Because this test is with just one sample, whether
    # enable_joint_realignment is True or False doesn't make a difference.
    # NOTE: When creating this test, I deliberately change the behavior of
    # enable_joint_realignment==False and confirm that this test can fail,
    # if the outputs are different when we alter enable_joint_realignment.
    self.assertDeepVariantExamplesEqual(examples1, examples2)

  @parameterized.parameters(
      dict(
          mode='calling',
          max_reads_per_partition=1500,
          expected_len_examples1=84,
          expected_len_examples2=17,
      ),
      dict(
          mode='calling',
          max_reads_per_partition=0,
          expected_len_examples1=84,
          expected_len_examples2=17,
      ),
      # For 'candidate_sweep' mode, it won't create examples, so we just aim
      # to run it through without errors.
      dict(
          mode='candidate_sweep',
          max_reads_per_partition=1500,
          expected_len_examples1=None,
          expected_len_examples2=None,
      ),
      dict(
          mode='candidate_sweep',
          max_reads_per_partition=0,
          expected_len_examples1=None,
          expected_len_examples2=None,
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_with_max_reads_for_dynamic_bases_per_region(
      self,
      mode,
      max_reads_per_partition,
      expected_len_examples1,
      expected_len_examples2,
  ):
    num_shards = 1
    FLAGS.reads = testdata.CHR20_BAM
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded(
            'test_max_reads_per_partition_and_bases.ex.tfrecord.gz', num_shards
        )
    )
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.partition_size = 1000
    FLAGS.mode = mode
    FLAGS.gvcf_gq_binsize = 5
    FLAGS.task = 0
    FLAGS.max_reads_per_partition = max_reads_per_partition

    if mode == 'candidate_sweep':
      FLAGS.candidate_positions = test_utils.test_tmpfile(
          _sharded(
              'test_max_reads_per_partition_and_bases.candidate_positions',
              num_shards,
          )
      )

    options = make_examples.default_options(add_flags=True)
    # We need to overwrite bam_fname for USE_CRAM test since Golden Set
    # generated from BAM file. BAM filename is stored in candidates. If we
    # don't overwrite default_options variants won't match and test fail.
    options.bam_fname = 'NA12878_S1.chr20.10_10p1mb.bam'
    make_examples_core.make_examples_runner(options)
    if expected_len_examples1 is not None:
      examples1 = self.verify_examples(
          FLAGS.examples, region, options, verify_labels=False
      )
      self.assertLen(examples1, expected_len_examples1)

    # Now, this is the main part of the test. I want to test the behavior after
    # I set max_reads_for_dynamic_bases_per_region.
    FLAGS.max_reads_for_dynamic_bases_per_region = 1
    options = make_examples.default_options(add_flags=True)
    options.bam_fname = 'NA12878_S1.chr20.10_10p1mb.bam'
    make_examples_core.make_examples_runner(options)
    if expected_len_examples2 is not None:
      examples2 = self.verify_examples(
          FLAGS.examples, region, options, verify_labels=False
      )
      self.assertLen(examples2, expected_len_examples2)

  # Golden sets are created with learning/genomics/internal/create_golden.sh
  @parameterized.parameters(
      # All tests are run with fast_pass_aligner enabled. There are no
      # golden sets version for ssw realigner.
      dict(mode='calling', num_shards=0),
      dict(mode='calling', num_shards=3),
      dict(mode='candidate_sweep', num_shards=0),
      dict(mode='candidate_sweep', num_shards=3),
      dict(
          mode='training', num_shards=0, labeler_algorithm='haplotype_labeler'
      ),
      dict(
          mode='training', num_shards=3, labeler_algorithm='haplotype_labeler'
      ),
      dict(
          mode='training', num_shards=0, labeler_algorithm='positional_labeler'
      ),
      dict(
          mode='training', num_shards=3, labeler_algorithm='positional_labeler'
      ),
      # The following tests are for CRAM input:
      dict(
          mode='calling', num_shards=0, test_condition=TestConditions.USE_CRAM
      ),
      dict(
          mode='training',
          num_shards=0,
          test_condition=TestConditions.USE_CRAM,
          labeler_algorithm='haplotype_labeler',
      ),
      # The following tests are for multiple BAM inputs:
      dict(
          mode='calling',
          num_shards=0,
          test_condition=TestConditions.USE_MULTI_BAMS,
      ),
      dict(
          mode='training',
          num_shards=0,
          test_condition=TestConditions.USE_MULTI_BAMS,
          labeler_algorithm='haplotype_labeler',
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_end2end(
      self,
      mode,
      num_shards,
      test_condition=TestConditions.USE_BAM,
      labeler_algorithm=None,
      use_fast_pass_aligner=True,
  ):
    self.assertIn(mode, {'calling', 'training', 'candidate_sweep'})
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.write_run_info = True
    FLAGS.ref = testdata.CHR20_FASTA
    if test_condition == TestConditions.USE_BAM:
      FLAGS.reads = testdata.CHR20_BAM
    elif test_condition == TestConditions.USE_CRAM:
      FLAGS.reads = testdata.CHR20_CRAM
    elif test_condition == TestConditions.USE_MULTI_BAMS:
      FLAGS.reads = ','.join(
          [testdata.CHR20_BAM_FIRST_HALF, testdata.CHR20_BAM_SECOND_HALF]
      )

    FLAGS.candidates = test_utils.test_tmpfile(
        _sharded('vsc.tfrecord', num_shards)
    )
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('examples.tfrecord', num_shards)
    )
    if mode == 'candidate_sweep':
      FLAGS.candidate_positions = test_utils.test_tmpfile(
          _sharded('candidate_positions', num_shards)
      )
      candidate_positions = test_utils.test_tmpfile(
          _sharded('candidate_positions', num_shards)
      )
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.partition_size = 1000
    FLAGS.mode = mode
    FLAGS.gvcf_gq_binsize = 5
    FLAGS.use_fast_pass_aligner = use_fast_pass_aligner
    if labeler_algorithm is not None:
      FLAGS.labeler_algorithm = labeler_algorithm

    if mode == 'calling':
      FLAGS.gvcf = test_utils.test_tmpfile(
          _sharded('gvcf.tfrecord', num_shards)
      )
    else:
      FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
      FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED

    if mode == 'candidate_sweep':
      golden_candidate_positions = _sharded(
          testdata.GOLDEN_CANDIDATE_POSITIONS, num_shards
      )

    for task_id in range(max(num_shards, 1)):
      FLAGS.task = task_id
      options = make_examples.default_options(add_flags=True)
      # We need to overwrite bam_fname for USE_CRAM test since Golden Set
      # generated from BAM file. BAM filename is stored in candidates. If we
      # don't overwrite default_options variants won't match and test fail.
      options.bam_fname = 'NA12878_S1.chr20.10_10p1mb.bam'
      make_examples_core.make_examples_runner(options)

      # Check that our run_info proto contains the basic fields we'd expect:
      # (a) our options are written to the run_info.options field.
      run_info = make_examples_core.read_make_examples_run_info(
          options.run_info_filename
      )
      self.assertEqual(run_info.options, options)
      # (b) run_info.resource_metrics is present and contains our hostname.
      self.assertTrue(run_info.HasField('resource_metrics'))
      self.assertEqual(run_info.resource_metrics.host_name, platform.node())

      # For candidate_sweep mode we verify that candidate positions match
      # golden set exactly.
      if mode == 'candidate_sweep':
        _, candidates_path = sharded_file_utils.resolve_filespecs(
            task_id, candidate_positions
        )
        _, gold_candidates_path = sharded_file_utils.resolve_filespecs(
            task_id, golden_candidate_positions
        )
        self.verify_candidate_positions(candidates_path, gold_candidates_path)

    # In candidate_sweep mode the test stops here.
    if mode == 'candidate_sweep':
      return

    # Test that our candidates are reasonable, calling specific helper functions
    # to check lots of properties of the output.
    candidates = sorted(
        tfrecord.read_tfrecords(
            FLAGS.candidates, proto=deepvariant_pb2.DeepVariantCall
        ),
        key=lambda c: variant_utils.variant_range_tuple(c.variant),
    )
    self.verify_deepvariant_calls(candidates, options)
    self.verify_variants(
        [call.variant for call in candidates], region, options, is_gvcf=False
    )

    # Verify that the variants in the examples are all good.
    examples = self.verify_examples(
        FLAGS.examples, region, options, verify_labels=mode == 'training'
    )
    example_variants = [dv_utils.example_variant(ex) for ex in examples]
    self.verify_variants(example_variants, region, options, is_gvcf=False)

    # Verify the integrity of the examples and then check that they match our
    # golden labeled examples. Note we expect the order for both training and
    # calling modes to produce deterministic order because we fix the random
    # seed.
    if mode in ['calling', 'candidate_sweep']:
      golden_file = _sharded(testdata.GOLDEN_CALLING_EXAMPLES, num_shards)
    else:
      golden_file = _sharded(testdata.GOLDEN_TRAINING_EXAMPLES, num_shards)

    compression_type = 'GZIP' if 'tfrecord.gz' in golden_file else None
    self.assertDeepVariantExamplesEqual(
        examples,
        list(
            tfrecord.read_tfrecords(
                golden_file, compression_type=compression_type
            )
        ),
    )

    if mode == 'calling':
      nist_reader = vcf.VcfReader(testdata.TRUTH_VARIANTS_VCF)
      nist_variants = list(nist_reader.query(region))
      self.verify_nist_concordance(example_variants, nist_variants)

      # Check the quality of our generated gvcf file.
      gvcfs = variant_utils.sorted_variants(
          tfrecord.read_tfrecords(FLAGS.gvcf, proto=variants_pb2.Variant)
      )
      self.verify_variants(gvcfs, region, options, is_gvcf=True)
      self.verify_contiguity(gvcfs, region)
      gvcf_golden_file = _sharded(
          testdata.GOLDEN_POSTPROCESS_GVCF_INPUT, num_shards
      )
      compression_type = 'GZIP' if 'tfrecord.gz' in gvcf_golden_file else None
      expected_gvcfs = list(
          tfrecord.read_tfrecords(
              gvcf_golden_file,
              proto=variants_pb2.Variant,
              compression_type=compression_type,
          )
      )
      # Despite the name, assertCountEqual checks that all elements match.
      self.assertCountEqual(gvcfs, expected_gvcfs)

    if (
        mode == 'training'
        and num_shards == 0
        and labeler_algorithm != 'positional_labeler'
    ):
      # The positional labeler doesn't track metrics, so don't try to read them
      # in when that's the mode.
      self.assertEqual(
          make_examples_core.read_make_examples_run_info(
              testdata.GOLDEN_MAKE_EXAMPLES_RUN_INFO
          ).labeling_metrics,
          run_info.labeling_metrics,
      )

  @flagsaver.flagsaver
  def test_make_examples_end2end_failed_on_mismatched_multi_bam(self):
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')

    FLAGS.write_run_info = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = ','.join([testdata.CHR20_BAM, testdata.NOCHR20_BAM])
    FLAGS.candidates = test_utils.test_tmpfile(
        _sharded('mismatched_multi_bam.vsc.tfrecord')
    )
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('mismatched_multi_bam.examples.tfrecord')
    )
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.partition_size = 1000
    FLAGS.mode = 'calling'
    FLAGS.gvcf_gq_binsize = 5
    options = make_examples.default_options(add_flags=True)
    # This shows an example of what the error message looks like:
    # TODO: OpError exception not propagated.
    with self.assertRaisesRegex(
        ValueError,
        (
            'NOT_FOUND: Unknown reference_name '
            'reference_name:[ \t]*"chr20" start: 9999999 end: 10000999'
        ),
    ):
      make_examples_core.make_examples_runner(options)

  @flagsaver.flagsaver
  def test_make_examples_end2end_failed_on_cram(self):
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')

    FLAGS.use_ref_for_cram = False
    FLAGS.write_run_info = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_CRAM
    FLAGS.candidates = test_utils.test_tmpfile(_sharded('failed.vsc.tfrecord'))
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('failed.examples.tfrecord')
    )
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.partition_size = 1000
    FLAGS.mode = 'calling'
    FLAGS.gvcf_gq_binsize = 5
    options = make_examples.default_options(add_flags=True)
    with self.assertRaisesRegex(ValueError, 'Failed to parse BAM/CRAM file.'):
      make_examples_core.make_examples_runner(options)

  # Golden sets are created with learning/genomics/internal/create_golden.sh
  @flagsaver.flagsaver
  def test_make_examples_training_end2end_with_customized_classes_labeler(self):
    FLAGS.labeler_algorithm = 'customized_classes_labeler'
    FLAGS.customized_classes_labeler_classes_list = 'ref,class1,class2'
    FLAGS.customized_classes_labeler_info_field_name = 'type'
    region = ranges.parse_literal('chr20:10,000,000-10,004,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.candidates = test_utils.test_tmpfile(_sharded('vsc.tfrecord'))
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.partition_size = 1000
    FLAGS.mode = 'training'
    FLAGS.gvcf_gq_binsize = 5
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF_WITH_TYPES
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    golden_file = _sharded(testdata.CUSTOMIZED_CLASSES_GOLDEN_TRAINING_EXAMPLES)
    # Verify that the variants in the examples are all good.
    examples = self.verify_examples(
        FLAGS.examples, region, options, verify_labels=True
    )
    self.assertDeepVariantExamplesEqual(
        examples, list(tfrecord.read_tfrecords(golden_file))
    )

  @flagsaver.flagsaver
  def test_make_examples_end2end_confirm_downsample_fraction_used(self):
    def _get_examples(downsample_fraction=None):
      if downsample_fraction is not None:
        FLAGS.downsample_fraction = downsample_fraction
      options = make_examples.default_options(add_flags=True)
      make_examples_core.make_examples_runner(options)
      examples = self.verify_examples(
          FLAGS.examples, region, options, verify_labels=False
      )
      return examples

    region = ranges.parse_literal('chr20:10,000,000-10,004,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.mode = 'calling'
    examples1 = _get_examples()
    examples2 = _get_examples(0.01)
    self.assertLess(len(examples2), len(examples1))

  @flagsaver.flagsaver
  def test_make_examples_end2end_confirm_vsc_min_fraction_used(self):
    """Set very low vsc_max_fraction_{snps,indels} and confirm they're used."""
    region = ranges.parse_literal('chr20:10,000,000-10,004,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('confirm_vsc_min.examples.tfrecord')
    )
    FLAGS.mode = 'calling'
    # Setting min_fractions to larger than 100% to confirm that this will end
    # up removing all examples.
    FLAGS.vsc_min_fraction_snps = 1.1
    FLAGS.vsc_min_fraction_indels = 1.1
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    examples = self.verify_examples(
        FLAGS.examples, region, options, verify_labels=False
    )
    self.assertEmpty(examples)

  # Golden sets are created with learning/genomics/internal/create_golden.sh
  @parameterized.parameters(
      dict(mode='calling'),
      dict(mode='training'),
  )
  @flagsaver.flagsaver
  def test_make_examples_end2end_vcf_candidate_importer(self, mode):
    FLAGS.variant_caller = 'vcf_candidate_importer'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.candidates = test_utils.test_tmpfile(
        _sharded('vcf_candidate_importer.{}.tfrecord'.format(mode))
    )
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('vcf_candidate_importer.examples.{}.tfrecord'.format(mode))
    )
    FLAGS.mode = mode

    if mode == 'calling':
      golden_file = _sharded(
          testdata.GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES
      )
      FLAGS.proposed_variants = testdata.VCF_CANDIDATE_IMPORTER_VARIANTS
      # Adding the following flags to match how the testdata was created.
      FLAGS.regions = 'chr20:59,777,000-60,000,000'
      FLAGS.realign_reads = False
    else:
      golden_file = _sharded(
          testdata.GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES
      )
      FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    # Verify that the variants in the examples are all good.
    examples = self.verify_examples(
        FLAGS.examples, None, options, verify_labels=mode == 'training'
    )
    self.assertDeepVariantExamplesEqual(
        examples, list(tfrecord.read_tfrecords(golden_file))
    )
    self.assertEqual(
        decode_example(examples[0])['image/shape'],
        [100, 221, dv_constants.PILEUP_NUM_CHANNELS],
    )

  @flagsaver.flagsaver
  def test_make_examples_training_vcf_candidate_importer_regions(self):
    """Confirms confident_regions is used in vcf_candidate_importer training."""

    def _get_examples(use_confident_regions=False):
      # `flag_name` can be either 'confident_regions' or 'regions'. Both should
      # be used to constrain the set of candidates generated, and as a result
      # generating the same examples.
      bed_path = test_utils.test_tmpfile('vcf_candidate_importer.bed')
      with gfile.Open(bed_path, 'w') as fout:
        fout.write('\t'.join(['chr20', '10000000', '10001000']) + '\n')
      if use_confident_regions:
        FLAGS.confident_regions = bed_path
        FLAGS.regions = ''
      else:
        FLAGS.confident_regions = ''
        FLAGS.regions = bed_path

      FLAGS.examples = test_utils.test_tmpfile(
          _sharded('vcf_candidate_importer.tfrecord')
      )
      FLAGS.mode = 'training'
      FLAGS.reads = testdata.CHR20_BAM
      FLAGS.ref = testdata.CHR20_FASTA
      FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
      FLAGS.variant_caller = 'vcf_candidate_importer'

      options = make_examples.default_options(add_flags=True)
      make_examples_core.make_examples_runner(options)
      # Verify that the variants in the examples are all good.
      examples = self.verify_examples(
          FLAGS.examples, None, options, verify_labels=False
      )
      return examples

    examples_with_regions = _get_examples(use_confident_regions=False)
    examples_with_confident_regions = _get_examples(use_confident_regions=True)
    self.assertNotEmpty(examples_with_regions)
    self.assertDeepVariantExamplesEqual(
        examples_with_regions, examples_with_confident_regions
    )

  # Golden sets are created with learning/genomics/internal/create_golden.sh
  @parameterized.parameters(
      dict(
          alt_align='rows',
          expected_shape=[300, 221, dv_constants.PILEUP_NUM_CHANNELS],
      ),
      dict(
          alt_align='diff_channels',
          expected_shape=[100, 221, dv_constants.PILEUP_NUM_CHANNELS + 2],
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_training_end2end_with_alt_aligned_pileup(
      self, alt_align, expected_shape
  ):
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.candidates = test_utils.test_tmpfile(_sharded('vsc.tfrecord'))
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.partition_size = 1000
    FLAGS.mode = 'training'
    FLAGS.gvcf_gq_binsize = 5
    FLAGS.alt_aligned_pileup = alt_align  # This is the only input change.
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    options = make_examples.default_options(add_flags=True)
    # Run make_examples with the flags above.
    make_examples_core.make_examples_runner(options)

    # Check the output for shape and against the golden file.
    if alt_align == 'rows':
      golden_file = _sharded(testdata.ALT_ALIGNED_ROWS_EXAMPLES)
    elif alt_align == 'diff_channels':
      golden_file = _sharded(testdata.ALT_ALIGNED_DIFF_CHANNELS_EXAMPLES)
    else:
      raise ValueError(
          "Golden data doesn't exist for this alt_align option: {}".format(
              alt_align
          )
      )
    # Verify that the variants in the examples are all good.
    examples = self.verify_examples(
        FLAGS.examples, region, options, verify_labels=True
    )
    self.assertDeepVariantExamplesEqual(
        examples, list(tfrecord.read_tfrecords(golden_file))
    )
    # Pileup image should have 3 rows of height 100, so resulting height is 300.
    self.assertEqual(decode_example(examples[0])['image/shape'], expected_shape)

  @flagsaver.flagsaver
  def test_make_examples_runtime_by_region(self):
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.mode = 'calling'
    num_shards = 4
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('examples.tfrecord', num_shards)
    )
    # Use same number of shards for profiling files as examples.
    output_prefix = test_utils.test_tmpfile('runtime_profile')
    FLAGS.runtime_by_region = output_prefix + '@{}'.format(num_shards)
    FLAGS.task = 2
    # Run make_examples with those FLAGS.
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    # Sharded output ending in @4 becomes -00002-of-00004 for task 2.
    expected_output_path = output_prefix + '-0000{}-of-00004'.format(FLAGS.task)
    expected_columns = [
        'region',
        'get reads',
        'find candidates',
        'make pileup images',
        'write outputs',
        'num reads',
        'num candidates',
        'num examples',
    ]

    with gfile.Open(expected_output_path, 'r') as fin:
      header = fin.readline()
      column_names = header.strip().split('\t')
      self.assertEqual(expected_columns, column_names)
      non_header_lines = fin.readlines()
      self.assertLen(non_header_lines, 3)
      one_row = non_header_lines[0].strip().split('\t')
      self.assertEqual(len(one_row), len(column_names))
      self.assertGreater(int(one_row[5]), 0, msg='num reads > 0')
      self.assertGreater(int(one_row[6]), 0, msg='num candidates > 0')
      self.assertGreater(int(one_row[7]), 0, msg='num examples > 0')

  @parameterized.parameters(
      dict(select_types=None, expected_count=78),
      dict(select_types='all', expected_count=78),
      dict(select_types='snps', expected_count=64),
      dict(select_types='indels', expected_count=11),
      dict(select_types='snps indels', expected_count=75),
      dict(select_types='multi-allelics', expected_count=3),
      dict(select_types=None, keep_legacy_behavior=True, expected_count=78),
      dict(select_types='all', keep_legacy_behavior=True, expected_count=78),
      dict(select_types='snps', keep_legacy_behavior=True, expected_count=64),
      dict(select_types='indels', keep_legacy_behavior=True, expected_count=11),
      dict(
          select_types='snps indels',
          keep_legacy_behavior=True,
          expected_count=75,
      ),
      dict(
          select_types='multi-allelics',
          keep_legacy_behavior=True,
          expected_count=3,
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_with_variant_selection(
      self, select_types, expected_count, keep_legacy_behavior=False
  ):
    if select_types is not None:
      FLAGS.select_variant_types = select_types
    region = ranges.parse_literal('chr20:10,000,000-10,010,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.candidates = test_utils.test_tmpfile(_sharded('vsc.tfrecord'))
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.partition_size = 1000
    FLAGS.mode = 'calling'
    FLAGS.keep_legacy_allele_counter_behavior = keep_legacy_behavior
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)

    candidates = list(tfrecord.read_tfrecords(FLAGS.candidates))
    self.assertLen(candidates, expected_count)

  @flagsaver.flagsaver
  def test_make_examples_with_allele_frequency_error_dup_chr(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.GRCH38_FASTA
    FLAGS.reads = testdata.GRCH38_CHR20_AND_21_BAM
    num_shards = 1
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('examples.tfrecord', num_shards)
    )
    region = ranges.parse_literal('chr20:61001-62000')
    FLAGS.use_allele_frequency = True
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.population_vcfs = ' '.join(
        [testdata.AF_VCF_CHR20_21_WILDCARD, testdata.AF_VCF_CHR20]
    )
    options = make_examples.default_options(add_flags=True)
    with self.assertRaisesRegex(
        ValueError, 'Variants on chr20 are included in multiple VCFs'
    ):
      # Run make_examples with the flags above.
      make_examples_core.make_examples_runner(options)

  @parameterized.parameters(
      dict(mode='one vcf'),
      dict(mode='two vcfs'),
      dict(mode='two vcfs with wildcard'),
  )
  @flagsaver.flagsaver
  def test_make_examples_with_allele_frequency(self, mode):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.GRCH38_FASTA
    FLAGS.reads = testdata.GRCH38_CHR20_AND_21_BAM
    num_shards = 1
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('examples.tfrecord', num_shards)
    )
    region = ranges.parse_literal('chr20:61001-62000')
    FLAGS.use_allele_frequency = True
    FLAGS.regions = [ranges.to_literal(region)]
    if mode == 'one vcf':
      FLAGS.population_vcfs = testdata.AF_VCF_CHR20_AND_21
    elif mode == 'two vcfs':
      FLAGS.population_vcfs = ' '.join(
          [testdata.AF_VCF_CHR20, testdata.AF_VCF_CHR21]
      )
    elif mode == 'two vcfs with wildcard':
      FLAGS.population_vcfs = testdata.AF_VCF_CHR20_21_WILDCARD
    else:
      raise ValueError('Invalid mode for parameterized test.')
    options = make_examples.default_options(add_flags=True)
    # Run make_examples with the flags above.
    make_examples_core.make_examples_runner(options)

    # Verify that the variants in the examples are all good.
    examples = self.verify_examples(
        FLAGS.examples, region, options, verify_labels=False
    )

    # Pileup images should have one extra channel.
    self.assertEqual(
        [100, 221, dv_constants.PILEUP_NUM_CHANNELS + 1],
        decode_example(examples[0])['image/shape'],
    )

    # Test there is something in the added channel.
    # Values capture whether each loci has been seen in the observed examples.
    population_matched_loci = {
        'chr20:61539_A': False,
        'chr20:61634_G': False,
        'chr20:61644_G': False,
    }

    for example in examples:
      locus_id = vis.locus_id_from_variant(vis.variant_from_example(example))
      if locus_id in population_matched_loci.keys():
        channels = vis.channels_from_example(example)
        self.assertGreater(
            np.sum(channels[dv_constants.PILEUP_NUM_CHANNELS]),
            0,
            msg='There should be something in the %s-th channel for variant %s'
            % (dv_constants.PILEUP_NUM_CHANNELS + 1, locus_id),
        )
        population_matched_loci[locus_id] = True
    self.assertTrue(
        all(population_matched_loci.values()),
        msg='Check that all 3 sample loci appeared in the examples.',
    )
    # Check against the golden file (same for all modes).
    golden_file = _sharded(testdata.GOLDEN_ALLELE_FREQUENCY_EXAMPLES)
    examples_from_golden = list(tfrecord.read_tfrecords(golden_file))
    self.assertDeepVariantExamplesEqual(examples_from_golden, examples)


  def verify_nist_concordance(self, candidates, nist_variants):
    # Tests that we call almost all of the real variants (according to NIST's
    # Genome in a Bottle callset for NA12878) in our candidate callset.
    # Tests that we don't have an enormous number of FP calls. We should have
    # no more than 5x (arbitrary) more candidate calls than real calls. If we
    # have more it's likely due to some major pipeline problem.
    self.assertLess(len(candidates), 5 * len(nist_variants))
    tp_count = 0
    for nist_variant in nist_variants:
      if self.assertVariantIsPresent(nist_variant, candidates):
        tp_count = tp_count + 1

    self.assertGreater(
        tp_count / len(nist_variants),
        0.983,
        'Recall must be greater than 0.983. TP={}, Truth variants={}'.format(
            tp_count, len(nist_variants)
        ),
    )

  def assertDeepVariantExamplesEqual(self, actual, expected):
    """Asserts that actual and expected tf.Examples from DeepVariant are equal.

    Args:
      actual: iterable of tf.Examples from DeepVariant. DeepVariant examples
        that we want to check.
      expected: iterable of tf.Examples. Expected results for actual.
    """
    self.assertEqual(len(actual), len(expected))
    for i in range(len(actual)):
      actual_example = decode_example(actual[i])
      expected_example = decode_example(expected[i])
      self.assertEqual(actual_example.keys(), expected_example.keys())
      for key in actual_example:
        self.assertEqual(
            actual_example[key], expected_example[key], 'Failed on %s' % key
        )

  def assertVariantIsPresent(self, to_find, variants):
    def variant_key(v):
      return (v.reference_bases, v.start, v.end)

    # Finds a call in our actual call set for each NIST variant, asserting
    # that we found exactly one.
    matches = [
        variant
        for variant in variants
        if variant_key(to_find) == variant_key(variant)
    ]
    if not matches:
      return False

    # Verify that every alt allele appears in the call (but the call might)
    # have more than just those calls.
    for alt in to_find.alternate_bases:
      if alt not in matches[0].alternate_bases:
        return False

    return True

  def verify_candidate_positions(
      self, candidate_positions_paths, candidate_positions_golden_set
  ):
    with epath.Path(candidate_positions_golden_set).open('rb') as my_file:
      positions_golden = np.frombuffer(my_file.read(), dtype=np.int32)
    with epath.Path(candidate_positions_paths).open('rb') as my_file:
      positions = np.frombuffer(my_file.read(), dtype=np.int32)
    logging.warning(
        '%d positions created, %d positions in golden file',
        len(positions),
        len(positions_golden),
    )
    self.assertCountEqual(positions, positions_golden)

  def verify_variants(self, variants, region, options, is_gvcf):
    # Verifies simple properties of the Variant protos in variants. For example,
    # checks that the reference_name() is our expected chromosome. The flag
    # is_gvcf determines how we check the VariantCall field of each variant,
    # enforcing expectations for gVCF records if true or variant calls if false.
    for variant in variants:
      if region:
        self.assertEqual(variant.reference_name, region.reference_name)
        self.assertGreaterEqual(variant.start, region.start)
        self.assertLessEqual(variant.start, region.end)
      self.assertNotEqual(variant.reference_bases, '')
      self.assertNotEmpty(variant.alternate_bases)
      self.assertLen(variant.calls, 1)

      call = variant_utils.only_call(variant)
      self.assertEqual(
          call.call_set_name,
          options.sample_options[0].variant_caller_options.sample_name,
      )
      if is_gvcf:
        # GVCF records should have 0/0 or ./. (un-called) genotypes as they are
        # reference sites, have genotype likelihoods and a GQ value.
        self.assertIn(list(call.genotype), [[0, 0], [-1, -1]])
        self.assertLen(call.genotype_likelihood, 3)
        self.assertGreaterEqual(variantcall_utils.get_gq(call), 0)

  def verify_contiguity(self, contiguous_variants, region):
    """Verifies region is fully covered by gvcf records."""
    # We expect that the intervals cover every base, so the first variant should
    # be at our interval start and the last one should end at our interval end.
    self.assertNotEmpty(contiguous_variants)
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
      # span over a larger interval when it's a deletion and we still produce
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
            options.sample_options[0].variant_caller_options.min_count_snps,
        )

  def sanity_check_example_info_json(self, example, examples_filename, task_id):
    """Checks `example_info.json` w/ examples_filename has the right info."""
    example_info_json = dv_utils.get_example_info_json_filename(
        examples_filename, task_id
    )
    example_info = json.load(gfile.GFile(example_info_json, 'r'))
    self.assertIn('shape', example_info)
    self.assertEqual(
        example_info['shape'], dv_utils.example_image_shape(example)
    )
    self.assertIn('channels', example_info)
    self.assertLen(example_info['channels'], example_info['shape'][2])

  def verify_examples(self, examples_filename, region, options, verify_labels):
    # Do some simple structural checks on the tf.Examples in the file.
    expected_features = [
        'variant/encoded',
        'locus',
        'image/encoded',
        'alt_allele_indices/encoded',
    ]
    if verify_labels:
      expected_features += ['label']

    compression_type = 'GZIP' if 'tfrecord.gz' in examples_filename else None

    examples = list(
        tfrecord.read_tfrecords(
            examples_filename,
            compression_type=compression_type,
        )
    )
    for example in examples:
      for label_feature in expected_features:
        self.assertIn(label_feature, example.features.feature)
      # pylint: disable=g-explicit-length-test
      self.assertNotEmpty(dv_utils.example_alt_alleles_indices(example))
    # Check that the variants in the examples are good.
    variants = [dv_utils.example_variant(x) for x in examples]
    self.verify_variants(variants, region, options, is_gvcf=False)

    if examples:
      self.sanity_check_example_info_json(
          examples[0], examples_filename, options.task_id
      )

    return examples


class DefaultOptionsTest(parameterized.TestCase):

  @flagsaver.flagsaver
  def test_keep_duplicates(self):
    FLAGS.keep_duplicates = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.keep_duplicates, True
    )

  @flagsaver.flagsaver
  def test_keep_supplementary_alignments(self):
    FLAGS.keep_supplementary_alignments = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.keep_supplementary_alignments,
        True,
    )

  @flagsaver.flagsaver
  def test_keep_secondary_alignments(self):
    FLAGS.keep_secondary_alignments = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.keep_secondary_alignments, True
    )

  @flagsaver.flagsaver
  def test_min_base_quality(self):
    FLAGS.min_base_quality = 5
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(options.pic_options.read_requirements.min_base_quality, 5)

  @flagsaver.flagsaver
  def test_min_mapping_quality(self):
    FLAGS.min_mapping_quality = 15
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.CHR20_BAM
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.min_mapping_quality, 15
    )

  @flagsaver.flagsaver
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
        options.sample_options[
            0
        ].variant_caller_options.fraction_reference_sites_to_emit,
        0.3,
    )

  @flagsaver.flagsaver
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
        options.sample_options[
            0
        ].variant_caller_options.fraction_reference_sites_to_emit,
        0.0,
    )

  @flagsaver.flagsaver
  def test_invalid_sequencing_type(self):
    FLAGS.mode = 'training'
    FLAGS.sequencing_type = 'wGs'
    with self.assertRaises(ValueError):
      make_examples.default_options(add_flags=True)

  @parameterized.parameters(
      ({'examples': ('foo', 'foo')},),
      ({'examples': ('foo', 'foo'), 'gvcf': ('bar', 'bar')},),
      ({'examples': ('foo@10', 'foo-00000-of-00010')},),
      ({'task': (0, 0), 'examples': ('foo@10', 'foo-00000-of-00010')},),
      ({'task': (1, 1), 'examples': ('foo@10', 'foo-00001-of-00010')},),
      (
          {
              'task': (1, 1),
              'examples': ('foo@10', 'foo-00001-of-00010'),
              'gvcf': ('bar@10', 'bar-00001-of-00010'),
          },
      ),
      (
          {
              'task': (1, 1),
              'examples': ('foo@10', 'foo-00001-of-00010'),
              'gvcf': ('bar@10', 'bar-00001-of-00010'),
              'candidates': ('baz@10', 'baz-00001-of-00010'),
          },
      ),
  )
  @flagsaver.flagsaver
  def test_sharded_outputs1(self, settings):
    # Set all of the requested flag values.
    for name, (flag_val, _) in settings.items():
      setattr(FLAGS, name, flag_val)

    FLAGS.mode = 'training'
    FLAGS.reads = ''
    FLAGS.ref = ''
    options = make_examples.default_options(add_flags=True)

    # Check all of the flags.
    for name, option_val in [
        ('examples', options.examples_filename),
        ('candidates', options.candidates_filename),
        ('gvcf', options.gvcf_filename),
    ]:
      expected = settings[name][1] if name in settings else ''
      self.assertEqual(expected, option_val)

  @flagsaver.flagsaver
  def test_add_supporting_other_alt_color(self):
    FLAGS.mode = 'training'
    FLAGS.gvcf = ''
    FLAGS.reads = ''
    FLAGS.ref = ''
    FLAGS.examples = ''
    FLAGS.add_supporting_other_alt_color = True
    options = make_examples.default_options(add_flags=True)
    self.assertAlmostEqual(
        options.pic_options.other_allele_supporting_read_alpha, 0.3
    )
    self.assertAlmostEqual(
        options.pic_options.allele_unsupporting_read_alpha, 0.6
    )


class MainTest(parameterized.TestCase):

  def test_catches_bad_argv(self):
    with (
        mock.patch.object(logging, 'error') as mock_logging,
        mock.patch.object(sys, 'exit') as mock_exit,
    ):
      make_examples.main(['make_examples.py', 'extra_arg'])
    mock_logging.assert_called_once_with(
        'Command line parsing failure: make_examples does not accept '
        'positional arguments but some are present on the command line: '
        "\"['make_examples.py', 'extra_arg']\"."
    )
    mock_exit.assert_called_once_with(errno.ENOENT)

  @flagsaver.flagsaver
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

    with (
        mock.patch.object(logging, 'error') as mock_logging,
        mock.patch.object(sys, 'exit') as mock_exit,
    ):
      make_examples.main(['make_examples.py'])
    mock_logging.assert_called_once_with(
        'confident_regions is required when in training mode.'
    )
    mock_exit.assert_called_once_with(errno.ENOENT)


if __name__ == '__main__':
  absltest.main()
