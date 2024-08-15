# Copyright 2020 Google LLC.
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
"""Tests for deeptrio.make_examples."""

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

from deeptrio import make_examples
from deeptrio import testdata
from deepvariant import dv_constants
from deepvariant import dv_utils
from deepvariant import make_examples_core
from deepvariant.protos import deepvariant_pb2
from tensorflow.python.platform import gfile
from third_party.nucleus.io import fasta
from third_party.nucleus.io import sharded_file_utils
from third_party.nucleus.io import tfrecord
from third_party.nucleus.io import vcf
from third_party.nucleus.protos import reference_pb2
from third_party.nucleus.protos import variants_pb2
from third_party.nucleus.testing import test_utils
from third_party.nucleus.util import ranges
from third_party.nucleus.util import variant_utils
from third_party.nucleus.util import variantcall_utils
from third_party.nucleus.util import vcf_constants

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
    'denovo_label': dv_utils.example_denovo_label,
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


def _make_contigs(specs):
  """Makes ContigInfo protos from specs.

  Args:
    specs: A list of 2- or 3-tuples. All tuples should be of the same length. If
      2-element, these should be the name and length in basepairs of each
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


def _sharded(basename, num_shards=None):
  if num_shards:
    return basename + '@' + str(num_shards)
  else:
    return basename


class MakeExamplesEnd2EndTest(parameterized.TestCase):

  # Golden sets are created with
  # learning/genomics/internal/create_golden_deep_trio.sh
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
  )
  @flagsaver.flagsaver
  def test_make_examples_end2end(
      self, mode, num_shards, labeler_algorithm=None, use_fast_pass_aligner=True
  ):
    self.assertIn(mode, {'calling', 'training', 'candidate_sweep'})
    region = ranges.parse_literal('20:10,000,000-10,010,000')
    FLAGS.write_run_info = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.candidates = test_utils.test_tmpfile(
        _sharded('vsc.tfrecord', num_shards)
    )
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('examples.tfrecord', num_shards)
    )
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    child_examples = test_utils.test_tmpfile(
        _sharded('examples_child.tfrecord', num_shards)
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
      child_gvcf = test_utils.test_tmpfile(
          _sharded('gvcf_child.tfrecord', num_shards)
      )
      child_candidates = test_utils.test_tmpfile(
          _sharded('vsc_child.tfrecord', num_shards)
      )
    else:
      FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
      FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
      child_candidates = test_utils.test_tmpfile(
          _sharded('vsc.tfrecord', num_shards)
      )

    if mode == 'candidate_sweep':
      golden_candidate_positions = _sharded(
          testdata.GOLDEN_CANDIDATE_POSITIONS, num_shards
      )
    for task_id in range(max(num_shards, 1)):
      FLAGS.task = task_id
      options = make_examples.default_options(add_flags=True)
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
            child_candidates,
            proto=deepvariant_pb2.DeepVariantCall,
            compression_type='GZIP',
        ),
        key=lambda c: variant_utils.variant_range_tuple(c.variant),
    )
    self.verify_deepvariant_calls(candidates, options)
    self.verify_variants(
        [call.variant for call in candidates], region, options, is_gvcf=False
    )

    # Verify that the variants in the examples are all good.
    if mode == 'calling':
      examples = self.verify_examples(
          child_examples,
          region,
          options,
          verify_labels=False,
          examples_filename=FLAGS.examples,
      )
    if mode == 'training':
      examples = self.verify_examples(
          FLAGS.examples, region, options, verify_labels=True
      )
    example_variants = [dv_utils.example_variant(ex) for ex in examples]
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
        examples,
        list(tfrecord.read_tfrecords(golden_file, compression_type='GZIP')),
    )

    if mode == 'calling':
      nist_reader = vcf.VcfReader(testdata.TRUTH_VARIANTS_VCF)
      nist_variants = list(nist_reader.query(region))
      self.verify_nist_concordance(example_variants, nist_variants)

      # Check the quality of our generated gvcf file.
      gvcfs = variant_utils.sorted_variants(
          tfrecord.read_tfrecords(
              child_gvcf, proto=variants_pb2.Variant, compression_type='GZIP'
          )
      )
      self.verify_variants(gvcfs, region, options, is_gvcf=True)
      self.verify_contiguity(gvcfs, region)
      gvcf_golden_file = _sharded(
          testdata.GOLDEN_POSTPROCESS_GVCF_INPUT, num_shards
      )
      expected_gvcfs = list(
          tfrecord.read_tfrecords(
              gvcf_golden_file,
              proto=variants_pb2.Variant,
              compression_type='GZIP',
          )
      )
      # Despite its name, assertCountEqual checks that all items are equal.
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

  @parameterized.parameters(
      dict(
          denovo_test=False,
          expected_denovo_variants=0,
      ),
      dict(
          denovo_test=True,
          expected_denovo_variants=3,
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_ont_end2end(
      self,
      denovo_test: bool,
      expected_denovo_variants: int,
  ):
    """Test end to end for long ONT reads with phasing enabled.

    Args:
      denovo_test: If true, denovo parameters will be set.
      expected_denovo_variants: Total number of denovo examples expected.

    This test runs ONT end to end and compares the output with the golden
    output. This test is introduced because previously in training mode the
    non training sample would not be phased. So this now tests to make sure
    all of the training examples are phased correctly.
    """
    region = ranges.parse_literal('chr20:5050000-5075000')
    FLAGS.write_run_info = True
    FLAGS.ref = testdata.GRCH38_CHR0_FASTA
    FLAGS.reads = testdata.ONT_HG002_BAM
    FLAGS.reads_parent1 = testdata.ONT_HG003_BAM
    FLAGS.reads_parent2 = testdata.ONT_HG004_BAM
    FLAGS.confident_regions = testdata.HG002_HIGH_CONFIDENCE_BED
    FLAGS.truth_variants = testdata.HG002_HIGH_CONFIDENCE_VCF
    FLAGS.sample_name = 'HG002'
    FLAGS.sample_name_to_train = 'HG002'
    FLAGS.sample_name_parent1 = 'HG003'
    FLAGS.sample_name_parent2 = 'HG004'
    FLAGS.alt_aligned_pileup = 'diff_channels'
    FLAGS.min_mapping_quality = 1
    FLAGS.mode = 'training'
    FLAGS.parse_sam_aux_fields = True
    FLAGS.partition_size = 25000
    FLAGS.phase_reads = True
    FLAGS.pileup_image_height_child = 100
    FLAGS.pileup_image_height_parent = 100
    FLAGS.pileup_image_width = 199
    FLAGS.realign_reads = False
    FLAGS.skip_parent_calling = True
    FLAGS.sort_by_haplotypes = True
    FLAGS.track_ref_reads = True
    FLAGS.vsc_min_fraction_indels = 0.12
    FLAGS.vsc_min_fraction_snps = 0.1
    num_shards = 0
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('examples.tfrecord', num_shards)
    )
    FLAGS.channel_list = ','.join(
        dv_constants.PILEUP_DEFAULT_CHANNELS + ['haplotype']
    )
    FLAGS.regions = [ranges.to_literal(region)]
    golden_file = _sharded(testdata.GOLDEN_ONT_MAKE_EXAMPLES_OUTPUT, num_shards)
    FLAGS.denovo_regions = None
    if denovo_test:
      # If denovo test is enabled, then set the parameters for denovo testing.
      golden_file = _sharded(
          testdata.GOLDEN_ONT_DENOVO_MAKE_EXAMPLES_OUTPUT, num_shards
      )
      FLAGS.write_run_info = True
      FLAGS.denovo_regions = testdata.HG002_DENOVO_BED

    for task_id in range(max(num_shards, 1)):
      FLAGS.task = task_id
      options = make_examples.default_options(add_flags=True)
      make_examples_core.make_examples_runner(options)

      examples = self.verify_examples(
          FLAGS.examples, region, options, verify_labels=True
      )

      self.assertDeepVariantExamplesEqual(
          examples,
          list(tfrecord.read_tfrecords(golden_file, compression_type='GZIP')),
      )
      if denovo_test:
        # Check total number of denovo examples.
        total_denovo = sum(
            [
                1
                for example in examples
                if dv_utils.example_denovo_label(example)
            ]
        )
        self.assertEqual(
            total_denovo,
            expected_denovo_variants,
            msg='ONT denovo golden test: denovo variants count.',
        )
        # Read the runinfo file
        runinfo = make_examples_core.read_make_examples_run_info(
            FLAGS.examples + '.run_info.pbtxt'
        )
        golden_runinfo = make_examples_core.read_make_examples_run_info(
            testdata.GOLDEN_ONT_DENOVO_MAKE_EXAMPLES_OUTPUT + '.run_info.pbtxt'
        )
        self.assertEqual(
            runinfo.stats.num_examples,
            golden_runinfo.stats.num_examples,
            msg='ONT denovo golden test: Run info comparison num_examples.',
        )
        self.assertEqual(
            runinfo.stats.num_denovo,
            golden_runinfo.stats.num_denovo,
            msg='ONT denovo golden test: Run info comparison num_denovo.',
        )
        self.assertEqual(
            runinfo.stats.num_nondenovo,
            golden_runinfo.stats.num_nondenovo,
            msg='ONT denovo golden test: Run info comparison num_nondenovo.',
        )

  # Golden sets are created with learning/genomics/internal/create_golden.sh
  @flagsaver.flagsaver
  def test_make_examples_training_end2end_with_customized_classes_labeler(self):
    FLAGS.labeler_algorithm = 'customized_classes_labeler'
    FLAGS.customized_classes_labeler_classes_list = 'ref,class1,class2'
    FLAGS.customized_classes_labeler_info_field_name = 'type'
    region = ranges.parse_literal('20:10,000,000-10,004,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.candidates = test_utils.test_tmpfile(_sharded('vsc.tfrecord'))
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
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
        examples,
        list(tfrecord.read_tfrecords(golden_file, compression_type='GZIP')),
    )

  # Golden sets are created with learning/genomics/internal/create_golden.sh
  @flagsaver.flagsaver
  def test_make_examples_training_end2end_with_alt_aligned_pileup(self):
    region = ranges.parse_literal('20:10,000,000-10,010,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.candidates = test_utils.test_tmpfile(_sharded('vsc.tfrecord'))
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.partition_size = 1000
    FLAGS.mode = 'training'
    FLAGS.gvcf_gq_binsize = 5

    # The following 4 lines are added.
    FLAGS.alt_aligned_pileup = 'diff_channels'
    FLAGS.pileup_image_height_child = 60
    FLAGS.pileup_image_height_parent = 40
    FLAGS.pileup_image_width = 199

    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    golden_file = _sharded(testdata.ALT_ALIGNED_PILEUP_GOLDEN_TRAINING_EXAMPLES)
    # Verify that the variants in the examples are all good.
    examples = self.verify_examples(
        FLAGS.examples, region, options, verify_labels=True
    )
    self.assertDeepVariantExamplesEqual(
        examples,
        list(tfrecord.read_tfrecords(golden_file, compression_type='GZIP')),
    )
    # Pileup image should now have 8 channels.
    # Height should be 60 + 40 * 2 = 140.
    self.assertEqual(decode_example(examples[0])['image/shape'], [140, 199, 8])

  @flagsaver.flagsaver
  def test_make_examples_compare_realignment_modes(self):
    def _run_with_realignment_mode(enable_joint_realignment, name):
      FLAGS.enable_joint_realignment = enable_joint_realignment
      region = ranges.parse_literal('20:10,000,000-10,010,000')
      FLAGS.ref = testdata.CHR20_FASTA
      FLAGS.reads = testdata.HG001_CHR20_BAM
      FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
      FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
      FLAGS.sample_name = 'child'
      FLAGS.sample_name_to_train = 'child'
      FLAGS.sample_name_parent1 = 'parent1'
      FLAGS.sample_name_parent2 = 'parent2'
      FLAGS.candidates = test_utils.test_tmpfile(f'{name}.vsc.tfrecord')
      FLAGS.examples = test_utils.test_tmpfile(f'{name}.examples.tfrecord')
      FLAGS.channel_list = ','.join(
          dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE
      )
      child_examples = test_utils.test_tmpfile(
          f'{name}_child.examples.tfrecord'
      )
      FLAGS.regions = [ranges.to_literal(region)]
      FLAGS.partition_size = 1000
      FLAGS.mode = 'calling'
      FLAGS.gvcf = test_utils.test_tmpfile(f'{name}.gvcf.tfrecord')
      # child_gvcf = test_utils.test_tmpfile(f'{name}.gvcf_child.tfrecord')
      # child_candidates = test_utils.test_tmpfile(f'{name}.vsc_child.tfrecord')
      options = make_examples.default_options(add_flags=True)
      make_examples_core.make_examples_runner(options)

      examples = self.verify_examples(
          child_examples,
          region,
          options,
          verify_labels=False,
          examples_filename=FLAGS.examples,
      )
      return examples

    examples1 = _run_with_realignment_mode(False, 'ex1')
    examples2 = _run_with_realignment_mode(True, 'ex2')
    self.assertNotEmpty(examples1)
    self.assertNotEmpty(examples2)
    # The assumption is just that these two lists of examples should be
    # different. In this case, it happens to be that we got different numbers
    # of examples:
    self.assertNotEmpty(examples1)
    self.assertDeepVariantExamplesNotEqual(examples1, examples2)

  @parameterized.parameters(
      dict(select_types=None, expected_count=79),
      dict(select_types='all', expected_count=79),
      dict(select_types='snps', expected_count=64),
      dict(select_types='indels', expected_count=12),
      dict(select_types='snps indels', expected_count=76),
      dict(select_types='multi-allelics', expected_count=3),
      dict(select_types=None, keep_legacy_behavior=True, expected_count=79),
      dict(select_types='all', keep_legacy_behavior=True, expected_count=79),
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
          expected_count=4,
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_with_variant_selection(
      self, select_types, expected_count, keep_legacy_behavior=False
  ):
    if select_types is not None:
      FLAGS.select_variant_types = select_types
    region = ranges.parse_literal('20:10,000,000-10,010,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.candidates = test_utils.test_tmpfile(_sharded('vsc.tfrecord'))
    child_candidates = test_utils.test_tmpfile(_sharded('vsc_child.tfrecord'))
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    FLAGS.partition_size = 1000
    FLAGS.mode = 'calling'
    FLAGS.keep_legacy_allele_counter_behavior = keep_legacy_behavior

    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)

    candidates = list(
        tfrecord.read_tfrecords(child_candidates, compression_type='GZIP')
    )
    self.assertLen(candidates, expected_count)

  @parameterized.parameters(
      dict(
          mode='calling', which_parent='parent1', sample_name_to_train='child'
      ),
      dict(
          mode='calling', which_parent='parent2', sample_name_to_train='child'
      ),
      dict(
          mode='training', which_parent='parent1', sample_name_to_train='child'
      ),
      dict(
          mode='training', which_parent='parent2', sample_name_to_train='child'
      ),
      dict(
          mode='calling', which_parent='parent1', sample_name_to_train='parent1'
      ),
      dict(
          mode='training',
          which_parent='parent1',
          sample_name_to_train='parent1',
      ),
      # Training on parent2 in a duo is not supported (with a clear error
      # message).
  )
  @flagsaver.flagsaver
  def test_make_examples_training_end2end_duos(
      self, mode, which_parent, sample_name_to_train
  ):
    region = ranges.parse_literal('20:10,000,000-10,010,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.examples = test_utils.test_tmpfile(_sharded('examples.tfrecord'))
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    FLAGS.partition_size = 1000

    FLAGS.mode = mode
    if mode == 'training':
      FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
      FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED

    if which_parent == 'parent1':
      FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
      FLAGS.sample_name_parent1 = 'parent1'
    elif which_parent == 'parent2':
      FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
      FLAGS.sample_name_parent2 = 'parent2'
    else:
      raise ValueError('Invalid `which_parent` value in test case.')
    FLAGS.sample_name_to_train = sample_name_to_train

    # This is only a simple test that it runs without errors.
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)

  @parameterized.parameters(
      dict(mode='calling'),
      dict(mode='training'),
  )
  @flagsaver.flagsaver
  def test_make_examples_end2end_vcf_candidate_importer(self, mode):
    FLAGS.variant_caller = 'vcf_candidate_importer'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.pileup_image_height_parent = 40
    FLAGS.pileup_image_height_child = 60
    FLAGS.candidates = test_utils.test_tmpfile(
        _sharded('vcf_candidate_importer.candidates.{}.tfrecord'.format(mode))
    )
    FLAGS.examples = test_utils.test_tmpfile(
        _sharded('vcf_candidate_importer.examples.{}.tfrecord'.format(mode))
    )
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    FLAGS.mode = mode
    FLAGS.regions = '20:10,000,000-10,010,000'

    if mode == 'calling':
      golden_file = _sharded(
          testdata.GOLDEN_VCF_CANDIDATE_IMPORTER_CALLING_EXAMPLES_CHILD
      )
      path_to_output_examples = test_utils.test_tmpfile(
          _sharded(
              'vcf_candidate_importer_child.examples.{}.tfrecord'.format(mode)
          )
      )
      FLAGS.proposed_variants_child = testdata.TRUTH_VARIANTS_VCF
      FLAGS.proposed_variants_parent1 = testdata.TRUTH_VARIANTS_VCF
      FLAGS.proposed_variants_parent2 = testdata.TRUTH_VARIANTS_VCF
    else:
      golden_file = _sharded(
          testdata.GOLDEN_VCF_CANDIDATE_IMPORTER_TRAINING_EXAMPLES
      )
      path_to_output_examples = test_utils.test_tmpfile(
          _sharded('vcf_candidate_importer.examples.{}.tfrecord'.format(mode))
      )
      FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
      FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED

    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    # Verify that the variants in the examples are all good.
    output_examples_to_compare = self.verify_examples(
        path_to_output_examples,
        None,
        options,
        verify_labels=mode == 'training',
        examples_filename=FLAGS.examples,
    )
    self.assertDeepVariantExamplesEqual(
        output_examples_to_compare,
        list(tfrecord.read_tfrecords(golden_file, compression_type='GZIP')),
    )

  @parameterized.parameters(
      dict(
          max_reads_per_partition=1500,
          expected_len_examples1=88,
          expected_len_examples2=32,
      ),
      dict(
          max_reads_per_partition=8,
          expected_len_examples1=34,
          expected_len_examples2=30,
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_with_max_reads_for_dynamic_bases_per_region(
      self,
      max_reads_per_partition,
      expected_len_examples1,
      expected_len_examples2,
  ):
    region = ranges.parse_literal('20:10,000,000-10,010,000')
    FLAGS.regions = [ranges.to_literal(region)]
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.examples = test_utils.test_tmpfile(_sharded('ex.tfrecord'))
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    child_examples = test_utils.test_tmpfile(_sharded('ex_child.tfrecord'))
    FLAGS.partition_size = 1000
    FLAGS.mode = 'calling'
    FLAGS.max_reads_per_partition = max_reads_per_partition

    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    examples1 = self.verify_examples(
        child_examples,
        region,
        options,
        verify_labels=False,
        examples_filename=FLAGS.examples,
    )
    self.assertLen(examples1, expected_len_examples1)
    # Now, this is the main part of the test. I want to test the behavior after
    # I set max_reads_for_dynamic_bases_per_region.
    FLAGS.max_reads_for_dynamic_bases_per_region = 1
    options = make_examples.default_options(add_flags=True)
    make_examples_core.make_examples_runner(options)
    examples2 = self.verify_examples(
        child_examples,
        region,
        options,
        verify_labels=False,
        examples_filename=FLAGS.examples,
    )
    self.assertLen(examples2, expected_len_examples2)

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
        0.9705,
        'Recall must be greater than 0.9705. TP={}, Truth variants={}'.format(
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
      self.assertEqual(decode_example(actual[i]), decode_example(expected[i]))

  def assertDeepVariantExamplesNotEqual(self, actual, expected):
    """Asserts that actual and expected tf.Examples are not equal.

    Args:
      actual: iterable of tf.Examples from DeepVariant. DeepVariant examples
        that we want to check.
      expected: iterable of tf.Examples. Expected results for actual.
    """
    pass_not_equal_check = False
    if len(actual) != len(expected):
      logging.warning(
          (
              'In assertDeepVariantExamplesNotEqual: '
              'actual(%d) and expected(%d) has different lengths'
          ),
          len(actual),
          len(expected),
      )
      pass_not_equal_check = True
    min_size = min(len(actual), len(expected))
    for i in range(min_size):
      if decode_example(actual[i]) != decode_example(expected[i]):
        logging.warning(
            (
                'assertDeepVariantExamplesNotEqual: '
                'actual example[%d] and expected example[%d] '
                'are different'
            ),
            i,
            i,
        )
        pass_not_equal_check = True
    self.assertTrue(
        pass_not_equal_check,
        (
            'assertDeepVariantExamplesNotEqual failed - '
            'actual and expected examples are identical.'
        ),
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
          options.sample_options[1].variant_caller_options.sample_name,
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
            options.sample_options[1].variant_caller_options.min_count_snps,
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

  def verify_examples(
      self,
      path_to_output_examples,
      region,
      options,
      verify_labels,
      examples_filename=None,
  ):
    # Do some simple structural checks on the tf.Examples in the file.
    expected_features = [
        'variant/encoded',
        'locus',
        'image/encoded',
        'alt_allele_indices/encoded',
    ]
    if verify_labels:
      expected_features += ['label']

    examples = list(
        tfrecord.read_tfrecords(
            path_to_output_examples, compression_type='GZIP'
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

    # In DeepTrio, path_to_output_examples can be pointing to the ones with
    # the suffixes (such as _child). In that case, we pass in the original
    # examples path to the `examples_filename` arg.
    # If `examples_filename` arg, directly use `path_to_output_examples`.
    if examples:
      if examples_filename is None:
        examples_filename = path_to_output_examples
      self.sanity_check_example_info_json(
          examples[0], examples_filename, options.task_id
      )
    return examples


class MakeExamplesUnitTest(parameterized.TestCase):

  def test_read_write_run_info(self):
    def _read_lines(path):
      with open(path) as fin:
        return list(fin.readlines())

    golden_actual = make_examples_core.read_make_examples_run_info(
        testdata.GOLDEN_MAKE_EXAMPLES_RUN_INFO
    )
    # We don't really want to inject too much knowledge about the golden right
    # here, so we only use a minimal test that (a) the run_info_filename is
    # a non-empty string and (b) the number of candidates sites in the labeling
    # metrics field is greater than 0. Any reasonable golden output will have at
    # least one candidate variant, and the reader should have filled in the
    # value.
    self.assertNotEmpty(golden_actual.options.run_info_filename)
    self.assertEqual(
        golden_actual.labeling_metrics.n_candidate_variant_sites,
        testdata.N_GOLDEN_TRAINING_EXAMPLES,
    )

    # Check that reading + writing the data produces the same lines:
    tmp_output = test_utils.test_tmpfile('written_run_info.pbtxt')
    make_examples_core.write_make_examples_run_info(golden_actual, tmp_output)
    print('*' * 100)
    print(_read_lines(tmp_output))
    print('*' * 100)
    self.assertEqual(
        _read_lines(testdata.GOLDEN_MAKE_EXAMPLES_RUN_INFO),
        _read_lines(tmp_output),
    )

  @flagsaver.flagsaver
  def test_keep_duplicates(self):
    FLAGS.keep_duplicates = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.keep_duplicates, True
    )

  @flagsaver.flagsaver
  def test_keep_supplementary_alignments(self):
    FLAGS.keep_supplementary_alignments = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.keep_supplementary_alignments,
        True,
    )

  @flagsaver.flagsaver
  def test_keep_secondary_alignments(self):
    FLAGS.keep_secondary_alignments = True
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.keep_secondary_alignments, True
    )

  @flagsaver.flagsaver
  def test_min_base_quality(self):
    FLAGS.min_base_quality = 5
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(options.pic_options.read_requirements.min_base_quality, 5)

  @flagsaver.flagsaver
  def test_min_mapping_quality(self):
    FLAGS.min_mapping_quality = 15
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    options = make_examples.default_options(add_flags=True)
    self.assertEqual(
        options.pic_options.read_requirements.min_mapping_quality, 15
    )

  @flagsaver.flagsaver
  def test_default_options_with_training_random_emit_ref_sites(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)

    FLAGS.training_random_emit_ref_sites = 0.3
    options = make_examples.default_options(add_flags=True)
    self.assertAlmostEqual(
        options.sample_options[
            1
        ].variant_caller_options.fraction_reference_sites_to_emit,
        0.3,
    )

  @flagsaver.flagsaver
  def test_default_options_without_training_random_emit_ref_sites(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)

    options = make_examples.default_options(add_flags=True)
    # In proto3, there is no way to check presence of scalar field:
    # redacted
    # As an approximation, we directly check that the value should be exactly 0.
    self.assertEqual(
        options.sample_options[
            1
        ].variant_caller_options.fraction_reference_sites_to_emit,
        0.0,
    )

  @flagsaver.flagsaver
  def test_confident_regions(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.truth_variants = testdata.TRUTH_VARIANTS_VCF
    FLAGS.confident_regions = testdata.CONFIDENT_REGIONS_BED
    FLAGS.mode = 'training'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)

    options = make_examples.default_options(add_flags=True)
    confident_regions = make_examples_core.read_confident_regions(options)

    # Our expected intervals, inlined from CONFIDENT_REGIONS_BED.
    expected = _from_literals_list([
        '20:10000847-10002407',
        '20:10002521-10004171',
        '20:10004274-10004964',
        '20:10004995-10006386',
        '20:10006410-10007800',
        '20:10007825-10008018',
        '20:10008044-10008079',
        '20:10008101-10008707',
        '20:10008809-10008897',
        '20:10009003-10009791',
        '20:10009934-10010531',
    ])
    # Our confident regions should be exactly those found in the BED file.
    self.assertCountEqual(expected, list(confident_regions))

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
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
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
    region = ranges.parse_literal('20:10,000,000-10,010,000')
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.candidates = test_utils.test_tmpfile('vsc.tfrecord')
    FLAGS.examples = test_utils.test_tmpfile('examples.tfrecord')
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
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

  @flagsaver.flagsaver
  def test_regions_and_exclude_regions_flags_with_trio_options(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.regions = '20:10,000,000-11,000,000'
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    FLAGS.exclude_regions = '20:10,010,000-10,100,000'

    options = make_examples.default_options(add_flags=True)
    _, _, regions_from_options = (
        make_examples_core.processing_regions_from_options(options)
    )
    self.assertCountEqual(
        list(ranges.RangeSet(regions_from_options)),
        _from_literals_list(
            ['20:10,000,000-10,009,999', '20:10,100,001-11,000,000']
        ),
    )

  @flagsaver.flagsaver
  def test_incorrect_empty_regions_with_trio_options(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    # Deliberately incorrect contig name.
    FLAGS.regions = 'xxx20:10,000,000-11,000,000'
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)

    options = make_examples.default_options(add_flags=True)
    with self.assertRaisesRegex(ValueError, 'The regions to call is empty.'):
      make_examples_core.processing_regions_from_options(options)


class RegionProcessorTest(parameterized.TestCase):

  def setUp(self):
    super(RegionProcessorTest, self).setUp()
    self.region = ranges.parse_literal('20:10,000,000-10,000,100')

    FLAGS.reads = ''
    self.options = make_examples.default_options(add_flags=False)
    self.options.reference_filename = testdata.CHR20_FASTA
    self.options.truth_variants_filename = testdata.TRUTH_VARIANTS_VCF
    self.options.mode = deepvariant_pb2.MakeExamplesOptions.TRAINING

    self.ref_reader = fasta.IndexedFastaReader(self.options.reference_filename)
    self.default_shape = [5, 5, 7]
    self.processor = make_examples_core.RegionProcessor(self.options)
    self.mock_init = self.add_mock('_initialize')
    for sample in self.processor.samples:
      sample.in_memory_sam_reader = mock.Mock()

  def add_mock(self, name, retval='dontadd', side_effect='dontadd'):
    patcher = mock.patch.object(self.processor, name, autospec=True)
    self.addCleanup(patcher.stop)
    mocked = patcher.start()
    if retval != 'dontadd':
      mocked.return_value = retval
    if side_effect != 'dontadd':
      mocked.side_effect = side_effect
    return mocked

  @parameterized.parameters([
      deepvariant_pb2.MakeExamplesOptions.TRAINING,
      deepvariant_pb2.MakeExamplesOptions.CALLING,
  ])
  def test_process_keeps_ordering_of_candidates_and_examples(self, mode):
    self.processor.options.mode = mode

    r1, r2 = mock.Mock(), mock.Mock()
    c1, c2 = mock.Mock(), mock.Mock()
    self.add_mock('region_reads_norealign', retval=[r1, r2])
    self.add_mock('candidates_in_region', retval=({'child': [c1, c2]}, {}))
    candidates_dict, gvcfs_dict, runtimes = self.processor.process(self.region)
    self.assertEqual({'child': [c1, c2]}, candidates_dict)
    self.assertEqual({}, gvcfs_dict)
    self.assertIsInstance(runtimes, dict)

    in_memory_sam_reader = self.processor.samples[1].in_memory_sam_reader
    in_memory_sam_reader.replace_reads.assert_called_once_with([r1, r2])

  @flagsaver.flagsaver
  def test_use_original_quality_scores_without_parse_sam_aux_fields(self):
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)
    FLAGS.use_original_quality_scores = True
    FLAGS.parse_sam_aux_fields = False

    with self.assertRaisesRegex(
        Exception,
        (
            'If --use_original_quality_scores is set then '
            '--parse_sam_aux_fields must be set too.'
        ),
    ):
      make_examples.default_options(add_flags=True)

  @parameterized.parameters(
      dict(height_parent=10, height_child=9),
      dict(height_parent=9, height_child=10),
      dict(height_parent=150, height_child=101),
      dict(height_parent=101, height_child=170),
  )
  @flagsaver.flagsaver
  def test_image_heights(self, height_parent, height_child):
    FLAGS.pileup_image_height_parent = height_parent
    FLAGS.pileup_image_height_child = height_child
    FLAGS.mode = 'calling'
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads = testdata.HG001_CHR20_BAM
    FLAGS.reads_parent1 = testdata.NA12891_CHR20_BAM
    FLAGS.reads_parent2 = testdata.NA12892_CHR20_BAM
    FLAGS.sample_name = 'child'
    FLAGS.sample_name_to_train = 'child'
    FLAGS.sample_name_parent1 = 'parent1'
    FLAGS.sample_name_parent2 = 'parent2'
    FLAGS.examples = 'examples.tfrecord'
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_CHANNELS_WITH_INSERT_SIZE)

    options = make_examples.default_options(add_flags=True)
    with self.assertRaisesRegex(
        Exception, 'Total pileup image heights must be between 75-362.'
    ):
      make_examples.check_options_are_valid(options)


if __name__ == '__main__':
  absltest.main()
