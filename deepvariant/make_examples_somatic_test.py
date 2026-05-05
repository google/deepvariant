# Copyright 2023 Google LLC.
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
import os

from absl import flags
from absl import logging
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from deepvariant import dv_constants
from deepvariant import make_examples_core
from deepvariant import make_examples_somatic
from deepvariant import testdata
from third_party.nucleus.testing import test_utils


FLAGS = flags.FLAGS


def setUpModule():
  logging.set_verbosity(logging.FATAL)
  testdata.init()


class MakeExamplesSomaticEnd2EndTest(parameterized.TestCase):

  @flagsaver.flagsaver
  def test_options_and_sample_names(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )
    self.assertLen(options.sample_options, 2)
    normal_sample_options = options.sample_options[
        make_examples_somatic.NORMAL_SAMPLE_INDEX
    ]
    tumor_sample_options = options.sample_options[1]
    self.assertEqual(normal_sample_options.name, 'NORMAL')
    self.assertEqual(tumor_sample_options.name, 'TUMOR')

  @flagsaver.flagsaver
  def test_make_examples_somatic_end2end_check_calling_examples_suffixes(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = test_utils.test_tmpfile('TEST_SUFFIX.tfrecord.gz')
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    with_normal_suffix = test_utils.test_tmpfile(
        'TEST_SUFFIX_normal.tfrecord.gz'
    )
    with_tumor_suffix = test_utils.test_tmpfile('TEST_SUFFIX_tumor.tfrecord.gz')
    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )
    FLAGS.regions = 'chr20:10,000,000-10,010,000'
    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )
    make_examples_core.make_examples_runner(options)
    # This shows that tumor examples in calling mode are generated without the
    # _tumor suffixes, because it's the only samples generating output.
    self.assertFalse(os.path.exists(with_normal_suffix))
    self.assertFalse(os.path.exists(with_tumor_suffix))
    self.assertTrue(os.path.exists(FLAGS.examples))

  @flagsaver.flagsaver
  def test_tumor_only_flag_options(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.reads_normal = None
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.sample_name_normal = None
    FLAGS.mode = 'calling'
    FLAGS.examples = test_utils.test_tmpfile('TEST_SUFFIX.tfrecord.gz')
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    result = make_examples_somatic.tumor_normal_samples_from_flags(FLAGS)
    self.assertEqual(result[0][0].order, [0])

  def test_tumor_normal_flag_options(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = test_utils.test_tmpfile('TEST_SUFFIX.tfrecord.gz')
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    result = make_examples_somatic.tumor_normal_samples_from_flags(FLAGS)
    self.assertEqual(result[0][1].order, [0, 1])

  @flagsaver.flagsaver
  def test_small_model_path(self):
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.trained_small_model_path = '/path/to/small_model'

    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )

    tumor_sample_options = options.sample_options[1]
    normal_sample_options = options.sample_options[
        make_examples_somatic.NORMAL_SAMPLE_INDEX
    ]

    self.assertEqual(
        tumor_sample_options.small_model_path, '/path/to/small_model'
    )
    self.assertEmpty(normal_sample_options.small_model_path)


class PhaseTumorReadsTrackRefReadsTest(parameterized.TestCase):
  """Tests that --phase_tumor_reads correctly sets track_ref_reads.

  Background:
    There are TWO separate track_ref_reads settings:
    1. AlleleCounterOptions.track_ref_reads (global): Controls whether the
       allele counter tracks reference-supporting reads at each position.
    2. VariantCallerOptions.track_ref_reads (per-sample): Controls whether
       the variant caller populates ref_support in DeepVariantCall protos.

    When --phase_tumor_reads is set, the code explicitly sets BOTH:
    - options.allele_counter_options.track_ref_reads = True  (global)
    - tumor sample_options.variant_caller_options.track_ref_reads = True

    Without the per-sample override, the variant caller would NOT include
    ref_support reads in its output because the global --track_ref_reads
    flag is NOT turned on by --phase_tumor_reads.
  """

  @flagsaver.flagsaver
  def test_phase_tumor_reads_sets_tumor_variant_caller_track_ref_reads(self):
    """The per-sample override is needed because the global flag stays False."""
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.phase_tumor_reads = True

    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )

    normal_opts = options.sample_options[
        make_examples_somatic.NORMAL_SAMPLE_INDEX
    ]
    tumor_opts = options.sample_options[1]

    # Global allele counter track_ref_reads is set to True.
    self.assertTrue(options.allele_counter_options.track_ref_reads)

    # Tumor variant caller MUST have track_ref_reads=True so it populates
    # ref_support in DeepVariantCall protos. Without the explicit override
    # in default_options(), this would be False because the global
    # --track_ref_reads flag was never set.
    self.assertTrue(
        tumor_opts.variant_caller_options.track_ref_reads,
        'Tumor variant_caller_options.track_ref_reads should be True '
        'when --phase_tumor_reads is set. Without this, the variant '
        'caller would discard reference-supporting read information '
        'even though the allele counter is tracking them.',
    )

    # Normal variant caller does NOT need track_ref_reads since we are
    # only phasing tumor reads.
    self.assertFalse(
        normal_opts.variant_caller_options.track_ref_reads,
        'Normal variant_caller_options.track_ref_reads should remain '
        'False when --phase_tumor_reads is set.',
    )

  @flagsaver.flagsaver
  def test_without_phase_tumor_reads_no_track_ref_reads(self):
    """Baseline: without --phase_tumor_reads, neither sample tracks refs."""
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.phase_tumor_reads = False

    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )

    normal_opts = options.sample_options[
        make_examples_somatic.NORMAL_SAMPLE_INDEX
    ]
    tumor_opts = options.sample_options[1]

    # Without --phase_tumor_reads, nothing changes.
    self.assertFalse(options.allele_counter_options.track_ref_reads)
    self.assertFalse(tumor_opts.variant_caller_options.track_ref_reads)
    self.assertFalse(normal_opts.variant_caller_options.track_ref_reads)

  @flagsaver.flagsaver
  def test_phase_tumor_reads_also_sets_phase_reads_and_skip_phasing(self):
    """--phase_tumor_reads enables phasing globally but skips normal."""
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.phase_tumor_reads = True

    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )

    normal_opts = options.sample_options[
        make_examples_somatic.NORMAL_SAMPLE_INDEX
    ]
    tumor_opts = options.sample_options[1]

    # phase_reads is implicitly enabled.
    self.assertTrue(options.phase_reads)

    # Normal sample should skip phasing; tumor should not.
    self.assertTrue(normal_opts.skip_phasing)
    self.assertFalse(tumor_opts.skip_phasing)

  @flagsaver.flagsaver
  def test_phase_reads_and_phase_tumor_reads_mutual_exclusion(self):
    """Setting both --phase_reads and --phase_tumor_reads should error."""
    FLAGS.ref = testdata.CHR20_FASTA
    FLAGS.reads_normal = testdata.CHR20_BAM
    FLAGS.reads_tumor = testdata.CHR20_BAM
    FLAGS.sample_name_normal = 'NORMAL'
    FLAGS.sample_name_tumor = 'TUMOR'
    FLAGS.mode = 'calling'
    FLAGS.examples = ''
    FLAGS.channel_list = ','.join(dv_constants.PILEUP_DEFAULT_CHANNELS)
    FLAGS.phase_tumor_reads = True
    FLAGS.phase_reads = True
    FLAGS.track_ref_reads = True

    options = make_examples_somatic.default_options(
        main_sample_index=1, add_flags=True
    )

    with self.assertRaisesRegex(
        Exception, 'Cannot use both --phase_reads and --phase_tumor_reads'
    ):
      make_examples_somatic.check_options_are_valid(
          options, main_sample_index=1
      )


if __name__ == '__main__':
  absltest.main()
