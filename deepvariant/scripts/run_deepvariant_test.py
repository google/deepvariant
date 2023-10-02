# Copyright 2019 Google LLC.
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
"""Tests for deepvariant .run_deepvariant."""

import io
from unittest import mock

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from deepvariant.scripts import run_deepvariant

FLAGS = flags.FLAGS


class RunDeepvariantTest(parameterized.TestCase):

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters(
      (model_type, use_slim_model)
      for model_type in ['WGS', 'WES', 'HYBRID_PACBIO_ILLUMINA']
      for use_slim_model in [False, True]
  )
  # pylint: enable=g-complex-comprehension
  @flagsaver.flagsaver
  def test_basic_commands(self, model_type, use_slim_model):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.use_slim_model = use_slim_model
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    extra_channel = ''
    if model_type == 'WGS' or model_type == 'WES':
      extra_channel = '--channels "insert_size" '

    extra_args_plus_gvcf = (
        '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
    )

    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '%s'
        '%s'
        '--task {}' % (extra_channel, extra_args_plus_gvcf),
    )
    call_variants_bin = (
        'call_variants_slim' if use_slim_model else 'call_variants'
    )
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/{} '
        '--outfile '
        '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '--checkpoint "/opt/models/{}"'.format(
            call_variants_bin, model_type.lower()
        ),
    )
    self.assertEqual(
        commands[2][0],
        (
            'time /opt/deepvariant/bin/postprocess_variants '
            '--ref "your_ref" '
            '--infile '
            '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
            '--outfile "your_vcf" '
            '--gvcf_outfile "your_gvcf" '
            '--nonvariant_site_tfrecord_path '
            '"/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
        ),
    )

  @flagsaver.flagsaver
  def test_pacbio_basic_commands(self):
    FLAGS.model_type = 'PACBIO'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.regions = None
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    # --gvcf is added in the middle because the flags are sorted
    # alphabetically. PacBio flags here is now mixed with some others.
    hp_args = (
        '--max_reads_per_partition "600" '
        + '--min_mapping_quality "1" '
        + '--parse_sam_aux_fields '
        + '--partition_size "25000" '
        + '--phase_reads '
        + '--pileup_image_width "199" '
        + '--norealign_reads '
        + '--sort_by_haplotypes '
        + '--track_ref_reads'
    )
    extra_args_plus_gvcf = (
        '--add_hp_channel '
        '--alt_aligned_pileup "diff_channels" '
        '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
        '{} --vsc_min_fraction_indels "0.12" '.format(hp_args)
    )
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '%s'
        '--task {}' % extra_args_plus_gvcf,
    )
    self.assertEqual(
        commands[1][0],
        (
            'time /opt/deepvariant/bin/call_variants --outfile'
            ' "/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
            ' --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz"'
            ' --checkpoint "/opt/models/pacbio"'
        ),
    )
    self.assertEqual(
        commands[2][0],
        (
            'time /opt/deepvariant/bin/postprocess_variants '
            '--ref "your_ref" '
            '--infile '
            '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
            '--outfile "your_vcf" '
            '--gvcf_outfile "your_gvcf" '
            '--nonvariant_site_tfrecord_path '
            '"/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
        ),
    )

  @flagsaver.flagsaver
  def test_ont_basic_commands(self):
    FLAGS.model_type = 'ONT_R104'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.regions = None
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    # --gvcf is added in the middle because the flags are sorted
    # alphabetically. ONT flags here is now mixed with some others.
    hp_args = (
        '--max_reads_per_partition "600" '
        + '--min_mapping_quality "5" '
        + '--parse_sam_aux_fields '
        + '--partition_size "25000" '
        + '--phase_reads '
        + '--pileup_image_width "199" '
        + '--norealign_reads '
        + '--sort_by_haplotypes '
        + '--track_ref_reads'
    )
    extra_args_plus_gvcf = (
        '--add_hp_channel '
        '--alt_aligned_pileup "diff_channels" '
        '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
        '{} --vsc_min_fraction_indels "0.12" --vsc_min_fraction_snps "0.08" '
        .format(hp_args)
    )
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '%s'
        '--task {}' % extra_args_plus_gvcf,
    )
    self.assertEqual(
        commands[1][0],
        (
            'time /opt/deepvariant/bin/call_variants --outfile'
            ' "/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
            ' --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz"'
            ' --checkpoint "/opt/models/ont_r104"'
        ),
    )
    self.assertEqual(
        commands[2][0],
        (
            'time /opt/deepvariant/bin/postprocess_variants '
            '--ref "your_ref" '
            '--infile '
            '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
            '--outfile "your_vcf" '
            '--gvcf_outfile "your_gvcf" '
            '--nonvariant_site_tfrecord_path '
            '"/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
        ),
    )

  @parameterized.parameters(None, 'your_sample_name')
  @flagsaver.flagsaver
  def test_sample_name_command(self, sample_name):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = None
    FLAGS.num_shards = 64
    FLAGS.regions = None
    FLAGS.sample_name = sample_name
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    extra_sample_name_flag = ''
    if FLAGS.sample_name:
      extra_sample_name_flag = ' --sample_name "your_sample_name"'

    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '--channels "insert_size"'
        '%s'
        ' --task {}' % extra_sample_name_flag,
    )
    self.assertEqual(
        commands[1][0],
        (
            'time /opt/deepvariant/bin/call_variants --outfile'
            ' "/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
            ' --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz"'
            ' --checkpoint "/opt/models/wgs"'
        ),
    )
    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/postprocess_variants '
        '--ref "your_ref" '
        '--infile '
        '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
        '--outfile "your_vcf"'
        '%s' % extra_sample_name_flag,
    )

  @parameterized.parameters(
      ('keep_secondary_alignments=true', '--keep_secondary_alignments'),
      ('keep_secondary_alignments=false', '--nokeep_secondary_alignments'),
      (
          'keep_secondary_alignments=true,keep_supplementary_alignments=true',
          '--keep_secondary_alignments --keep_supplementary_alignments',
      ),
      (
          'use_ref_for_cram=true,keep_secondary_alignments=true,'
          + 'keep_supplementary_alignments=false',
          '--keep_secondary_alignments --nokeep_supplementary_alignments '
          + '--use_ref_for_cram',
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_extra_args_boolean(
      self, make_examples_extra_args, expected_args
  ):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.make_examples_extra_args = make_examples_extra_args
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '--channels "insert_size" '
        '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
        '%s '
        '--task {}' % expected_args,
    )

  @parameterized.parameters(
      (
          None,
          '--add_hp_channel '
          + '--alt_aligned_pileup "diff_channels" '
          + '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
          + '--max_reads_per_partition "600" '
          + '--min_mapping_quality "1" '
          + '--parse_sam_aux_fields '
          + '--partition_size "25000" '
          + '--phase_reads '
          + '--pileup_image_width "199" '
          + '--norealign_reads '
          + '--sort_by_haplotypes '
          + '--track_ref_reads '
          + '--vsc_min_fraction_indels "0.12" ',
          '',
      ),
      (
          'alt_aligned_pileup="rows",vsc_min_fraction_indels=0.03',
          '--add_hp_channel '
          + '--alt_aligned_pileup "rows" '
          + '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
          + '--max_reads_per_partition "600" '
          + '--min_mapping_quality "1" '
          + '--parse_sam_aux_fields '
          + '--partition_size "25000" '
          + '--phase_reads '
          + '--pileup_image_width "199" '
          + '--norealign_reads '
          + '--sort_by_haplotypes '
          + '--track_ref_reads '
          + '--vsc_min_fraction_indels "0.03" ',
          '\nWarning: --alt_aligned_pileup is previously set to diff_channels, '
          + 'now to "rows".\n'
          + '\nWarning: --vsc_min_fraction_indels is previously set to 0.12, '
          + 'now to 0.03.\n',
      ),
  )
  @flagsaver.flagsaver
  def test_pacbio_args_overwrite(
      self, make_examples_extra_args, expected_args, expected_stdout
  ):
    """Confirms that adding extra flags can overwrite the default from mode."""
    FLAGS.model_type = 'PACBIO'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.regions = None
    FLAGS.make_examples_extra_args = make_examples_extra_args
    with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
      commands = run_deepvariant.create_all_commands_and_logfiles(
          '/tmp/deepvariant_tmp_output'
      )
      self.assertEqual(mock_stdout.getvalue(), expected_stdout)
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '%s'
        '--task {}' % expected_args,
    )

  @parameterized.parameters(
      (
          None,
          '--add_hp_channel '
          + '--alt_aligned_pileup "diff_channels" '
          + '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
          + '--max_reads_per_partition "600" '
          + '--min_mapping_quality "5" '
          + '--parse_sam_aux_fields '
          + '--partition_size "25000" '
          + '--phase_reads '
          + '--pileup_image_width "199" '
          + '--norealign_reads '
          + '--sort_by_haplotypes '
          + '--track_ref_reads '
          + '--vsc_min_fraction_indels "0.12" '
          + '--vsc_min_fraction_snps "0.08" ',
          '',
      ),
      (
          'alt_aligned_pileup="rows",vsc_min_fraction_indels=0.03',
          '--add_hp_channel '
          + '--alt_aligned_pileup "rows" '
          + '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
          + '--max_reads_per_partition "600" '
          + '--min_mapping_quality "5" '
          + '--parse_sam_aux_fields '
          + '--partition_size "25000" '
          + '--phase_reads '
          + '--pileup_image_width "199" '
          + '--norealign_reads '
          + '--sort_by_haplotypes '
          + '--track_ref_reads '
          + '--vsc_min_fraction_indels "0.03" '
          + '--vsc_min_fraction_snps "0.08" ',
          '\nWarning: --alt_aligned_pileup is previously set to diff_channels, '
          + 'now to "rows".\n'
          + '\nWarning: --vsc_min_fraction_indels is previously set to 0.12, '
          + 'now to 0.03.\n',
      ),
  )
  @flagsaver.flagsaver
  def test_ont_args_overwrite(
      self, make_examples_extra_args, expected_args, expected_stdout
  ):
    """Confirms that adding extra flags can overwrite the default from mode."""
    FLAGS.model_type = 'ONT_R104'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.regions = None
    FLAGS.make_examples_extra_args = make_examples_extra_args
    with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
      commands = run_deepvariant.create_all_commands_and_logfiles(
          '/tmp/deepvariant_tmp_output'
      )
      self.assertEqual(mock_stdout.getvalue(), expected_stdout)
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '%s'
        '--task {}' % expected_args,
    )

  @flagsaver.flagsaver
  def test_logging_dir(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/deepvariant_tmp_output/LOGDIR'
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        commands[0][0],
        (
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples --mode calling --ref'
            ' "your_ref" --reads "your_bam" --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz"'
            ' --channels "insert_size" --gvcf'
            ' "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" --task {}'
        ),
    )

  @flagsaver.flagsaver
  def test_with_runtime_report(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.report_title = 'custom title'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/deepvariant_tmp_output/LOGDIR'
    FLAGS.runtime_report = True
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertLen(commands, 5)
    self.assertEqual(
        commands[4][0],
        (
            'time /opt/deepvariant/bin/runtime_by_region_vis '
            '--input "/tmp/deepvariant_tmp_output/LOGDIR/'
            'make_examples_runtime_by_region/make_examples_runtime@64.tsv" '
            '--title "custom title" '
            '--output "/tmp/deepvariant_tmp_output/LOGDIR/'
            'make_examples_runtime_by_region_report.html"'
        ),
    )

  @flagsaver.flagsaver
  def test_without_runtime_report(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/deepvariant_tmp_output/LOGDIR'
    FLAGS.runtime_report = False
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertLen(commands, 4)  # No runtime report command.
    # Check that make_examples isn't asked to output runtimes:
    self.assertNotIn('--runtime_by_region', commands[0][0])

  @flagsaver.flagsaver
  def test_with_vcf_stats_report(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.report_title = 'custom title'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/deepvariant_tmp_output/LOGDIR'
    FLAGS.vcf_stats_report = True
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertLen(commands, 4)
    self.assertEqual(
        commands[3][0],
        (
            'time /opt/deepvariant/bin/vcf_stats_report --input_vcf "your_vcf" '
            '--outfile_base "your_vcf" --title "custom title"'
        ),
    )
    self.assertEqual(
        commands[3][1],
        '/tmp/deepvariant_tmp_output/LOGDIR/vcf_stats_report.log',
    )

  @flagsaver.flagsaver
  def test_without_vcf_stats_report(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/deepvariant_tmp_output/LOGDIR'
    FLAGS.vcf_stats_report = False
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertLen(commands, 3)  # No vcf stats report command.

  @parameterized.parameters(
      ('chr1:20-30', '--regions "chr1:20-30"'),
      ('chr1:20-30 chr2:100-200', '--regions "chr1:20-30 chr2:100-200"'),
      ("'chr1:20-30 chr2:100-200'", "--regions 'chr1:20-30 chr2:100-200'"),
  )
  @flagsaver.flagsaver
  def test_make_examples_regions(self, regions, expected_args):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.regions = regions
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/make_examples --mode calling '
        '--ref "your_ref" --reads "your_bam" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '--channels "insert_size" '
        '--gvcf "/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
        '%s '
        '--task {}' % expected_args,
    )

  @flagsaver.flagsaver
  def test_make_examples_extra_args_invalid(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.make_examples_extra_args = 'keep_secondary_alignments'
    with self.assertRaisesRegex(ValueError, 'not enough values to unpack'):
      _ = run_deepvariant.create_all_commands_and_logfiles(
          '/tmp/deepvariant_tmp_output'
      )

  @flagsaver.flagsaver
  def test_early_exit_if_use_openvino(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.call_variants_extra_args = 'use_openvino=true'
    with self.assertRaisesRegex(
        RuntimeError,
        (
            'OpenVINO is not installed by default in DeepVariant Docker images.'
            ' Please rerun without use_openvino flag.'
        ),
    ):
      _ = run_deepvariant.create_all_commands_and_logfiles(
          '/tmp/deepvariant_tmp_output'
      )

  @parameterized.parameters(
      ('batch_size=1024', '--batch_size "1024"'),
      (
          (
              'batch_size=4096,config_string="gpu_options:'
              ' {per_process_gpu_memory_fraction: 0.5}"'
          ),
          (
              '--batch_size "4096" --config_string "gpu_options:'
              ' {per_process_gpu_memory_fraction: 0.5}"'
          ),
      ),
  )
  @flagsaver.flagsaver
  def test_call_variants_extra_args(
      self, call_variants_extra_args, expected_args
  ):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.call_variants_extra_args = call_variants_extra_args
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/call_variants '
        '--outfile '
        '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
        '--examples "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '--checkpoint "/opt/models/wgs" '
        '%s' % expected_args,
    )

  @parameterized.parameters(
      ('qual_filter=3.0', '--qual_filter "3.0"'),
  )
  @flagsaver.flagsaver
  def test_postprocess_variants_extra_args(
      self, postprocess_variants_extra_args, expected_args
  ):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.postprocess_variants_extra_args = postprocess_variants_extra_args
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/postprocess_variants '
        '--ref "your_ref" '
        '--infile '
        '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
        '--outfile "your_vcf" '
        '--gvcf_outfile "your_gvcf" '
        '--nonvariant_site_tfrecord_path '
        '"/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
        '%s' % expected_args,
    )

  @flagsaver.flagsaver
  def test_postprocess_variants_haploid_flags(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.sample_name = 'your_sample'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.haploid_contigs = 'chrY'
    FLAGS.par_regions_bed = 'foo.bed'
    commands = run_deepvariant.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/postprocess_variants '
        '--ref "your_ref" '
        '--infile '
        '"/tmp/deepvariant_tmp_output/call_variants_output.tfrecord.gz" '
        '--outfile "your_vcf" '
        '--gvcf_outfile "your_gvcf" '
        '--haploid_contigs "chrY" '
        '--nonvariant_site_tfrecord_path '
        '"/tmp/deepvariant_tmp_output/gvcf.tfrecord@64.gz" '
        '--par_regions_bed "foo.bed" '
        '--sample_name "your_sample"',
    )


if __name__ == '__main__':
  absltest.main()
