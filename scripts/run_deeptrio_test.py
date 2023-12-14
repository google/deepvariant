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
"""Tests for deepvariant .run_deeptrio."""

import io
from unittest import mock

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized
from deepvariant.opensource_only.scripts import run_deeptrio

FLAGS = flags.FLAGS


# pylint: disable=line-too-long
class RunDeeptrioTest(parameterized.TestCase):

  def setUp(self):
    super().setUp()
    self._saved_flags = flagsaver.save_flag_values()

  def tearDown(self):
    super().tearDown()
    flagsaver.restore_flag_values(self._saved_flags)

  def _create_all_commands_and_check_stdout(self, expected_stdout=None):
    with mock.patch('sys.stdout', new_callable=io.StringIO) as mock_stdout:
      commands, postprocess_cmds, report_commands = (
          run_deeptrio.create_all_commands('/tmp/deeptrio_tmp_output')
      )
      # Confirm that these basic commands don't have extra messages printed out
      # to stdout.
      if expected_stdout is None:
        self.assertEmpty(mock_stdout.getvalue())
      else:
        self.assertEqual(mock_stdout.getvalue(), expected_stdout)
    return commands, postprocess_cmds, report_commands

  @parameterized.parameters('WGS', 'WES', 'PACBIO')
  @flagsaver.flagsaver
  def test_call_variants_postprocess_variants_commands(self, model_type):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    commands, postprocess_cmds, report_commands = (
        self._create_all_commands_and_check_stdout()
    )
    # Because PACBIO model will always have use_candidate_partition on,
    # so there will be one extra make_examples command.
    call_variants_commands_start_index = 1
    if model_type == 'PACBIO':
      call_variants_commands_start_index = 2
    self.assertEqual(
        commands[call_variants_commands_start_index],
        'time /opt/deepvariant/bin/call_variants --outfile'
        ' "/tmp/deeptrio_tmp_output/call_variants_output_child.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deeptrio_tmp_output/make_examples_child.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/deeptrio/{}/child"'.format(
            model_type.lower()
        ),
    )
    self.assertEqual(
        commands[call_variants_commands_start_index + 1],
        'time /opt/deepvariant/bin/call_variants --outfile'
        ' "/tmp/deeptrio_tmp_output/call_variants_output_parent1.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deeptrio_tmp_output/make_examples_parent1.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/deeptrio/{}/parent"'.format(
            model_type.lower()
        ),
    )
    self.assertEqual(
        commands[call_variants_commands_start_index + 2],
        'time /opt/deepvariant/bin/call_variants --outfile'
        ' "/tmp/deeptrio_tmp_output/call_variants_output_parent2.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deeptrio_tmp_output/make_examples_parent2.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/deeptrio/{}/parent"'.format(
            model_type.lower()
        ),
    )
    self.assertEqual(
        postprocess_cmds[0],
        (
            'time /opt/deepvariant/bin/postprocess_variants '
            '--ref "your_ref" '
            '--infile '
            '"/tmp/deeptrio_tmp_output/call_variants_output_child.tfrecord.gz" '
            '--outfile "your_vcf_child" '
            '--cpus 0 '
            '--nonvariant_site_tfrecord_path '
            '"/tmp/deeptrio_tmp_output/gvcf_child.tfrecord@64.gz" '
            '--gvcf_outfile "your_gvcf_child"'
        ),
    )
    self.assertEqual(
        postprocess_cmds[1],
        (
            'time /opt/deepvariant/bin/postprocess_variants --ref "your_ref"'
            ' --infile'
            ' "/tmp/deeptrio_tmp_output/call_variants_output_parent1.tfrecord.gz"'
            ' --outfile "your_vcf_parent1" --cpus 0'
            ' --nonvariant_site_tfrecord_path'
            ' "/tmp/deeptrio_tmp_output/gvcf_parent1.tfrecord@64.gz"'
            ' --gvcf_outfile "your_gvcf_parent1"'
        ),
    )
    self.assertEqual(
        postprocess_cmds[2],
        (
            'time /opt/deepvariant/bin/postprocess_variants --ref "your_ref"'
            ' --infile'
            ' "/tmp/deeptrio_tmp_output/call_variants_output_parent2.tfrecord.gz"'
            ' --outfile "your_vcf_parent2" --cpus 0'
            ' --nonvariant_site_tfrecord_path'
            ' "/tmp/deeptrio_tmp_output/gvcf_parent2.tfrecord@64.gz"'
            ' --gvcf_outfile "your_gvcf_parent2"'
        ),
    )
    # Because PACBIO model will always have use_candidate_partition on,
    # so there will be one extra make_examples command.
    if model_type == 'PACBIO':
      self.assertLen(commands, 5)
    else:
      self.assertLen(commands, 4)
    self.assertLen(postprocess_cmds, 3)
    self.assertLen(report_commands, 3)

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters(
      (model_type, use_slim_model)
      for model_type in ['WGS', 'WES', 'PACBIO']
      for use_slim_model in [False, True]
  )
  # pylint: enable=g-complex-comprehension
  @flagsaver.flagsaver
  def test_duo_call_variants_postprocess_variants_commands(
      self, model_type, use_slim_model
  ):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.use_slim_model = use_slim_model
    FLAGS.num_shards = 64
    commands, postprocess_cmds, report_commands = (
        self._create_all_commands_and_check_stdout()
    )

    call_variants_bin = (
        'call_variants_slim' if use_slim_model else 'call_variants'
    )
    # Because PACBIO model will always have use_candidate_partition on,
    # so there will be one extra make_examples command.
    call_variants_commands_start_index = 1
    if model_type == 'PACBIO':
      call_variants_commands_start_index = 2
    self.assertEqual(
        commands[call_variants_commands_start_index],
        f'time /opt/deepvariant/bin/{call_variants_bin} --outfile'
        ' "/tmp/deeptrio_tmp_output/call_variants_output_child.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deeptrio_tmp_output/make_examples_child.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/deeptrio/{}/child"'.format(
            model_type.lower()
        ),
    )
    self.assertEqual(
        commands[call_variants_commands_start_index + 1],
        f'time /opt/deepvariant/bin/{call_variants_bin} --outfile'
        ' "/tmp/deeptrio_tmp_output/call_variants_output_parent1.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deeptrio_tmp_output/make_examples_parent1.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/deeptrio/{}/parent"'.format(
            model_type.lower()
        ),
    )
    self.assertEqual(
        postprocess_cmds[0],
        (
            'time /opt/deepvariant/bin/postprocess_variants '
            '--ref "your_ref" '
            '--infile '
            '"/tmp/deeptrio_tmp_output/call_variants_output_child.tfrecord.gz" '
            '--outfile "your_vcf_child" '
            '--cpus 0 '
            '--nonvariant_site_tfrecord_path '
            '"/tmp/deeptrio_tmp_output/gvcf_child.tfrecord@64.gz" '
            '--gvcf_outfile "your_gvcf_child"'
        ),
    )
    self.assertEqual(
        postprocess_cmds[1],
        (
            'time /opt/deepvariant/bin/postprocess_variants --ref "your_ref"'
            ' --infile'
            ' "/tmp/deeptrio_tmp_output/call_variants_output_parent1.tfrecord.gz"'
            ' --outfile "your_vcf_parent1" --cpus 0'
            ' --nonvariant_site_tfrecord_path'
            ' "/tmp/deeptrio_tmp_output/gvcf_parent1.tfrecord@64.gz"'
            ' --gvcf_outfile "your_gvcf_parent1"'
        ),
    )
    # pylint: disable=g-generic-assert
    # Because PACBIO model will always have use_candidate_partition on,
    # so there will be one extra make_examples command.
    if model_type == 'PACBIO':
      self.assertLen(commands, 4)
    else:
      self.assertLen(commands, 3)
    self.assertLen(postprocess_cmds, 2)
    self.assertLen(report_commands, 2)

  @parameterized.parameters(
      (
          'WGS',
          False,
          '--channels "insert_size" '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--pileup_image_height_child "60" '
          + '--pileup_image_height_parent "40" ',
      ),
      (
          'WES',
          False,
          '--channels "insert_size" '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--pileup_image_height_child "100" '
          + '--pileup_image_height_parent "100" ',
      ),
      (
          # For pacbio candidate paritionning is turned on by default.
          # Currently, there is no way to disable it.
          'PACBIO',
          True,
          '--add_hp_channel --alt_aligned_pileup "diff_channels" '
          + '--candidate_positions '
          + '"/tmp/deeptrio_tmp_output/candidate_positions@64" '
          + '--discard_non_dna_regions '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--max_reads_for_dynamic_bases_per_region "200" '
          + '--min_mapping_quality "1" '
          + '--parse_sam_aux_fields --partition_size "25000" --phase_reads '
          + '--pileup_image_height_child "60" '
          + '--pileup_image_height_parent "40" --pileup_image_width "199" '
          + '--norealign_reads --sort_by_haplotypes '
          + '--track_ref_reads '
          + '--vsc_min_fraction_indels "0.12" ',
      ),
      (
          'WGS',
          True,
          '--candidate_positions '
          + '"/tmp/deeptrio_tmp_output/candidate_positions@64" '
          + '--channels "insert_size" '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--pileup_image_height_child "60" --pileup_image_height_parent'
          ' "40" ',
      ),
      (
          'WES',
          True,
          '--candidate_positions '
          + '"/tmp/deeptrio_tmp_output/candidate_positions@64" '
          + '--channels "insert_size" '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--pileup_image_height_child "100" --pileup_image_height_parent'
          ' "100" ',
      ),
  )
  @flagsaver.flagsaver
  def test_make_examples_commands_with_types(
      self, model_type, use_candidate_partition, extra_args_plus_gvcf
  ):
    FLAGS.model_type = model_type
    FLAGS.use_candidate_partition = use_candidate_partition
    FLAGS.ref = 'your_ref'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    if model_type == 'PACBIO':
      use_candidate_partition = True
    make_examples_command_index = 1 if use_candidate_partition else 0
    commands, _, _ = self._create_all_commands_and_check_stdout()
    self.assertEqual(
        commands[make_examples_command_index],
        'time seq 0 63 '
        '| parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/deeptrio/make_examples '
        '--mode calling '
        '--ref "your_ref" '
        '--reads_parent1 "your_bam_parent1" '
        '--reads_parent2 "your_bam_parent2" '
        '--reads "your_bam_child" '
        '--examples "/tmp/deeptrio_tmp_output/make_examples.tfrecord@64.gz" '
        '--sample_name "your_sample_child" '
        '--sample_name_parent1 "your_sample_parent1" '
        '--sample_name_parent2 "your_sample_parent2" '
        '%s'
        '--task {}' % extra_args_plus_gvcf,
    )

  @parameterized.parameters((
      'PACBIO',
      # make_examples command with candidat_sweep mode
      'time seq 0 63 '
      + '| parallel -q --halt 2 --line-buffer '
      + '/opt/deepvariant/bin/deeptrio/make_examples '
      + '--mode candidate_sweep '
      + '--ref "your_ref" '
      + '--reads_parent1 "your_bam_parent1" '
      + '--reads_parent2 "your_bam_parent2" '
      + '--reads "your_bam_child" '
      + '--examples "/tmp/deeptrio_tmp_output/make_examples.tfrecord@64.gz" '
      + '--sample_name "your_sample_child" '
      + '--sample_name_parent1 "your_sample_parent1" '
      + '--sample_name_parent2 "your_sample_parent2" '
      + '--add_hp_channel --alt_aligned_pileup "diff_channels" '
      + '--candidate_positions '
      + '"/tmp/deeptrio_tmp_output/candidate_positions@64" '
      + '--discard_non_dna_regions '
      + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
      + '--max_reads_for_dynamic_bases_per_region "200" '
      + '--min_mapping_quality "1" '
      + '--parse_sam_aux_fields --partition_size "10000" --phase_reads '
      + '--pileup_image_height_child'
      ' "60" --pileup_image_height_parent "40" --pileup_image_width'
      ' "199" '
      + '--norealign_reads --sort_by_haplotypes '
      + '--track_ref_reads '
      + '--vsc_min_fraction_indels "0.12" '
      + '--task {}',
      # # make_examples command with call mode
      'time seq 0 63 '
      + '| parallel -q --halt 2 --line-buffer '
      + '/opt/deepvariant/bin/deeptrio/make_examples '
      + '--mode calling '
      + '--ref "your_ref" '
      + '--reads_parent1 "your_bam_parent1" '
      + '--reads_parent2 "your_bam_parent2" '
      + '--reads "your_bam_child" '
      + '--examples "/tmp/deeptrio_tmp_output/make_examples.tfrecord@64.gz" '
      + '--sample_name "your_sample_child" '
      + '--sample_name_parent1 "your_sample_parent1" '
      + '--sample_name_parent2 "your_sample_parent2" '
      + '--add_hp_channel --alt_aligned_pileup "diff_channels" '
      + '--candidate_positions '
      + '"/tmp/deeptrio_tmp_output/candidate_positions@64" '
      + '--discard_non_dna_regions '
      + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
      + '--max_reads_for_dynamic_bases_per_region "200" '
      + '--min_mapping_quality "1" '
      + '--parse_sam_aux_fields --partition_size "25000" --phase_reads '
      + '--pileup_image_height_child "60" '
      + '--pileup_image_height_parent "40" --pileup_image_width "199" '
      + '--norealign_reads --sort_by_haplotypes '
      + '--track_ref_reads '
      + '--vsc_min_fraction_indels "0.12" '
      + '--task {}',
  ))
  @flagsaver.flagsaver
  def test_make_examples_commands_with_candidate_partition(
      self, model_type, extra_args_1, extra_args_2
  ):
    FLAGS.model_type = model_type
    FLAGS.use_candidate_partition = True
    FLAGS.ref = 'your_ref'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    commands, _, _ = self._create_all_commands_and_check_stdout()
    self.assertEqual(commands[0], extra_args_1)
    self.assertEqual(commands[1], extra_args_2)

  @parameterized.parameters(
      (
          'WGS',
          '--channels "insert_size" '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--pileup_image_height_child "60" '
          + '--pileup_image_height_parent "40" ',
      ),
      (
          'WES',
          '--channels "insert_size" '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--pileup_image_height_child "100" '
          + '--pileup_image_height_parent "100" ',
      ),
      (
          'PACBIO',
          '--add_hp_channel --alt_aligned_pileup "diff_channels" '
          + '--candidate_positions '
          + '"/tmp/deeptrio_tmp_output/candidate_positions@64" '
          + '--discard_non_dna_regions '
          + '--gvcf "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz" '
          + '--max_reads_for_dynamic_bases_per_region "200" '
          + '--min_mapping_quality "1" '
          + '--parse_sam_aux_fields --partition_size "25000" --phase_reads '
          + '--pileup_image_height_child "60" '
          + '--pileup_image_height_parent "40" --pileup_image_width "199" '
          + '--norealign_reads --sort_by_haplotypes '
          + '--track_ref_reads '
          + '--vsc_min_fraction_indels "0.12" ',
      ),
  )
  @flagsaver.flagsaver
  def test_duo_make_examples_commands_with_types(
      self, model_type, extra_args_plus_gvcf
  ):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.num_shards = 64
    commands, _, _ = self._create_all_commands_and_check_stdout()
    use_candidate_partition = False
    if model_type == 'PACBIO':
      use_candidate_partition = True
    make_examples_command_index = 1 if use_candidate_partition else 0
    self.assertEqual(
        commands[make_examples_command_index],
        'time seq 0 63 '
        '| parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/deeptrio/make_examples '
        '--mode calling '
        '--ref "your_ref" '
        '--reads_parent1 "your_bam_parent1" '
        '--reads "your_bam_child" '
        '--examples "/tmp/deeptrio_tmp_output/make_examples.tfrecord@64.gz" '
        '--sample_name "your_sample_child" '
        '--sample_name_parent1 "your_sample_parent1" '
        '%s'
        '--task {}' % extra_args_plus_gvcf,
    )

  @parameterized.parameters(
      (
          None,
          (
              '--add_hp_channel --alt_aligned_pileup "diff_channels"'
              ' --candidate_positions'
              ' "/tmp/deeptrio_tmp_output/candidate_positions@64"'
              ' --discard_non_dna_regions --gvcf'
              ' "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz"'
              ' --max_reads_for_dynamic_bases_per_region "200"'
              ' --min_mapping_quality "1" --parse_sam_aux_fields'
              ' --partition_size "10000" --phase_reads'
              ' --pileup_image_height_child "60" --pileup_image_height_parent'
              ' "40" --pileup_image_width "199" --norealign_reads'
              ' --sort_by_haplotypes --track_ref_reads'
              ' --vsc_min_fraction_indels "0.12" '
          ),
          None,
      ),
      (
          'alt_aligned_pileup="rows",vsc_min_fraction_indels=0.03',
          (
              '--add_hp_channel --alt_aligned_pileup "rows"'
              ' --candidate_positions'
              ' "/tmp/deeptrio_tmp_output/candidate_positions@64"'
              ' --discard_non_dna_regions --gvcf'
              ' "/tmp/deeptrio_tmp_output/gvcf.tfrecord@64.gz"'
              ' --max_reads_for_dynamic_bases_per_region "200"'
              ' --min_mapping_quality "1" --parse_sam_aux_fields'
              ' --partition_size "10000" --phase_reads'
              ' --pileup_image_height_child "60" --pileup_image_height_parent'
              ' "40" --pileup_image_width "199" --norealign_reads'
              ' --sort_by_haplotypes --track_ref_reads'
              ' --vsc_min_fraction_indels "0.03" '
          ),
          # Because PacBio uses candidate_sweep, make_examples got run twice.
          (
              '\nWarning: --alt_aligned_pileup is previously set to'
              ' diff_channels, now to "rows".\n\nWarning:'
              ' --vsc_min_fraction_indels is previously set to 0.12, now to'
              ' 0.03.\n'
              '\nWarning: --alt_aligned_pileup is previously set to'
              ' diff_channels, now to "rows".\n\nWarning:'
              ' --vsc_min_fraction_indels is previously set to 0.12, now to'
              ' 0.03.\n'
          ),
      ),
  )
  @flagsaver.flagsaver
  def test_pacbio_args_overwrite(
      self, make_examples_extra_args, expected_args, expected_stdout
  ):
    """Confirms that adding extra flags can overwrite the default from mode."""
    FLAGS.model_type = 'PACBIO'
    FLAGS.ref = 'your_ref'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    FLAGS.regions = None
    FLAGS.make_examples_extra_args = make_examples_extra_args
    commands, _, _ = self._create_all_commands_and_check_stdout(expected_stdout)
    self.assertEqual(
        commands[0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/deeptrio/make_examples --mode candidate_sweep '
        '--ref "your_ref" --reads_parent1 "your_bam_parent1" '
        '--reads_parent2 "your_bam_parent2" '
        '--reads "your_bam_child" '
        '--examples "/tmp/deeptrio_tmp_output/make_examples.tfrecord@64.gz" '
        '--sample_name "your_sample_child" '
        '--sample_name_parent1 "your_sample_parent1" '
        '--sample_name_parent2 "your_sample_parent2" '
        '%s'
        '--task {}' % expected_args,
    )

  @parameterized.parameters(
      (
          'chr1:20-30',
          '--channels "insert_size" --pileup_image_height_child "60" '
          + '--pileup_image_height_parent "40" '
          + '--regions "chr1:20-30"',
      ),
      (
          'chr1:20-30 chr2:100-200',
          '--channels "insert_size" '
          + '--pileup_image_height_child "60" --pileup_image_height_parent'
          ' "40" '
          + '--regions "chr1:20-30 chr2:100-200"',
      ),
      (
          "'chr1:20-30 chr2:100-200'",
          '--channels "insert_size" '
          + '--pileup_image_height_child "60" --pileup_image_height_parent'
          ' "40" '
          + "--regions 'chr1:20-30 chr2:100-200'",
      ),
  )
  def test_make_examples_regions(self, regions, expected_args):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.num_shards = 64
    FLAGS.regions = regions
    commands, _, _ = self._create_all_commands_and_check_stdout()

    self.assertEqual(
        commands[0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer '
        '/opt/deepvariant/bin/deeptrio/make_examples --mode calling '
        '--ref "your_ref" --reads_parent1 "your_bam_parent1" '
        '--reads_parent2 "your_bam_parent2" '
        '--reads "your_bam_child" '
        '--examples "/tmp/deeptrio_tmp_output/make_examples.tfrecord@64.gz" '
        '--sample_name "your_sample_child" '
        '--sample_name_parent1 "your_sample_parent1" '
        '--sample_name_parent2 "your_sample_parent2" '
        '%s '
        '--task {}' % expected_args,
    )

  @flagsaver.flagsaver
  def test_make_examples_extra_args_invalid(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    FLAGS.make_examples_extra_args = 'keep_secondary_alignments'
    with self.assertRaisesRegex(ValueError, 'not enough values to unpack'):
      _, _, _ = run_deeptrio.create_all_commands('/tmp/deeptrio_tmp_output')

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
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    FLAGS.call_variants_extra_args = call_variants_extra_args
    commands, _, _ = self._create_all_commands_and_check_stdout()

    self.assertEqual(
        commands[1],
        'time /opt/deepvariant/bin/call_variants --outfile'
        ' "/tmp/deeptrio_tmp_output/call_variants_output_child.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deeptrio_tmp_output/make_examples_child.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/deeptrio/wgs/child" %s' % expected_args,
    )

  @parameterized.parameters(
      ('qual_filter=3.0', '--qual_filter "3.0"'),
  )
  @flagsaver.flagsaver
  def test_postprocess_variants_child_extra_args(
      self, postprocess_variants_child_extra_args, expected_args
  ):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    FLAGS.postprocess_variants_child_extra_args = (
        postprocess_variants_child_extra_args
    )
    _, commands_post_process, _ = self._create_all_commands_and_check_stdout()

    self.assertEqual(
        commands_post_process[0],
        'time /opt/deepvariant/bin/postprocess_variants '
        '--ref "your_ref" '
        '--infile '
        '"/tmp/deeptrio_tmp_output/call_variants_output_child.tfrecord.gz" '
        '--outfile "your_vcf_child" '
        '--cpus 0 '
        '--nonvariant_site_tfrecord_path '
        '"/tmp/deeptrio_tmp_output/gvcf_child.tfrecord@64.gz" '
        '--gvcf_outfile "your_gvcf_child" '
        '%s' % expected_args,
    )
    self.assertEqual(
        commands_post_process[1],
        'time /opt/deepvariant/bin/postprocess_variants '
        '--ref "your_ref" '
        '--infile '
        '"/tmp/deeptrio_tmp_output/call_variants_output_parent1.tfrecord.gz" '
        '--outfile "your_vcf_parent1" '
        '--cpus 0 '
        '--nonvariant_site_tfrecord_path '
        '"/tmp/deeptrio_tmp_output/gvcf_parent1.tfrecord@64.gz" '
        '--gvcf_outfile "your_gvcf_parent1"',
    )

  def test_all_report_commands(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.sample_name_child = 'your_sample_child'
    FLAGS.sample_name_parent1 = 'your_sample_parent1'
    FLAGS.sample_name_parent2 = 'your_sample_parent2'
    FLAGS.reads_child = 'your_bam_child'
    FLAGS.reads_parent1 = 'your_bam_parent1'
    FLAGS.reads_parent2 = 'your_bam_parent2'
    FLAGS.output_vcf_child = 'your_vcf_child'
    FLAGS.output_vcf_parent1 = 'your_vcf_parent1'
    FLAGS.output_vcf_parent2 = 'your_vcf_parent2'
    FLAGS.output_gvcf_child = 'your_gvcf_child'
    FLAGS.output_gvcf_parent1 = 'your_gvcf_parent1'
    FLAGS.output_gvcf_parent2 = 'your_gvcf_parent2'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/deeptrio_tmp_output/LOGDIR'
    FLAGS.runtime_report = True
    FLAGS.vcf_stats_report = True

    _, _, report_commands = self._create_all_commands_and_check_stdout()

    self.assertEqual(
        report_commands[0],
        (
            'time /opt/deepvariant/bin/vcf_stats_report --input_vcf'
            ' "your_vcf_child" --outfile_base "your_vcf_child"'
        ),
    )
    self.assertEqual(
        report_commands[1],
        (
            'time /opt/deepvariant/bin/vcf_stats_report --input_vcf'
            ' "your_vcf_parent1" --outfile_base "your_vcf_parent1"'
        ),
    )
    self.assertEqual(
        report_commands[2],
        (
            'time /opt/deepvariant/bin/vcf_stats_report --input_vcf'
            ' "your_vcf_parent2" --outfile_base "your_vcf_parent2"'
        ),
    )
    self.assertEqual(
        report_commands[3],
        (
            'time /opt/deepvariant/bin/runtime_by_region_vis --input'
            ' "/tmp/deeptrio_tmp_output/LOGDIR/make_examples_runtime_by_region/make_examples_runtime@64.tsv"'
            ' --title "DeepTrio" --output'
            ' "/tmp/deeptrio_tmp_output/LOGDIR/make_examples_runtime_by_region_report.html"'
        ),
    )


if __name__ == '__main__':
  absltest.main()
