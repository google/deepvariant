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
"""Tests for deepvariant .run_deepsomatic."""


from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from deepvariant.scripts import run_deepsomatic

FLAGS = flags.FLAGS


class RunDeepSomaticTest(parameterized.TestCase):

  @parameterized.parameters(
      ('WGS', True),
      ('WGS', False),
  )
  @flagsaver.flagsaver
  def test_basic_commands(self, model_type, use_keras_model):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.use_keras_model = use_keras_model
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )

    extra_args_plus_gvcf = (
        '--gvcf "/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz"'
    )

    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_somatic --mode calling --ref'
        ' "your_ref" --reads_tumor "your_tumor_bam" --reads_normal'
        ' "your_normal_bam" --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --channels "insert_size"'
        ' %s'
        ' --vsc_max_fraction_indels_for_non_target_sample "0.5"'
        ' --vsc_max_fraction_snps_for_non_target_sample "0.5"'
        ' --vsc_min_fraction_indels "0.05" --vsc_min_fraction_snps "0.029"'
        ' --task {}' % (extra_args_plus_gvcf),
    )
    call_variants_bin = (
        'call_variants_keras' if use_keras_model else 'call_variants'
    )
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/{} --outfile'
        ' "/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/{}/model.ckpt"'.format(
            call_variants_bin, model_type.lower()
        ),
    )
    self.assertEqual(
        commands[2][0],
        (
            'time /opt/deepvariant/bin/postprocess_variants '
            '--ref "your_ref" '
            '--infile '
            '"/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz" '
            '--outfile "your_vcf" '
            '--process_somatic=true '
            '--gvcf_outfile "your_gvcf" '
            '--nonvariant_site_tfrecord_path '
            '"/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz"'
        ),
    )

  @parameterized.parameters(
      ('WGS', True, True),
      ('WGS', False, False),
  )
  @flagsaver.flagsaver
  def test_process_somatic(self, model_type, use_keras_model, process_somatic):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.use_keras_model = use_keras_model
    FLAGS.process_somatic = process_somatic
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )

    extra_args_plus_gvcf = (
        '--gvcf "/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz"'
    )

    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_somatic --mode calling --ref'
        ' "your_ref" --reads_tumor "your_tumor_bam" --reads_normal'
        ' "your_normal_bam" --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --channels "insert_size"'
        ' %s'
        ' --vsc_max_fraction_indels_for_non_target_sample "0.5"'
        ' --vsc_max_fraction_snps_for_non_target_sample "0.5"'
        ' --vsc_min_fraction_indels "0.05" --vsc_min_fraction_snps "0.029"'
        ' --task {}' % (extra_args_plus_gvcf),
    )
    call_variants_bin = (
        'call_variants_keras' if use_keras_model else 'call_variants'
    )
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/{} --outfile'
        ' "/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/{}/model.ckpt"'.format(
            call_variants_bin, model_type.lower()
        ),
    )
    if process_somatic:
      self.assertEqual(
          commands[2][0],
          (
              'time /opt/deepvariant/bin/postprocess_variants '
              '--ref "your_ref" '
              '--infile '
              '"/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz" '
              '--outfile "your_vcf" '
              '--process_somatic=true '
              '--gvcf_outfile "your_gvcf" '
              '--nonvariant_site_tfrecord_path '
              '"/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz"'
          ),
      )
    else:
      self.assertEqual(
          commands[2][0],
          (
              'time /opt/deepvariant/bin/postprocess_variants '
              '--ref "your_ref" '
              '--infile '
              '"/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz" '
              '--outfile "your_vcf" '
              '--noprocess_somatic '
              '--gvcf_outfile "your_gvcf" '
              '--nonvariant_site_tfrecord_path '
              '"/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz"'
          ),
      )

  @parameterized.parameters((None, None),
                            ('your_tumor_sample_name',
                             'your_normal_sample_name'),
                            ('just_tumor_sample_name', None))
  @flagsaver.flagsaver
  def test_sample_name_command(self, sample_name_tumor, sample_name_normal):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = None
    FLAGS.num_shards = 64
    FLAGS.regions = None
    FLAGS.sample_name_tumor = sample_name_tumor
    FLAGS.sample_name_normal = sample_name_normal
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )

    extra_sample_name_flag = ''
    if sample_name_normal is not None:
      extra_sample_name_flag += (
          f' --sample_name_normal "{sample_name_normal}"'
      )
    if sample_name_tumor is not None:
      extra_sample_name_flag += (
          f' --sample_name_tumor "{sample_name_tumor}"'
      )

    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_somatic --mode calling --ref'
        ' "your_ref" --reads_tumor "your_tumor_bam" --reads_normal'
        ' "your_normal_bam" --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --channels "insert_size"'
        '%s'
        ' --vsc_max_fraction_indels_for_non_target_sample "0.5"'
        ' --vsc_max_fraction_snps_for_non_target_sample "0.5"'
        ' --vsc_min_fraction_indels "0.05" --vsc_min_fraction_snps "0.029"'
        ' --task {}' % extra_sample_name_flag,
    )
    self.assertEqual(
        commands[1][0],
        (
            'time /opt/deepvariant/bin/call_variants --outfile'
            ' "/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz"'
            ' --examples'
            ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
            ' --checkpoint "/opt/models/wgs/model.ckpt"'
        ),
    )
    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/postprocess_variants '
        '--ref "your_ref" '
        '--infile '
        '"/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz" '
        '--outfile "your_vcf" '
        '--process_somatic=true',
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
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.make_examples_extra_args = make_examples_extra_args
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_somatic --mode calling --ref'
        ' "your_ref" --reads_tumor "your_tumor_bam" --reads_normal'
        ' "your_normal_bam" --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --channels "insert_size"'
        ' --gvcf'
        ' "/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz" %s'
        ' --vsc_max_fraction_indels_for_non_target_sample "0.5"'
        ' --vsc_max_fraction_snps_for_non_target_sample "0.5"'
        ' --vsc_min_fraction_indels "0.05" --vsc_min_fraction_snps "0.029"'
        ' --task {}' % expected_args,
    )

  @flagsaver.flagsaver
  def test_logging_dir(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/deepsomatic_tmp_output/LOGDIR'
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )

    self.assertEqual(
        commands[0][0],
        (
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples_somatic --mode calling --ref'
            ' "your_ref" --reads_tumor "your_tumor_bam" --reads_normal'
            ' "your_normal_bam" --examples'
            ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
            ' --channels "insert_size" --gvcf'
            ' "/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz"'
            ' --vsc_max_fraction_indels_for_non_target_sample "0.5"'
            ' --vsc_max_fraction_snps_for_non_target_sample "0.5"'
            ' --vsc_min_fraction_indels "0.05" --vsc_min_fraction_snps "0.029"'
            ' --task {}'
        ),
    )

  @parameterized.parameters(
      ('chr1:20-30', '--regions "chr1:20-30"'),
      ('chr1:20-30 chr2:100-200', '--regions "chr1:20-30 chr2:100-200"'),
      ("'chr1:20-30 chr2:100-200'", "--regions 'chr1:20-30 chr2:100-200'"),
  )
  @flagsaver.flagsaver
  def test_make_examples_regions(self, regions, expected_args):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.regions = regions
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )

    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_somatic --mode calling --ref'
        ' "your_ref" --reads_tumor "your_tumor_bam" --reads_normal'
        ' "your_normal_bam" --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --channels "insert_size"'
        ' --gvcf'
        ' "/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz" %s'
        ' --vsc_max_fraction_indels_for_non_target_sample "0.5"'
        ' --vsc_max_fraction_snps_for_non_target_sample "0.5"'
        ' --vsc_min_fraction_indels "0.05" --vsc_min_fraction_snps "0.029"'
        ' --task {}' % expected_args,
    )

  @flagsaver.flagsaver

  @flagsaver.flagsaver
  def test_make_examples_extra_args_invalid(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.make_examples_extra_args = 'keep_secondary_alignments'
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    with self.assertRaisesRegex(ValueError, 'not enough values to unpack'):
      _ = run_deepsomatic.create_all_commands_and_logfiles(
          '/tmp/deepsomatic_tmp_output', used_in_test=True
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
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.call_variants_extra_args = call_variants_extra_args
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )

    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/call_variants --outfile'
        ' "/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz"'
        ' --examples'
        ' "/tmp/deepsomatic_tmp_output/make_examples_somatic.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/wgs/model.ckpt" %s' % expected_args,
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
    FLAGS.reads_tumor = 'your_tumor_bam'
    FLAGS.reads_normal = 'your_normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.postprocess_variants_extra_args = postprocess_variants_extra_args
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    commands = run_deepsomatic.create_all_commands_and_logfiles(
        '/tmp/deepsomatic_tmp_output', used_in_test=True
    )

    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/postprocess_variants '
        '--ref "your_ref" '
        '--infile '
        '"/tmp/deepsomatic_tmp_output/call_variants_output.tfrecord.gz" '
        '--outfile "your_vcf" '
        '--process_somatic=true '
        '--gvcf_outfile "your_gvcf" '
        '--nonvariant_site_tfrecord_path '
        '"/tmp/deepsomatic_tmp_output/gvcf.tfrecord@64.gz" '
        '%s' % expected_args,
    )


if __name__ == '__main__':
  absltest.main()
