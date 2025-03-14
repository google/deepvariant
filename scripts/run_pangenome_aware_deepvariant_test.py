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
"""Tests for deepvariant .run_pangenome_aware_deepvariant."""

from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from deepvariant.opensource_only.scripts import run_pangenome_aware_deepvariant

FLAGS = flags.FLAGS


class RunPangenomeAwareDeepVariantTest(parameterized.TestCase):

  @parameterized.parameters(('WGS',))
  @flagsaver.flagsaver
  def test_basic_commands(self, model_type):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.disable_small_model = False
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )

    extra_args_plus_gvcf = (
        '--gvcf'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
    )

    # pyformat: disable
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_pangenome_aware_dv --mode calling'
        ' --ref "your_ref" --reads "your_bam" --pangenome "your_pangenome_bam"'
        ' --examples'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/wgs/model.ckpt"'
        ' --call_small_model_examples %s'
        ' --keep_legacy_allele_counter_behavior'
        ' --keep_only_window_spanning_haplotypes'
        ' --keep_supplementary_alignments'
        ' --min_mapping_quality "0" --normalize_reads'
        ' --sample_name_pangenome "hprc_v1.1"'
        ' --small_model_indel_gq_threshold "30"'
        ' --small_model_snp_gq_threshold "25"'
        ' --small_model_vaf_context_window_size "51"'
        ' --sort_by_haplotypes'
        ' --trained_small_model_path "/opt/smallmodels/wgs"'
        ' --trim_reads_for_pileup'
        ' --task {}'
        % (extra_args_plus_gvcf),
    )
    # pyformat: enable
    call_variants_bin = 'call_variants'
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/{} --outfile'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
        ' --examples'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/{}/model.ckpt"'.format(
            call_variants_bin, model_type.lower()
        ),
    )
    # pyformat: disable
    self.assertEqual(
        commands[2][0],
        (
            'time /opt/deepvariant/bin/postprocess_variants --ref "your_ref"'
            ' --infile'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
            ' --outfile "your_vcf"'
            ' --cpus "64"'
            ' --small_model_cvo_records'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv_call_variant_outputs.tfrecord@64.gz"'
            ' --gvcf_outfile "your_gvcf"'
            ' --nonvariant_site_tfrecord_path'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
        ),
    )
    # pyformat: enable

  @parameterized.parameters(
      # If I set pangenome_sample_name to None, I get:
      # `absl.flags._exceptions.IllegalFlagValueError: Unexpected None value for
      #  flag sample_name_pangenome`
      # which I couldn't figured why yet.
      # So, setting it to '' for now.
      (None, ''),
      (None, 'your_pangenome_sample_name'),
      ('your_sample_name', 'your_pangenome_sample_name'),
      ('your_sample_name', ''),
  )
  @flagsaver.flagsaver
  def test_sample_name_command(self, sample_name_reads, sample_name_pangenome):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = None
    FLAGS.num_shards = 64
    FLAGS.regions = None
    FLAGS.sample_name_reads = sample_name_reads
    FLAGS.sample_name_pangenome = sample_name_pangenome
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.disable_small_model = False
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )

    extra_sample_name_flag = ''
    if sample_name_pangenome is not None:
      extra_sample_name_flag += (
          f' --sample_name_pangenome "{sample_name_pangenome}"'
      )
    if sample_name_reads is not None:
      extra_sample_name_flag += f' --sample_name_reads "{sample_name_reads}"'

    # pyformat: disable
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_pangenome_aware_dv --mode calling'
        ' --ref "your_ref" --reads "your_bam" --pangenome "your_pangenome_bam"'
        ' --examples'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/wgs/model.ckpt"'
        ' --call_small_model_examples'
        ' --keep_legacy_allele_counter_behavior'
        ' --keep_only_window_spanning_haplotypes'
        ' --keep_supplementary_alignments'
        ' --min_mapping_quality "0" --normalize_reads'
        '%s'
        ' --small_model_indel_gq_threshold "30"'
        ' --small_model_snp_gq_threshold "25"'
        ' --small_model_vaf_context_window_size "51"'
        ' --sort_by_haplotypes'
        ' --trained_small_model_path "/opt/smallmodels/wgs"'
        ' --trim_reads_for_pileup --task {}' % extra_sample_name_flag,
    )
    self.assertEqual(
        commands[1][0],
        (
            'time /opt/deepvariant/bin/call_variants --outfile'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
            ' --examples'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
            ' --checkpoint "/opt/models/wgs/model.ckpt"'
        ),
    )
    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/postprocess_variants --ref "your_ref"'
        ' --infile'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
        ' --outfile "your_vcf"'
        ' --cpus "64"'
        ' --small_model_cvo_records'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv_call_variant_outputs.tfrecord@64.gz"',
    )
    # pyformat: enable

  @parameterized.parameters(
      # pyformat: disable
      (
          'keep_secondary_alignments=true',
          ' --keep_only_window_spanning_haplotypes'
          + ' --keep_secondary_alignments'
          + ' --keep_supplementary_alignments'
          + ' --min_mapping_quality "0" --normalize_reads'
          + ' --sample_name_pangenome "hprc_v1.1"'
          + ' --small_model_indel_gq_threshold "30"'
          + ' --small_model_snp_gq_threshold "25"'
          + ' --small_model_vaf_context_window_size "51"'
          + ' --sort_by_haplotypes'
          + ' --trained_small_model_path "/opt/smallmodels/wgs"'
          + ' --trim_reads_for_pileup',
      ),
      (
          'keep_secondary_alignments=false',
          ' --keep_only_window_spanning_haplotypes'
          + ' --nokeep_secondary_alignments'
          + ' --keep_supplementary_alignments'
          + ' --min_mapping_quality "0" --normalize_reads'
          + ' --sample_name_pangenome "hprc_v1.1"'
          + ' --small_model_indel_gq_threshold "30"'
          + ' --small_model_snp_gq_threshold "25"'
          + ' --small_model_vaf_context_window_size "51"'
          + ' --sort_by_haplotypes'
          + ' --trained_small_model_path "/opt/smallmodels/wgs"'
          + ' --trim_reads_for_pileup',
      ),
      (
          'keep_secondary_alignments=true,keep_supplementary_alignments=true',
          ' --keep_only_window_spanning_haplotypes'
          + ' --keep_secondary_alignments'
          + ' --keep_supplementary_alignments'
          + ' --min_mapping_quality "0" --normalize_reads'
          + ' --sample_name_pangenome "hprc_v1.1"'
          + ' --small_model_indel_gq_threshold "30"'
          + ' --small_model_snp_gq_threshold "25"'
          + ' --small_model_vaf_context_window_size "51"'
          + ' --sort_by_haplotypes'
          + ' --trained_small_model_path "/opt/smallmodels/wgs"'
          + ' --trim_reads_for_pileup',
      ),
      (
          'use_ref_for_cram=true,keep_secondary_alignments=true,'
          + 'keep_supplementary_alignments=false',
          ' --keep_only_window_spanning_haplotypes'
          + ' --keep_secondary_alignments'
          + ' --nokeep_supplementary_alignments'
          + ' --min_mapping_quality "0" --normalize_reads'
          + ' --sample_name_pangenome "hprc_v1.1"'
          + ' --small_model_indel_gq_threshold "30"'
          + ' --small_model_snp_gq_threshold "25"'
          + ' --small_model_vaf_context_window_size "51"'
          + ' --sort_by_haplotypes'
          + ' --trained_small_model_path "/opt/smallmodels/wgs"'
          + ' --trim_reads_for_pileup'
          + ' --use_ref_for_cram',
      ),
      # pyformat: enable
  )
  @flagsaver.flagsaver
  def test_make_examples_extra_args_boolean(
      self, make_examples_extra_args, full_expected_args
  ):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.make_examples_extra_args = make_examples_extra_args
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.disable_small_model = False
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )
    # pyformat: disable
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_pangenome_aware_dv --mode calling'
        ' --ref "your_ref" --reads "your_bam" --pangenome "your_pangenome_bam"'
        ' --examples'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/wgs/model.ckpt"'
        ' --call_small_model_examples'
        ' --gvcf'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
        ' --keep_legacy_allele_counter_behavior%s'
        ' --task {}' % full_expected_args,
    )
    # pyformat: enable

  @flagsaver.flagsaver
  def test_logging_dir(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.logging_dir = '/tmp/pangenome_aware_deepvariant_tmp_output/LOGDIR'
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.disable_small_model = False
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )

    # pyformat: disable
    self.assertEqual(
        commands[0][0],
        (
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples_pangenome_aware_dv --mode'
            ' calling --ref "your_ref" --reads "your_bam" --pangenome'
            ' "your_pangenome_bam" --examples'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
            ' --checkpoint "/opt/models/wgs/model.ckpt"'
            ' --call_small_model_examples'
            ' --gvcf'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
            ' --keep_legacy_allele_counter_behavior'
            ' --keep_only_window_spanning_haplotypes'
            ' --keep_supplementary_alignments'
            ' --min_mapping_quality "0" --normalize_reads'
            ' --sample_name_pangenome "hprc_v1.1"'
            ' --small_model_indel_gq_threshold "30"'
            ' --small_model_snp_gq_threshold "25"'
            ' --small_model_vaf_context_window_size "51"'
            ' --sort_by_haplotypes'
            ' --trained_small_model_path "/opt/smallmodels/wgs"'
            ' --trim_reads_for_pileup'
            ' --task {}'
        ),
    )
    # pyformat: enable

  @flagsaver.flagsaver
  def test_gbz_shared_memory_name(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    # Having --pangenome end with .gbz is important to test that the
    # load_gbz_into_shared_memory command is generated correctly.
    FLAGS.pangenome = 'your_pangenome_bam.gbz'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.num_shards = 64
    FLAGS.gbz_shared_memory_name = 'NEW_SHARED_MEMORY_NAME'
    FLAGS.disable_small_model = False
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )
    # Because --pangenome ends with .gbz, we expect a
    # load_gbz_into_shared_memory command to be generated.
    self.assertLen(commands, 4)
    self.assertEqual(
        commands[0][0],
        (
            'time /opt/deepvariant/bin/load_gbz_into_shared_memory'
            ' --pangenome_gbz "your_pangenome_bam.gbz"'
            ' --shared_memory_name "NEW_SHARED_MEMORY_NAME"'
            ' --shared_memory_size_gb 12'
            ' --num_shards "64"'
        ),
    )
    # pyformat: disable
    self.assertEqual(
        commands[1][0],
        (
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples_pangenome_aware_dv --mode'
            ' calling --ref "your_ref" --reads "your_bam" --pangenome'
            ' "your_pangenome_bam.gbz" --examples'
            ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
            ' --checkpoint "/opt/models/pangenome_aware_deepvariant/wgs"'
            ' --call_small_model_examples'
            ' --gbz_shared_memory_name "NEW_SHARED_MEMORY_NAME"'
            ' --keep_legacy_allele_counter_behavior'
            ' --keep_only_window_spanning_haplotypes'
            ' --keep_supplementary_alignments'
            ' --min_mapping_quality "0" --normalize_reads'
            ' --sample_name_pangenome "hprc_v1.1"'
            ' --small_model_indel_gq_threshold "30"'
            ' --small_model_snp_gq_threshold "25"'
            ' --small_model_vaf_context_window_size "51"'
            ' --sort_by_haplotypes'
            ' --trained_small_model_path "/opt/smallmodels/wgs"'
            ' --trim_reads_for_pileup'
            ' --use_loaded_gbz_shared_memory'
            ' --task {}'
        ),
    )
    # pyformat: enable

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
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.regions = regions
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.disable_small_model = False
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )

    # pyformat: disable
    self.assertEqual(
        commands[0][0],
        'time seq 0 63 | parallel -q --halt 2 --line-buffer'
        ' /opt/deepvariant/bin/make_examples_pangenome_aware_dv --mode calling'
        ' --ref "your_ref" --reads "your_bam" --pangenome "your_pangenome_bam"'
        ' --examples'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
        ' --checkpoint "/opt/models/wgs/model.ckpt"'
        ' --call_small_model_examples'
        ' --gvcf'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/gvcf.tfrecord@64.gz"'
        ' --keep_legacy_allele_counter_behavior'
        ' --keep_only_window_spanning_haplotypes'
        ' --keep_supplementary_alignments'
        ' --min_mapping_quality "0" --normalize_reads %s'
        ' --sample_name_pangenome "hprc_v1.1"'
        ' --small_model_indel_gq_threshold "30"'
        ' --small_model_snp_gq_threshold "25"'
        ' --small_model_vaf_context_window_size "51"'
        ' --sort_by_haplotypes'
        ' --trained_small_model_path "/opt/smallmodels/wgs"'
        ' --trim_reads_for_pileup'
        ' --task {}'
        % expected_args,
    )
    # pyformat: enable

  @flagsaver.flagsaver
  @flagsaver.flagsaver
  def test_make_examples_extra_args_invalid(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.make_examples_extra_args = 'keep_secondary_alignments'
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.disable_small_model = False
    with self.assertRaisesRegex(ValueError, 'not enough values to unpack'):
      _ = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
          '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
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
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.call_variants_extra_args = call_variants_extra_args
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.disable_small_model = False
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )

    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/call_variants --outfile'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
        ' --examples'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv.tfrecord@64.gz"'
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
    FLAGS.reads = 'your_bam'
    FLAGS.pangenome = 'your_pangenome_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.output_gvcf = 'your_gvcf'
    FLAGS.num_shards = 64
    FLAGS.postprocess_variants_extra_args = postprocess_variants_extra_args
    FLAGS.customized_model = '/opt/models/wgs/model.ckpt'
    FLAGS.disable_small_model = False
    commands = run_pangenome_aware_deepvariant.create_all_commands_and_logfiles(
        '/tmp/pangenome_aware_deepvariant_tmp_output', used_in_test=True
    )

    # pyformat: disable
    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/postprocess_variants --ref "your_ref"'
        ' --infile'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/call_variants_output.tfrecord.gz"'
        ' --outfile "your_vcf"'
        ' --cpus "64"'
        ' --small_model_cvo_records'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/make_examples_pangenome_aware_dv_call_variant_outputs.tfrecord@64.gz"'
        ' --gvcf_outfile "your_gvcf"'
        ' --nonvariant_site_tfrecord_path'
        ' "/tmp/pangenome_aware_deepvariant_tmp_output/gvcf.tfrecord@64.gz" %s'
        % expected_args,
    )
    # pyformat: enable


if __name__ == '__main__':
  absltest.main()
