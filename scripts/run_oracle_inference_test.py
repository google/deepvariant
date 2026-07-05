# Copyright 2025 Google LLC.
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
from absl import flags
from absl.testing import absltest
from absl.testing import flagsaver
from absl.testing import parameterized

from deepvariant.opensource_only.scripts import run_oracle_inference

FLAGS = flags.FLAGS


class RunOracleInferenceTest(parameterized.TestCase):

  @flagsaver.flagsaver
  def test_basic_commands_wgs(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.num_shards = 64
    FLAGS.truth_variants = 'your_truth.vcf'
    FLAGS.confident_regions = 'your_conf.bed'
    FLAGS.labeler_algorithm = 'HAPLOTYPE_LABELER'
    commands = run_oracle_inference.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        first=commands[0][0],
        second=(
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples --mode training --ref'
            ' "your_ref" --reads "your_bam" --labeler_algorithm'
            ' "HAPLOTYPE_LABELER" --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz"'
            ' --channel_list "BASE_CHANNELS" --max_reads_per_partition 1500'
            ' --partition_size "1000" --confident_regions "your_conf.bed"'
            ' --truth_variants "your_truth.vcf" --task {}'
        ),
    )
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/labeled_examples_to_vcf '
        '--ref "your_ref" --examples'
        ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '--output_vcf "your_vcf"',
    )

  # pylint: disable=g-complex-comprehension
  @parameterized.parameters(model_type for model_type in ['PACBIO', 'ONT_R104'])
  # pylint: enable=g-complex-comprehension
  @flagsaver.flagsaver
  def test_basic_commands_long_reads(self, model_type):
    FLAGS.model_type = model_type
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.num_shards = 64
    FLAGS.truth_variants = 'your_truth.vcf'
    FLAGS.confident_regions = 'your_conf.bed'
    FLAGS.labeler_algorithm = 'HAPLOTYPE_LABELER'
    commands = run_oracle_inference.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        first=commands[0][0],
        second=(
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples --mode training --ref'
            ' "your_ref" --reads "your_bam" --labeler_algorithm'
            ' "HAPLOTYPE_LABELER" --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz"'
            ' --channel_list "BASE_CHANNELS" --max_reads_per_partition 1500'
            ' --partition_size "25000" --confident_regions "your_conf.bed"'
            ' --truth_variants "your_truth.vcf" --task {}'
        ),
    )
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/labeled_examples_to_vcf '
        '--ref "your_ref" --examples'
        ' "/tmp/deepvariant_tmp_output/make_examples.tfrecord@64.gz" '
        '--output_vcf "your_vcf"',
    )

  @flagsaver.flagsaver
  def test_basic_commands_pangenome(self):
    FLAGS.model_type = 'WGS'
    FLAGS.ref = 'your_ref'
    FLAGS.reads = 'your_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.num_shards = 64
    FLAGS.truth_variants = 'your_truth.vcf'
    FLAGS.confident_regions = 'your_conf.bed'
    FLAGS.labeler_algorithm = 'HAPLOTYPE_LABELER'
    FLAGS.pangenome = 'your_pangenome.gbz'
    FLAGS.gbz_shared_memory_size_gb = 20
    FLAGS.channel_list = 'read_base,base_quality,mapping_quality'
    commands = run_oracle_inference.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    # First command should be load_gbz_into_shared_memory.
    self.assertIn(
        '/opt/deepvariant/bin/load_gbz_into_shared_memory',
        commands[0][0],
    )
    self.assertIn('your_pangenome.gbz', commands[0][0])
    self.assertIn('--shared_memory_size_gb 20', commands[0][0])

    # Second command should be make_examples_pangenome_aware_dv.
    self.assertIn(
        '/opt/deepvariant/bin/make_examples_pangenome_aware_dv',
        commands[1][0],
    )
    self.assertIn('--mode training', commands[1][0])
    self.assertIn('--pangenome "your_pangenome.gbz"', commands[1][0])
    self.assertIn('--use_loaded_gbz_shared_memory', commands[1][0])
    self.assertIn('make_examples_pangenome.tfrecord@64.gz', commands[1][0])
    self.assertIn(
        '--channel_list "read_base,base_quality,mapping_quality"',
        commands[1][0],
    )

    # Third command should be labeled_examples_to_vcf.
    self.assertEqual(
        commands[2][0],
        'time /opt/deepvariant/bin/labeled_examples_to_vcf '
        '--ref "your_ref" --examples'
        ' "/tmp/deepvariant_tmp_output/make_examples_pangenome.tfrecord@64.gz" '
        '--output_vcf "your_vcf"',
    )

  @flagsaver.flagsaver
  def test_basic_commands_somatic(self):
    FLAGS.model_type = 'WES'
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'tumor_bam'
    FLAGS.reads_normal = 'normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.num_shards = 64
    FLAGS.truth_variants = 'your_truth.vcf'
    FLAGS.confident_regions = 'your_conf.bed'
    FLAGS.labeler_algorithm = 'HAPLOTYPE_LABELER'
    commands = run_oracle_inference.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        first=commands[0][0],
        second=(
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples_somatic --mode training --ref'
            ' "your_ref" --reads_tumor "tumor_bam" --reads_normal "normal_bam"'
            ' --labeler_algorithm "HAPLOTYPE_LABELER" --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples_somatic.tfrecord@64.gz"'
            ' --channel_list "BASE_CHANNELS" --max_reads_per_partition 1500'
            ' --partition_size "1000" --confident_regions "your_conf.bed"'
            ' --truth_variants "your_truth.vcf" --task {}'
        ),
    )
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/labeled_examples_to_vcf '
        '--ref "your_ref" --examples'
        ' "/tmp/deepvariant_tmp_output/make_examples_somatic.tfrecord@64.gz" '
        '--output_vcf "your_vcf"',
    )

  @flagsaver.flagsaver
  def test_basic_commands_somatic_with_sample_name(self):
    FLAGS.model_type = 'WES'
    FLAGS.ref = 'your_ref'
    FLAGS.reads_tumor = 'tumor_bam'
    FLAGS.reads_normal = 'normal_bam'
    FLAGS.output_vcf = 'your_vcf'
    FLAGS.num_shards = 64
    FLAGS.truth_variants = 'your_truth.vcf'
    FLAGS.confident_regions = 'your_conf.bed'
    FLAGS.labeler_algorithm = 'HAPLOTYPE_LABELER'
    FLAGS.sample_name_tumor = 'tumor_sample'
    FLAGS.sample_name_normal = 'normal_sample'
    commands = run_oracle_inference.create_all_commands_and_logfiles(
        '/tmp/deepvariant_tmp_output'
    )

    self.assertEqual(
        first=commands[0][0],
        second=(
            'time seq 0 63 | parallel -q --halt 2 --line-buffer'
            ' /opt/deepvariant/bin/make_examples_somatic --mode training --ref'
            ' "your_ref" --reads_tumor "tumor_bam" --reads_normal "normal_bam"'
            ' --labeler_algorithm "HAPLOTYPE_LABELER" --examples'
            ' "/tmp/deepvariant_tmp_output/make_examples_somatic.tfrecord@64.gz"'
            ' --channel_list "BASE_CHANNELS" --max_reads_per_partition 1500'
            ' --partition_size "1000" --confident_regions "your_conf.bed"'
            ' --sample_name_normal "normal_sample" --sample_name_tumor'
            ' "tumor_sample" --truth_variants "your_truth.vcf" --task {}'
        ),
    )
    self.assertEqual(
        commands[1][0],
        'time /opt/deepvariant/bin/labeled_examples_to_vcf '
        '--ref "your_ref" --examples'
        ' "/tmp/deepvariant_tmp_output/make_examples_somatic.tfrecord@64.gz" '
        '--output_vcf "your_vcf" --sample_name "tumor_sample"',
    )


if __name__ == '__main__':
  absltest.main()
