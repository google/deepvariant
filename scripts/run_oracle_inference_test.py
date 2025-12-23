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


if __name__ == '__main__':
  absltest.main()
