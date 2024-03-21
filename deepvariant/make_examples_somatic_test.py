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
"""Tests for make_examples_somatic."""

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


if __name__ == '__main__':
  absltest.main()
