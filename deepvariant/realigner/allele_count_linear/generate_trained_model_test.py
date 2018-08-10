# Copyright 2018 Google Inc.
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
"""Tests for :generate_trained_model."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tempfile
from absl import flags
from absl.testing import absltest
from absl.testing import parameterized
import tensorflow as tf

from deepvariant import testdata
from deepvariant.protos import realigner_pb2
from deepvariant.realigner.allele_count_linear import generate_trained_model
from google.protobuf import text_format

FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


class GenerateTrainedModelEnd2EndTest(parameterized.TestCase):

  def test_runs_reproducible(self):
    """Makes sure that generate_trained_model returns the expected proto."""
    output_model_proto = tempfile.NamedTemporaryFile(
        mode='w', dir=FLAGS.test_tmpdir, delete=False)
    fname = output_model_proto.name

    generate_trained_model.generate_trained_model_runner(
        truth_variants=testdata.TRUTH_VARIANTS_VCF,
        reads=testdata.CHR20_BAM,
        ref=testdata.CHR20_FASTA,
        output_model_proto=fname,
        output_model_pckl=None,
        exclude_contig=None,
        from_contig='chr20',
        random_seed=42,
        indel_weight=1)

    with tf.gfile.GFile(fname) as f:
      window_selector_model = text_format.Parse(
          f.read(), realigner_pb2.WindowSelectorModel())
      # Hardcoded value obtained from a "golden" run.
      expected = realigner_pb2.WindowSelectorModel(
          model_type=realigner_pb2.WindowSelectorModel.ALLELE_COUNT_LINEAR,
          allele_count_linear_model=(
              realigner_pb2.WindowSelectorModel.AlleleCountLinearModel(
                  bias=0.0259438883513,
                  coeff_soft_clip=0.00196795910597,
                  coeff_substitution=-0.545202672482,
                  coeff_insertion=0.267441004515,
                  coeff_deletion=0.211069211364,
                  coeff_reference=0.191676750779,
                  decision_boundary=3.0)))
      self.assertEqual(window_selector_model, expected)


if __name__ == '__main__':
  absltest.main()
