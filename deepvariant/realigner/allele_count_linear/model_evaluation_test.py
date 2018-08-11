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
"""Tests for :model_evaluation."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import tempfile
from absl import flags
from absl.testing import absltest
from absl.testing import parameterized
import pandas as pd

from deepvariant import testdata
from deepvariant.realigner.allele_count_linear import model_evaluation

FLAGS = flags.FLAGS


def setUpModule():
  testdata.init()


class ModelEvaluationEnd2EndTest(parameterized.TestCase):

  def test_runs_reproducible(self):
    """Makes sure that model_evaluation returns the expected values."""
    output_report_csv = tempfile.NamedTemporaryFile(
        mode='w', dir=FLAGS.test_tmpdir, delete=False)
    fname = output_report_csv.name

    model_evaluation.model_evaluation_runner(
        truth_variants=testdata.TRUTH_VARIANTS_VCF,
        reads=testdata.CHR20_BAM,
        ref=testdata.CHR20_FASTA,
        input_model_pckl=testdata.WS_ALLELE_COUNT_LINEAR_MODEL_PCKL,
        eval_region='chr20:10000000-10010000',
        output_report_csv=fname)

    actual_report = pd.read_csv(fname).to_dict('records')
    # Hardcoded value obtained from a "golden" run.
    expected = [{
        'threshold': -5.0,
        'recall': 1.0,
        'precision': 0.0007035175879396985
    }, {
        'threshold': -4.5,
        'recall': 1.0,
        'precision': 0.000720016457519029
    }, {
        'threshold': -4.0,
        'recall': 1.0,
        'precision': 0.00078405017921146967
    }, {
        'threshold': -3.5,
        'recall': 1.0,
        'precision': 0.0011762728953117123
    }, {
        'threshold': -3.0,
        'recall': 1.0,
        'precision': 0.003116651825467498
    }, {
        'threshold': -2.5,
        'recall': 1.0,
        'precision': 0.011217948717948718
    }, {
        'threshold': -2.0,
        'recall': 1.0,
        'precision': 0.04929577464788732
    }, {
        'threshold': -1.5,
        'recall': 1.0,
        'precision': 0.08045977011494253
    }, {
        'threshold': -1.0,
        'recall': 1.0,
        'precision': 0.09210526315789473
    }, {
        'threshold': -0.5,
        'recall': 1.0,
        'precision': 0.11666666666666667
    }, {
        'threshold': 0.0,
        'recall': 1.0,
        'precision': 0.1794871794871795
    }, {
        'threshold': 0.5,
        'recall': 1.0,
        'precision': 0.2413793103448276
    }, {
        'threshold': 1.0,
        'recall': 1.0,
        'precision': 0.3684210526315789
    }, {
        'threshold': 1.5,
        'recall': 1.0,
        'precision': 0.3888888888888889
    }, {
        'threshold': 2.0,
        'recall': 1.0,
        'precision': 0.3888888888888889
    }, {
        'threshold': 2.5,
        'recall': 1.0,
        'precision': 0.3888888888888889
    }, {
        'threshold': 3.0,
        'recall': 1.0,
        'precision': 0.4666666666666667
    }, {
        'threshold': 3.5,
        'recall': 1.0,
        'precision': 0.4666666666666667
    }, {
        'threshold': 4.0,
        'recall': 1.0,
        'precision': 0.5384615384615384
    }, {
        'threshold': 4.5,
        'recall': 1.0,
        'precision': 0.5384615384615384
    }, {
        'threshold': 5.0,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 5.5,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 6.0,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 6.5,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 7.0,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 7.5,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 8.0,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 8.5,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 9.0,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }, {
        'threshold': 9.5,
        'recall': 1.0,
        'precision': 0.5833333333333334
    }]
    self.assertEqual(len(actual_report), len(expected))
    for a_record, e_record in zip(actual_report, expected):
      self.assertEqual(a_record['threshold'], e_record['threshold'])
      self.assertEqual(a_record['recall'], e_record['recall'])
      self.assertAlmostEqual(
          a_record['precision'], e_record['precision'], places=6)


if __name__ == '__main__':
  absltest.main()
