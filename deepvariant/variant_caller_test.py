# Copyright 2017 Google Inc.
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
"""Tests for deepvariant .variant_labeler."""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function



from absl.testing import absltest
from absl.testing import parameterized
import mock
import numpy as np
import numpy.testing as npt

from deepvariant import test_utils
from deepvariant import variant_caller
from deepvariant.core import variantutils
from deepvariant.protos import deepvariant_pb2


def setUpModule():
  test_utils.init()


def _reference_model_options(p_error, max_gq):
  return deepvariant_pb2.VariantCallerOptions(
      sample_name='UNKNOWN',
      p_error=p_error,
      max_gq=max_gq,
      gq_resolution=1,
      ploidy=2)


class VariantCallerTests(parameterized.TestCase):

  def make_test_caller(self, p_error, max_gq):
    options = _reference_model_options(p_error, max_gq)
    return variant_caller.VariantCaller(options, use_cache_table=False)

  def fake_allele_counter(self, start_pos, counts):
    allele_counter = mock.Mock()
    allele_counter.summary_counts.return_value = [
        deepvariant_pb2.AlleleCountSummary(
            ref_supporting_read_count=n_ref,
            total_read_count=n_ref + n_alt,
            ref_base=ref,
            reference_name='chr1',
            position=start_pos + i)
        for i, (n_alt, n_ref, ref) in enumerate(counts)
    ]
    return allele_counter

  # R code to produce the testdata expectation table.
  # expected <- function(n_ref, n_alt, perr, max_gq = 100) {
  #   p_ref <- dbinom(n_alt, n_ref, perr)
  #   p_het <- dbinom(n_alt, n_ref, 0.5)
  #   p_alt <- dbinom(n_ref - n_alt, n_ref, perr)
  #   raw <- c(p_ref, p_het, p_alt)
  #   norm <- raw / sum(raw)
  #   gq = min(floor(-10 * log10(1 - norm[1])), max_gq)
  #   likelihoods = paste(sprintf("%.6f", log10(norm)), collapse=", ")
  #   likelihoods = paste("[", likelihoods, "]", sep="")
  #   result = paste(n_ref, n_alt, perr, 100, 1, likelihoods, gq, sep=", ")
  #   cat(paste("[", result, "],\n", sep=""))
  # }
  #
  # for (n in c(10, 20)) {
  #  for (k in seq(0, n)) {
  #     expected(n, k, 0.01)
  #   }
  # }
  #
  # for (perr in c(0.1, 0.01, 0.001, 0.0001)) {
  #   expected(10, 0, perr)
  #   expected(10, 1, perr)
  # }
  #
  # for (n_ref in c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 1000, 10000)) {
  #   expected(n_ref, 0, 0.01)
  # }
  @parameterized.parameters(
      # No coverage case.
      [0, 0, 0.01, 100, [-0.477121, -0.477121, -0.477121], 1],
      # Test systematically values of n and k.
      [10, 0, 0.01, 100, [-0.000469, -2.967121, -19.956821], 29],
      [10, 1, 0.01, 100, [-0.044109, -1.015126, -16.009190], 10],
      [10, 2, 0.01, 100, [-1.063830, -0.039211, -13.037641], 0],
      [10, 3, 0.01, 100, [-3.020668, -0.000414, -11.003209], 0],
      [10, 4, 0.01, 100, [-5.015893, -0.000004, -9.007163], 0],
      [10, 5, 0.01, 100, [-7.011524, -0.000000, -7.011524], 0],
      [10, 6, 0.01, 100, [-9.007163, -0.000004, -5.015893], 0],
      [10, 7, 0.01, 100, [-11.003209, -0.000414, -3.020668], 0],
      [10, 8, 0.01, 100, [-13.037641, -0.039211, -1.063830], 0],
      [10, 9, 0.01, 100, [-16.009190, -1.015126, -0.044109], 0],
      [10, 10, 0.01, 100, [-19.956821, -2.967121, -0.000469], 0],
      [20, 0, 0.01, 100, [-0.000001, -5.933304, -39.912704], 59],
      [20, 1, 0.01, 100, [-0.000050, -3.937719, -35.921484], 39],
      [20, 2, 0.01, 100, [-0.004935, -1.946968, -31.935098], 19],
      [20, 3, 0.01, 100, [-0.328657, -0.275056, -28.267550], 2],
      [20, 4, 0.01, 100, [-2.053097, -0.003860, -26.000720], 0],
      [20, 5, 0.01, 100, [-4.044911, -0.000039, -24.001263], 0],
      [20, 6, 0.01, 100, [-6.040508, -0.000000, -22.005589], 0],
      [20, 7, 0.01, 100, [-8.036143, -0.000000, -20.009954], 0],
      [20, 8, 0.01, 100, [-10.031778, -0.000000, -18.014319], 0],
      [20, 9, 0.01, 100, [-12.027413, -0.000000, -16.018683], 0],
      [20, 10, 0.01, 100, [-14.023048, -0.000000, -14.023048], 0],
      [20, 11, 0.01, 100, [-16.018683, -0.000000, -12.027413], 0],
      [20, 12, 0.01, 100, [-18.014319, -0.000000, -10.031778], 0],
      [20, 13, 0.01, 100, [-20.009954, -0.000000, -8.036143], 0],
      [20, 14, 0.01, 100, [-22.005589, -0.000000, -6.040508], 0],
      [20, 15, 0.01, 100, [-24.001263, -0.000039, -4.044911], 0],
      [20, 16, 0.01, 100, [-26.000720, -0.003860, -2.053097], 0],
      [20, 17, 0.01, 100, [-28.267550, -0.275056, -0.328657], 0],
      [20, 18, 0.01, 100, [-31.935098, -1.946968, -0.004935], 0],
      [20, 19, 0.01, 100, [-35.921484, -3.937719, -0.000050], 0],
      [20, 20, 0.01, 100, [-39.912704, -5.933304, -0.000001], 0],
      # Testing different values of p_error.
      [10, 0, 0.1, 100, [-0.001215, -2.553940, -9.543640], 25],
      [10, 1, 0.1, 100, [-0.010811, -1.609294, -7.644752], 16],
      [10, 0, 0.01, 100, [-0.000469, -2.967121, -19.956821], 29],
      [10, 1, 0.01, 100, [-0.044109, -1.015126, -16.009190], 10],
      [10, 0, 0.001, 100, [-0.000428, -3.006383, -29.996083], 30],
      [10, 1, 0.001, 100, [-0.297847, -0.304236, -24.294371], 3],
      [10, 0, 1e-04, 100, [-0.000424, -3.010290, -39.999990], 30],
      [10, 1, 1e-04, 100, [-1.032394, -0.042303, -33.032046], 0],
      # Test scaling of calculation with more coverage, hitting max_gq.
      [10, 0, 0.01, 100, [-0.000469, -2.967121, -19.956821], 29],
      [20, 0, 0.01, 100, [-0.000001, -5.933304, -39.912704], 59],
      [30, 0, 0.01, 100, [-0.000000, -8.899956, -59.869056], 88],
      [40, 0, 0.01, 100, [-0.000000, -11.866608, -79.825408], 100],
      [50, 0, 0.01, 100, [-0.000000, -14.833260, -99.781760], 100],
      [60, 0, 0.01, 100, [0.000000, -17.799911, -119.738112], 100],
      [70, 0, 0.01, 100, [0.000000, -20.766563, -139.694464], 100],
      [80, 0, 0.01, 100, [0.000000, -23.733215, -159.650816], 100],
      [90, 0, 0.01, 100, [0.000000, -26.699867, -179.607168], 100],
      [100, 0, 0.01, 100, [0.000000, -29.666519, -199.563519], 100],
  )
  def test_ref_calc(self, total_n, alt_n, p_error, max_gq, expected_likelihoods,
                    expected_gq):
    caller = self.make_test_caller(p_error, max_gq)
    gq, likelihoods = caller.reference_confidence(total_n - alt_n, total_n)
    npt.assert_allclose(expected_likelihoods, likelihoods, atol=1e-6)
    self.assertEqual(expected_gq, gq)

  @parameterized.parameters(
      # Values below max_allowed_reads are returned without modification.
      [0, 10, 100, (0, 10)],
      [5, 10, 100, (5, 10)],
      [10, 10, 100, (10, 10)],
      [10, 100, 100, (10, 100)],
      [100, 100, 100, (100, 100)],

      # Checks that the rescaling works when n_total_reads > max_allowed.
      [0, 200, 100, (0, 100)],
      [0, 200, 100, (0, 100)],
      [0, 1000, 100, (0, 100)],
      [0, 10000, 100, (0, 100)],
      [1, 200, 100, (1, 100)],
      [1, 1000, 100, (1, 100)],
      [1, 10000, 100, (1, 100)],
      [1, 100000, 100, (1, 100)],
      [2, 200, 100, (1, 100)],
      [3, 200, 100, (2, 100)],
      [4, 200, 100, (2, 100)],
      [10, 200, 100, (5, 100)],
      [50, 200, 100, (25, 100)],
      [100, 200, 100, (50, 100)],
      [200, 200, 100, (100, 100)],
      # I saw a bug at runtime, and the testcase makes sure we scale values of
      # n_ref_reads close to n_total_reads appropriately.
      [99, 100, 100, (99, 100)],
  )
  def test_rescale_read_counts(self, n_ref, n_total, max_allowed_reads,
                               expected):
    actual = variant_caller._rescale_read_counts_if_necessary(
        n_ref, n_total, max_allowed_reads)
    self.assertEqual(actual, expected)

  @parameterized.parameters((n_ref, n_alt_fraction)
                            for n_ref in [1000, 10000, 100000, 1000000]
                            for n_alt_fraction in [0.0, 0.01, 0.02])
  def test_handles_large_reference_counts(self, n_ref, n_alt_fraction):
    """Tests that we don't blow up when the coverage gets really high."""
    caller = self.make_test_caller(0.01, 100)
    n_alt = int(n_alt_fraction * n_ref)
    gq, likelihoods = caller._calc_reference_confidence(n_ref, n_ref + n_alt)
    self.assertTrue(
        np.isfinite(likelihoods).all(),
        'Non-finite likelihoods {}'.format(likelihoods))
    self.assertEqual(100, gq)

  @parameterized.parameters(*variant_caller.CANONICAL_DNA_BASES)
  def test_gvcf_basic(self, ref):
    options = _reference_model_options(0.01, 100)
    caller = variant_caller.VariantCaller(options, use_cache_table=False)
    allele_counter = self.fake_allele_counter(100, [(0, 0, ref)])
    gvcfs = list(caller.make_gvcfs(allele_counter.summary_counts()))
    self.assertLen(gvcfs, 1)
    self.assertGVCF(
        gvcfs[0],
        ref=ref,
        gq=1.0,
        start=100,
        end=101,
        chrom='chr1',
        gls=[-0.47712125472] * 3,
        sample_name=options.sample_name)

  @parameterized.parameters('N', 'R', 'W', 'B')
  def test_gvcf_basic_skips_iupac_ref_base(self, ref):
    options = _reference_model_options(0.01, 100)
    caller = variant_caller.VariantCaller(options, use_cache_table=False)
    allele_counter = self.fake_allele_counter(100, [(0, 0, ref)])
    self.assertEmpty(list(caller.make_gvcfs(allele_counter.summary_counts())))

  @parameterized.parameters('X', '>', '!')
  def test_gvcf_basic_raises_with_bad_ref_base(self, ref):
    options = _reference_model_options(0.01, 100)
    caller = variant_caller.VariantCaller(options, use_cache_table=False)
    allele_counter = self.fake_allele_counter(100, [(0, 0, ref)])
    with self.assertRaisesRegexp(ValueError,
                                 'Invalid reference base={}'.format(ref)):
      list(caller.make_gvcfs(allele_counter.summary_counts()))

  @parameterized.parameters(True, False)
  def test_calls_from_allele_counts(self, include_gvcfs):
    # Our test AlleleCounts are 5 positions:
    #
    # 10: A ref [no reads]
    # 11: G/C variant
    # 12: G ref [no reads]
    # 13: G ref [no reads]
    # 14: T/C variant
    #
    # The ref sites have no reads for ref or any alt simply because it
    # simplifies comparing them with the expected variant genotype likelihoods.
    # We aren't testing the correctness of the gvcf calculation here (that's
    # elsewhere) but rather focusing here on the separation of variants from
    # gvcf records, and the automatic merging of the gvcf blocks.
    allele_counter = self.fake_allele_counter(10, [
        (0, 0, 'A'),
        (10, 10, 'G'),
        (0, 0, 'G'),
        (0, 0, 'G'),
        (10, 10, 'T'),
    ])
    fake_candidates = [
        deepvariant_pb2.DeepVariantCall(
            variant=test_utils.make_variant(alleles=['G', 'C'], start=11)),
        deepvariant_pb2.DeepVariantCall(
            variant=test_utils.make_variant(alleles=['T', 'C'], start=14)),
    ]

    caller = self.make_test_caller(0.01, 100)
    with mock.patch.object(caller, 'cpp_variant_caller') as mock_cpp:
      mock_cpp.calls_from_allele_counter.return_value = fake_candidates
      candidates, gvcfs = caller.calls_from_allele_counter(
          allele_counter, include_gvcfs)

    mock_cpp.calls_from_allele_counter.assert_called_once_with(allele_counter)
    self.assertEqual(candidates, fake_candidates)

    # We expect our gvcfs to occur at the 10 position and that 12 and 13 have
    # been merged into a 2 bp block, if enabled. Otherwise should be empty.
    if include_gvcfs:
      self.assertLen(gvcfs, 4)
      # Expected diploid genotype likelihoods when there's no coverage. The
      # chance of having each genotype is 1/3, in log10 space.
      flat_gls = np.log10([1.0 / 3] * 3)
      self.assertGVCF(gvcfs[0], ref='A', start=10, end=11, gq=1, gls=flat_gls)
      self.assertGVCF(
          gvcfs[1],
          ref='G',
          start=11,
          end=12,
          gq=0,
          gls=np.array([-14.0230482368, -8.32667268469e-15, -14.0230482368]))
      self.assertGVCF(gvcfs[2], ref='G', start=12, end=14, gq=1, gls=flat_gls)
    else:
      self.assertEmpty(gvcfs)

  def assertGVCF(self,
                 gvcf,
                 ref,
                 gq,
                 start,
                 end,
                 chrom='chr1',
                 gls=None,
                 sample_name=None):
    if chrom:
      self.assertEqual(gvcf.reference_name, chrom)
    self.assertNotEmpty(gvcf.reference_name)
    self.assertEqual(gvcf.reference_bases, ref)
    self.assertEqual(gvcf.alternate_bases, ['<*>'])
    self.assertEqual(gvcf.start, start)
    self.assertEqual(gvcf.end, end if end else start + 1)
    self.assertEqual(len(gvcf.calls), 1)
    self.assertEqual(variantutils.genotype_quality(gvcf), gq)
    self.assertNotEmpty(gvcf.calls[0].genotype_likelihood)
    if gls is not None:
      npt.assert_allclose(list(gvcf.calls[0].genotype_likelihood), gls)
    if sample_name:
      self.assertEqual(gvcf.calls[0].call_set_name, sample_name)

  @parameterized.parameters(
      # Check some basics.
      ([(0, 0, 'A')], [dict(start=1, end=2, ref='A', gq=1)]),
      # Two equal records are merged, and the reference base is the first one.
      ([(0, 0, 'A'), (0, 0, 'C')], [dict(start=1, end=3, ref='A', gq=1)]),
      ([(0, 0, 'C'), (0, 0, 'A')], [dict(start=1, end=3, ref='C', gq=1)]),
      # Three equal records are merged into a single block.
      ([(0, 0, 'A'), (0, 0, 'C'),
        (0, 0, 'T')], [dict(start=1, end=4, ref='A', gq=1)]),
      # We don't merge together different GQ value blocks:
      ([(0, 0, 'A'), (0, 100, 'C')], [
          dict(start=1, end=2, ref='A', gq=1),
          dict(start=2, end=3, ref='C', gq=100),
      ]),
      ([(0, 100, 'A'), (0, 0, 'C')], [
          dict(start=1, end=2, ref='A', gq=100),
          dict(start=2, end=3, ref='C', gq=1),
      ]),
      ([(0, 0, 'A'), (0, 20, 'C'), (0, 100, 'T')], [
          dict(start=1, end=2, ref='A', gq=1),
          dict(start=2, end=3, ref='C', gq=59),
          dict(start=3, end=4, ref='T', gq=100),
      ]),
  )
  def test_make_gvcfs(self, counts, expecteds):
    allele_counts = self.fake_allele_counter(1, counts).summary_counts()
    caller = self.make_test_caller(0.01, 100)
    gvcfs = list(caller.make_gvcfs(allele_counts))

    self.assertLen(gvcfs, len(expecteds))
    for actual, expected in zip(gvcfs, expecteds):
      self.assertGVCF(actual, **expected)


_CACHE_COVERAGE = 20  # Outside class so we can refer to it in @Parameters.


class VariantCallerCacheTests(parameterized.TestCase):

  @classmethod
  def setUpClass(cls):
    options = _reference_model_options(0.1, 50)
    cls.raw_caller = variant_caller.VariantCaller(
        options, use_cache_table=False)
    cls.cache_caller = variant_caller.VariantCaller(
        options, use_cache_table=True, max_cache_coverage=_CACHE_COVERAGE)

  @parameterized.parameters((n_alt, n_total)
                            for n_total in range(_CACHE_COVERAGE + 1)
                            for n_alt in range(n_total + 1))
  def test_caching(self, n_alt, n_total):
    # Note that we only expect the gq and gls to be close if we are not
    # rescaling the counts, so we are only looping over values that should be
    # cached. In practice the cache is set to values sufficiently large that
    # these differences don't matter, but for this test we are limiting the
    # cache size to a small value in _CACHE_COVERAGE so we can test that the
    # cache lookups are correct.
    raw_gq, raw_gls = self.raw_caller.reference_confidence(n_alt, n_total)
    cache_gq, cache_gls = self.cache_caller.reference_confidence(n_alt, n_total)
    self.assertEqual(raw_gq, cache_gq)
    npt.assert_allclose(raw_gls, cache_gls)


if __name__ == '__main__':
  absltest.main()
