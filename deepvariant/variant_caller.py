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
"""A VariantCaller producing DeepVariantCall and gVCF records.

This module provides the primary interface for calling candidate variants using
for the AlleleCounts in an AlleleCounter by wrapping the low-level C++ code and
adding a nicer API and functions to compute gVCF records as well.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

import itertools
import math as py_math
import operator


import numpy as np

from deepvariant.core import math
from deepvariant.core import variantutils
from deepvariant.core.genomics import variants_pb2
from deepvariant.python import variant_calling

# Reference bases with genotype calls must be one of these four values.
CANONICAL_DNA_BASES = frozenset('ACGT')

# Possible DNA base codes seen in a reference genome.
EXTENDED_IUPAC_CODES = frozenset('ACGTRYSWKMBDHVN')


def _rescale_read_counts_if_necessary(n_ref_reads, n_total_reads,
                                      max_allowed_reads):
  """Ensures that n_total_reads <= max_allowed_reads, rescaling if necessary.

  This function ensures that n_total_reads <= max_allowed_reads. If
  n_total_reads is <= max_allowed_reads, n_ref_reads and n_total_reads are just
  returned. However, if n_total_reads > max_allowed_reads, then n_ref_reads and
  n_total_reads are rescaled to new values n_ref_reads' and n_total_reads' so
  that n_total_reads' == max_allowed_reads and n_ref_reads' / n_total_reads' ~
  n_ref_reads / n_total_reads.

  Args:
    n_ref_reads: int. Number of reference supporting reads.
    n_total_reads: int. Total number of reads.
    max_allowed_reads: int. The maximum value allowed for n_total after
      rescaling, if necessary.

  Returns:
    New values for n_ref_reads and n_total_reads.
  """
  if n_total_reads > max_allowed_reads:
    ratio = n_ref_reads / (1.0 * n_total_reads)
    n_ref_reads = int(py_math.ceil(ratio * max_allowed_reads))
    n_total_reads = max_allowed_reads
  return n_ref_reads, n_total_reads


class VariantCaller(object):
  """Call variants and gvcf records from an AlleleCounter."""

  def __init__(self, options, use_cache_table=True, max_cache_coverage=100):
    self.options = options
    self.cpp_variant_caller = variant_calling.VariantCaller(self.options)

    self.max_cache_coverage = max_cache_coverage
    if use_cache_table:
      self.table = [[
          self._calc_reference_confidence(n_ref, n_total)
          for n_ref in range(n_total + 1)
      ]
                    for n_total in range(self.max_cache_coverage + 1)]
    else:
      self.table = None

  def reference_confidence(self, n_ref, n_total):
    """Computes the confidence that a site in the genome has no variation.

    Computes this confidence using only the counts of the number of reads
    supporting the reference allele and the total number of reads at this site.

    See: https:www.broadinstitute.org/gatk/guide/article?id=4017 for
    background.

    Computes the reference confidence for site allele_count.

    Examines the number of reference supporting and alternate supporting reads
    in allele_count and estimates the genotype likelihoods and confidence that
    this site is homozygous reference. These values are written into the first
    VariantCall record of variant, into the repeated field genotype_likelihood
    and the map field GQ.

    The genotype likelihoods are computed against any possible alternative
    allele, the so-called <*> allele, which boils down to a model that looks
    like:

     log10_p_ref = Log10PBinomial(non_ref_n, total_n, p_error_);
     log10_p_het = Log10PBinomial(ref_n, total_n, 0.5);
     log10_p_hom_alt = Log10PBinomial(ref_n, total_n, p_error_);

    ref_n is the number of reference supporting reads and non_ref_n is the sum
    of any reads supporting any alternate alleles. Non-informative reads are
    excluded from the calculation.

    and written in as the normalized log10 values so that:

     sum(10^genotype_likelihoods) = 1

    The GQ, according to the VCF specification, is the conditional genotype
    quality, encoded as a phred quality -10 * log10 p(genotype call is wrong,
    conditioned on the site's being variant, as an integer. See:
    https:samtools.github.io/hts-specs/VCFv4.3.pdf
    We are calculating the GQ not for the best genotype, but the GQ of the 0/0
    genotype, regardless of the likelihoods.
    1 = pRR + pRA + pAA; [R is reference, A=<*> is any alternative alternative]
    GQ of 0/0 = -10 * log10(pRA + pAA) [prob that any other differen genotype]
              = -10 * log10(1 - pRR) [substitution from the previous equation]
    Here we don't have pRR directly, but rather log10(pRR).

    Args:
      n_ref: int >= 0 and <= n_total: The number of reads supporting the
        reference allele.
      n_total: int >= 0 and >= n_ref: The number of reads supporting any allele
        at this site.

    Returns:
      A tuple of two values. The first is an integer value for the GQ (genotype
      quality) and the second is an array-like of the log10 probabilities for
      each of the three genotype configurations.
    """
    if self.table is None:
      return self._calc_reference_confidence(n_ref, n_total)
    else:
      ref_index, total_index = _rescale_read_counts_if_necessary(
          n_ref, n_total, self.max_cache_coverage)
      return self.table[total_index][ref_index]

  def _calc_reference_confidence(self, n_ref, n_total):
    """Performs the calculation described in reference_confidence()."""
    if n_ref < 0:
      raise ValueError('n_ref={} must be >= 0'.format(n_ref))
    if n_total < n_ref:
      raise ValueError('n_total={} must be >= n_ref={}'.format(n_total, n_ref))
    if self.options.ploidy != 2:
      raise ValueError('ploidy={} but we only support ploidy=2'.format(
          self.options.ploidy))

    if n_total == 0:
      # No coverage case - all likelihoods are log10 of 1/3, 1/3, 1/3.
      log10_probs = math.normalize_log10_probs([-1.0, -1.0, -1.0])
    else:
      n_alts = n_total - n_ref
      log10_p_ref = math.log10_binomial(n_alts, n_total, self.options.p_error)
      log10_p_het = math.log10_binomial(n_ref, n_total,
                                        1.0 / self.options.ploidy)
      log10_p_hom_alt = math.log10_binomial(n_ref, n_total,
                                            self.options.p_error)
      log10_probs = math.normalize_log10_probs(
          [log10_p_ref, log10_p_het, log10_p_hom_alt])

    gq = math.log10_ptrue_to_phred(log10_probs[0], self.options.max_gq)
    gq = int(min(np.floor(gq), self.options.max_gq))
    return gq, log10_probs

  def make_gvcfs(self, allele_count_summaries):
    """Primary interface function for computing gVCF confidence at a site.

    Looks at the counts in the provided list of AlleleCountSummary protos and
    returns properly-formatted Variant protos containing gVCF reference
    blocks for all sites in allele_count_summaries. The returned Variant has
    reference_name, start, end are set and contains a single VariantCall in the
    calls field with call_set_name of options.sample_name, genotypes set to 0/0
    (diploid reference), and a GQ value bound in the info field appropriate to
    the data in allele_count.

    The provided allele count must have either a canonical DNA sequence base (
    A, C, G, T) or be "N".

    Args:
      allele_count_summaries: iterable of AlleleCountSummary protos in
        coordinate-sorted order. Each proto is used to get the read counts for
        reference and alternate alleles, the reference position, and reference
        base.

    Yields:
      learning.genomics.deepvariant.core.genomics.Variant proto in
      coordinate-sorted order containing gVCF records.
    """

    def with_gq_and_likelihoods(summary_counts):
      """Returns summary_counts along with GQ and genotype likelihoods.

      If the reference base is not in CANONICAL_DNA_BASES, both GQ and genotype
      likelihoods are set to None.

      Args:
        summary_counts: A single AlleleCountSummary.

      Returns:
        A tuple of summary_counts, GQ, and genotype likelihoods for
        summary_counts where GQ and genotype_likelihood are calculated by
        self.reference_confidence.

      Raises:
        ValueError: The reference base is not a valid DNA or IUPAC base.
      """
      if summary_counts.ref_base not in CANONICAL_DNA_BASES:
        if summary_counts.ref_base in EXTENDED_IUPAC_CODES:
          # Skip calculating gq and likelihoods, since this is an ambiguous
          # reference base.
          gq, likelihoods = None, None
        else:
          raise ValueError('Invalid reference base={} found during gvcf '
                           'calculation'.format(summary_counts.ref_base))
      else:
        n_ref = summary_counts.ref_supporting_read_count
        n_total = summary_counts.total_read_count
        gq, likelihoods = self.reference_confidence(n_ref, n_total)
      return summary_counts, gq, likelihoods

    # Combines contiguous, compatible single-bp blocks into larger gVCF blocks,
    # respecting non-reference variants interspersed among them. Yields each
    # combined gVCF Variant proto, in order. Compatible right now means that the
    # blocks to be merged have the same non-None GQ value.
    for key, combinable in itertools.groupby(
        (with_gq_and_likelihoods(sc) for sc in allele_count_summaries),
        key=operator.itemgetter(1)):
      if key is None:
        # A None key indicates that a non-DNA reference base was encountered, so
        # skip this group.
        continue
      combinable = list(combinable)
      summary_counts, gq, likelihoods = combinable[0]
      call = variants_pb2.VariantCall(
          call_set_name=self.options.sample_name,
          genotype=[0, 0],
          genotype_likelihood=likelihoods)
      variantutils.set_variantcall_gq(call, gq)
      yield variants_pb2.Variant(
          reference_name=summary_counts.reference_name,
          reference_bases=summary_counts.ref_base,
          alternate_bases=[variantutils.GVCF_ALT_ALLELE],
          start=summary_counts.position,
          end=combinable[-1][0].position + 1,
          calls=[call])

  def calls_from_allele_counter(self, allele_counter, include_gvcfs):
    """Gets variant calls and gvcf records for all sites in allele_counter.

    Args:
      allele_counter: AlleleCounter object holding the allele counts we will use
        to find candidate variants and create gvcf records.
      include_gvcfs: boolean. If True, we will compute gVCF records for all of
        the AlleleCounts in AlleleCounter.

    Returns:
      Two values. The first is a list of DeepVariantCall protos containing our
      candidate variants. The second is a list of gVCF blocks in Variant proto
      format, if include_gvcfs is True. If False, an empty list is returned.
    """
    candidates = self.cpp_variant_caller.calls_from_allele_counter(
        allele_counter)

    gvcfs = []
    if include_gvcfs:
      gvcfs = list(self.make_gvcfs(allele_counter.summary_counts()))
    return candidates, gvcfs
