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
"""A VerySensitiveCaller producing DeepVariantCall and gVCF records.

This module provides the primary interface for calling candidate variants using
for the AlleleCounts in an AlleleCounter by wrapping the low-level C++ code and
adding a nicer API and functions to compute gVCF records as well.
"""

from typing import Dict, Sequence

from deepvariant import variant_caller
from deepvariant.protos import deepvariant_pb2
from deepvariant.python import allelecounter


class VerySensitiveCaller(variant_caller.VariantCaller):
  """Call variants and gvcf records from an AlleleCounter."""

  def __init__(self, options, use_cache_table=True, max_cache_coverage=100):
    super(VerySensitiveCaller, self).__init__(
        options=options,
        use_cache_table=use_cache_table,
        max_cache_coverage=max_cache_coverage)

  def get_candidates(
      self, allele_counters: Dict[str, allelecounter.AlleleCounter],
      sample_name: str) -> Sequence[deepvariant_pb2.DeepVariantCall]:
    # AlleleCounter doesn't have copy constructor therefore we pass
    # allele counts only (which is a member of AlleleCounter).
    allele_counts = {
        sample_id: allele_counters[sample_id].counts()
        for sample_id in allele_counters
    }
    return self.cpp_variant_caller.calls_from_allele_counts(
        allele_counts, sample_name)

  def get_candidate_positions(
      self, allele_counters: Dict[str, allelecounter.AlleleCounter],
      sample_name: str):
    allele_counts = {
        sample_id: allele_counters[sample_id].counts()
        for sample_id in allele_counters
    }
    return self.cpp_variant_caller.call_positions_from_allele_counts(
        allele_counts, sample_name)
