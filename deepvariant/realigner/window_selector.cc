/*
 * Copyright 2017 Google LLC.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
#include "deepvariant/realigner/window_selector.h"

#include <algorithm>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/protos/realigner.pb.h"
#include "absl/log/check.h"
#include "absl/log/log.h"

namespace learning {
namespace genomics {
namespace deepvariant {

// Helper function to update counts from `start` (inclusive) to `end`
// (exclusive) by `by`. This function also provides some tolerance to invalid
// start and end values: if `start < 0`, a value of `start = 0` is used instead.
// If `end > counts->size()` then `end = counts->size()` will be used instead.
// This simplifies the call site where the bounding of start/end can be
// overloaded to this function instead of duplicating it at all call sites.
template <class T>
void UpdateCounts(T by, int start, int end, std::vector<T>* counts) {
  CHECK_LE(start, end) << "Start should be <= end";

  start = std::max(start, 0);
  end = std::min(end, static_cast<int>(counts->size()));
  for (int i = start; i < end; i++) {
    (*counts)[i] += by;
  }
}

bool AlleleFilter(
    const Allele& allele, int total_count,
    const WindowSelectorOptions& config) {
  if (allele.type() == AlleleType::REFERENCE) {
    return false;
  }
  // No alleles with read support less than 2 are used. The same filter is used
  // for candidate creaetion later. However, candidates filter is more strict
  // and also filters out low quality alleles and allele fractions.
  if (allele.count() < config.min_allele_support()) {
    return false;
  }
  return true;
}


std::vector<Allele> SelectAltAlleles(
    const AlleleCount& allele_count,
    const WindowSelectorOptions& config
    ) {
  const std::vector<Allele> target_sample_alleles =
      SumAlleleCounts(allele_count);
  const int target_samples_total_count =
      TotalAlleleCounts(allele_count);

  std::vector<Allele> alt_alleles;
  for (const auto& allele : target_sample_alleles) {
    if (AlleleFilter(allele, target_samples_total_count, config)) {
      alt_alleles.push_back(allele);
    }
  }
  return alt_alleles;
}

std::vector<int> VariantReadsWindowSelectorCandidates(
    const AlleleCounter& allele_counter, const WindowSelectorOptions& config) {
  // We start with a vector of 0s, one for each position in allele_counter.
  std::vector<int> window_counts(allele_counter.Counts().size(), 0);

  // Now loop over all of the counts. Alleles are created for each position
  // in the allele_counter. window_counts are increamented for each allele by
  // the number of reads supporting the allele.
  const std::vector<AlleleCount>& counts = allele_counter.Counts();
  for (int i = 0; i < counts.size(); ++i) {
    auto alt_alleles = SelectAltAlleles(counts[i], config);
    for (const auto& allele : alt_alleles) {
      int start = 0;
      int end = 0;
      switch (allele.type()) {
        case SUBSTITUTION:
          start = i; end = i + 1;
          UpdateCounts(allele.count(), start, end, &window_counts);
          break;
        case SOFT_CLIP:
        case INSERTION:
          start = i + 1 - (allele.bases().length() - 1);
          end = i + allele.bases().length();
          UpdateCounts(allele.count(), start, end, &window_counts);
          break;
        case DELETION:
          start = i + 1;
          end = i + allele.bases().length();
          UpdateCounts(allele.count(), start, end, &window_counts);
          break;
        case REFERENCE:
          // We don't update our counts for reference positions.
          break;
        case UNSPECIFIED:
        default:
          LOG(FATAL) << "Saw an Allele " << allele.DebugString()
                     << " with an unexpected type " << allele.type()
                     << " in AlleleCount " << counts[i].DebugString()
                     << " which should never happen.";
      }
    }
  }

  return window_counts;
}

std::vector<float> AlleleCountLinearWindowSelectorCandidates(
    const AlleleCounter& allele_counter,
    const WindowSelectorModel::AlleleCountLinearModel& config) {
  std::vector<float> window_scores(allele_counter.Counts().size(),
                                   config.bias());

  const std::vector<AlleleCount>& counts = allele_counter.Counts();
  for (int i = 0; i < counts.size(); ++i) {
    UpdateCounts(
        counts[i].ref_supporting_read_count() * config.coeff_reference(), i,
        i + 1, &window_scores);

    for (const auto& entry : counts[i].read_alleles()) {
      const Allele& allele = entry.second;

      int start, end;
      switch (allele.type()) {
        case SUBSTITUTION:
          start = i;
          end = i + 1;
          UpdateCounts(allele.count() * config.coeff_substitution(), start,
                       end, &window_scores);
          break;
        case SOFT_CLIP:
          start = i + 1 - (allele.bases().length() - 1);
          end = i + allele.bases().length();
          UpdateCounts(allele.count() * config.coeff_soft_clip(), start, end,
                       &window_scores);
          break;
        case INSERTION:
          start = i + 1 - (allele.bases().length() - 1);
          end = i + allele.bases().length();
          UpdateCounts(allele.count() * config.coeff_insertion(), start, end,
                       &window_scores);
          break;
        case DELETION:
          start = i + 1;
          end = i + allele.bases().length();
          UpdateCounts(allele.count() * config.coeff_deletion(), start, end,
                       &window_scores);
          break;
        case REFERENCE:
          start = i;
          end = i + 1;
          UpdateCounts(allele.count() * config.coeff_reference(), start, end,
                       &window_scores);
          break;
        case UNSPECIFIED:
        default:
          LOG(FATAL) << "Saw an Allele " << allele.DebugString()
                     << " with an unexpected type " << allele.type()
                     << " in AlleleCount " << counts[i].DebugString()
                     << " which should never happen.";
      }
    }
  }

  return window_scores;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

