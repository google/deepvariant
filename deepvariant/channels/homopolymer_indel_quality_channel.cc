/*
 * Copyright 2025 Google LLC.
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

#include "deepvariant/channels/homopolymer_indel_quality_channel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/channels/channel_utils.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

HomopolymerInDelQualityChannel::HomopolymerInDelQualityChannel(
    int width, const PileupImageOptions& options)
    : Channel(width, options) {}

// Explanation of what TP means in Ultima Genomics data:
// Given read R of size N,
// Cram will have a QUAL vector of size N, and a TP vector of size N.
//
// TP[i] refers to the direction and magnitude of the error which is encoded by
// QUAL[i]
//
// For example:
// If read has sequence ...TAAAAAG...,
// decided_homopolymer_size for the A homopolymer is len(AAAAA) = 5.
// Let i be an index in the read in the range of this homopolymer: [start, end].
// When TP[i] = 1, it means QUAL[i] will contain the PHRED score for the
// homopolymer to be 5 + 1 = 6.
// When TP[i] = -1, it means QUAL[i] will contain the PHRED score for the
// homopolymer to be 5 - 1 = 4.
// When TP[i] = 2, it means QUAL[i] will contain the PHRED score for the
// homopolymer to be 5 + 2 = 7.
std::vector<int8_t> HomopolymerInDelQualityChannel::GetTPValues(
    const Read& read) {
  std::vector<int8_t> int_tps(read.aligned_sequence().size());
  if (!read.info().contains("tp")) {
    // Return a vector of zeros if tp tag is not present.
    return int_tps;
  }
  const auto& tps = read.info().at("tp").values();

  if (tps.empty()) {
    return int_tps;
  }

  for (int i = 0; i < tps.size() && i < int_tps.size(); i++) {
    int_tps[i] = tps[i].int_value();
  }
  return int_tps;
}

std::vector<std::uint8_t> HomopolymerInDelQualityChannel::HomoPolymerWeighted(
    const Read& read) {
  // Generates a vector reflecting the number of repeats observed
  std::vector<std::uint8_t> homopolymer_weighted(read.aligned_sequence().size(),
                                                 1);
  const std::string& seq = read.aligned_sequence();

  if (seq.empty()) {
    return homopolymer_weighted;
  }

  int i = 0;
  while (i < seq.size()) {
    int hmer_length = 1;
    char current_base = seq[i];
    int j = i + 1;

    // Count consecutive identical bases
    while (j < seq.size() && seq[j] == current_base) {
      hmer_length++;
      j++;
    }

    // Fill all positions in this homopolymer with its length
    for (int k = i; k < j; k++) {
      // Maximum value of uint8_t is 255.
      homopolymer_weighted[k] = std::min(
          hmer_length, static_cast<int>(std::numeric_limits<uint8_t>::max()));
    }

    i = j;
  }

  return homopolymer_weighted;
}

std::vector<std::uint8_t>
HomopolymerInDelQualityChannel::HomoPolymerInDelQuality(const Read& read,
                                                        bool is_deletion) {
  // Decode base quality and tp tags to generate a vector indicating the
  // probability of not having an hmer deletion/insertion for each homopolymer
  // in the read. The probabilities are encoded in phred-scores.
  // The parameter is_deletion determines the direction of the computed error
  // rate per homopolymer.
  std::vector<std::uint8_t> hmer_directed_qualities(
      read.aligned_sequence().size(),
      channels::internal::BaseQualityColor(kMaxQScore));

  auto seq = read.aligned_sequence();
  auto hmer_lengths = HomoPolymerWeighted(read);
  auto tps = GetTPValues(read);

  // If tp tag is not present, return default quality values
  if (tps.empty() || tps.size() != seq.size()) {
    return hmer_directed_qualities;
  }

  int i = 0;
  while (i < hmer_lengths.size()) {
    int hmer_length = hmer_lengths[i];
    float hmer_directed_error_prob = 0;

    // Iterate quality encodings for hmer [i], and sum them up to upward
    // direction or downward direction quality PHRED scores.
    for (int j = 0; j < hmer_length; j++) {
      if (tps[i + j] == 0) {
        continue;
      }
      bool is_deletion_err = tps[i + j] < 0;
      if (is_deletion_err == is_deletion) {
        std::uint8_t encoded_hmer_qual = read.aligned_quality()[i + j];
        float error_prob = std::pow(10, (encoded_hmer_qual / -10.0));
        hmer_directed_error_prob += error_prob;
      }
    }

    int hmer_directed_quality =
        hmer_directed_error_prob == 0
            ? kMaxQScore
            : static_cast<int>(-10 * std::log10(hmer_directed_error_prob));
    // Clamp to valid range
    if (hmer_directed_quality > kMaxQScore) {
      hmer_directed_quality = kMaxQScore;
    }

    for (int j = 0; j < hmer_length; j++) {
      hmer_directed_qualities[i + j] =
          channels::internal::BaseQualityColor(hmer_directed_quality);
    }
    i += hmer_length;
  }
  return hmer_directed_qualities;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
