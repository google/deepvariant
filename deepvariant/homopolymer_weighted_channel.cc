/*
 * Copyright 2024 Google LLC.
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

#include "deepvariant/homopolymer_weighted_channel.h"

#include <cstdint>
#include <string>
#include <vector>

namespace learning {
namespace genomics {
namespace deepvariant {
void HomopolymerWeightedChannel::FillReadLevelData(
    const Read& read, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles,
    std::vector<unsigned char>& read_level_data) {
  std::vector<std::uint8_t> homopolymer_weighted = HomopolymerWeighted(read);
  read_level_data =
      ScaleColorVector(homopolymer_weighted, kMaxHomopolymerWeighted);
}
void HomopolymerWeightedChannel::FillRefData(
    const std::string& ref_bases, std::vector<unsigned char>& ref_data) {
  Read refRead;
  refRead.set_aligned_sequence(ref_bases);
  std::vector<std::uint8_t> homopolymer_weighted = HomopolymerWeighted(refRead);
  ref_data = ScaleColorVector(homopolymer_weighted, kMaxHomopolymerWeighted);
}

std::vector<std::uint8_t> HomopolymerWeightedChannel::HomopolymerWeighted(
    const Read& read) {
  // Generates a vector reflecting the number of repeats observed.
  // ATCGGGAA
  // 11133322
  std::vector<std::uint8_t> homopolymer(read.aligned_sequence().size());
  auto seq = read.aligned_sequence();
  homopolymer[0] = 1;
  int current_weight = 1;
  for (int i = 1; i <= seq.size(); i++) {
    if (seq[i] == seq[i - 1]) {
      current_weight += 1;
    } else {
      for (int cw = current_weight; cw >= 1; cw--) {
        homopolymer[i - cw] = current_weight;
      }
      current_weight = 1;
    }
  }
  return homopolymer;
}

// Scales an input vector to pixel range 0-254
std::vector<std::uint8_t> HomopolymerWeightedChannel::ScaleColorVector(
    std::vector<std::uint8_t>& channel_values, float max_val) {
  for (int i = 0; i < channel_values.size(); i++) {
    int value = channel_values[i];
    if (static_cast<float>(value) > max_val) {
      value = max_val;
    }
    channel_values[i] = static_cast<int>(kMaxPixelValueAsFloat *
                                         (static_cast<float>(value) / max_val));
  }
  return channel_values;
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
