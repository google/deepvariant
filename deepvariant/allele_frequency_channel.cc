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

#include "deepvariant/allele_frequency_channel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "absl/types/span.h"

namespace learning {
namespace genomics {
namespace deepvariant {

void AlleleFrequencyChannel::FillReadLevelData(
    const Read& read, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles,
    std::vector<unsigned char>& read_level_data) {
  const float allele_frequency =
      ReadAlleleFrequency(dv_call, read, alt_alleles);
  read_level_data = std::vector<unsigned char>(
      1, static_cast<std::uint8_t>(AlleleFrequencyColor(allele_frequency)));
}
void AlleleFrequencyChannel::FillRefData(const std::string& ref_bases,
                                         std::vector<unsigned char>& ref_data) {
  int allele_frequency_color = AlleleFrequencyColor(0);
  ref_data = std::vector<unsigned char>(
      width_, static_cast<std::uint8_t>(allele_frequency_color));
}

// Get allele frequency color for a read.
// Convert a frequency value in float to color intensity (int) and normalize.
unsigned char AlleleFrequencyChannel::AlleleFrequencyColor(
    float allele_frequency) {
  if (allele_frequency <= options_.min_non_zero_allele_frequency()) {
    return 0;
  } else {
    float log10_af = log10(allele_frequency);
    float log10_min = log10(options_.min_non_zero_allele_frequency());
    return ((log10_min - log10_af) / log10_min) *
           static_cast<int>(kMaxPixelValueAsFloat);
  }
}

// Get the allele frequency of the alt allele that is carried by a read.
float AlleleFrequencyChannel::ReadAlleleFrequency(
    const DeepVariantCall& dv_call, const Read& read,
    absl::Span<const std::string> alt_alleles) {
  std::string key =
      (read.fragment_name() + "/" + std::to_string(read.read_number()));

  // Iterate over all alts, not just alt_alleles.
  for (const std::string& alt_allele : dv_call.variant().alternate_bases()) {
    const auto& allele_support = dv_call.allele_support();
    auto it_read = allele_support.find(alt_allele);

    if (it_read != allele_support.end()) {
      const auto& supp_read_names = it_read->second.read_names();
      for (const std::string& read_name : supp_read_names) {
        const bool alt_in_alt_alleles =
            std::find(alt_alleles.begin(), alt_alleles.end(), alt_allele) !=
            alt_alleles.end();
        // If the read supports an alt we are currently considering, return the
        // associated allele frequency.
        if (read_name == key && alt_in_alt_alleles) {
          auto it = dv_call.allele_frequency().find(alt_allele);
          if (it != dv_call.allele_frequency().end())
            return it->second;
          else
            return 0;
        }
      }
    }
  }
  // If cannot find the matching variant, set the frequency to 0.
  return 0;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
