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

#include "deepvariant/channels/read_supports_variant_channel.h"

#include <algorithm>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/log/check.h"
#include "absl/types/span.h"
namespace learning {
namespace genomics {
namespace deepvariant {

ReadSupportsVariantChannel::ReadSupportsVariantChannel(
    int width,
    const learning::genomics::deepvariant::PileupImageOptions& options)
    : Channel(width, options) {
  supports_variant_color_ = std::nullopt;
}

void ReadSupportsVariantChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  if (!supports_variant_color_.has_value()) {
    int read_supports_alt = ReadSupportsAlt(dv_call, read, alt_alleles);
    supports_variant_color_ = std::optional<unsigned char>{
        static_cast<unsigned char>(SupportsAltColor(read_supports_alt))};
  }
  data[col] = supports_variant_color_.value();
}

void ReadSupportsVariantChannel::FillRefBase(
    std::vector<unsigned char>& ref_data, int col, char ref_base,
    const std::string& ref_bases) {
  ref_data[col] = SupportsAltColor(0);
}

// Does this read support ref, one of the alternative alleles, or an allele we
// aren't considering?
int ReadSupportsVariantChannel::ReadSupportsAlt(
    const DeepVariantCall& dv_call, const Read& read,
    absl::Span<const std::string> alt_alleles) {
  std::string key =
      (read.fragment_name() + "/" + std::to_string(read.read_number()));

  // Iterate over all alts, not just alt_alleles.
  for (const std::string& alt_allele : dv_call.variant().alternate_bases()) {
    const auto& allele_support = dv_call.allele_support();
    const bool alt_allele_present_in_call =
        allele_support.find(alt_allele) != allele_support.cend();

    if (alt_allele_present_in_call) {
      const auto& supp_read_names = allele_support.at(alt_allele).read_names();
      for (const std::string& read_name : supp_read_names) {
        const bool alt_in_alt_alleles =
            std::find(alt_alleles.begin(), alt_alleles.end(), alt_allele) !=
            alt_alleles.end();
        // Read can support an alt we are currently considering (1), a different
        // alt not present in alt_alleles (2), or ref (0).
        if (read_name == key && alt_in_alt_alleles) {
          return 1;
        } else if (read_name == key && !alt_in_alt_alleles) {
          return 2;
        }
      }
    }
  }
  return 0;
}

int ReadSupportsVariantChannel::SupportsAltColor(int read_supports_alt) const {
  float alpha;
  if (read_supports_alt == 0) {
    alpha = options_.allele_unsupporting_read_alpha();
  } else if (read_supports_alt == 1) {
    alpha = options_.allele_supporting_read_alpha();
  } else {
    CHECK_EQ(read_supports_alt, 2) << "read_supports_alt can only be 0/1/2.";
    alpha = options_.other_allele_supporting_read_alpha();
  }
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
