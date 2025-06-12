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

// This code implements the ReadSupportsVariantFuzzyChannel, a new image channel
// for DeepVariant designed to improve indel calling. The channel highlights
// reads that provide partial or "fuzzy" support for a candidate alternate
// allele.
// This is achieved by identifying reads that support an indel of a slightly
// different size (e.g., +/- 1-2 bases) but share the same haplotype as the main
// candidate variant. By doing so, it provides the model with evidence from a
// wider set of reads that belong to the same phase.
// This feature relies critically on phasing information (comparing a read's HP
// tag to a variant's PS tag) and is therefore only suitable for use with
// long-read sequencing data where robust phasing is available.

#include "deepvariant/channels/read_supports_variant_fuzzy_channel.h"

#include <algorithm>
#include <cstdlib>
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

// These values are chosen to make the color close to the one that is used for
// reads fully supporting alt alleles (which is 1.0).
float kReadSupportAltWithinOneBase = 0.90;
float kReadSupportAltWithinTwoBases = 0.80;
float kReadSupportAltWithinThreeBases = 0.70;

ReadSupportsVariantFuzzyChannel::ReadSupportsVariantFuzzyChannel(
    int width,
    const learning::genomics::deepvariant::PileupImageOptions& options)
    : Channel(width, options) {
  supports_variant_color_ = std::nullopt;
}

void ReadSupportsVariantFuzzyChannel::FillReadBase(
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

void ReadSupportsVariantFuzzyChannel::FillRefBase(
    std::vector<unsigned char>& ref_data, int col, char ref_base,
    const std::string& ref_bases) {
  ref_data[col] = SupportsAltColor(0);
}


// The ReadSupportsAlt method is the core logic of the fuzzy channel. For a
// given read, it determines how it supports the candidate variant and returns
// an integer code to represent that support level.
// Its categorization is as follows:
// * Exact Support (returns 1): The read perfectly supports one of the primary
//   alternate alleles that the image is being generated for.
// * Fuzzy Support (returns 10 or 9): The read supports a different indel that
//   is not the primary alternate, but shares the same haplotype phase (by
//   comparing the read's HP tag with the variant's ALT_PS info). The return
//   value indicates the closeness of the indel size (10 for a 1bp difference,
//   9 for a 2bp difference).
// * Reference Support (returns 0): The read does not support any of the
//   alternate alleles and is treated as supporting the reference.
// * Other Alt Support (returns 2): The read supports a different alternate
//   allele that does not qualify for fuzzy matching.
int ReadSupportsVariantFuzzyChannel::ReadSupportsAlt(
     const DeepVariantCall& dv_call, const Read& read,
     absl::Span<const std::string> alt_alleles) {
  std::string key =
       (read.fragment_name() + "/" + std::to_string(read.read_number()));

  const int num_alt_alleles = dv_call.variant().alternate_bases().size();
  // Store haplotag value for each alt allele in the candidate.
  std::vector<int> alt_allele_phases(num_alt_alleles, 0);
  if (dv_call.variant().info().contains("ALT_PS")) {
    const auto& alt_ps = dv_call.variant().info().at("ALT_PS").values();
    int allele_number = 0;
    for (int alt_allele_index = 0;
        alt_allele_index <  dv_call.variant().alternate_bases().size();
        alt_allele_index++) {
      if (alt_ps.size() > alt_allele_index+1) {
        alt_allele_phases[allele_number] =
            alt_ps[alt_allele_index+1].int_value();
      }
      allele_number++;
    }
  }

  // Candidate may have many alt alleles, but pileup image is created with only
  // one or two alt alleles. Here we try to find alt alleles from the candidate
  // that are close to the alt alleles in the pileup image.
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
        // Find allele from alt_alleles that has the same phase as alt_allele.
        for (int image_alt_allele_index = 0;
            image_alt_allele_index < alt_alleles.size();
            image_alt_allele_index++) {
          int image_alt_allele_global_index = 0;
          for (const auto& candidate_alt_allele :
                   dv_call.variant().alternate_bases()) {
            if (candidate_alt_allele == alt_alleles[image_alt_allele_index]) {
              break;
            }
            image_alt_allele_global_index++;
          }
          CHECK_LT(image_alt_allele_global_index,
                   dv_call.variant().alternate_bases().size());
          // Find the phase assigned to this read.
          int hp_value = 0;
          if (read.info().contains("HP")) {
            const auto& hp_values = read.info().at("HP").values();
            if (!hp_values.empty()) {
              hp_value = hp_values[0].int_value();
            }
          }
          if (alt_allele_phases[image_alt_allele_global_index] == hp_value
              && hp_value != 0) {
            // If read supports an alt that is close to alt_allele and has the
            // same phase.
            if (std::abs((int)alt_alleles[image_alt_allele_index].size() -
                         (int)alt_allele.size()) == 1) {
              return 10;
            }
            if (std::abs((int)alt_alleles[image_alt_allele_index].size() -
                         (int)alt_allele.size()) == 2) {
              return 9;
            }
          }
        }
        return 2;
        }
      }
    }
  }
  return 0;
}

int ReadSupportsVariantFuzzyChannel::SupportsAltColor(
    int read_supports_alt) const {
  float alpha;
  if (read_supports_alt == 0) {
    alpha = options_.allele_unsupporting_read_alpha();
  } else if (read_supports_alt == 1) {
    alpha = options_.allele_supporting_read_alpha();
  } else if (read_supports_alt == 10) {
    alpha = kReadSupportAltWithinOneBase;
  } else if (read_supports_alt == 9) {
    alpha = kReadSupportAltWithinTwoBases;
  } else if (read_supports_alt == 8) {
    alpha = kReadSupportAltWithinThreeBases;
  } else {
    CHECK_EQ(read_supports_alt, 2) << "read_supports_alt can only be 0/1/2.";
    alpha = options_.other_allele_supporting_read_alpha();
  }
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
