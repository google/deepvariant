
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

#include "deepvariant/channels/allele_sample_probability_channel.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning::genomics::deepvariant {

using learning::genomics::deepvariant::DeepVariantCall;
using nucleus::genomics::v1::Read;

void AlleleSampleProbabilityChannel::FillReadBase(
    std::vector<unsigned char>& data, int col, char read_base, char ref_base,
    int base_quality, const Read& read, int read_index,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  int total_reads = 0;
  int total_reads_supporting_allele = 0;
  const std::string read_key =
      (read.fragment_name() + "/" + std::to_string(read.read_number()));

  bool found_allele_group = false;
  for (const auto& [allele, reads] : dv_call.allele_support()) {
    total_reads += reads.read_names_size();
    for (const auto& read_name : reads.read_names()) {
      if (read_name == read_key) {
        total_reads_supporting_allele = reads.read_names_size();
        found_allele_group = true;
        break;
      }
    }
    if (found_allele_group) {
      break;
    }
  }

  if (!found_allele_group) {
    total_reads_supporting_allele = dv_call.ref_support().size();
  }

  total_reads += dv_call.ref_support().size();
  data[col] = ScaleColor(total_reads_supporting_allele, total_reads);
}

void AlleleSampleProbabilityChannel::FillRefBase(
    std::vector<unsigned char>& ref_data, int col, char ref_base,
    const std::string& ref_bases) {
  ref_data[col] = 0;
}

std::uint8_t AlleleSampleProbabilityChannel::ScaleColor(int value,
                                                        float max_val) const {
  if (max_val == 0) {
    return 0;
  }
  float value_as_float = static_cast<float>(value);
  value_as_float = std::clamp<float>(value_as_float, 0.0f, max_val);

  // probability in (0, 1]
  double probability = value_as_float / max_val;
  // We care more about capturing low probabilities, so we take a square root.
  double scaled_probability = std::sqrt(probability);
  // To return zero, the original probability must be below 1/kMaxPixelValue^2.
  return static_cast<int>(kMaxPixelValueAsFloat * scaled_probability);
}
}  // namespace learning::genomics::deepvariant
