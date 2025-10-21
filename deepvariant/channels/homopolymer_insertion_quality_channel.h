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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_CHANNELS_HOMOPOLYMER_INSERTION_QUALITY_CHANNEL_H_
#define LEARNING_GENOMICS_DEEPVARIANT_CHANNELS_HOMOPOLYMER_INSERTION_QUALITY_CHANNEL_H_

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using learning::genomics::deepvariant::DeepVariantCall;
using nucleus::genomics::v1::Read;

class HomopolymerInsertionQualityChannel : public Channel {
 public:
  HomopolymerInsertionQualityChannel(
      int width,
      const learning::genomics::deepvariant::PileupImageOptions& options);

  void FillReadBase(std::vector<unsigned char>& data, int col, char read_base,
                    char ref_base, int base_quality, const Read& read,
                    int read_index, const DeepVariantCall& dv_call,
                    const std::vector<std::string>& alt_alleles) override;

  void FillRefBase(std::vector<unsigned char>& ref_data, int col, char ref_base,
                   const std::string& ref_bases) override;

  // Public for testing
  std::vector<std::uint8_t> HomoPolymerInDelQuality(const Read& read,
                                                    bool is_deletion);

 private:
  // Helper functions for reading and processing tags from reads
  std::vector<int8_t> GetTPValues(const Read& read);
  std::vector<std::uint8_t> HomoPolymerWeighted(const Read& read);
  std::uint8_t BaseQualityColor(int base_qual);

  // Scales an input vector to pixel range 0-254
  std::vector<std::uint8_t> ScaleColorVector(
      std::vector<std::uint8_t>& channel_values, float max_val);

  static const constexpr int kMaxQScore = 93;
  static const constexpr int kMaxHomoPolymerWeighted = 30;

  std::optional<std::vector<unsigned char>>
      homopolymer_insertion_quality_vector_;
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_CHANNELS_HOMOPOLYMER_INSERTION_QUALITY_CHANNEL_H_
