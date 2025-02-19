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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_CHANNELS_GAP_COMPRESSED_IDENTITY_CHANNEL_H_
#define LEARNING_GENOMICS_DEEPVARIANT_CHANNELS_GAP_COMPRESSED_IDENTITY_CHANNEL_H_

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/channels/channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using learning::genomics::deepvariant::DeepVariantCall;
using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Read;

class GapCompressedIdentityChannel : public Channel {
 public:
  GapCompressedIdentityChannel(
      int width,
      const learning::genomics::deepvariant::PileupImageOptions& options);

  void FillReadBase(std::vector<unsigned char>& data, int col, char read_base,
                    char ref_base, int base_quality, const Read& read,
                    int read_index, const DeepVariantCall& dv_call,
                    const std::vector<std::string>& alt_alleles) override;

  void FillRefBase(std::vector<unsigned char>& ref_data, int col, char ref_base,
                   const std::string& ref_bases) override;

  // Gap Compressed Identity: Ins/Del treated as individual events.
  // Public for testing.
  int GapCompressedIdentity(const Read& read);

 private:
  // Scales an input value to pixel range 0-254.
  std::uint8_t ScaleColor(int value, float max_val) const;

  std::optional<unsigned char> gap_compressed_identity_color_;

  static const constexpr int kMaxIdentity = 100;
};
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
#endif  // LEARNING_GENOMICS_DEEPVARIANT_CHANNELS_GAP_COMPRESSED_IDENTITY_CHANNEL_H_
