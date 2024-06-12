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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_CHANNEL_H_
#define LEARNING_GENOMICS_DEEPVARIANT_CHANNEL_H_

#include <string>
#include <vector>

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

// Class for channels that only depend on information at the
// granularity of an entire read
class Channel {
 public:
  Channel(int width,
          const learning::genomics::deepvariant::PileupImageOptions& options);
  virtual ~Channel();

  // Fills the provided read_level_data reference. The length of read_level_data
  // is width_, unless all values are the same, in which case we use a vector of
  // size 1 to save memory.
  virtual void FillReadLevelData(
      const Read& read, const DeepVariantCall& dv_call,
      const std::vector<std::string>& alt_alleles,
      std::vector<unsigned char>& read_level_data) = 0;

  // Fills the provided ref_data reference.
  // The length of ref_data_ is always width_.
  virtual void FillRefData(const std::string& ref_bases,
                           std::vector<unsigned char>& ref_data) = 0;

  int width_;

 protected:
  learning::genomics::deepvariant::PileupImageOptions options_;

  // The maximum value a pixel can have as a float. We use the 254.0
  // value as originally set in DeepVariant v1. This means our pixel
  // values can go from 0 to 254. Which, when converted to an int,
  // gives us 255 or 256 possible pixel values.
  static const constexpr float kMaxPixelValueAsFloat = 254;

  // TODO: make this value configurable as a flag
  static const constexpr float kMaxFragmentLength = 1000;
};
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
#endif  // LEARNING_GENOMICS_DEEPVARIANT_CHANNEL_H_
