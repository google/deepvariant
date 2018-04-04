/*
 * Copyright 2017 Google Inc.
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

#ifndef LEARNING_GENOMICS_DEEPVARIANT_PILEUP_IMAGE_NATIVE_H_
#define LEARNING_GENOMICS_DEEPVARIANT_PILEUP_IMAGE_NATIVE_H_

#include <memory>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "tensorflow/core/platform/types.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using tensorflow::string;

struct ImageRow {
  std::vector<unsigned char> base;
  std::vector<unsigned char> base_quality;
  std::vector<unsigned char> mapping_quality;
  std::vector<unsigned char> on_positive_strand;
  std::vector<unsigned char> supports_alt;
  std::vector<unsigned char> matches_ref;

  int Width() const;
  explicit ImageRow(int width);
};

class PileupImageEncoderNative {
 public:
  // Essential API methods.
  explicit PileupImageEncoderNative(const PileupImageOptions& options);

  // Encode one read into a row of pixels for our image.
  std::unique_ptr<ImageRow> EncodeRead(
      const learning::genomics::deepvariant::DeepVariantCall& dv_call,
      const string& ref_bases, const nucleus::genomics::v1::Read& read,
      int image_start_pos, const std::vector<string>& alt_alleles);

  // Encode the reference bases into a single row of pixels.
  std::unique_ptr<ImageRow> EncodeReference(const string& ref_bases);

 public:
  // Get the pixel color (int) for a base.
  int BaseColor(char base) const;
  // Overload of the above provided for CLIF.
  int BaseColor(const string& base) const;
  // Get the strand pixel color (int) for a positive strand observation.
  int StrandColor(bool on_positive_strand) const;
  // Get the pixel color (int) for a read that supports an alt allele.
  int SupportsAltColor(bool read_supports_alt) const;
  // Get the pixel color (int) for a read that matches ref.
  int MatchesRefColor(bool base_matches_ref) const;
  // Get the pixel color (int) for a base read quality.
  int BaseQualityColor(int base_qual) const;
  // Get the pixel color (int) for a base mapping quality.
  int MappingQualityColor(int mapping_qual) const;

 private:
  const PileupImageOptions options_;
};


}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning



#endif  // LEARNING_GENOMICS_DEEPVARIANT_PILEUP_IMAGE_NATIVE_H_
