/*
 * Copyright 2017 Google LLC.
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
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using std::string;

constexpr int NUM_CHANNELS = 6;
constexpr int NUM_SEQ_TYPE = PileupImageOptions::SequencingType_ARRAYSIZE;

template <class T>
std::vector<T> ToVector(const google::protobuf::RepeatedPtrField<T> container) {
  return std::vector<T>(std::make_move_iterator(container.begin()),
                        std::make_move_iterator(container.end()));
}

struct ImageRow {
  int width;
  int num_channels;
  std::vector<DeepVariantChannelEnum> channel_enums;
  std::vector<std::vector<unsigned char>> channel_data;

  int Width() const;
  explicit ImageRow(int width, int num_channels);
};

class PileupImageEncoderNative {
 public:
  // All channel enums that will be processed by this PileupImageEncoder
  std::vector<DeepVariantChannelEnum> channel_enums_;

  // Essential API methods.
  explicit PileupImageEncoderNative(const PileupImageOptions& options);

  // Return all the channels in this image, as a list of enums.
  std::vector<DeepVariantChannelEnum> AllChannelsEnum(
      const std::string& alt_aligned_representation);

  // Create read pileup image section for one sample.
  std::vector<std::unique_ptr<ImageRow>> BuildPileupForOneSample(
      const DeepVariantCall& dv_call, const string& ref_bases,
      const std::vector<const ::nucleus::genomics::v1::Read *>& reads,
      int image_start_pos, const std::vector<std::string>& alt_alleles,
      const SampleOptions& sample_options);

  // Simple wrapper around BuildPileupForOneSample that allows us to efficiently
  // pass large protobufs in from Python. Simply unwraps the ConstProtoPtr
  // objects and calls BuildPileupForOneSample().
  std::vector<std::unique_ptr<ImageRow>> BuildPileupForOneSamplePython(
      const nucleus::ConstProtoPtr<
          const learning::genomics::deepvariant::DeepVariantCall>&
          wrapped_dv_call,
      const string& ref_bases,
      const std::vector<
          nucleus::ConstProtoPtr<const ::nucleus::genomics::v1::Read>>&
          wrapped_reads,
      int image_start_pos, const std::vector<std::string>& alt_alleles,
      const SampleOptions& sample_options) {
    std::vector<const nucleus::genomics::v1::Read *> reads;
    reads.reserve(wrapped_reads.size());
    for (const auto& wrapped_read : wrapped_reads) {
      reads.emplace_back(wrapped_read.p_);
    }
    return BuildPileupForOneSample(*(wrapped_dv_call.p_), ref_bases,
                                   reads, image_start_pos, alt_alleles,
                                   sample_options);
  }

  // Encode one read into a row of pixels for our image.
  std::unique_ptr<ImageRow> EncodeRead(
      const learning::genomics::deepvariant::DeepVariantCall& dv_call,
      const string& ref_bases, const nucleus::genomics::v1::Read& read,
      int image_start_pos, const std::vector<std::string>& alt_alleles);

  // Simple wrapper around EncodeRead that allows us to efficiently pass large
  // protobufs in from Python. Simply unwraps the ConstProtoPtr objects and
  // calls EncodeRead().
  std::unique_ptr<ImageRow> EncodeReadPython(
      const nucleus::ConstProtoPtr<
          const learning::genomics::deepvariant::DeepVariantCall>&
          wrapped_dv_call,
      const string& ref_bases,
      const nucleus::ConstProtoPtr<const ::nucleus::genomics::v1::Read>&
          wrapped_read,
      int image_start_pos, const std::vector<std::string>& alt_alleles) {
    return EncodeRead(*(wrapped_dv_call.p_), ref_bases, *(wrapped_read.p_),
                      image_start_pos, alt_alleles);
  }

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
  int SupportsAltColor(int read_supports_alt) const;
  // Get the pixel color (int) for a read with an allele frequency.
  int AlleleFrequencyColor(float allele_frequency) const;
  // Get the pixel color (int) for a read that matches ref.
  int MatchesRefColor(bool base_matches_ref) const;
  // Get the pixel color (int) for a base read quality.
  int BaseQualityColor(int base_qual) const;
  // Get the pixel color (int) for a base mapping quality.
  int MappingQualityColor(int mapping_qual) const;

 private:
  int GetHapIndex(const nucleus::genomics::v1::Read& read);

  const PileupImageOptions options_;
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_PILEUP_IMAGE_NATIVE_H_
