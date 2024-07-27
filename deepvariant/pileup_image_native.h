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
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/string_view.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using std::string;

constexpr int NUM_CHANNELS = 6;
constexpr int NUM_SEQ_TYPE = PileupImageOptions::SequencingType_ARRAYSIZE;
constexpr unsigned char kChannelValue255 = 255;
constexpr unsigned char kChannelValue200 = 200;
// Different ways alt aligned reads can be expressed.
// This enum matches flag values in make_example_options.
// It helps to avoid string comparison in performance critical parts.
enum AltAlignedPileup {
  kNone = 0,
  kBaseChannels = 1,
  kDiffChannels = 2,
  kRows = 3,
};

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
  bool operator==(const ImageRow& other) const {
    if (channel_enums != other.channel_enums) {
      LOG(INFO) << "ImageRow channel_enums mismatch";
    }
    if (channel_data != other.channel_data) {
      LOG(INFO) << "ImageRow channel_data mismatch";
    }
    return width == other.width && num_channels == other.num_channels &&
           channel_enums == other.channel_enums &&
           channel_data == other.channel_data;
  }
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
      const std::vector<const ::nucleus::genomics::v1::Read*>& reads,
      int image_start_pos, const std::vector<std::string>& alt_alleles,
      const SampleOptions& sample_options, float mean_coverage,
      // Contains original alignment positions for trimmed reads. This array has
      // the same order as reads.
      const std::vector<int64_t>* alignment_positions = nullptr,
      absl::flat_hash_set<DeepVariantChannelEnum> channels_enum_to_blank = {});

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
    return BuildPileupForOneSample(*(wrapped_dv_call.p_), ref_bases, reads,
                                   image_start_pos, alt_alleles, sample_options,
                                   /*mean_coverage=*/0.0);
  }

  // Encode one read into a row of pixels for our image.
  std::unique_ptr<ImageRow> EncodeRead(
      const learning::genomics::deepvariant::DeepVariantCall& dv_call,
      absl::string_view ref_bases, const nucleus::genomics::v1::Read& read,
      int image_start_pos, const std::vector<std::string>& alt_alleles,
      absl::flat_hash_set<DeepVariantChannelEnum> channels_enum_to_blank = {});

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

 private:
  int GetHapIndex(const nucleus::genomics::v1::Read& read);

  const PileupImageOptions options_;
};

// Converting a vector<ImageRow> into a 3-d vector where channels are stored in
// the lowest dimension. Special processing is required for alt aligned
// channels.
// * --alt_aligned_pileup=diff_channels Add channel 5 (base differs
// from ref) of both alts as channels
// * --alt_aligned_pileup=base_channels Add channel 0 (bases ATCG) of both alts
// as channels.
// * --alt_aligned_pileup=rows alt aligned data is stored as extra rows.
//
// TODO: Consolidate streaming examples and normal examples so
// that we don't need to use the template anymore.
template <class T>
void FillPileupArray(
    absl::Span<const std::unique_ptr<ImageRow>> image,
    absl::Span<const std::vector<std::unique_ptr<ImageRow>>> alt_image,
    AltAlignedPileup alt_aligned_representation, T* pileup_array,
    int buffer_size, int buffer_pos = 0) {
  // TODO Although channels 0 and 5 should not change we need to set
  // channel number programmatically.
  int alt_channel_index = 0;
  switch (alt_aligned_representation) {
    case AltAlignedPileup::kDiffChannels:
      alt_channel_index = 5;
      break;
    case AltAlignedPileup::kBaseChannels:
      alt_channel_index = 0;
      break;
    default:
      alt_channel_index = 0;
  }
  for (int row = 0; row < image.size(); row++) {
    for (int column = 0; column < image[row]->Width(); column++) {
      if (!image[row]->channel_data.empty()) {
        // Lower dimension is a channel data. Here we iterate all channels to
        // fill one position of the pileup image.
        for (int channel = 0; channel < image[row]->channel_data.size();
             channel++) {
          CHECK_LT(buffer_pos, buffer_size);
          (*pileup_array)[buffer_pos] =
              image[row]->channel_data[channel][column];
          buffer_pos++;
        }
        // Fill alt aligned channels if needed.
        if (alt_aligned_representation == AltAlignedPileup::kBaseChannels ||
            alt_aligned_representation == AltAlignedPileup::kDiffChannels) {
          CHECK_EQ(alt_image.size(), 2);
          unsigned char alt_aligned_channel_1_value = 0;
          if (!alt_image[0].empty() &&
              !alt_image[0][row]->channel_data.empty()) {
            alt_aligned_channel_1_value =
                alt_image[0][row]->channel_data[alt_channel_index][column];
          }
          CHECK_LT(buffer_pos, buffer_size);
          // Fill alt aligned channel 1.
          (*pileup_array)[buffer_pos] = alt_aligned_channel_1_value;
          buffer_pos++;
          // Fill alt aligned channel 2.
          if (alt_image[1].empty() || alt_image[1][row]->channel_data.empty()) {
            // Fill with alt aligned channel 1 if alt2 is empty.
            CHECK_LT(buffer_pos, buffer_size);
            (*pileup_array)[buffer_pos] = alt_aligned_channel_1_value;
            buffer_pos++;
          } else {
            CHECK_LT(buffer_pos, buffer_size);
            (*pileup_array)[buffer_pos] =
                alt_image[1][row]->channel_data[alt_channel_index][column];
            buffer_pos++;
          }
        }  // if need_alt_alignment
      }  // if !channel_data.empty()
    }  // for row->Width
  }  // for row

  // Fill alt aligned channels as rows if AltAlignedPileup::kRows
  if (alt_aligned_representation == AltAlignedPileup::kRows) {
    for (const auto& one_alt_image : alt_image) {
      if (one_alt_image.empty()) {
        CHECK_LT(buffer_pos, buffer_size);
        auto pos_offset_to_end =
            image.size() * image[0]->Width() * image[0]->channel_data.size();
        CHECK_LE(buffer_pos + pos_offset_to_end, buffer_size);
        auto* const fill_from = &(*pileup_array)[buffer_pos];
        std::fill(fill_from, fill_from + pos_offset_to_end,
                  static_cast<unsigned char>(0));
        buffer_pos +=
            image.size() * image[0]->Width() * image[0]->channel_data.size();
        continue;
      }
      for (int row = 0; row < one_alt_image.size(); row++) {
        for (int column = 0; column < one_alt_image[row]->Width(); column++) {
          if (!one_alt_image[row]->channel_data.empty()) {
            // Lower dimension is a channel data. Here we iterate all channels
            // to fill one position of the pileup image.
            for (int channel = 0;
                 channel < one_alt_image[row]->channel_data.size(); channel++) {
              CHECK_LT(buffer_pos, buffer_size);
              (*pileup_array)[buffer_pos] =
                  one_alt_image[row]->channel_data[channel][column];
              buffer_pos++;
            }
          }
        }
      }
    }
  }
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_PILEUP_IMAGE_NATIVE_H_
