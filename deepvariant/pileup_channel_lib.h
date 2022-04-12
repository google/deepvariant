/*
 * Copyright 2021 Google LLC.
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
#ifndef LEARNING_GENOMICS_DEEPVARIANT_PILEUP_CHANNEL_LIB_H_
#define LEARNING_GENOMICS_DEEPVARIANT_PILEUP_CHANNEL_LIB_H_

#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "tensorflow/core/platform/logging.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Read;
using tensorflow::uint8;

//--------------//
// Opt Channels //
//--------------//

static const auto& ch_read_mapping_percent = "read_mapping_percent";
static const auto& ch_avg_base_quality = "avg_base_quality";
static const auto& ch_identity = "identity";
static const auto& ch_gap_compressed_identity = "gap_compressed_identity";
static const auto& ch_gc_content = "gc_content";
static const auto& ch_is_homopolymer = "is_homopolymer";
static const auto& ch_homopolymer_weighted = "homopolymer_weighted";
static const auto& ch_blank = "blank";
static const auto& ch_insert_size = "insert_size";

//-------//
// Utils //
//-------//

// The maximum value a pixel can have as a float. We use the 254.0
// value as originally set in DeepVariant v1. This means our pixel
// values can go from 0 to 254. Which, when converted to an int,
// gives us 255 or 256 possible pixel values.
const float kMaxPixelValueAsFloat = 254.0;
// The maximum value that we will consider for fragment length.
// TODO: make this value configurable as a flag
const float MaxFragmentLength = 1000;

// Scales an input value to pixel range 0-254.
inline uint8 ScaleColor(int value, float max_val) {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}

// Scales an input vector to pixel range 0-254
inline std::vector<uint8> ScaleColorVector(std::vector<uint8>& channel_values,
                                           float max_val) {
  for (int i = 0; i < channel_values.size(); i++) {
    int value = channel_values[i];
    if (static_cast<float>(value) > max_val) {
      value = max_val;
    }
    channel_values[i] = static_cast<int>(kMaxPixelValueAsFloat *
                                         (static_cast<float>(value) / max_val));
  }
  return channel_values;
}

//-----------------------//
// Experimental Channels //
//-----------------------//

// Read Mapping Percent: Calculates percentage of bases mapped to reference.
inline int ReadMappingPercent(const Read& read) {
  int match_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        break;
      default:
        break;
    }
  }
  float mapping_percent = (static_cast<float>(match_len) /
                           static_cast<float>(read.aligned_sequence().size())) *
                          100;
  return static_cast<int>(mapping_percent);
}

// Average Base Quality: Averages base quality over length of read.
inline int AvgBaseQuality(const Read& read) {
  int base_qual_sum = 0;
  for (const auto& base_qual : read.aligned_quality()) {
    base_qual_sum += base_qual;
    // Base qualities range between 0 and 93
    if (base_qual < 0 || base_qual > 93) {
      LOG(FATAL) << "Encountered base quality outside of bounds (0,93):"
                 << base_qual << ", read=" << read.DebugString();
    }
  }
  float avg_base_qual = (static_cast<float>(base_qual_sum) /
                         static_cast<float>(read.aligned_quality().size()));
  return static_cast<int>(avg_base_qual);
}

// Identity: Similar to mapping percent but with a slightly different def.
inline int Identity(const Read& read) {
  int match_len = 0;
  int total_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        total_len += op_len;
        break;
      case CigarUnit::SEQUENCE_MISMATCH:
        total_len += op_len;
        break;
      case CigarUnit::INSERT:
        total_len += op_len;
        break;
      case CigarUnit::DELETE:
        total_len += op_len;
        break;
      default:
        break;
    }
  }
  float mapping_percent = (static_cast<float>(match_len) /
                           static_cast<float>(read.aligned_sequence().size())) *
                          100;
  return static_cast<int>(mapping_percent);
}

// Gap Compressed Identity: Ins/Del treated as individual events.
inline int GapCompressedIdentity(const Read& read) {
  // Calculates percentage of the read mapped to the reference
  int match_len = 0;
  int gap_compressed_len = 0;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        gap_compressed_len += op_len;
        break;
      case CigarUnit::SEQUENCE_MISMATCH:
        gap_compressed_len += op_len;
        break;
      case CigarUnit::INSERT:
        // Add a single event for insertion.
        gap_compressed_len += 1;
        break;
      case CigarUnit::DELETE:
        // Add a single event for a deletion.
        gap_compressed_len += 1;
        break;
      default:
        break;
    }
  }
  float gap_compressed_identity = static_cast<float>(match_len) /
                                  static_cast<float>(gap_compressed_len) * 100;
  return static_cast<int>(gap_compressed_identity);
}

inline int GcContent(const Read& read) {
  int gc_count{};

  for (const auto& base : read.aligned_sequence()) {
    if (base == 'G' || base == 'C') {
      gc_count += 1;
    }
  }

  return static_cast<int>((static_cast<float>(gc_count) /
                           static_cast<float>(read.aligned_sequence().size())) *
                          100);
}

inline std::vector<uint8> IsHomoPolymer(const Read& read) {
  // Generates a vector indicating homopolymers of 3 or more.
  // ATCGGGAG
  // 00011100
  std::vector<uint8> homopolymer(read.aligned_sequence().size());
  auto seq = read.aligned_sequence();
  for (int i = 2; i < seq.size(); i++) {
    if (seq[i] == seq[i - 1] && seq[i - 1] == seq[i - 2]) {
      homopolymer[i] = 1;
      homopolymer[i - 1] = 1;
      homopolymer[i - 2] = 1;
    }
  }
  return homopolymer;
}

inline std::vector<uint8> HomoPolymerWeighted(const Read& read) {
  // Generates a vector reflecting the number of repeats observed.
  // ATCGGGAA
  // 11133322
  std::vector<uint8> homopolymer(read.aligned_sequence().size());
  auto seq = read.aligned_sequence();
  homopolymer[0] = 1;
  int current_weight = 1;
  for (int i = 1; i <= seq.size(); i++) {
    if (seq[i] == seq[i - 1]) {
      current_weight += 1;
    } else {
      for (int cw = current_weight; cw >= 1; cw--) {
        homopolymer[i - cw] = current_weight;
      }
      current_weight = 1;
    }
  }
  return homopolymer;
}

inline std::vector<uint8> Blank(const Read& read) {
  // Used to return a blank channel.
  std::vector<uint8> blank(read.aligned_sequence().size(), 0);
  return blank;
}

inline bool channel_exists(std::vector<std::string>& channels,
                           const std::string& channel_name) {
  if (std::find(channels.begin(), channels.end(), channel_name) !=
      channels.end()) {
    return true;
  }
  return false;
}

// normalizes a Read's `fragment_length` to a pixel value
inline int normalizeFragmentLength(const Read& read) {
  int fragment_length = std::abs(read.fragment_length());
  if (static_cast<float>(fragment_length) > MaxFragmentLength) {
    fragment_length = static_cast<int>(MaxFragmentLength);
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
    (static_cast<float>(fragment_length) / MaxFragmentLength));
}

inline std::vector<uint8> ReadInsertSize(
    const Read& read) {
  // Generates a vector reflecting the fragment length of the read
  std::vector<uint8> reads_with_insert_size(
      read.aligned_sequence().size(), normalizeFragmentLength(read));
  return reads_with_insert_size;
}


inline DeepVariantChannelEnum ChannelStrToEnum(const std::string& channel) {
  if (channel == ch_read_mapping_percent)
    return DeepVariantChannelEnum::CH_READ_MAPPING_PERCENT;
  if (channel == ch_avg_base_quality)
    return DeepVariantChannelEnum::CH_AVG_BASE_QUALITY;
  if (channel == ch_identity) return DeepVariantChannelEnum::CH_IDENTITY;
  if (channel == ch_gap_compressed_identity)
    return DeepVariantChannelEnum::CH_GAP_COMPRESSED_IDENTITY;
  if (channel == ch_gc_content) return DeepVariantChannelEnum::CH_GC_CONTENT;
  if (channel == ch_is_homopolymer)
    return DeepVariantChannelEnum::CH_IS_HOMEOPOLYMER;
  if (channel == ch_homopolymer_weighted)
    return DeepVariantChannelEnum::CH_HOMEOPOLYMER_WEIGHTED;
  if (channel == ch_blank) return DeepVariantChannelEnum::CH_BLANK;
  if (channel == ch_insert_size)
    return DeepVariantChannelEnum::CH_INSERT_SIZE;
  CHECK(false) << "Channel '" << channel << "' should have a corresponding "
      << "enum in DeepVariantChannelEnum.";
}


//-------------------//
// Channels Accessor //
//-------------------//

// Max values for scaling
const int MaxMappingPercent = 100;
const int MaxAvgBaseQuality = 93;
const int MaxIdentity = 100;
const int MaxGapCompressedIdentity = 100;
const int MaxGcContent = 100;
const int MaxIsHomoPolymer = 1;
const int MaxHomoPolymerWeighted = 30;

class OptChannels {
 public:
  std::map<std::string, std::vector<unsigned char>> data_;
  std::map<std::string, std::vector<unsigned char>> ref_data_;
  void CalculateChannels(const std::vector<std::string>& channels,
                         const Read& read) {
    // Calculates values for each channel
    for (const std::string& channel : channels) {
      if (channel == ch_read_mapping_percent) {
        data_[channel].assign(
            {ScaleColor(ReadMappingPercent(read), MaxMappingPercent)});
      } else if (channel == ch_avg_base_quality) {
        data_[channel].assign(
            {ScaleColor(AvgBaseQuality(read), MaxAvgBaseQuality)});
      } else if (channel == ch_identity) {
        data_[channel].assign({ScaleColor(Identity(read), MaxIdentity)});
      } else if (channel == ch_gap_compressed_identity) {
        data_[channel].assign(
            {ScaleColor(GapCompressedIdentity(read), MaxIdentity)});
      } else if (channel == ch_gc_content) {
        data_[channel].assign({ScaleColor(GcContent(read), MaxGcContent)});
      } else if (channel == ch_is_homopolymer) {
        std::vector<uint8> is_homopolymer = IsHomoPolymer(read);
        data_[channel] = ScaleColorVector(is_homopolymer, MaxIsHomoPolymer);
      } else if (channel == ch_homopolymer_weighted) {
        std::vector<uint8> homopolymer_weighted = HomoPolymerWeighted(read);
        data_[channel] =
            ScaleColorVector(homopolymer_weighted, MaxHomoPolymerWeighted);
      } else if (channel == ch_blank) {
        data_[channel] = Blank(read);
      } else if (channel == ch_insert_size) {
        data_[channel] = ReadInsertSize(read);
      }
    }
  }

  inline uint8 GetChannelData(const std::string& channel, int pos) {
    // Returns values for each channel
    if (data_[channel].size() == 1) {
      return data_[channel][0];
    } else {
      return data_[channel][pos];
    }
  }

  void CalculateRefRows(const std::vector<std::string>& channels,
                        const std::string& ref_bases) {
    // Calculates reference row values for each channel
    // Create a fake read to represent reference bases.
    Read refRead;
    for (const std::string& channel : channels) {
      if (channel == ch_read_mapping_percent) {
        ref_data_[channel].assign({static_cast<uint8>(kMaxPixelValueAsFloat)});
      } else if (channel == ch_avg_base_quality) {
        ref_data_[channel].assign({static_cast<uint8>(kMaxPixelValueAsFloat)});
      } else if (channel == ch_identity) {
        ref_data_[channel].assign({static_cast<uint8>(kMaxPixelValueAsFloat)});
      } else if (channel == ch_gap_compressed_identity) {
        ref_data_[channel].assign({static_cast<uint8>(kMaxPixelValueAsFloat)});
      } else if (channel == ch_insert_size) {
        ref_data_[channel].assign({static_cast<uint8>(kMaxPixelValueAsFloat)});
      } else if (channel == ch_gc_content) {
        refRead.set_aligned_sequence(ref_bases);
        ref_data_[channel].assign(
            {ScaleColor(GcContent(refRead), MaxGcContent)});
      } else if (channel == ch_is_homopolymer) {
        refRead.set_aligned_sequence(ref_bases);
        std::vector<uint8> is_homopolymer = IsHomoPolymer(refRead);
        ref_data_[channel] = ScaleColorVector(is_homopolymer, MaxIsHomoPolymer);
      } else if (channel == ch_homopolymer_weighted) {
        refRead.set_aligned_sequence(ref_bases);
        std::vector<uint8> homopolymer_weighted = HomoPolymerWeighted(refRead);
        ref_data_[channel] =
            ScaleColorVector(homopolymer_weighted, MaxHomoPolymerWeighted);
      } else {
        ref_data_[channel].assign({0});
      }
    }
  }

  inline uint8 GetRefRows(const std::string& channel, int col) {
    // Returns first value if size 1 else return specific column.
    // Note that ref_data is indexed by col and not pos.
    if (ref_data_[channel].size() == 1) {
      return ref_data_[channel][0];
    } else {
      return ref_data_[channel][col];
    }
  }
};

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_PILEUP_CHANNEL_LIB_H_
