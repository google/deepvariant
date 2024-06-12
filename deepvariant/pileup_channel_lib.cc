

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
#include "deepvariant/pileup_channel_lib.h"

#include <math.h>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "deepvariant/allele_frequency_channel.h"
#include "deepvariant/avg_base_quality_channel.h"
#include "deepvariant/blank_channel.h"
#include "deepvariant/channel.h"
#include "deepvariant/gap_compressed_identity_channel.h"
#include "deepvariant/gc_content_channel.h"
#include "deepvariant/haplotype_tag_channel.h"
#include "deepvariant/homopolymer_weighted_channel.h"
#include "deepvariant/identity_channel.h"
#include "deepvariant/insert_size_channel.h"
#include "deepvariant/is_homopolymer_channel.h"
#include "deepvariant/mapping_quality_channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/read_mapping_percent_channel.h"
#include "deepvariant/read_supports_variant_channel.h"
#include "deepvariant/strand_channel.h"
#include "absl/container/flat_hash_set.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Read;
using std::vector;

using learning::genomics::deepvariant::DeepVariantCall;

namespace learning {
namespace genomics {
namespace deepvariant {

bool Channels::CalculateChannels(
    const std::vector<DeepVariantChannelEnum>& channel_enums, const Read& read,
    absl::string_view ref_bases, const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles, int image_start_pos,
    const absl::flat_hash_set<DeepVariantChannelEnum> channels_enum_to_blank) {
  int maxEnumValue = getMaxEnumValue(channel_enums);
  channel_enum_to_index_ = std::vector<int>(maxEnumValue + 1);
  int currIndex = 0;
  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    channel_enum_to_index_[channel_enum] = currIndex;
    currIndex++;
  }

  /*--------------------------------------
  Calculate read-level channels
  ---------------------------------------*/
  data_ = std::vector<std::vector<unsigned char>>(channel_enums.size());
  read_level_data_ =
      std::vector<std::vector<unsigned char>>(channel_enums.size());

  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    int index = channel_enum_to_index_[channel_enum];
    data_[index] = std::vector<unsigned char>(ref_bases.size(), 0);
    if (!isBaseLevelChannel(channel_enum) &&
        !channels_enum_to_blank.contains(channel_enum)) {
      bool ok =
          CalculateReadLevelData(channel_enum, read, dv_call, alt_alleles);
      if (!ok) return false;
    }
  }
  /*--------------------------------------
  Calculate base-level channels
  ---------------------------------------*/

  // Handler for each component of the CIGAR string, as subdivided
  // according the rules below.
  // Side effect: draws in img_row
  // Return value: true on normal exit; false if we determine that we
  // have a low quality base at the call position (in which case we
  // should return null) from EncodeRead.
  std::function<bool(int, int, const CigarUnit::Operation&)>
      action_per_cigar_unit = [&, ref_bases](
                                  int ref_i, int read_i,
                                  const CigarUnit::Operation& cigar_op) {
        char read_base = 0;
        if (cigar_op == CigarUnit::INSERT) {
          read_base = options_.indel_anchoring_base_char()[0];
        } else if (cigar_op == CigarUnit::DELETE) {
          ref_i -= 1;  // Adjust anchor base on reference
          read_base = options_.indel_anchoring_base_char()[0];
        } else if (cigar_op == CigarUnit::ALIGNMENT_MATCH ||
                   cigar_op == CigarUnit::SEQUENCE_MATCH ||
                   cigar_op == CigarUnit::SEQUENCE_MISMATCH) {
          read_base = read.aligned_sequence()[read_i];
        }

        size_t col = ref_i - image_start_pos;
        if (read_base && 0 <= col && col < ref_bases.size()) {
          int base_quality = read.aligned_quality(read_i);
          // Bail out if we found this read had a low-quality base at the
          // call site.
          if (ref_i == dv_call.variant().start() &&
              base_quality < options_.read_requirements().min_base_quality()) {
            return false;
          }

          // Calculate base level values for channels
          for (int i = 0; i < channel_enums.size(); ++i) {
            DeepVariantChannelEnum channel_enum = channel_enums[i];
            int index = channel_enum_to_index_[channel_enum];
            if (isBaseLevelChannel(channel_enum) &&
                !channels_enum_to_blank.contains(channel_enum)) {
              if (channel_enum == DeepVariantChannelEnum::CH_READ_BASE) {
                data_[index][col] = BaseColor(read_base, options_);
              } else if (channel_enum ==
                         DeepVariantChannelEnum::CH_BASE_QUALITY) {
                data_[index][col] =
                    ScaleColor(base_quality, options_.base_quality_cap());
              } else if (channel_enum ==
                         DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF) {
                bool matches_ref = (read_base == ref_bases[col]);
                data_[index][col] = MatchesRefColor(matches_ref, options_);
              }
            }
          }

          // Fill in base level value for read level channels from
          // previously calculated read level values
          for (int i = 0; i < channel_enums.size(); ++i) {
            DeepVariantChannelEnum channel_enum = channel_enums[i];
            int index = channel_enum_to_index_[channel_enum];
            if (!isBaseLevelChannel(channel_enum) &&
                !channels_enum_to_blank.contains(channel_enum)) {
              // Channels that assign the same value for all columns use a
              // size 1 vector.
              if (read_level_data_[index].size() == 1) {
                data_[index][col] = read_level_data_[index][0];
              } else {
                data_[index][col] = read_level_data_[index][read_i];
              }
            }
          }
        }
        return true;
      };
  return CalculateBaseLevelData(read, action_per_cigar_unit);
}

// Calculate values for channels that only depend on information at the
// granularity of an entire read.
bool Channels::CalculateReadLevelData(
    const DeepVariantChannelEnum channel_enum, const Read& read,
    const DeepVariantCall& dv_call,
    const std::vector<std::string>& alt_alleles) {
  int index = channel_enum_to_index_[channel_enum];

  const int mapping_quality = read.alignment().mapping_quality();
  const int min_mapping_quality =
      options_.read_requirements().min_mapping_quality();
  // Bail early if this read's mapping quality is too low.
  if (mapping_quality < min_mapping_quality) {
    return false;
  }

  std::unique_ptr<Channel> readLevelChannel = Channels::ChannelEnumToObject(
      channel_enum, read.aligned_sequence().size(), options_);
  readLevelChannel->FillReadLevelData(read, dv_call, alt_alleles,
                                      read_level_data_[index]);
  return true;
}

// Calculate values for channels that depend on information at the
// granularity of bases within the read and/or reference sequence.
bool Channels::CalculateBaseLevelData(
    const Read& read, std::function<bool(int, int, const CigarUnit::Operation&)>
                          action_per_cigar_unit) {
  // In the following, we iterate over alignment information for each
  // base of read creating an association between a reference index,
  // read index and cigar operation and storing in a vector for subsequent
  // base level channel value calculation(s).
  //
  // The handling of each cigar element type is given below, assuming
  // it has length n.
  //
  // ALIGNMENT_MATCH, SEQUENCE_MATCH, SEQUENCE_MISMATCH:
  //   Provide a segment ref_i, read_i for each of the n bases in the
  //   operator, where ref_i is the position on the genome where this
  //   base aligns.
  //
  // INSERT, CLIP_SOFT:
  //   Provides a single ref_i, read_i segment regardless of n. ref_i
  //   is set to the preceding base of the insertion; i.e., the anchor
  //   base. Beware that ref_i could be -1 if the insertion is aligned
  //   to the first base of a contig.  read_i points to the first base
  //   of the insertion. So if our cigar is 1M2I1M for a read starting
  //   at S, we'd see first (S, 0, '1M'), followed by one (S, 1,
  //   '2I'), and then (S + 1, 3, '1M').
  //
  // DELETE, SKIP:
  //   Provides a single ref_i, read_i segment regardless of n. ref_i
  //   is set to the first base of the deletion, just like in an
  //   ALIGNMENT_MATCH. read_i points to the previous base in the
  //   read, as there's no actual read sequence associated with a
  //   deletion. Beware that read_i could be -1 if the deletion is the
  //   first cigar of the read. So if our cigar is 1M2D1M for a read
  //   starting at S, we'd see first (S, 0, '1M'), followed by one (S
  //   + 1, 0, '2D'), and then (S + 3, 1, '1M').
  //
  // CLIP_HARD, PAD:
  //   These operators are ignored by as they don't impact the
  //   alignment of the read w.r.t. the reference.
  //
  // Any other CIGAR op:
  //   Fatal error, at present; later we should fail with a status encoding.

  int ref_i = read.alignment().position().position();
  int read_i = 0;
  bool ok = true;
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();

    switch (op) {
      case CigarUnit::ALIGNMENT_MATCH:
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::SEQUENCE_MISMATCH:
        // Alignment op.
        for (int i = 0; i < op_len; i++) {
          ok = ok && action_per_cigar_unit(ref_i, read_i, op);
          ref_i++;
          read_i++;
        }
        break;
      case CigarUnit::INSERT:
      case CigarUnit::CLIP_SOFT:
        // Insert op.
        ok = action_per_cigar_unit(ref_i - 1, read_i, op);
        read_i += op_len;
        break;
      case CigarUnit::DELETE:
      case CigarUnit::SKIP:
        // Delete op.
        ok = action_per_cigar_unit(ref_i, read_i - 1, op);
        ref_i += op_len;
        break;
      case CigarUnit::CLIP_HARD:
      case CigarUnit::PAD:
        // Ignored ops.  Do nothing.
        break;
      default:
        LOG(FATAL) << "Unrecognized CIGAR op";
    }

    if (!ok) {
      return false;
    }
  }

  return true;
}

std::uint8_t Channels::GetChannelData(const std::string& channel, int col) {
  DeepVariantChannelEnum channel_enum = ChannelStrToEnum(channel);
  int index = channel_enum_to_index_[channel_enum];
  return data_[index][col];
}

void Channels::CalculateRefRows(
    const std::vector<DeepVariantChannelEnum>& channel_enums,
    const std::string& ref_bases) {
  int maxEnumValue = getMaxEnumValue(channel_enums);
  channel_enum_to_index_ = std::vector<int>(maxEnumValue + 1);
  int currIndex = 0;
  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    channel_enum_to_index_[channel_enum] = currIndex;
    currIndex++;
  }
  ref_data_ = std::vector<std::vector<unsigned char>>(channel_enums.size());

  for (const DeepVariantChannelEnum channel_enum : channel_enums) {
    int index = channel_enum_to_index_[channel_enum];

    std::unique_ptr<Channel> readLevelChannel =
        Channels::ChannelEnumToObject(channel_enum, ref_bases.size(), options_);
    if (readLevelChannel != nullptr) {
      readLevelChannel->FillRefData(ref_bases, ref_data_[index]);
      // Base-level channels
    } else if (channel_enum == DeepVariantChannelEnum::CH_READ_BASE) {
      ref_data_[index] = BaseColorVector(ref_bases, options_);
    } else if (channel_enum == DeepVariantChannelEnum::CH_BASE_QUALITY) {
      int ref_qual = options_.reference_base_quality();
      ref_data_[index] = std::vector<unsigned char>(
          ref_bases.size(), ScaleColor(ref_qual, options_.base_quality_cap()));
    } else if (channel_enum ==
               DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF) {
      int ref = MatchesRefColor(true, options_);
      ref_data_[index] = std::vector<unsigned char>(
          ref_bases.size(), static_cast<std::uint8_t>(ref));
    } else {
      ref_data_[index] = std::vector<unsigned char>(ref_bases.size(), 0);
    }
  }
}

std::uint8_t Channels::GetRefRows(DeepVariantChannelEnum channel_enum,
                                  int col) {
  // Returns first value if size 1 else return specific column.
  // Note that ref_data is indexed by col and not pos.
  int index = channel_enum_to_index_[channel_enum];
  return ref_data_[index][col];
}

int Channels::getMaxEnumValue(
    const std::vector<DeepVariantChannelEnum>& channel_enums) {
  // Get the max channel enum value we need to support, used to size vectors
  // returns -1 if channel_enums is empty
  int maxEnumValue = -1;
  for (int i = 0; i < channel_enums.size(); ++i) {
    if (channel_enums[i] > maxEnumValue) {
      maxEnumValue = channel_enums[i];
    }
  }
  return maxEnumValue;
}
// Scales an input value to pixel range 0-254.
std::uint8_t Channels::ScaleColor(int value, float max_val) {
  if (static_cast<float>(value) > max_val) {
    value = max_val;
  }
  return static_cast<int>(kMaxPixelValueAsFloat *
                          (static_cast<float>(value) / max_val));
}

// Scales an input vector to pixel range 0-254
std::vector<std::uint8_t> Channels::ScaleColorVector(
    std::vector<std::uint8_t>& channel_values, float max_val) {
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
//---------------//
// Base Channels //
//---------------//
int Channels::BaseColor(char base, const PileupImageOptions& options) {
  switch (base) {
    case 'A':
      return (options.base_color_offset_a_and_g() +
              options.base_color_stride() * 3);
    case 'G':
      return (options.base_color_offset_a_and_g() +
              options.base_color_stride() * 2);
    case 'T':
      return (options.base_color_offset_t_and_c() +
              options.base_color_stride() * 1);
    case 'C':
      return (options.base_color_offset_t_and_c() +
              options.base_color_stride() * 0);
    default:
      return 0;
  }
}

std::vector<std::uint8_t> Channels::BaseColorVector(
    const std::string& bases, const PileupImageOptions& options) {
  std::vector<std::uint8_t> base_colors;
  base_colors.reserve(bases.size());
  for (const char base : bases) {
    int color = BaseColor(base, options);
    base_colors.push_back(color);
  }
  return base_colors;
}

// Returns a value based on whether the current read base matched the
// reference base it was compared to.
int Channels::MatchesRefColor(bool base_matches_ref,
                              const PileupImageOptions& options) {
  float alpha = (base_matches_ref ? options.reference_matching_read_alpha()
                                  : options.reference_mismatching_read_alpha());
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

//-----------------------//
// Experimental Channels //
//-----------------------//
constexpr bool Channels::isBaseLevelChannel(
    DeepVariantChannelEnum channelEnum) {
  switch (channelEnum) {
    case DeepVariantChannelEnum::CH_READ_BASE:
    case DeepVariantChannelEnum::CH_BASE_QUALITY:
    case DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF:
      return true;
    default:
      return false;
  }
}

// Given a channel enum for a read-level channel instantiates a corresponding
// channel object, containing the methods to populate that channel.
std::unique_ptr<Channel> Channels::ChannelEnumToObject(
    DeepVariantChannelEnum channel_enum, int width,
    const learning::genomics::deepvariant::PileupImageOptions& options) {
  if (isBaseLevelChannel(channel_enum)) {
    return nullptr;
  }
  switch (channel_enum) {
    case DeepVariantChannelEnum::CH_MAPPING_QUALITY:
      return std::unique_ptr<Channel>(
          new MappingQualityChannel(width, options));
    case DeepVariantChannelEnum::CH_STRAND:
      return std::unique_ptr<Channel>(new StrandChannel(width, options));
    case DeepVariantChannelEnum::CH_READ_SUPPORTS_VARIANT:
      return std::unique_ptr<Channel>(
          new ReadSupportsVariantChannel(width, options));
    case DeepVariantChannelEnum::CH_READ_MAPPING_PERCENT:
      return std::unique_ptr<Channel>(
          new ReadMappingPercentChannel(width, options));
    case DeepVariantChannelEnum::CH_HAPLOTYPE_TAG:
      return std::unique_ptr<Channel>(new HaplotypeTagChannel(width, options));
    case DeepVariantChannelEnum::CH_ALLELE_FREQUENCY:
      return std::unique_ptr<Channel>(
          new AlleleFrequencyChannel(width, options));
    case DeepVariantChannelEnum::CH_AVG_BASE_QUALITY:
      return std::unique_ptr<Channel>(
          new AvgBaseQualityChannel(width, options));
    case DeepVariantChannelEnum::CH_IDENTITY:
      return std::unique_ptr<Channel>(new IdentityChannel(width, options));
    case DeepVariantChannelEnum::CH_GAP_COMPRESSED_IDENTITY:
      return std::unique_ptr<Channel>(
          new GapCompressedIdentityChannel(width, options));
    case DeepVariantChannelEnum::CH_GC_CONTENT:
      return std::unique_ptr<Channel>(new GcContentChannel(width, options));
    case DeepVariantChannelEnum::CH_IS_HOMOPOLYMER:
      return std::unique_ptr<Channel>(new IsHomopolymerChannel(width, options));
    case DeepVariantChannelEnum::CH_HOMOPOLYMER_WEIGHTED:
      return std::unique_ptr<Channel>(
          new HomopolymerWeightedChannel(width, options));
    case DeepVariantChannelEnum::CH_BLANK:
      return std::unique_ptr<Channel>(new BlankChannel(width, options));
    case DeepVariantChannelEnum::CH_INSERT_SIZE:
      return std::unique_ptr<Channel>(new InsertSizeChannel(width, options));
    default:
      return nullptr;
  }
}

DeepVariantChannelEnum Channels::ChannelStrToEnum(const std::string& channel) {
  if (channel == ch_read_base) return DeepVariantChannelEnum::CH_READ_BASE;
  if (channel == ch_base_quality)
    return DeepVariantChannelEnum::CH_BASE_QUALITY;
  if (channel == ch_mapping_quality)
    return DeepVariantChannelEnum::CH_MAPPING_QUALITY;
  if (channel == ch_strand) return DeepVariantChannelEnum::CH_STRAND;
  if (channel == ch_read_supports_variant)
    return DeepVariantChannelEnum::CH_READ_SUPPORTS_VARIANT;
  if (channel == ch_base_differs_from_ref)
    return DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF;
  if (channel == ch_read_mapping_percent)
    return DeepVariantChannelEnum::CH_READ_MAPPING_PERCENT;
  if (channel == ch_haplotype_tag)
    return DeepVariantChannelEnum::CH_HAPLOTYPE_TAG;
  if (channel == ch_allele_frequency)
    return DeepVariantChannelEnum::CH_ALLELE_FREQUENCY;
  if (channel == ch_diff_channels_alternate_allele_1)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  if (channel == ch_diff_channels_alternate_allele_2)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  if (channel == ch_avg_base_quality)
    return DeepVariantChannelEnum::CH_AVG_BASE_QUALITY;
  if (channel == ch_identity) return DeepVariantChannelEnum::CH_IDENTITY;
  if (channel == ch_gap_compressed_identity)
    return DeepVariantChannelEnum::CH_GAP_COMPRESSED_IDENTITY;
  if (channel == ch_gc_content) return DeepVariantChannelEnum::CH_GC_CONTENT;
  if (channel == ch_is_homopolymer)
    return DeepVariantChannelEnum::CH_IS_HOMOPOLYMER;
  if (channel == ch_homopolymer_weighted)
    return DeepVariantChannelEnum::CH_HOMOPOLYMER_WEIGHTED;
  if (channel == ch_blank) return DeepVariantChannelEnum::CH_BLANK;
  if (channel == ch_insert_size) return DeepVariantChannelEnum::CH_INSERT_SIZE;
  if (channel == ch_base_channels_alternate_allele_1)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  if (channel == ch_base_channels_alternate_allele_2)
    return DeepVariantChannelEnum::CH_UNSPECIFIED;
  CHECK(false) << "Channel '" << channel << "' should have a corresponding "
               << "enum in DeepVariantChannelEnum.";
}

bool Channels::channel_exists(std::vector<std::string>& channels,
                              absl::string_view channel_name) {
  if (std::find(channels.begin(), channels.end(), channel_name) !=
      channels.end()) {
    return true;
  }
  return false;
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
