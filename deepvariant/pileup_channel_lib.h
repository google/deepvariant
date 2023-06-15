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

#include <algorithm>
#include <cstdlib>
#include <functional>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/container/btree_set.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "tensorflow/core/platform/logging.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::genomics::v1::CigarUnit;
using nucleus::genomics::v1::Read;
using tensorflow::uint8;

//------------------------//
// Default Channels Names //
//------------------------//
static const auto& ch_read_base = "read_base";
static const auto& ch_base_quality = "base_quality";
static const auto& ch_mapping_quality = "mapping_quality";
static const auto& ch_strand = "strand";
static const auto& ch_read_supports_variant = "read_supports_variant";
static const auto& ch_base_differs_from_ref = "base_differs_from_ref";

//--------------------//
// Opt Channels Names //
//--------------------//

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

//---------------//
// Base Channels //
//---------------//

inline int BaseColor(char base, const PileupImageOptions& options) {
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

inline std::vector<uint8> BaseColorVector(const std::string& bases,
                                          const PileupImageOptions& options) {
  std::vector<uint8> base_colors;
  base_colors.reserve(bases.size());
  for (const char base : bases) {
    int color = BaseColor(base, options);
    base_colors.push_back(color);
  }
  return base_colors;
}

inline int StrandColor(bool on_positive_strand,
                       const PileupImageOptions& options) {
  return (on_positive_strand ? options.positive_strand_color()
                             : options.negative_strand_color());
}

inline int SupportsAltColor(int read_supports_alt,
                            const PileupImageOptions& options) {
  float alpha;
  if (read_supports_alt == 0) {
    alpha = options.allele_unsupporting_read_alpha();
  } else if (read_supports_alt == 1) {
    alpha = options.allele_supporting_read_alpha();
  } else {
    CHECK_EQ(read_supports_alt, 2) << "read_supports_alt can only be 0/1/2.";
    alpha = options.other_allele_supporting_read_alpha();
  }
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
}

// Does this read support ref, one of the alternative alleles, or an allele we
// aren't considering?
inline int ReadSupportsAlt(const DeepVariantCall& dv_call, const Read& read,
                           const std::vector<std::string>& alt_alleles) {
  std::string key =
      (read.fragment_name() + "/" + std::to_string(read.read_number()));

  // Iterate over all alts, not just alt_alleles.
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
          return 2;
        }
      }
    }
  }
  return 0;
}

// Returns a value based on whether the current read base matched the
// reference base it was compared to.
inline int MatchesRefColor(bool base_matches_ref,
                           const PileupImageOptions& options) {
  float alpha = (base_matches_ref ? options.reference_matching_read_alpha()
                                  : options.reference_mismatching_read_alpha());
  return static_cast<int>(kMaxPixelValueAsFloat * alpha);
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
  for (const auto& cigar_elt : read.alignment().cigar()) {
    const CigarUnit::Operation& op = cigar_elt.operation();
    int op_len = cigar_elt.operation_length();
    switch (op) {
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::ALIGNMENT_MATCH:
        match_len += op_len;
        break;
      case CigarUnit::SEQUENCE_MISMATCH:
        break;
      case CigarUnit::INSERT:
        break;
      case CigarUnit::DELETE:
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
  if (channel == ch_read_base) return DeepVariantChannelEnum::CH_READ_BASE;
  if (channel == ch_base_quality)
    return DeepVariantChannelEnum::CH_BASE_QUALITY;
  if (channel == ch_mapping_quality)
    return DeepVariantChannelEnum::CH_MAPPING_QUALITY;
  if (channel == ch_strand)
    return DeepVariantChannelEnum::CH_STRAND;
  if (channel == ch_read_supports_variant)
    return DeepVariantChannelEnum::CH_READ_SUPPORTS_VARIANT;
  if (channel == ch_base_differs_from_ref)
    return DeepVariantChannelEnum::CH_BASE_DIFFERS_FROM_REF;
  if (channel == ch_read_mapping_percent)
    return DeepVariantChannelEnum::CH_READ_MAPPING_PERCENT;
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
  const PileupImageOptions& options_;
  std::map<std::string, std::vector<unsigned char>> ref_data_;

  std::map<std::string, std::vector<unsigned char>> read_level_data_;
  std::map<std::string, std::vector<unsigned char>> data_;
  const std::set<std::string> base_level_channels_set_ = {
    ch_read_base,
    ch_base_quality,
    ch_base_differs_from_ref,
  };

  bool CalculateChannels(const std::vector<std::string>& channels,
                         const Read& read, const std::string& ref_bases,
                         const DeepVariantCall& dv_call,
                         const std::vector<std::string>& alt_alleles,
                         int image_start_pos) {
    absl::btree_set<std::string> included_base_level_channels;

    /*--------------------------------------
    Calculate read-level channels
    ---------------------------------------*/
    for (const std::string& channel : channels) {
      // Instantiate each channel data row
      data_[channel] = std::vector<unsigned char>(ref_bases.size(), 0);

      // If we are looking at a base level channel we will fill that data out
      // later. For read level channels we can calculate values now
      if (base_level_channels_set_.find(channel) !=
                  base_level_channels_set_.end()) {
        included_base_level_channels.insert(channel);
      } else {
        bool ok = CalculateReadLevelData(channel, read, dv_call, alt_alleles);
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
        action_per_cigar_unit =
          [&](int ref_i, int read_i, const CigarUnit::Operation& cigar_op) {
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
                  base_quality <
                    options_.read_requirements().min_base_quality()) {
                return false;
              }

              // Calculate base level values for channels
              for (const std::string& channel : included_base_level_channels)
                {
                if (channel == ch_read_base) {
                  data_[channel][col] = BaseColor(read_base, options_);
                } else if (channel == ch_base_quality) {
                  data_[channel][col] =
                      ScaleColor(base_quality, options_.base_quality_cap());
                } else if (channel == ch_base_differs_from_ref) {
                  bool matches_ref = (read_base == ref_bases[col]);
                  data_[channel][col] =
                      MatchesRefColor(matches_ref, options_);
                }
              }

              // Fill in base level value for read level channels from
              // previously calculated read level values
              for (
                std::map<std::string, std::vector<unsigned char>>::iterator
                          iter = read_level_data_.begin();
                    iter != read_level_data_.end(); ++iter)
              {
                std::string channel =  iter->first;
                if (iter->second.size() == 1) {
                  data_[channel][col] = iter->second[0];
                } else {
                  data_[channel][col] = iter->second[read_i];
                }
                }
            }
            return true;
          };

    return CalculateBaseLevelData(read, action_per_cigar_unit);
  }

  // Calculate values for channels that only depend on information at the
  // granularity of an entire read.
  bool CalculateReadLevelData(const std::string& channel, const Read& read,
                              const DeepVariantCall& dv_call,
                              const std::vector<std::string>& alt_alleles) {
    if (channel == ch_mapping_quality) {
      const int mapping_quality = read.alignment().mapping_quality();
      const int min_mapping_quality =
          options_.read_requirements().min_mapping_quality();
      // Bail early if this read's mapping quality is too low.
      if (mapping_quality < min_mapping_quality) {
        return false;
      }
      read_level_data_[channel].assign(
          {ScaleColor(mapping_quality, options_.mapping_quality_cap())});
    } else if (channel == ch_strand) {
      const bool is_forward_strand =
          !read.alignment().position().reverse_strand();
      read_level_data_[channel].assign(
          {static_cast<uint8>(StrandColor(is_forward_strand, options_))});
    } else if (channel == ch_read_supports_variant) {
      int supports_alt = ReadSupportsAlt(dv_call, read, alt_alleles);
      read_level_data_[channel].assign(
          {static_cast<uint8>(SupportsAltColor(supports_alt, options_))});
    } else if (channel == ch_read_mapping_percent) {
      read_level_data_[channel].assign(
          {ScaleColor(ReadMappingPercent(read), MaxMappingPercent)});
    } else if (channel == ch_avg_base_quality) {
      read_level_data_[channel].assign(
          {ScaleColor(AvgBaseQuality(read), MaxAvgBaseQuality)});
    } else if (channel == ch_identity) {
      read_level_data_[channel].assign({ScaleColor(Identity(read),
                                                    MaxIdentity)});
    } else if (channel == ch_gap_compressed_identity) {
      read_level_data_[channel].assign(
          {ScaleColor(GapCompressedIdentity(read), MaxIdentity)});
    } else if (channel == ch_gc_content) {
      read_level_data_[channel].assign({ScaleColor(GcContent(read),
                                                    MaxGcContent)});
    } else if (channel == ch_is_homopolymer) {
      std::vector<uint8> is_homopolymer = IsHomoPolymer(read);
      read_level_data_[channel] = ScaleColorVector(is_homopolymer,
                                                    MaxIsHomoPolymer);
    } else if (channel == ch_homopolymer_weighted) {
      std::vector<uint8> homopolymer_weighted = HomoPolymerWeighted(read);
      read_level_data_[channel] =
          ScaleColorVector(homopolymer_weighted, MaxHomoPolymerWeighted);
    } else if (channel == ch_blank) {
      read_level_data_[channel] = Blank(read);
    } else if (channel == ch_insert_size) {
      read_level_data_[channel] = ReadInsertSize(read);
    }

    return true;
  }

  // Calculate values for channels that depend on information at the
  // granularity of bases within the read and/or reference sequence.
  bool CalculateBaseLevelData(const Read& read,
                              std::function<
                                  bool(int, int, const CigarUnit::Operation&)>
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
          ok = action_per_cigar_unit(ref_i - 1 , read_i, op);
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

      if (!ok) {return false;}
    }

    return true;
  }

  inline uint8 GetChannelData(const std::string& channel, int col) {
    return data_[channel][col];
  }

  void CalculateRefRows(const std::vector<std::string>& channels,
                        const std::string& ref_bases) {
    // Calculates reference row values for each channel
    // Create a fake read to represent reference bases.
    Read refRead;
    for (const std::string& channel : channels) {
      if (channel == ch_read_base) {
        ref_data_[channel] = BaseColorVector(ref_bases, options_);
      } else if (channel == ch_base_quality) {
        int ref_qual = options_.reference_base_quality();
        ref_data_[channel].assign(
            {ScaleColor(ref_qual, options_.base_quality_cap())});
      } else if (channel == ch_mapping_quality) {
        int ref_qual = options_.reference_base_quality();
        ref_data_[channel].assign(
            {ScaleColor(ref_qual, options_.base_quality_cap())});
      } else if (channel == ch_strand) {
        int strand = StrandColor(true, options_);
        ref_data_[channel].assign({static_cast<uint8>(strand)});
      } else if (channel == ch_read_supports_variant) {
        int alt = SupportsAltColor(0, options_);
        ref_data_[channel].assign({static_cast<uint8>(alt)});
      } else if (channel == ch_base_differs_from_ref) {
        int ref = MatchesRefColor(true, options_);
        ref_data_[channel].assign({static_cast<uint8>(ref)});
      } else if (channel == ch_read_mapping_percent) {
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
