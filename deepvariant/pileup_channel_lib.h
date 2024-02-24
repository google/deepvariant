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

#include <cstdint>
#include <functional>
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

using nucleus::genomics::v1::Read;

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

static const auto& ch_haplotype_tag = "haplotype";
static const auto& ch_allele_frequency = "allele_frequency";
static const auto& ch_diff_channels_alternate_allele_1 =
    "diff_channels_alternate_allele_1";
static const auto& ch_diff_channels_alternate_allele_2 =
    "diff_channels_alternate_allele_2";
static const auto& ch_read_mapping_percent = "read_mapping_percent";
static const auto& ch_avg_base_quality = "avg_base_quality";
static const auto& ch_identity = "identity";
static const auto& ch_gap_compressed_identity = "gap_compressed_identity";
static const auto& ch_gc_content = "gc_content";
static const auto& ch_is_homopolymer = "is_homopolymer";
static const auto& ch_homopolymer_weighted = "homopolymer_weighted";
static const auto& ch_blank = "blank";
static const auto& ch_insert_size = "insert_size";
static const auto& ch_base_channels_alternate_allele_1 =
    "base_channels_alternate_allele_1";
static const auto& ch_base_channels_alternate_allele_2 =
    "base_channels_alternate_allele_2";

//-------------------//
// Channels Accessor //
//-------------------//

// Max values for scaling
static const constexpr int MaxMappingPercent = 100;
static const constexpr int MaxAvgBaseQuality = 93;
static const constexpr int MaxIdentity = 100;
static const constexpr int MaxGapCompressedIdentity = 100;
static const constexpr int MaxGcContent = 100;
static const constexpr int MaxIsHomoPolymer = 1;
static const constexpr int MaxHomoPolymerWeighted = 30;

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

class Channels {
 public:  // public only for tests
  int getMaxEnumValue(const std::vector<DeepVariantChannelEnum>& channel_enums);

  constexpr bool isBaseLevelChannel(DeepVariantChannelEnum channelEnum);

  // Scales an input value to pixel range 0-254.
  std::uint8_t ScaleColor(int value, float max_val);

  // Scales an input vector to pixel range 0-254
  std::vector<std::uint8_t> ScaleColorVector(
      std::vector<std::uint8_t>& channel_values, float max_val);

  //---------------//
  // Base Channels //
  //---------------//
  int BaseColor(char base, const PileupImageOptions& options);
  std::vector<std::uint8_t> BaseColorVector(const std::string& bases,
                                            const PileupImageOptions& options);

  int StrandColor(bool on_positive_strand, const PileupImageOptions& options);

  int SupportsAltColor(int read_supports_alt,
                       const PileupImageOptions& options);

  // Does this read support ref, one of the alternative alleles, or an allele we
  // aren't considering?
  int ReadSupportsAlt(const DeepVariantCall& dv_call,
                      const nucleus::genomics::v1::Read& read,
                      const std::vector<std::string>& alt_alleles);

  // Returns a value based on whether the current read base matched the
  // reference base it was compared to.
  int MatchesRefColor(bool base_matches_ref, const PileupImageOptions& options);

  //-----------------------//
  // Experimental Channels //
  //-----------------------//
  std::vector<std::uint8_t> ReadInsertSize(
      const nucleus::genomics::v1::Read& read);

  // normalizes a Read's `fragment_length` to a pixel value
  int normalizeFragmentLength(const Read& read);

  bool channel_exists(std::vector<std::string>& channels,
                      const std::string& channel_name);

  std::vector<std::uint8_t> Blank(const Read& read);
  std::vector<std::uint8_t> HomoPolymerWeighted(
      const nucleus::genomics::v1::Read& read);
  std::vector<std::uint8_t> IsHomoPolymer(
      const nucleus::genomics::v1::Read& read);
  int GcContent(const nucleus::genomics::v1::Read& read);
  // Gap Compressed Identity: Ins/Del treated as individual events.
  int GapCompressedIdentity(const nucleus::genomics::v1::Read& read);
  // Identity: Similar to mapping percent but with a slightly different def.
  int Identity(const nucleus::genomics::v1::Read& read);
  // Average Base Quality: Averages base quality over length of read.
  int AvgBaseQuality(const nucleus::genomics::v1::Read& read);
  // Read Mapping Percent: Calculates percentage of bases mapped to reference.
  int ReadMappingPercent(const nucleus::genomics::v1::Read& read);

 public:
  static DeepVariantChannelEnum ChannelStrToEnum(const std::string& channel);
  static int GetHPValueForHPChannel(const nucleus::genomics::v1::Read& read,
                                    int hp_tag_for_assembly_polishing);
  // Get allele frequency color for a read.
  // Convert a frequency value in float to color intensity (int) and normalize.
  static unsigned char AlleleFrequencyColor(float allele_frequency,
                                            const PileupImageOptions& options);

  // Get the allele frequency of the alt allele that is carried by a read.
  static float ReadAlleleFrequency(const DeepVariantCall& dv_call,
                                   const nucleus::genomics::v1::Read& read,
                                   const std::vector<std::string>& alt_alleles);

  const PileupImageOptions& options_;
  std::vector<std::vector<unsigned char>> ref_data_;

  std::vector<std::vector<unsigned char>> read_level_data_;
  std::vector<std::vector<unsigned char>> data_;
  std::vector<int> channel_enum_to_index_;

  bool CalculateChannels(
      const std::vector<DeepVariantChannelEnum>& channel_enums,
      const nucleus::genomics::v1::Read& read, const std::string& ref_bases,
      const DeepVariantCall& dv_call,
      const std::vector<std::string>& alt_alleles, int image_start_pos);

  // Calculate values for channels that only depend on information at the
  // granularity of an entire read.
  bool CalculateReadLevelData(const DeepVariantChannelEnum channel_enum,
                              const nucleus::genomics::v1::Read& read,
                              const DeepVariantCall& dv_call,
                              const std::vector<std::string>& alt_alleles);

  // Calculate values for channels that depend on information at the
  // granularity of bases within the read and/or reference sequence.
  bool CalculateBaseLevelData(
      const Read& read,
      std::function<bool(int, int,
                         const nucleus::genomics::v1::CigarUnit::Operation&)>
          action_per_cigar_unit);

  std::uint8_t GetChannelData(const std::string& channel, int col);

  void CalculateRefRows(
      const std::vector<DeepVariantChannelEnum>& channel_enums,
      const std::string& ref_bases);
  std::uint8_t GetRefRows(DeepVariantChannelEnum channel_enum, int col);
};
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

#endif  // LEARNING_GENOMICS_DEEPVARIANT_PILEUP_CHANNEL_LIB_H_
