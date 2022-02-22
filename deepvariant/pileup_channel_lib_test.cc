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

#include "deepvariant/pileup_channel_lib.h"

#include <vector>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"

using nucleus::genomics::v1::Read;

namespace learning {
namespace genomics {
namespace deepvariant {

TEST(ScaleColor, BasicCase) {
  EXPECT_EQ(ScaleColor(50, 100), 127);
  EXPECT_EQ(ScaleColor(127, kMaxPixelValueAsFloat), 127);
  // Beyond max scales to max.
  EXPECT_EQ(ScaleColor(500, kMaxPixelValueAsFloat), kMaxPixelValueAsFloat);
}

TEST(ScaleColorVector, BasicCase) {
  std::vector<uint8> test_vector{5, 10, 25};
  std::vector<uint8> expect_vector{25, 50, 127};
  EXPECT_EQ(ScaleColorVector(test_vector, 50), expect_vector);
}

TEST(ScaleColorVectorLarge, OverMaxCase) {
  std::vector<uint8> test_vector;
  test_vector.resize(500);
  for (int i = 0; i < test_vector.size(); i++) {
    test_vector[i] = i + 1;
  }
  test_vector = ScaleColorVector(test_vector, kMaxPixelValueAsFloat);
  uint8 j = 0;
  for (auto& i : test_vector) {
    j++;
    if (i < kMaxPixelValueAsFloat) {
      EXPECT_EQ(i, j);
    } else {
      EXPECT_EQ(i, static_cast<int>(kMaxPixelValueAsFloat));
    }
  }
}

TEST(ReadMappingPercentTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"5M", "5D"});
  uint8 rmp = ReadMappingPercent(read);
  EXPECT_EQ(rmp, 50);
}

TEST(AvgBaseQualityTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"10M"});
  for (size_t i = 0; i < read.aligned_sequence().size(); ++i) {
    read.set_aligned_quality(i, i + 1);
  }
  uint8 rmp = AvgBaseQuality(read);
  EXPECT_EQ(rmp, 5);
}

TEST(AvgBaseQualityTest, BaseQualityTooHigh) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"10M"});
  for (size_t i = 0; i < read.aligned_sequence().size(); ++i) {
    read.set_aligned_quality(i, 100);
  }
  EXPECT_DEATH(AvgBaseQuality(read),
               "Encountered base quality outside of bounds");
}

TEST(IdentityTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"5M", "1I", "4M"});
  uint8 id = Identity(read);
  EXPECT_EQ(id, 90);
}

TEST(IdentityTest, PacBioStyleCigar) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"5=", "1X", "4="});
  uint8 id = Identity(read);
  EXPECT_EQ(id, 90);
}

TEST(GapCompressedIdentityTest, InsertionCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"3M", "4I", "3M"});
  uint8 id = GapCompressedIdentity(read);
  EXPECT_EQ(id, 85);
}

TEST(GapCompressedIdentityTest, DeletionCase) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"3M", "4D", "3M"});
  uint8 id = GapCompressedIdentity(read);
  EXPECT_EQ(id, 85);
}

TEST(GapCompressedIdentityTest, PacBioStyleCigar) {
  Read read =
      nucleus::MakeRead("chr1", 1, "AAAAATTTTT", {"3=", "2X", "2I", "3="});
  uint8 id = GapCompressedIdentity(read);
  EXPECT_EQ(id, 66);
}

TEST(GcContestTest, AllGc) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGGGCCCCC", {"10M"});
  uint8 gc = GcContent(read);
  EXPECT_EQ(gc, 100);
}

TEST(GcContestTest, HalfGc) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGGGTTTTT", {"10M"});
  uint8 gc = GcContent(read);
  EXPECT_EQ(gc, 50);
}

TEST(IsHomoPolymerTest, IsHomopolymerBeginning) {
  Read read = nucleus::MakeRead("chr1", 1, "GGGATAATA", {"9M"});
  std::vector<uint8> is_homopolymer = IsHomoPolymer(read);
  std::vector<uint8> expected{1, 1, 1, 0, 0, 0, 0, 0, 0};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(IsHomoPolymerTest, IsHomopolymerMiddle) {
  Read read = nucleus::MakeRead("chr1", 1, "ATTGGGTTA", {"9M"});
  std::vector<uint8> is_homopolymer = IsHomoPolymer(read);
  std::vector<uint8> expected{0, 0, 0, 1, 1, 1, 0, 0, 0};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(IsHomoPolymerTest, IsHomopolymerEnd) {
  Read read = nucleus::MakeRead("chr1", 1, "ATAATAGGG", {"9M"});
  std::vector<uint8> is_homopolymer = IsHomoPolymer(read);
  std::vector<uint8> expected{0, 0, 0, 0, 0, 0, 1, 1, 1};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(IsHomoPolymerTest, IsHomopolymerAll) {
  Read read = nucleus::MakeRead("chr1", 1, "AAAAAAAAA", {"9M"});
  std::vector<uint8> is_homopolymer = IsHomoPolymer(read);
  std::vector<uint8> expected{1, 1, 1, 1, 1, 1, 1, 1, 1};
  EXPECT_EQ(is_homopolymer, expected);
}

TEST(HomoPolymerWeightedTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  std::vector<uint8> w_homopolymer = HomoPolymerWeighted(read);
  std::vector<uint8> expected{1, 1, 2, 2, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 5};
  EXPECT_EQ(w_homopolymer, expected);
}

TEST(HomoPolymerWeightedTest, WeightedHomoPolymerMax) {
  std::string bases;
  bases.insert(0, 20, 'A');
  bases.insert(0, 10, 'G');
  LOG(INFO) << bases;
  Read read = nucleus::MakeRead("chr1", 1, bases, {"60M"});
  std::vector<uint8> w_homopolymer = HomoPolymerWeighted(read);
  std::vector<uint8> expected;
  expected.insert(expected.begin(), 20, 20);
  expected.insert(expected.begin(), 10, 10);
  EXPECT_EQ(w_homopolymer, expected);
}

TEST(ReadInsertSizeTest, BasicCase) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  read.set_fragment_length(22);
  std::vector<uint8> w_insert_size = ReadInsertSize(read);
  std::vector<uint8> expected{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(ReadInsertSizeTest, NegativeValueCase) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  read.set_fragment_length(-22);
  std::vector<uint8> w_insert_size = ReadInsertSize(read);
  std::vector<uint8> expected{5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(ReadInsertSizeTest, ExceedsLimit) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  read.set_fragment_length(1001);
  std::vector<uint8> w_insert_size = ReadInsertSize(read);
  std::vector<uint8> expected{254, 254, 254, 254, 254, 254, 254, 254,
                              254, 254, 254, 254, 254, 254, 254, 254};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(ReadInsertSizeTest, NoValue) {
  Read read = nucleus::MakeRead("chr1", 1, "GATTGGGCCCCAAAAA", {"15M"});
  std::vector<uint8> w_insert_size = ReadInsertSize(read);
  std::vector<uint8> expected{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
  EXPECT_EQ(w_insert_size, expected);
}

TEST(GetChannelDataTest, ReadData) {
  OptChannels channel_set{};
  std::vector<std::string> channels{ch_read_mapping_percent,
                                    ch_avg_base_quality,
                                    ch_identity,
                                    ch_gap_compressed_identity,
                                    ch_gc_content,
                                    ch_is_homopolymer,
                                    ch_homopolymer_weighted,
                                    ch_blank};
  Read read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTT", {"8M"});
  const int base_quality = 33;
  for (size_t i = 0; i < read.aligned_sequence().size(); ++i) {
    read.set_aligned_quality(i, base_quality);
  }
  channel_set.CalculateChannels(channels, read);
  EXPECT_EQ(channel_set.GetChannelData(ch_read_mapping_percent, 3), 203);
  EXPECT_EQ(channel_set.GetChannelData(ch_avg_base_quality, 3), 90);
  EXPECT_EQ(channel_set.GetChannelData(ch_identity, 9), 203);
  EXPECT_EQ(channel_set.GetChannelData(ch_gap_compressed_identity, 9),
            static_cast<uint8>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetChannelData(ch_gc_content, 3), 152);
  EXPECT_EQ(channel_set.GetChannelData(ch_is_homopolymer, 1),
            static_cast<uint8>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetChannelData(ch_is_homopolymer, 4), 0);
  EXPECT_EQ(channel_set.GetChannelData(ch_homopolymer_weighted, 1), 25);
  EXPECT_EQ(channel_set.GetChannelData(ch_homopolymer_weighted, 9), 33);
  EXPECT_EQ(channel_set.GetChannelData(ch_blank, 0), 0);
}

TEST(GetRefChannelDataTest, ReadData) {
  OptChannels channel_set{};
  std::vector<std::string> channels{ch_read_mapping_percent,
                                    ch_avg_base_quality,
                                    ch_identity,
                                    ch_gap_compressed_identity,
                                    ch_gc_content,
                                    ch_is_homopolymer,
                                    ch_homopolymer_weighted,
                                    ch_blank};
  Read ref_read = nucleus::MakeRead("chr1", 1, "GGGCGCTTTT", {"8M"});
  const int base_quality = 33;
  for (size_t i = 0; i < ref_read.aligned_sequence().size(); ++i) {
    ref_read.set_aligned_quality(i, base_quality);
  }
  channel_set.CalculateRefRows(channels, ref_read.aligned_sequence());
  EXPECT_EQ(channel_set.GetRefRows(ch_read_mapping_percent, 3),
            static_cast<uint8>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_avg_base_quality, 3),
            static_cast<uint8>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_identity, 9), 254);
  EXPECT_EQ(channel_set.GetRefRows(ch_gap_compressed_identity, 9),
            static_cast<uint8>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_gc_content, 3), 152);
  EXPECT_EQ(channel_set.GetRefRows(ch_is_homopolymer, 1),
            static_cast<uint8>(kMaxPixelValueAsFloat));
  EXPECT_EQ(channel_set.GetRefRows(ch_is_homopolymer, 4), 0);
  EXPECT_EQ(channel_set.GetRefRows(ch_homopolymer_weighted, 1), 25);
  EXPECT_EQ(channel_set.GetRefRows(ch_homopolymer_weighted, 9), 33);
  EXPECT_EQ(channel_set.GetRefRows(ch_blank, 0), 0);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
