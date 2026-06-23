/*
 * Copyright 2025 Google LLC.
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

// Tests for the Ultima Genomics homopolymer quality channels:
//   - HomopolymerInsertionQualityChannel (uses TP tag)
//   - HomopolymerDeletionQualityChannel  (uses TP tag)
//   - InterHomopolymerInsertionQualityChannel (uses T0 tag)
//
// Background:
//   Ultima Genomics data encodes per-base error information in BAM AUX tags:
//
//   TP tag (signed int array, one per read base):
//     Indicates the direction and magnitude of homopolymer length errors.
//     Given a homopolymer of decided length L:
//       tp[i] > 0: QUAL[i] encodes P(true length = L + tp[i])  (insertion)
//       tp[i] < 0: QUAL[i] encodes P(true length = L + tp[i])  (deletion)
//       tp[i] = 0: no error information at this position
//
//   T0 tag (ASCII Phred string, one char per read base):
//     Encodes inter-homopolymer insertion probability. Each char is
//     a Phred score + 33. For example, '5' = Q20, 'I' = Q40.

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

#include "deepvariant/channels/channel_utils.h"
#include "deepvariant/channels/homopolymer_deletion_quality_channel.h"
#include "deepvariant/channels/homopolymer_insertion_quality_channel.h"
#include "deepvariant/channels/inter_homopolymer_insertion_quality_channel.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "tensorflow/core/platform/test.h"
#include "absl/types/span.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace {

using ::nucleus::genomics::v1::Read;

// ---------------------------------------------------------------------------
// Test helpers
// ---------------------------------------------------------------------------

// Creates a Read with TP tag set as a repeated int array and custom qualities.
Read MakeReadWithTP(const std::string& seq, absl::Span<const int> tp_values,
                    absl::Span<const int> base_qualities) {
  Read read = nucleus::MakeRead("chr1", 100, seq,
                                {std::to_string(seq.size()) + "M"});
  for (int i = 0; i < base_qualities.size(); i++) {
    read.mutable_aligned_quality()->at(i) = base_qualities[i];
  }
  auto& tp_info = (*read.mutable_info())["tp"];
  for (int val : tp_values) {
    tp_info.add_values()->set_int_value(val);
  }
  return read;
}

// Creates a Read with T0 tag set as an ASCII-encoded Phred string.
Read MakeReadWithT0(const std::string& seq, const std::string& t0_string) {
  Read read = nucleus::MakeRead("chr1", 100, seq,
                                {std::to_string(seq.size()) + "M"});
  auto& t0_info = (*read.mutable_info())["t0"];
  t0_info.add_values()->set_string_value(t0_string);
  return read;
}

PileupImageOptions MakeOptions(int quality_cap = 40) {
  PileupImageOptions options;
  options.set_base_quality_cap(quality_cap);
  return options;
}

// ---------------------------------------------------------------------------
// HomoPolymerInDelQuality tests (public method on base class)
// ---------------------------------------------------------------------------
// HomoPolymerInDelQuality computes a per-base quality vector from TP and QUAL.
// It aggregates insertion or deletion error probabilities within each
// homopolymer run, then broadcasts the result to all positions in that run.

class HomopolymerInDelQualityTest : public ::testing::Test {
 protected:
  PileupImageOptions options_ = MakeOptions();
  // Test through HomopolymerInsertionQualityChannel since
  // HomoPolymerInDelQuality is public on the base class.
  HomopolymerInsertionQualityChannel channel_{10, options_};
};

TEST_F(HomopolymerInDelQualityTest, InsertionDirection_BasicCase) {
  // Read: AAAGGG (3xA hmer, 3xG hmer)
  // TP:   1, 0, 0, -1, 0, 0
  // QUAL: 20,30,30, 15,30,30
  //
  // For insertion (is_deletion=false):
  //   A-hmer: tp[0]=1 (insertion) → contributes QUAL[0]=20
  //   G-hmer: tp[3]=-1 (deletion) → ignored for insertion → max quality
  Read read = MakeReadWithTP("AAAGGG",
                              {1, 0, 0, -1, 0, 0},
                              {20, 30, 30, 15, 30, 30});
  auto ins_qual = channel_.HomoPolymerInDelQuality(read, /*is_deletion=*/false);
  ASSERT_EQ(ins_qual.size(), 6);

  // All positions within same hmer should be equal.
  EXPECT_EQ(ins_qual[0], ins_qual[1]);
  EXPECT_EQ(ins_qual[1], ins_qual[2]);
  EXPECT_EQ(ins_qual[3], ins_qual[4]);
  EXPECT_EQ(ins_qual[4], ins_qual[5]);

  // G-hmer should be at max quality (no insertion errors).
  uint8_t max_color =
      channels::internal::MaxQualityColor(options_.base_quality_cap());
  EXPECT_EQ(ins_qual[3], max_color);

  // A-hmer had an insertion error, so its quality should be lower.
  EXPECT_LT(ins_qual[0], max_color);
}

TEST_F(HomopolymerInDelQualityTest, DeletionDirection_BasicCase) {
  // Same read, now checking deletion direction.
  // For deletion (is_deletion=true):
  //   A-hmer: tp[0]=1 (insertion) → ignored → max quality
  //   G-hmer: tp[3]=-1 (deletion) → contributes QUAL[3]=15
  Read read = MakeReadWithTP("AAAGGG",
                              {1, 0, 0, -1, 0, 0},
                              {20, 30, 30, 15, 30, 30});
  auto del_qual = channel_.HomoPolymerInDelQuality(read, /*is_deletion=*/true);
  ASSERT_EQ(del_qual.size(), 6);

  uint8_t max_color =
      channels::internal::MaxQualityColor(options_.base_quality_cap());
  // A-hmer: max quality (no deletion errors).
  EXPECT_EQ(del_qual[0], max_color);

  // G-hmer: has a deletion error, should be lower.
  EXPECT_LT(del_qual[3], max_color);

  // Uniform within each hmer.
  EXPECT_EQ(del_qual[0], del_qual[1]);
  EXPECT_EQ(del_qual[3], del_qual[4]);
}

TEST_F(HomopolymerInDelQualityTest, NoTPTag_ReturnsMaxQuality) {
  Read read = nucleus::MakeRead("chr1", 100, "AAAGGG", {"6M"});
  auto qual = channel_.HomoPolymerInDelQuality(read, false);
  ASSERT_EQ(qual.size(), 6);
  uint8_t max_color =
      channels::internal::MaxQualityColor(options_.base_quality_cap());
  for (auto v : qual) {
    EXPECT_EQ(v, max_color);
  }
}

TEST_F(HomopolymerInDelQualityTest, MultipleErrorsSameHomopolymer) {
  // Read: AAAAA (5-base A hmer)
  // TP:   1, 2, 0, 0, 0  — two insertion errors
  // Both contribute to the total insertion error probability:
  //   P(total) = 10^(-QUAL[0]/10) + 10^(-QUAL[1]/10)
  // This should produce a lower quality than a single error.
  Read read_multi = MakeReadWithTP("AAAAA",
                                    {1, 2, 0, 0, 0},
                                    {20, 10, 30, 30, 30});
  Read read_single = MakeReadWithTP("AAAAA",
                                     {1, 0, 0, 0, 0},
                                     {20, 30, 30, 30, 30});

  auto qual_multi = channel_.HomoPolymerInDelQuality(read_multi, false);
  auto qual_single = channel_.HomoPolymerInDelQuality(read_single, false);

  // All positions equal within the hmer.
  for (int i = 1; i < 5; i++) {
    EXPECT_EQ(qual_multi[i], qual_multi[0]);
    EXPECT_EQ(qual_single[i], qual_single[0]);
  }

  // Two errors should yield lower quality than one.
  EXPECT_LT(qual_multi[0], qual_single[0]);
}

TEST_F(HomopolymerInDelQualityTest, SingleBaseHomopolymers) {
  // Each base is its own 1-base hmer: A, C, G, T
  // TP: 1, -1, 0, 1
  Read read = MakeReadWithTP("ACGT", {1, -1, 0, 1}, {20, 15, 30, 10});
  auto ins_qual = channel_.HomoPolymerInDelQuality(read, false);
  ASSERT_EQ(ins_qual.size(), 4);

  uint8_t max_color =
      channels::internal::MaxQualityColor(options_.base_quality_cap());
  // Position 0: tp=1 (insertion) → has insertion error
  EXPECT_LT(ins_qual[0], max_color);
  // Position 1: tp=-1 (deletion) → no insertion error → max
  EXPECT_EQ(ins_qual[1], max_color);
  // Position 2: tp=0 → no error → max
  EXPECT_EQ(ins_qual[2], max_color);
  // Position 3: tp=1 (insertion) → has insertion error
  EXPECT_LT(ins_qual[3], max_color);
}

// Regression test: A homopolymer run longer than 255 bases previously caused a
// crash. The old code stored homopolymer lengths in a uint8_t (max 255), so a
// 300-base run was seen as length 255. This caused the loop to only advance by
// 255, then re-process positions 255-299 as a new 45-base run, reading
// tps[255+j] out of bounds.
TEST_F(HomopolymerInDelQualityTest, LongHomopolymer_DoesNotCrash) {
  const int long_hmer_len = 300;  // Longer than uint8_t max (255).
  std::string seq(long_hmer_len, 'A');
  seq += "GG";  // Append a short second homopolymer.
  const size_t total_len = seq.size();

  // All TP = 0 (no errors). Any quality values work.
  std::vector<int> tps(total_len, 0);
  std::vector<int> quals(total_len, 30);

  Read read = MakeReadWithTP(seq, tps, quals);
  // This call would crash before the fix due to out-of-bounds read.
  auto result = channel_.HomoPolymerInDelQuality(read, /*is_deletion=*/false);

  ASSERT_EQ(result.size(), total_len);
  // With no errors (all TP=0), every position should be at max quality.
  uint8_t max_color =
      channels::internal::MaxQualityColor(options_.base_quality_cap());
  for (size_t i = 0; i < total_len; i++) {
    EXPECT_EQ(result[i], max_color) << "Mismatch at position " << i;
  }
}

// ---------------------------------------------------------------------------
// FillReadBase tests — alignment correctness
// ---------------------------------------------------------------------------
// These verify that FillReadBase maps read_index (position in the read's
// aligned_sequence) to the correct image column.

class FillReadBaseAlignmentTest : public ::testing::Test {
 protected:
  PileupImageOptions options_ = MakeOptions();
  DeepVariantCall dv_call_ = DeepVariantCall::default_instance();
  std::vector<std::string> alt_alleles_;
};

TEST_F(FillReadBaseAlignmentTest, InsertionChannel_CorrectColumnMapping) {
  HomopolymerInsertionQualityChannel channel(10, options_);
  Read read = MakeReadWithTP("AAAGGG",
                              {1, 0, 0, -1, 0, 0},
                              {20, 30, 30, 15, 30, 30});
  auto expected = channel.HomoPolymerInDelQuality(read, false);

  // Simulate CIGAR walking for a simple 6M alignment.
  std::vector<unsigned char> data(10, 0);
  for (int i = 0; i < 6; i++) {
    channel.FillReadBase(data, /*col=*/i,
                         read.aligned_sequence()[i], 'A', 30, read,
                         /*read_index=*/i, dv_call_, alt_alleles_);
  }

  for (int i = 0; i < 6; i++) {
    EXPECT_EQ(data[i], expected[i]) << "Mismatch at column " << i;
  }
}

TEST_F(FillReadBaseAlignmentTest, DeletionChannel_CorrectColumnMapping) {
  HomopolymerDeletionQualityChannel channel(10, options_);
  Read read = MakeReadWithTP("AAAGGG",
                              {1, 0, 0, -1, 0, 0},
                              {20, 30, 30, 15, 30, 30});
  auto expected = channel.HomoPolymerInDelQuality(read, true);

  std::vector<unsigned char> data(10, 0);
  for (int i = 0; i < 6; i++) {
    channel.FillReadBase(data, i, read.aligned_sequence()[i], 'A', 30, read,
                         i, dv_call_, alt_alleles_);
  }
  for (int i = 0; i < 6; i++) {
    EXPECT_EQ(data[i], expected[i]) << "Mismatch at column " << i;
  }
}

TEST_F(FillReadBaseAlignmentTest, InterHPChannel_CorrectColumnMapping) {
  InterHomopolymerInsertionQualityChannel channel(10, options_);
  // T0: '5' = Q20, 'I' = Q40
  Read read = MakeReadWithT0("AAAATTTT", "5555IIII");

  std::vector<unsigned char> data(10, 0);
  for (int i = 0; i < 8; i++) {
    channel.FillReadBase(data, i, read.aligned_sequence()[i], 'A', 30, read,
                         i, dv_call_, alt_alleles_);
  }

  // Verify expected color values.
  uint8_t color_q20 = channels::internal::BaseQualityColor(20, 40);
  uint8_t color_q40 = channels::internal::BaseQualityColor(40, 40);
  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(data[i], color_q20) << "Position " << i << " should be Q20";
  }
  for (int i = 4; i < 8; i++) {
    EXPECT_EQ(data[i], color_q40) << "Position " << i << " should be Q40";
  }
}

TEST_F(FillReadBaseAlignmentTest, OutOfBoundsReadIndex_ReturnsZero) {
  // During DELETE CIGAR ops, read_index can be -1 (if deletion is the
  // first cigar element). FillReadBase should safely return 0.
  HomopolymerInsertionQualityChannel channel(10, options_);
  Read read = MakeReadWithTP("AAA", {0, 0, 0}, {30, 30, 30});

  std::vector<unsigned char> data(10, 99);  // Fill with sentinel.
  channel.FillReadBase(data, 0, 'A', 'A', 30, read,
                       /*read_index=*/-1, dv_call_, alt_alleles_);
  EXPECT_EQ(data[0], 0);

  channel.FillReadBase(data, 1, 'A', 'A', 30, read,
                       /*read_index=*/100, dv_call_, alt_alleles_);
  EXPECT_EQ(data[1], 0);
}

TEST_F(FillReadBaseAlignmentTest, InterHP_MissingT0Tag_AllZero) {
  // When T0 tag is absent, all positions should be zero.
  InterHomopolymerInsertionQualityChannel channel(10, options_);
  Read read = nucleus::MakeRead("chr1", 100, "ACGT", {"4M"});

  std::vector<unsigned char> data(10, 99);
  for (int i = 0; i < 4; i++) {
    channel.FillReadBase(data, i, read.aligned_sequence()[i], 'A', 30, read,
                         i, dv_call_, alt_alleles_);
  }
  for (int i = 0; i < 4; i++) {
    EXPECT_EQ(data[i], 0) << "Position " << i << " should be 0 without T0";
  }
}

TEST_F(FillReadBaseAlignmentTest, InterHP_ShortT0_PartialData) {
  // T0 string shorter than read: first positions get values, rest get 0.
  InterHomopolymerInsertionQualityChannel channel(10, options_);
  Read read = MakeReadWithT0("AAAAGGGG", "55");

  std::vector<unsigned char> data(10, 99);
  for (int i = 0; i < 8; i++) {
    channel.FillReadBase(data, i, read.aligned_sequence()[i], 'A', 30, read,
                         i, dv_call_, alt_alleles_);
  }

  uint8_t color_q20 = channels::internal::BaseQualityColor(20, 40);
  EXPECT_EQ(data[0], color_q20);
  EXPECT_EQ(data[1], color_q20);
  for (int i = 2; i < 8; i++) {
    EXPECT_EQ(data[i], channels::internal::BaseQualityColor(0, 40))
        << "Position " << i << " should map to Q0 color";
  }
}

// ---------------------------------------------------------------------------
// channel_utils tests
// ---------------------------------------------------------------------------

TEST(BaseQualityColorTest, BasicScaling) {
  // BaseQualityColor(qual, max) = min(qual, max) * 254.0 / max
  EXPECT_EQ(channels::internal::BaseQualityColor(0, 40), 0);
  EXPECT_EQ(channels::internal::BaseQualityColor(20, 40), 127);
  EXPECT_EQ(channels::internal::BaseQualityColor(40, 40), 254);
}

TEST(BaseQualityColorTest, ClampedAboveCap) {
  EXPECT_EQ(channels::internal::BaseQualityColor(60, 40),
            channels::internal::BaseQualityColor(40, 40));
}

TEST(MaxQualityColorTest, EqualsBaseQualityColorAtMax) {
  EXPECT_EQ(channels::internal::MaxQualityColor(40),
            channels::internal::BaseQualityColor(40, 40));
}

}  // namespace
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
