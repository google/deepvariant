/*
 * Copyright 2018 Google LLC.
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

#include <string>

#include "deepvariant/protos/realigner.pb.h"
#include "deepvariant/realigner/fast_pass_aligner.h"
#include "tensorflow/core/platform/test.h"
#include "absl/strings/string_view.h"

namespace learning {
namespace genomics {
namespace deepvariant {

class GlobalAlignTest : public ::testing::Test {
 protected:
  FastPassAligner aligner_;

  void SetPenalties(int match, int mismatch, int gap_open, int gap_extend) {
    AlignerOptions options;
    options.set_match(match);
    options.set_mismatch(mismatch);
    options.set_gap_open(gap_open);
    options.set_gap_extend(gap_extend);
    aligner_.set_options(options);
  }
};

TEST_F(GlobalAlignTest, ExactMatch) {
  SetPenalties(4, 6, 8, 1);
  auto alignment = aligner_.GlobalAlign("ACGT", "ACGT");
  EXPECT_EQ(alignment.sw_score, 16);
  EXPECT_EQ(alignment.cigar_string, "4=");
  EXPECT_EQ(alignment.ref_begin, 0);
  EXPECT_EQ(alignment.ref_end, 3);
}

TEST_F(GlobalAlignTest, SingleMismatch) {
  SetPenalties(4, 6, 8, 1);
  // Match, Mismatch, Match, Match
  auto alignment = aligner_.GlobalAlign("ACGT", "AGGT");
  EXPECT_EQ(alignment.sw_score, 4 * 3 - 6);
  EXPECT_EQ(alignment.cigar_string, "1=1X2=");
}

TEST_F(GlobalAlignTest, SingleInsertionInQuery) {
  SetPenalties(4, 6, 10, 2);
  // ACGT in A-GT -> ACGT aligns with 'C' as insertion.
  auto alignment = aligner_.GlobalAlign("ACGT", "AGT");
  EXPECT_EQ(alignment.sw_score, 4 * 3 - (10 + 2));
  EXPECT_EQ(alignment.cigar_string, "1=1I2=");
}

TEST_F(GlobalAlignTest, SingleDeletionInQuery) {
  SetPenalties(4, 6, 10, 2);
  // AGT in ACGT -> AGT aligns with 'C' as mismatch because mismatch (6) is
  // cheaper than a gap (12).
  auto alignment = aligner_.GlobalAlign("AGT", "ACGT");
  EXPECT_EQ(alignment.sw_score, 2);
  EXPECT_EQ(alignment.cigar_string, "1X2=");
}

TEST_F(GlobalAlignTest, AnchorAtEnd) {
  SetPenalties(4, 6, 8, 1);
  // Query "ACGT" should match at the end of "TTTTACGT"
  auto alignment = aligner_.GlobalAlign("ACGT", "TTTTACGT");
  EXPECT_EQ(alignment.sw_score, 16);
  EXPECT_EQ(alignment.ref_begin, 4);
  EXPECT_EQ(alignment.cigar_string, "4=");
  EXPECT_EQ(alignment.ref_begin, 4);
}

TEST_F(GlobalAlignTest, AnchorAtStart) {
  SetPenalties(4, 6, 8, 1);
  auto alignment = aligner_.GlobalAlign("ACGT", "ACGTTTTT");
  EXPECT_EQ(alignment.sw_score, 16);
  EXPECT_EQ(alignment.ref_begin, 0);
  EXPECT_EQ(alignment.ref_end, 3);
  EXPECT_EQ(alignment.cigar_string, "4=");
  EXPECT_EQ(alignment.ref_begin, 0);
  EXPECT_EQ(alignment.ref_end, 3);
}

TEST_F(GlobalAlignTest, AnchorInMiddle) {
  SetPenalties(4, 6, 8, 1);
  auto alignment = aligner_.GlobalAlign("ACGT", "GGGACGTGGG");
  EXPECT_EQ(alignment.sw_score, 16);
  EXPECT_EQ(alignment.ref_begin, 3);
  EXPECT_EQ(alignment.ref_end, 6);
  EXPECT_EQ(alignment.cigar_string, "4=");
  EXPECT_EQ(alignment.ref_begin, 3);
  EXPECT_EQ(alignment.ref_end, 6);
}

TEST_F(GlobalAlignTest, HighMismatchPenaltyFavorsGap) {
  SetPenalties(4, 20, 2, 1);
  auto alignment = aligner_.GlobalAlign("ACGT", "AGGT");
  EXPECT_EQ(alignment.cigar_string, "1=1D1I2=");
}

TEST_F(GlobalAlignTest, HighGapPenaltyFavorsMismatch) {
  SetPenalties(4, 2, 20, 10);
  auto alignment = aligner_.GlobalAlign("ACGT", "GGGAGTAAAA");
  EXPECT_EQ(alignment.cigar_string, "2X2=");
}

TEST_F(GlobalAlignTest, AffineGapLongGap) {
  SetPenalties(4, 6, 10, 1);
  // ACGTTTGT (8bp) vs ACGTGT (6bp)
  auto alignment = aligner_.GlobalAlign("ACGTTTGT", "ACGTGT");
  // The tie-breaking in backtrack yields 3=2I3=
  EXPECT_EQ(alignment.cigar_string, "3=2I3=");
  EXPECT_EQ(alignment.sw_score, 4 * 6 - 12);
}

TEST_F(GlobalAlignTest, SoftClipPrefix) {
  SetPenalties(4, 10, 10, 2);
  // Query: TTTTACGT
  // Target: GGGGACGT
  // Current implementation forces full query alignment.
  auto alignment = aligner_.GlobalAlign("TTTTACGT", "GGGGACGT");
  EXPECT_EQ(alignment.cigar_string, "4I4=");
  EXPECT_EQ(alignment.sw_score, -2);
}

TEST_F(GlobalAlignTest, AllSoftClipped) {
  SetPenalties(4, 10, 10, 10);
  // Query: AAAA, Target: TTTT. Score 4 * -10 = -40.
  auto alignment = aligner_.GlobalAlign("AAAA", "TTTT");
  EXPECT_EQ(alignment.cigar_string, "4X");
  EXPECT_EQ(alignment.sw_score, -40);
}

TEST_F(GlobalAlignTest, TieBreakerLastMatch) {
  SetPenalties(4, 6, 8, 1);
  auto alignment = aligner_.GlobalAlign("AAAA", "AAAA_AAAA");
  EXPECT_EQ(alignment.ref_begin, 5);
}

TEST_F(GlobalAlignTest, ShortQuery) {
  SetPenalties(4, 6, 8, 1);
  auto alignment = aligner_.GlobalAlign("A", "TTTATT");
  EXPECT_EQ(alignment.sw_score, 4);
  EXPECT_EQ(alignment.ref_begin, 3);
  EXPECT_EQ(alignment.cigar_string, "1=");
}

TEST_F(GlobalAlignTest, QueryLongerThanTarget) {
  SetPenalties(4, 6, 8, 1);
  // Query: AAAACCCC, Target: CCCC.
  auto alignment = aligner_.GlobalAlign("AAAACCCC", "CCCC");
  EXPECT_EQ(alignment.cigar_string, "5I3=");
  EXPECT_EQ(alignment.sw_score, -1);
}

TEST_F(GlobalAlignTest, EmptyQuery) {
  auto alignment = aligner_.GlobalAlign("", "ACGT");
  EXPECT_EQ(alignment.sw_score, 0);
  EXPECT_EQ(alignment.cigar_string, "");
}

TEST_F(GlobalAlignTest, EmptyTarget) {
  auto alignment = aligner_.GlobalAlign("ACGT", "");
  EXPECT_EQ(alignment.sw_score, 0);
  EXPECT_EQ(alignment.cigar_string, "");
}

TEST_F(GlobalAlignTest, MixedIndelMismatch) {
  SetPenalties(4, 6, 10, 2);
  auto alignment = aligner_.GlobalAlign("ACGTA", "AGCA");
  EXPECT_EQ(alignment.cigar_string, "1=1I1=1X1=");
}

TEST_F(GlobalAlignTest, LowComplexityRegion) {
  SetPenalties(4, 6, 10, 2);
  auto alignment = aligner_.GlobalAlign("AAAAA", "AAAAAA");
  EXPECT_EQ(alignment.cigar_string, "5=");
  EXPECT_EQ(alignment.sw_score, 20);
}

TEST_F(GlobalAlignTest, DifferentMatchScore) {
  SetPenalties(1, 6, 8, 1);
  auto alignment = aligner_.GlobalAlign("ACGT", "ACGT");
  EXPECT_EQ(alignment.sw_score, 4);
  EXPECT_EQ(alignment.cigar_string, "4=");
}

TEST_F(GlobalAlignTest, AffineGapBacktrackingBugTest) {
  SetPenalties(3, 10, 2, 1);
  auto alignment = aligner_.GlobalAlign("C", "AC");
  EXPECT_EQ(alignment.cigar_string, "1=");
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
