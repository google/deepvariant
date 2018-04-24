/*
 * Copyright 2018 Google Inc.
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

#include <map>

#include "deepvariant/realigner/fast_pass_aligner.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace learning {
namespace genomics {
namespace deepvariant {

class FastPassAlignerTest : public ::testing::Test {
 protected:
  FastPassAligner aligner_;

  void SetUp() override {
    aligner_.set_reference(
        "ATCAAGGGAAAAAGTGCCCAGGGCCAAATATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
            "CTGAAGATATG");
  }
};

TEST_F(FastPassAlignerTest, ReadsIndexIntegrationTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});

  KmerIndexType expected_index = {
    {"AAA", {KmerOccurrence(ReadId(0), KmerOffset(0))}},
    {"AAC", {KmerOccurrence(ReadId(0), KmerOffset(1))}},
    {"ACC", {KmerOccurrence(ReadId(0), KmerOffset(2))}},
    {"CCC", {KmerOccurrence(ReadId(0), KmerOffset(3))}},
    {"CTC",
     {KmerOccurrence(ReadId(1), KmerOffset(0)),
      KmerOccurrence(ReadId(1), KmerOffset(2))}},
    {"TCT",
     {KmerOccurrence(ReadId(1), KmerOffset(1)),
      KmerOccurrence(ReadId(1), KmerOffset(3))}},
    {"TGA",
     {KmerOccurrence(ReadId(2), KmerOffset(0)),
      KmerOccurrence(ReadId(2), KmerOffset(5))}},
    {"GAG", {KmerOccurrence(ReadId(2), KmerOffset(1))}},
    {"AGC", {KmerOccurrence(ReadId(2), KmerOffset(2))}},
    {"GCT", {KmerOccurrence(ReadId(2), KmerOffset(3))}},
    {"CTG", {KmerOccurrence(ReadId(2), KmerOffset(4))}},
    {"GAA", {KmerOccurrence(ReadId(2), KmerOffset(6))}},
    {"AAG", {KmerOccurrence(ReadId(2), KmerOffset(7))}}};

  aligner_.set_kmer_size(3);
  aligner_.BuildIndex();

  EXPECT_EQ(aligner_.GetKmerIndex(), expected_index);
}

// Test haplotype consists of read1 and read3. We check expected haplotype
// score and all read alignments.
TEST_F(FastPassAlignerTest, FastAlignReadsToHaplotypeTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});
  aligner_.set_kmer_size(3);
  aligner_.BuildIndex();

  // Haplotype made from read 3 + TT + read 1
  string haplotype = "TGAGCTGAAGTTAAACCC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<string> aligner_reads = aligner_.get_reads();
  int expected_hap_score =
      aligner_reads[2].length() * aligner_.get_match_score() +
      aligner_reads[0].length() * aligner_.get_match_score();

  std::vector<ReadAlignment> expected_read_scores(aligner_reads.size());
  expected_read_scores[0] = ReadAlignment(
      12, "6=", aligner_reads[0].length() * aligner_.get_match_score());
  expected_read_scores[1] = ReadAlignment();
  expected_read_scores[2] = ReadAlignment(
      0, "10=", aligner_reads[2].length() * aligner_.get_match_score());

  aligner_.FastAlignReadsToHaplotype(haplotype, &haplotype_score, &read_scores);
  EXPECT_EQ(expected_hap_score, haplotype_score);
  EXPECT_THAT(read_scores,
              testing::UnorderedElementsAreArray(expected_read_scores));
}

// One of the reads only partially overlaps haplotype. In this case the read
// has to be skipped.
TEST_F(FastPassAlignerTest, FastAlignReadsToHaplotypePartialReadOverlapTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});
  aligner_.set_kmer_size(3);
  aligner_.BuildIndex();

  // Haplotype made from read 3 + TT + read first 4 bases of read 1
  // In this case we cannot align read as a while and skip it.
  string haplotype = "TGAGCTGAAGTTAAAC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<string> aligner_reads = aligner_.get_reads();
  int expected_hap_score =
      aligner_reads[2].length() * aligner_.get_match_score();

  std::vector<ReadAlignment> expected_read_scores(aligner_reads.size());
  expected_read_scores[0] = ReadAlignment();
  expected_read_scores[1] = ReadAlignment();
  expected_read_scores[2] = ReadAlignment(
      0, "10=", aligner_reads[2].length() * aligner_.get_match_score());

  aligner_.FastAlignReadsToHaplotype(haplotype, &haplotype_score, &read_scores);
  EXPECT_EQ(expected_hap_score, haplotype_score);
  EXPECT_THAT(read_scores,
              testing::UnorderedElementsAreArray(expected_read_scores));
}

// One read is aligned with one mismatch.
TEST_F(FastPassAlignerTest,
       FastAlignReadsToHaplotypeReadAlignedWithMismatchesTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});
  aligner_.set_kmer_size(3);
  aligner_.BuildIndex();

  // Haplotype made from read 3 with one mismatch + TT + read 1
  string haplotype = "TGAGCCGAAGTTAAACCC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<string> aligner_reads = aligner_.get_reads();
  int expected_hap_score =
      (aligner_reads[2].length() - 1) * aligner_.get_match_score()
      - 1 * aligner_.get_mismatch_penalty()
      + aligner_reads[0].length() * aligner_.get_match_score();

  std::vector<ReadAlignment> expected_read_scores(aligner_reads.size());
  expected_read_scores[0] = ReadAlignment(
      12, "6=", aligner_reads[0].length() * aligner_.get_match_score());
  expected_read_scores[1] = ReadAlignment();
  expected_read_scores[2] = ReadAlignment(
      0,
      "10=",
      (aligner_reads[2].length() - 1) * aligner_.get_match_score()
      - 1 * aligner_.get_mismatch_penalty());

  aligner_.FastAlignReadsToHaplotype(haplotype, &haplotype_score, &read_scores);
  EXPECT_EQ(expected_hap_score, haplotype_score);
  EXPECT_THAT(read_scores,
              testing::UnorderedElementsAreArray(expected_read_scores));
}

// One read is aligned with more than maximum allowed mismatches.
TEST_F(FastPassAlignerTest,
       FastAlignReadsToHaplotypeReadAlignedWithMoreThanMismatchesTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});
  aligner_.set_kmer_size(3);
  aligner_.set_max_num_of_mismatches(2);
  aligner_.BuildIndex();

  // Haplotype made from read 3 with 3 mismatches + TT + read 1
  // max_num_of_mismatches_ is set to 2
  string haplotype = "TTTGCCGAAGTTAAACCC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<string> aligner_reads = aligner_.get_reads();
  int expected_hap_score =
      aligner_reads[0].length() * aligner_.get_match_score();

  std::vector<ReadAlignment> expected_read_scores(aligner_reads.size());
  expected_read_scores[0] = ReadAlignment(
      12, "6=", aligner_reads[0].length() * aligner_.get_match_score());
  expected_read_scores[1] = ReadAlignment();
  expected_read_scores[2] = ReadAlignment();

  aligner_.FastAlignReadsToHaplotype(haplotype, &haplotype_score, &read_scores);
  EXPECT_EQ(expected_hap_score, haplotype_score);
  EXPECT_THAT(read_scores,
              testing::UnorderedElementsAreArray(expected_read_scores));
}

// This test is not intended to test SSW library. It is a sanity check that
// library can be called and results are as excepted.
TEST_F(FastPassAlignerTest, SswAlignerSanityCheck) {
  aligner_.InitSswLib();
  Filter filter;
  aligner_.SswSetReference("TTTGCCGAAGTTAAACCC");
  Alignment alignment = aligner_.SswAlign("GCCGAAGTTA");
  EXPECT_EQ(alignment.cigar_string, "10=");
  EXPECT_EQ(alignment.ref_begin, 3);
}

TEST_F(FastPassAlignerTest, AlignHaplotypesToReference_Test) {
  aligner_.InitSswLib();
  Filter filter;
  const string REF_SEQ = "AGAAGGTCCCTTTGCCGAAGTTAAACCCTTTCGCGC";
  aligner_.set_reference(REF_SEQ);
  std::vector<string> haplotypes = {
      "GTCCCTTTGCCGAAGTTAAACCCTTT",  // equals to reference
      "GTCCCTTTGCCGAGTTAAACCCTTT",   // has deletion
      "GTCCCTATGCCGAAGTTAAACCCTTT"   // has mismatch
  };

  HaplotypeReadsAlignment ha1(0, -1, std::vector<ReadAlignment>());
  ha1.cigar = "26=";
  ha1.cigar_ops = std::list<CigarOp>({CigarOp(
      nucleus::genomics::v1::CigarUnit_Operation_ALIGNMENT_MATCH, 26)});
  ha1.ref_pos = 5;

  HaplotypeReadsAlignment ha2(1, -1, std::vector<ReadAlignment>());
  ha2.cigar = "12=1D13=";
  ha2.cigar_ops = std::list<CigarOp>(
      {CigarOp(nucleus::genomics::v1::CigarUnit_Operation_ALIGNMENT_MATCH, 12),
       CigarOp(nucleus::genomics::v1::CigarUnit_Operation_DELETE, 1),
       CigarOp(nucleus::genomics::v1::CigarUnit_Operation_ALIGNMENT_MATCH,
               13)});
  ha2.ref_pos = 5;

  HaplotypeReadsAlignment ha3(2, -1, std::vector<ReadAlignment>());
  ha3.cigar = "6=1X19=";
  ha3.cigar_ops = std::list<CigarOp>(
      {CigarOp(nucleus::genomics::v1::CigarUnit_Operation_ALIGNMENT_MATCH, 6),
       CigarOp(nucleus::genomics::v1::CigarUnit_Operation_ALIGNMENT_MATCH, 1),
       CigarOp(nucleus::genomics::v1::CigarUnit_Operation_ALIGNMENT_MATCH,
               19)});
  ha3.ref_pos = 5;

  std::vector<HaplotypeReadsAlignment> expectedHaplotypeAlignments = {ha1, ha2,
                                                                      ha3};

  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  EXPECT_THAT(aligner_.GetReadToHaplotypeAlignments(),
              testing::UnorderedElementsAreArray(expectedHaplotypeAlignments));
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
