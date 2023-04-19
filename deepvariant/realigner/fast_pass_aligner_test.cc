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

#include "deepvariant/realigner/fast_pass_aligner.h"

#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "deepvariant/protos/realigner.pb.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/string_view.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "google/protobuf/text_format.h"

namespace learning {
namespace genomics {
namespace deepvariant {

class FastPassAlignerTest : public ::testing::Test {
 protected:
  FastPassAligner aligner_;
  static constexpr int kSomeScore = 100;

  void SetUp() override {
    // Test reference for most of the test. Some tests create a different
    // reference sequence.
    aligner_.set_reference(
        "ATCAAGGGAAAAAGTGCCCAGGGCCAAATATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
            "CTGAAGATATG");
  }

  // Load read protos from test data human readable file.
  std::vector<nucleus::genomics::v1::Read> LoadReadProtosFromFile(
      absl::string_view filename) {
    string file_path = nucleus::GetTestData(filename,
        "deepvariant/testdata/input");
    std::ifstream reads_stream(file_path);
    QCHECK(!reads_stream.fail()) << "Failed to read " << file_path;

    std::vector<nucleus::genomics::v1::Read> reads;
    string proto_str;
    string line;
    nucleus::genomics::v1::Read read;
    while (std::getline(reads_stream, line)) {
      if (line.empty()) {
        google::protobuf::TextFormat::ParseFromString(proto_str, &read);
        reads.push_back(read);
        proto_str.clear();
      } else {
        proto_str.append(line);
        proto_str.append("\n");
      }
    }
    reads_stream.close();
    return reads;
  }

  void LoadReferenceFromFile(absl::string_view filename, string* reference) {
    string file_path = nucleus::GetTestData(filename,
        "deepvariant/testdata/input");
    std::ifstream ref_stream(file_path);
    QCHECK(!ref_stream.fail()) << "Failed to read " << file_path;
    ref_stream >> *reference;
    ref_stream.close();
  }

  void LoadHaplotypesFromFile(absl::string_view filename,
                              std::vector<std::string>* haplotypes) {
    string file_path = nucleus::GetTestData(filename,
        "deepvariant/testdata/input");
    std::ifstream haps_stream(file_path);
    QCHECK(!haps_stream.fail()) << "Failed to read " << file_path;
    string line;
    while (std::getline(haps_stream, line)) {
      haplotypes->push_back(line);
    }
    haps_stream.close();
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

  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_.set_options(aligner_options);
  aligner_.BuildIndex();

  EXPECT_EQ(aligner_.GetKmerIndex(), expected_index);
}

// Test checks that we can handle reads shorter than kmer. In this test
// Kmer size is set to 4, and the first read has size 3.
TEST_F(FastPassAlignerTest, ReadsIndexIgnoreReadsShorterThanKmerTest) {
  aligner_.set_reads({"AAC", "TGAGCTG"});
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(4);
  aligner_.set_options(aligner_options);
  KmerIndexType expected_index = {
      {"TGAG", {KmerOccurrence(ReadId(1), KmerOffset(0))}},
      {"GAGC", {KmerOccurrence(ReadId(1), KmerOffset(1))}},
      {"AGCT", {KmerOccurrence(ReadId(1), KmerOffset(2))}},
      {"GCTG", {KmerOccurrence(ReadId(1), KmerOffset(3))}}
  };
  aligner_.BuildIndex();
  EXPECT_EQ(aligner_.GetKmerIndex(), expected_index);
}

// Test haplotype consists of read1 and read3. We check expected haplotype
// score and all read alignments.
TEST_F(FastPassAlignerTest, FastAlignReadsToHaplotypeTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_.set_options(aligner_options);
  aligner_.BuildIndex();

  // Haplotype made from read 3 + TT + read 1
  string haplotype = "TGAGCTGAAGTTAAACCC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<std::string> aligner_reads = aligner_.get_reads();
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
  aligner_.set_reads({"TGAGCTGAAGTT", "AAACCC", "AGTTAAAC"});
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_.set_options(aligner_options);
  aligner_.BuildIndex();

  // Haplotype made from read 3 + TT + read first 4 bases of read 1
  // In this case we cannot align read as a while and skip it.
  string haplotype = "TGAGCTGAAGTTAAAC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<std::string> aligner_reads = aligner_.get_reads();
  int expected_hap_score =
      aligner_reads[0].length() * aligner_.get_match_score() +
      aligner_reads[2].length() * aligner_.get_match_score();

  std::vector<ReadAlignment> expected_read_scores(aligner_reads.size());
  expected_read_scores[0] = ReadAlignment(
      0, "12=", aligner_reads[0].length() * aligner_.get_match_score());
  expected_read_scores[1] = ReadAlignment();
  expected_read_scores[2] = ReadAlignment(
      8, "8=", aligner_reads[2].length() * aligner_.get_match_score());

  aligner_.FastAlignReadsToHaplotype(haplotype, &haplotype_score, &read_scores);
  EXPECT_EQ(expected_hap_score, haplotype_score);
  EXPECT_THAT(read_scores,
              testing::UnorderedElementsAreArray(expected_read_scores));
}

// One read is aligned with one mismatch.
TEST_F(FastPassAlignerTest,
       FastAlignReadsToHaplotypeReadAlignedWithMismatchesTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_.set_options(aligner_options);
  aligner_.BuildIndex();

  // Haplotype made from read 3 with one mismatch + TT + read 1
  string haplotype = "TGAGCCGAAGTTAAACCC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<std::string> aligner_reads = aligner_.get_reads();
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

// Two reads are aligned with more than maximum allowed mismatches.
TEST_F(FastPassAlignerTest,
       FastAlignReadsToHaplotypeReadAlignedWithMoreThanMismatchesTest) {
  // First read spans the whole haplotype to make sure that we have coverage
  // over the entire haplotype (otherwise realignement stops as soon as it
  // hits the position with zero coverage).
  aligner_.set_reads({"TTTGCCGAAGTTAAACCC", "CTCTCT", "TGAGCTGAAG"});
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_options.set_max_num_of_mismatches(2);
  aligner_.set_options(aligner_options);
  aligner_.set_ref_prefix_len(0);
  aligner_.set_ref_suffix_len(0);
  aligner_.BuildIndex();

  string haplotype = "TTTGCCGAAGTTAAACCC";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<std::string> aligner_reads = aligner_.get_reads();
  int expected_hap_score =
      aligner_reads[0].length() * aligner_.get_match_score();

  std::vector<ReadAlignment> expected_read_scores(aligner_reads.size());
  expected_read_scores[0] = ReadAlignment(
      0, "18=", aligner_reads[0].length() * aligner_.get_match_score());
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
  std::vector<std::string> haplotypes = {
      "GTCCCTTTGCCGAAGTTAAACCCTTT",  // equals to reference
      "GTCCCTTTGCCGAGTTAAACCCTTT",   // has deletion
      "GTCCCTATGCCGAAGTTAAACCCTTT"   // has mismatch
  };

  HaplotypeReadsAlignment ha1(0, -1, std::vector<ReadAlignment>());
  ha1.cigar = "26=";
  ha1.cigar_ops = std::list<CigarOp>(
      {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 26)});
  ha1.ref_pos = 5;
  ha1.is_reference = true;

  HaplotypeReadsAlignment ha2(1, -1, std::vector<ReadAlignment>());
  ha2.cigar = "12=1D13=";
  ha2.cigar_ops = std::list<CigarOp>(
      {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 12),
       CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
       CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 13)});
  ha2.ref_pos = 5;
  ha2.is_reference = false;

  HaplotypeReadsAlignment ha3(2, -1, std::vector<ReadAlignment>());
  ha3.cigar = "6=1X19=";
  ha3.cigar_ops = std::list<CigarOp>(
      {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 6),
       CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 1),
       CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 19)});
  ha3.ref_pos = 5;
  ha3.is_reference = false;

  std::vector<HaplotypeReadsAlignment> expectedHaplotypeAlignments = {ha1, ha2,
                                                                      ha3};

  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  EXPECT_THAT(aligner_.GetReadToHaplotypeAlignments(),
              testing::UnorderedElementsAreArray(expectedHaplotypeAlignments));
}

TEST_F(FastPassAlignerTest, SetPositionsMapForNotStructuralAlignment_Test) {
  HaplotypeReadsAlignment haplotype_alignment;
  haplotype_alignment.cigar = "10=1X3=";  // Length = 24
  SetPositionsMap(24, &haplotype_alignment);
  std::vector<int> expected_positions_map(24, 0);
  EXPECT_THAT(haplotype_alignment.hap_to_ref_positions_map,
              testing::ElementsAreArray(expected_positions_map));
}

// Reference:  A  C  T  -  -  -  -  G  A
// Haplotype:  A  C  T  C  A  C  A  G  A
// Shifts:     0  0  0  0 -1 -2 -3 -4 -4
// If read to haplotype alignment has position = 3 then we can calculate
// read to reference position as pos_map[3] + 3 = 0 + 3 = 3
// for position = 4, pos_map[4] + 4 = -1 + 4 = 3
// etc.
TEST_F(FastPassAlignerTest, SetPositionsMapWithIns_Test) {
  HaplotypeReadsAlignment haplotype_alignment;
  haplotype_alignment.cigar = "3=4I2=";  // Length = 9
  SetPositionsMap(9, &haplotype_alignment);
  std::vector<int> expected_positions_map = {0, 0, 0, 0, -1, -2, -3, -4, -4};
  EXPECT_THAT(haplotype_alignment.hap_to_ref_positions_map,
              testing::ElementsAreArray(expected_positions_map));
}

TEST_F(FastPassAlignerTest, SetPositionsMapWithDel_Test) {
  HaplotypeReadsAlignment haplotype_alignment;
  haplotype_alignment.cigar = "3=4D2=";  // Length = 5
  SetPositionsMap(5, &haplotype_alignment);
  std::vector<int> expected_positions_map = {0, 0, 0, 4, 4};
  EXPECT_THAT(haplotype_alignment.hap_to_ref_positions_map,
              testing::ElementsAreArray(expected_positions_map));
}

// Reference:  A  C  T  A  C  C  T  G  T  -  -  G  A
// Haplotype:  A  C  T  -  -  -  -  G  T  C  A  G  A
// Shifts:     0  0  0              4  4  4  3  2  2
TEST_F(FastPassAlignerTest, SetPositionsMapWithDelAndIns_Test) {
  HaplotypeReadsAlignment haplotype_alignment;
  haplotype_alignment.cigar = "3=4D2=2I2=";  // Length = 9
  SetPositionsMap(9, &haplotype_alignment);
  std::vector<int> expected_positions_map = {0, 0, 0, 4, 4, 4, 3, 2, 2};
  EXPECT_THAT(haplotype_alignment.hap_to_ref_positions_map,
              testing::ElementsAreArray(expected_positions_map));
}

// Reference:  A  C  T  -  -  -  -  G  T  C  A  G  A
// Haplotype:  A  C  T  A  C  C  T  G  T  -  -  G  A
// Shifts:     0  0  0  0 -1 -2 -3 -4 -4       -2 -2
TEST_F(FastPassAlignerTest, SetPositionsMapWithInsAndDel_Test) {
  HaplotypeReadsAlignment haplotype_alignment;
  haplotype_alignment.cigar = "3=4I2=2D2=";  // Length = 11
  SetPositionsMap(11, &haplotype_alignment);
  std::vector<int> expected_positions_map = {
      0, 0, 0, 0, -1, -2, -3, -4, -4, -2, -2,
  };
  EXPECT_THAT(haplotype_alignment.hap_to_ref_positions_map,
              testing::ElementsAreArray(expected_positions_map));
}

// This test checks SswAlignReadsToHaplotypes.
// There are 2 haplotypes, and 5 reads. Some of the reads are better aligned
// to haplotype 1, some of the reads are better aligned to hap 2. This is
// reflected in the comments for each read.
// For example, last read has the best alignment to haplotype 1, but this
// alignment has a bad score (32). SswAlignReadsToHaplotypes is expected to
// not realign this read because it's score is lower than a threshold.
TEST_F(FastPassAlignerTest, SswAlignReadsToHaplotypes_Test) {
  aligner_.InitSswLib();
  std::vector<std::string> haplotypes = {
      // reference with 1 del
      "AAGTGCCCAGGGCCAAATGTTTTGGGTTTTGCAGGACAAAGTATGGTT",
      // reference with 1 sub
      "AAGTGCCCAGGGCCAAATATGCACAGGGTTTTGCAGGACAAAGTATGGTT"};

  // Read 1: exactly matches a subset of haplotype_1
  // Read 2: matches haplotype_2 with 2 mismatches
  // Read 3: subset of haploype_1 with 2 base del
  // Read 4: subset of haplotype_2 with one ins
  // Read 5: alignment score is less than threshold
  aligner_.set_reads({
                         "CAGGGCCAAATGTTT",         // "15=", 60
                         "GCCATATATGCACAGGGTTATG",  // "4=1X14=1X2=", 68
                         "TTGGGTTGCAGGACA",         // "5=2D10=", 51
                         "ACAGGGTTTTTTGCAGGACAA",   // "6=2I13=", 67
                         "TGTTGGGTTCAGCAGTTTT"      // "2S7=2X4=4S", 32
                     });
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_.set_options(aligner_options);
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  aligner_.SswAlignReadsToHaplotypes(40);
  std::vector<std::string> aligner_reads = aligner_.get_reads();
  std::vector<ReadAlignment> expected_read_alignments_for_hap1(
      aligner_reads.size());
  std::vector<ReadAlignment> expected_read_alignments_for_hap2(
      aligner_reads.size());
  expected_read_alignments_for_hap1[0] = ReadAlignment(7, "15=", 60);
  expected_read_alignments_for_hap1[1] = ReadAlignment();
  expected_read_alignments_for_hap1[2] = ReadAlignment(21, "5=2D10=", 51);
  expected_read_alignments_for_hap1[3] = ReadAlignment(23, "3S3=2I13=", 55);
  expected_read_alignments_for_hap1[4] = ReadAlignment();

  expected_read_alignments_for_hap2[0] = ReadAlignment(7, "11=4S", 44);
  expected_read_alignments_for_hap2[1] = ReadAlignment(11, "4=1X14=1X2=", 68);
  expected_read_alignments_for_hap2[2] = ReadAlignment(25, "2S3=2D10=", 43);
  expected_read_alignments_for_hap2[3] = ReadAlignment(22, "6=2I13=", 67);
  expected_read_alignments_for_hap2[4] = ReadAlignment();

  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();
  EXPECT_THAT(haplotype_alignments[0].read_alignment_scores,
              testing::ElementsAreArray(expected_read_alignments_for_hap1));
  EXPECT_THAT(haplotype_alignments[1].read_alignment_scores,
              testing::ElementsAreArray(expected_read_alignments_for_hap2));
}

// Haplotype to ref has one mismatch. Read matches haplotype exactly.
TEST_F(FastPassAlignerTest, CalculateReadToRefAlignment_MatchMismatch_Test) {
  aligner_.InitSswLib();
  std::vector<std::string> haplotypes = {
      // reference with 1 mismatch at position 5 (T->A).
      "TGTTTAGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"};
  // Reference is set in init.
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();

  // Read matches haplotype.
  aligner_.set_reads({
      "TGTTTAGGGTTTTGCAGGA",  // "19="
  });

  std::list<CigarOp> read_to_ref_cigar_ops;
  aligner_.CalculateReadToRefAlignment(0, ReadAlignment(7, "19=", kSomeScore),
                                       haplotype_alignments[0].cigar_ops,
                                       &read_to_ref_cigar_ops);

  std::list<CigarOp> expected_read_to_ref_cigar_ops = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 19)};

  EXPECT_THAT(read_to_ref_cigar_ops,
              testing::ElementsAreArray(expected_read_to_ref_cigar_ops));
}

// Haplotype alignment to ref starts before reference contig.
// Read matches haplotype.
TEST_F(FastPassAlignerTest,
       CalculateReadToRefAlignment_HaplotypeSoftClipped_Test) {
  aligner_.InitSswLib();
  aligner_.set_reference(
      "nnnnnnnnnnnTGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
      "CTGAAGATATG");
  std::vector<std::string> haplotypes = {
      // reference with 1 mismatch at position 5 (T->A).
      "GATCATGTTTAGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"};

  // Reference is set in init.
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();

  // Read matches haplotype exactly.
  aligner_.set_reads({
      "GATCATGTTTAGGGTTTT",  // "18="
  });

  std::list<CigarOp> read_to_ref_cigar_ops;
  aligner_.CalculateReadToRefAlignment(0, ReadAlignment(0, "19=", kSomeScore),
                                       haplotype_alignments[0].cigar_ops,
                                       &read_to_ref_cigar_ops);

  std::list<CigarOp> expected_read_to_ref_cigar_ops = {
      CigarOp(nucleus::genomics::v1::CigarUnit::CLIP_SOFT, 5),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 13)};

  EXPECT_THAT(read_to_ref_cigar_ops,
              testing::ElementsAreArray(expected_read_to_ref_cigar_ops));
}

struct CalculateReadToRefAlignmentTestCase {
  std::string test_case_name;
  std::vector<std::string> haplotypes;
  std::vector<std::string> reads;
  ReadAlignment read_alignment;
  std::list<CigarOp> expected_cigar;
};

TEST_F(FastPassAlignerTest, CalculateReadToRefAlignmentCommonTest) {
  aligner_.InitSswLib();
  aligner_.set_reference(
      "CTCTGTAATCGGATCATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
      "CTGAAGATATG");

  std::vector<CalculateReadToRefAlignmentTestCase> test_cases{
    {
      // This test verifies a correct merging of INS and SNP at the same
      // position. Specifically we verify that after merging INS and SNP we
      // correctly generate a merged insertion. More details can be found
      // in internal.
      "ins_snp_merge",
      {"CGGATCATGTTTTTTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 9=2I11=1D24=
      {"GATCATGATTTTTGGGTTTTCAG"},  // "7=1X15="
      // Read to haplotype has one 1 base ins at 7. Read aligns to haplotype
      // from position 2.
      // After trimming haplotype to reference cigar we get: 7=1I34=
      // Merging 7=2I11=1D24= and 7=1X15= we should get 7=2I11=1D3=
      ReadAlignment(2, "7=1X15=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 11),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3)
      }
    },

    // This test verifies a correct merging of two INSertions at the same
    // position. Specifically we verify that after merging two INSes we
    // correctly merge operations that follow. More details can be found in
    // internal.
    {
      "ins_ins_merge",
      {"CGGATCATGTTTTTTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 9=2I11=1D24=
      {"GATCATGTTTTTTTGGGTTTTCAG"},  // "7=1I16="
      // Read to haplotype has one 1 base ins at 7. Read aligns to haplotype
      // from position 2.
      // After trimming haplotype to reference cigar we get: 7=1I34=
      // Merging 7=2I11=1D24= and 7=1I16= we should get 7=3I11=1D3=
      ReadAlignment(2, "7=1I16=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 3),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 11),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3)
      }
    },

    // This test verifies a correct merging of two DELs at the same position.
    // Specifically we verify that after merging two dels we correctly merge
    // operations that follow. More details can be found in internal.
    {
      "del_del_merge",
      {"CGGATCATGTTTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 9=1D10=1D24=
      {"GATCATGTTGGGTTTTCAGGACAAA"},  // "7=1D18="
      // Read to haplotype has one del at 7. Read aligns to haplotype
      // from position 2.
      // After trimming haplotype to reference cigar we get: 7=1D34=
      // Merging 7=1D10=1D24= and 7=1D18= we should get 7=2D9=1D9=
      ReadAlignment(2, "7=1D18=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 2),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 9),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 9)
      }
    },

    // This test verifies a correct merging of DEL and INS at the same position.
    {
      "del_ins_merge",
      {"CGGATCATGTTTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 9=1D10=1D24=
      {"GATCATGTTTTTGGGTTTTCAGGACAAA"},  // "7=2I19="
      ReadAlignment(2, "7=2I19=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 11),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 9)
      }
    },

    // This test verifies a correct merging of DEL and INS at the same position.
    {
      "del_ins_merge2",
      {"CGGATCATGTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 9=3D8=1D24=
      {"GATCATGTTTGGGTTTTCAGGACAAA"},  // "7=2I17="
      ReadAlignment(2, "7=2I17=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 10),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 9)
      }
    },

    // This test verifies a correct merging of INS and DEL at the same position.
    // The read in this test is identical to the read in del_ins_merge test,
    // therefore read to ref alignments should be identical as well.
    {
      "ins_del_merge",
      {"CGGATCATGTTTTTTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 9=2I11=1D24=
      {"GATCATGTTTTTGGGTTTTCAGGACAAA"},  // "7=1D21="
      ReadAlignment(2, "7=1D21=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 11),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 9)
      }
    },

    // This test verifies a correct merging of INS and DEL at the same position.
    // The read in this test is identical to the read in del_ins_merge test,
    // therefore read to ref alignments should be identical as well.
    {
      "2ins_3del_merge",
      {"CGGATCATGTTTTTTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 9=2I11=1D24=
      {"GATCATGTTTGGGTTTTCAGGACAAA"},  // "7=3D19="
      ReadAlignment(2, "7=3D19=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 10),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 9)
      }
    },

    {
      "1ins_1del_back_to_back",
      {"CGGATCATGTTTTGGGTTTTCAGGACAAAGTATGGTTGAAACTG"},  // 20=1D24=
      {"GATCATGTTTTGGGTTTTCCAGGACAAA"},  // "18=1I9="
      ReadAlignment(2, "18=1I9=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 28)
      }
    },

    {
      "1ins_1del_consecutive",
      {"CGGATCATGTTTTGGGTTTTTTGCAGGACAAAGTATGGTTGAAACTG"},  // 16=2I29=
      {"GATCATGTTTTGGGTTTTGCAGGACAAA"},  // "16=2D12="
      ReadAlignment(2, "16=2D12=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 28)
      }
    },

    {
      "1del_1ins_consecutive2",
      {"CGGATCATGTTTTGGGTTTTGCGCAGGACAAAGTATGGTTGAAACTG"},  // 20=2I25=
      {"GATCATGTTTTGGGTTGCGCAGGACAAA"},  // "16=2D12="
      ReadAlignment(2, "16=2D12=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 28)
      }
    },

    // Haplotype alignment to ref has one del.
    // Read alignment to haplotype has one del.
    {
      "two_dels_different_positions",
      // reference with 1 del at 9.
      // haplotype to ref: 9=1D34=
      {"CGGATCATGTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"},
      {"GATCATGTTTGGTTTT"},  // "10=1D6="
      ReadAlignment(2, "10=1D6=", kSomeScore),
      {
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 6)
      }
    }
  };

  for (auto& test_case : test_cases) {
    aligner_.set_haplotypes(test_case.haplotypes);
    aligner_.AlignHaplotypesToReference();
    const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();
    aligner_.set_reads(test_case.reads);
    std::list<CigarOp> read_to_ref_cigar_ops;
    std::cout << "-----------------------------------------------" << std::endl;
    std::cout << "Testing " << test_case.test_case_name << " case" << std::endl;
    aligner_.CalculateReadToRefAlignment(
        0, test_case.read_alignment,
        haplotype_alignments[0].cigar_ops, &read_to_ref_cigar_ops);
    EXPECT_THAT(read_to_ref_cigar_ops,
                testing::ElementsAreArray(test_case.expected_cigar));
  }
}


// Haplotype alignment to ref has one del.
// Read alignment to happlotype has one del. When merged both DELs happen at
// the same position. The test verifies that DELs are properly merged in a
// single DEL.
TEST_F(FastPassAlignerTest, CalculateReadToRefAlignment_MergedDels_Test) {
  aligner_.InitSswLib();
  aligner_.set_reference(
      "CTCTGTAATCGGATCATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
      "CTGAAGATATG");
  std::vector<std::string> haplotypes = {
      // reference with 1 del at 9.
      // haplotype to ref: 9=1D34=
      "CGGATCATGTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"};

  // Reference is set in init.
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();

  // Read to haplotype has one del at the same position as in haplotype to
  // reference alignment.
  // Read aligns to haplotype from position 2.
  // After trimming haplotype to reference cigar we get: 7=1D34=
  // Merging 7=1D34= and 7=1D9= we should get 7=2D9=
  aligner_.set_reads({
      "GATCATGTTGGGTTTT",  // "7=1D9="
  });

  std::list<CigarOp> read_to_ref_cigar_ops;
  aligner_.CalculateReadToRefAlignment(
      0, ReadAlignment(2, "7=1D9=", kSomeScore),
      haplotype_alignments[0].cigar_ops, &read_to_ref_cigar_ops);

  std::list<CigarOp> expected_read_to_ref_cigar_ops = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
      CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 2),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 9)};
}

// Haplotype alignment to ref has one INS.
// Read alignment to happlotype has one INS.
// Insertions are located at different positions. Test verifies that two
// different INSertions are created  in the merged alignment.
TEST_F(FastPassAlignerTest, CalculateReadToRefAlignment_Ins_Test) {
  aligner_.InitSswLib();
  aligner_.set_reference(
      "CTCTGTAATCGGATCATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
      "CTGAAGATATG");
  std::vector<std::string> haplotypes = {
      // reference with 1 ins at 12 (AA).
      // haplotype to ref: 13=2I32=
      "CGGATCATGTTTTAAGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"};

  // Reference is set in init.
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();

  // Read to haplotype has one INS at 16 (CC). Read aligns to haplotype from
  // position 2.
  // After trimming haplotype to reference cigar we get: 11=2I32=
  // Merging 11=2I32= and 16=2I4= we should get 11=2I3=2I4=
  aligner_.set_reads({
      "GATCATGTTTTAAGGGCCTTTT",  // "16=2I4="
  });

  std::list<CigarOp> read_to_ref_cigar_ops;
  aligner_.CalculateReadToRefAlignment(
      0, ReadAlignment(2, "16=2I4=", kSomeScore),
      haplotype_alignments[0].cigar_ops, &read_to_ref_cigar_ops);

  std::list<CigarOp> expected_read_to_ref_cigar_ops = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 11),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 4),
  };

  EXPECT_THAT(read_to_ref_cigar_ops,
              testing::ElementsAreArray(expected_read_to_ref_cigar_ops));
}

// Haplotype alignment to ref has one INS.
// Read alignment to happlotype has one INS.
// Insertions are located at the same position. Test verifies that two
// different INSertions are merged into one.
TEST_F(FastPassAlignerTest, CalculateReadToRefAlignment_MergedIns_Test) {
  aligner_.InitSswLib();
  aligner_.set_reference(
      "CTCTGTAATCGGATCATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
      "CTGAAGATATG");
  std::vector<std::string> haplotypes = {
      // reference with 1 ins at 12 (AA).
      // haplotype to ref: 13=2I32=
      "CGGATCATGTTTTAAGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"};

  // Reference is set in init.
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();

  // Read to haplotype has one INS at 13 (CC). Read aligns to haplotype from
  // position 2.
  // After trimming haplotype to reference cigar we get: 11=2I32=
  // Merging 11=2I32= and 16=2I4= we should get 11=4I7=
  aligner_.set_reads({
      "GATCATGTTTTAAAAGGGTTTT",  // "13=2I7="
  });

  std::list<CigarOp> read_to_ref_cigar_ops;
  aligner_.CalculateReadToRefAlignment(
      0, ReadAlignment(2, "13=2I7=", kSomeScore),
      haplotype_alignments[0].cigar_ops, &read_to_ref_cigar_ops);

  std::list<CigarOp> expected_read_to_ref_cigar_ops = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 11),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 4),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7)};

  EXPECT_THAT(read_to_ref_cigar_ops,
              testing::ElementsAreArray(expected_read_to_ref_cigar_ops));
}

// Haplotype alignment to ref has one DEL.
// Read alignment to happlotype has one INS.
// Test verifies that this type of merge is handeled correctly.
TEST_F(FastPassAlignerTest, CalculateReadToRefAlignment_DelIns_Test) {
  aligner_.InitSswLib();
  aligner_.set_reference(
      "CTCTGTAATCGGATCATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
      "CTGAAGATATG");
  std::vector<std::string> haplotypes = {
      // reference with 1 del at 9 (missing T).
      // haplotype to ref: 9=1D34=
      "CGGATCATGTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"};

  // Reference is set in init.
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();

  // Read to haplotype has one del at the same position as in haplotype to
  // reference alignment.
  // Read aligns to haplotype from position 2.
  // After trimming haplotype to reference cigar we get: 7=1D34=
  // Merging 7=1D34= and 13=2I4= we should get 7=1D6=2I4=
  aligner_.set_reads({
      "GATCATGTTTGGGAATTTT",  // "13=2I4="
  });

  std::list<CigarOp> read_to_ref_cigar_ops;
  aligner_.CalculateReadToRefAlignment(
      0, ReadAlignment(2, "13=2I4=", kSomeScore),
      haplotype_alignments[0].cigar_ops, &read_to_ref_cigar_ops);

  std::list<CigarOp> expected_read_to_ref_cigar_ops = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
      CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 6),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 4),
  };
  EXPECT_THAT(read_to_ref_cigar_ops,
              testing::ElementsAreArray(expected_read_to_ref_cigar_ops));
}

// Haplotype alignment to ref has one INS.
// Read alignment to happlotype has one DEL.
// Test verifies that this type of merge is handeled correctly.
TEST_F(FastPassAlignerTest, CalculateReadToRefAlignment_InsDel_Test) {
  aligner_.InitSswLib();
  aligner_.set_reference(
      "CTCTGTAATCGGATCATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
      "CTGAAGATATG");
  std::vector<std::string> haplotypes = {
      // reference with 1 ins at 9 (AA).
      // haplotype to ref: 9=2I36=
      "CGGATCATGAATTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTG"};

  // Reference is set in init.
  aligner_.set_haplotypes(haplotypes);
  aligner_.AlignHaplotypesToReference();
  const auto& haplotype_alignments = aligner_.GetReadToHaplotypeAlignments();

  // Read to haplotype has one del at 16 (missing T)
  // Read aligns to haplotype from position 2.
  // After trimming haplotype to reference cigar we get: 7=2I36=
  // Merging 7=2I36= and 16=1D3= we should get 7=2I11=1D3=
  aligner_.set_reads({
      "GATCATGAATTTTGGGTTT",  // "16=1D3="
  });

  std::list<CigarOp> read_to_ref_cigar_ops;
  aligner_.CalculateReadToRefAlignment(
      0, ReadAlignment(2, "16=1D3=", kSomeScore),
      haplotype_alignments[0].cigar_ops, &read_to_ref_cigar_ops);

  std::list<CigarOp> expected_read_to_ref_cigar_ops = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 7),
      CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 1),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
  };
  EXPECT_THAT(read_to_ref_cigar_ops,
              testing::ElementsAreArray(expected_read_to_ref_cigar_ops));
}

TEST_F(FastPassAlignerTest, MergeCigarOp_emtyCigar_Test) {
  std::list<CigarOp> cigar;
  MergeCigarOp(CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
               10, &cigar);
  EXPECT_THAT(cigar,
              testing::ElementsAreArray({CigarOp(
                  nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3)}));
}

TEST_F(FastPassAlignerTest, MergeCigarOp_mergeDifferentOp_Test) {
  std::list<CigarOp> cigar = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 5)};
  MergeCigarOp(CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 2),
               10, &cigar);
  EXPECT_THAT(
      cigar,
      testing::ElementsAreArray(
          {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
           CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 5),
           CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 2)}));
}

TEST_F(FastPassAlignerTest, MergeCigarOp_mergeSameOp_Test) {
  std::list<CigarOp> cigar = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 5)};
  MergeCigarOp(CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2), 10,
               &cigar);
  EXPECT_THAT(
      cigar, testing::ElementsAreArray(
                 {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
                  CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 7)}));
}

TEST_F(FastPassAlignerTest, MergeCigarOp_alignedLengthOverflow_Test) {
  std::list<CigarOp> cigar = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
      CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 5)};
  MergeCigarOp(CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 20), 10,
               &cigar);
  EXPECT_THAT(
      cigar, testing::ElementsAreArray(
                 {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
                  CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 7)}));
}

// DEL does not count towards aligned length. This test verifies that DEL can
// be merged doesn't matter what it's length is.
TEST_F(FastPassAlignerTest, MergeCigarOp_alignedLengthOverflowDel_Test) {
  std::list<CigarOp> cigar = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
      CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 5)};
  MergeCigarOp(CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 20), 10,
               &cigar);
  EXPECT_THAT(
      cigar, testing::ElementsAreArray(
                 {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
                  CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 25)}));
}

// Try to merge an operation when alined read is already equals read length.
// The operation should not merge in this case.
TEST_F(FastPassAlignerTest, MergeCigarOp_alignedLengthOverflow2_Test) {
  std::list<CigarOp> cigar = {
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
      CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 5),
      CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 5)};
  MergeCigarOp(CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 20), 8,
               &cigar);
  EXPECT_THAT(
      cigar,
      testing::ElementsAreArray(
          {CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 3),
           CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 5),
           CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 5)}));
}

// Test that ssw_alignment_score_threshold does not go negative if similarity
// threshold is less than 0.5.
TEST_F(FastPassAlignerTest, CalculateSswAlignmentScoreThreshold_Test) {
  const int read_size = 10;
  AlignerOptions aligner_options;
  aligner_options.set_read_size(read_size);
  aligner_options.set_realignment_similarity_threshold(0.1);
  aligner_.set_options(aligner_options);
  aligner_.CalculateSswAlignmentScoreThreshold();
  EXPECT_GE(aligner_.get_ssw_alignment_score_threshold(), 0);
  EXPECT_LE(aligner_.get_ssw_alignment_score_threshold(),
            read_size * aligner_.get_match_score());
}

// This test verifies a situation described in internal. Sometime haplotype
// generated by DeBruijn graph cannot be recreated with reads. In that case we
// want to skip that haplotype and do not realign reads to it.
// If reads cannot be aligned to haplotype within ref_prefix or ref_suffix then
// it doesn't govern haplotype invalidation (this is tested with read 6)
TEST_F(FastPassAlignerTest, HaplotypeHasZeroCoverageOutsideInterval_Test) {
  // Read 1: exactly matches haplotype
  // Read 2: exactly matches haplotype
  // Read 3: matches haplotype with 2 mismatches
  // Read 4: matches haplotype with 1 mismatch
  // Read 5: matches haplotype exactly
  // Read 6: matches haplotype with 4 bases del
  //
  // 0         1         2         3         4         5         6         7
  // 012345678901234567890123456789012345678901234567890123456789012345678901
  //
  // ATCAAGGGAAAAAGTGCCCAGGGCCAAATATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAGCT
  // ATCAAGGGAAAAAGTGCCCA
  //                     GGGCCAAATATGTTTTG
  //                            ATATGTTATGGGTTATGCAGGA
  //                                GTTTTGGGTTTTGCAGGTCA
  //                                              AGGACAAAGTATGGTT
  //                                                  CAAAGTATGGTTG---TGAGCT
  aligner_.set_reads({
      "ATCAAGGGAAAAAGTGCCCA",    // "20=", 80
      "GGGCCAAATATGTTTTG",       // "17=", 68
      "ATATGTTATGGGTTATGCAGGA",  // "7=1X6=1X7=", 68
      "GTTTTGGGTTTTGCAGGTCA",    // "17=1X2=", 51
      "AGGACAAAGTATGGTT",        // "16=", 64
      "CAAAGTATGGTTGTGAGCT"      // "13=4D6=", 65
  });
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_.set_ref_prefix_len(11);
  aligner_.set_ref_suffix_len(11);
  aligner_options.set_max_num_of_mismatches(2);
  aligner_.set_options(aligner_options);
  aligner_.BuildIndex();

  // reference, with ref_prefix_len = 11, and ref_suffix_len = 11
  string haplotype =
      "ATCAAGGGAAAAAGTGCCCAGGGCCAAATATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAGC"
      "T";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  std::vector<std::string> aligner_reads = aligner_.get_reads();
  int expected_hap_score = 350;

  std::vector<ReadAlignment> expected_read_scores(aligner_reads.size());
  expected_read_scores[0] = ReadAlignment(0, "20=", 80);
  expected_read_scores[1] = ReadAlignment(20, "17=", 68);
  expected_read_scores[2] = ReadAlignment(27, "22=", 68);
  expected_read_scores[3] = ReadAlignment(31, "20=", 70);
  expected_read_scores[4] = ReadAlignment(45, "16=", 64);
  expected_read_scores[5] = ReadAlignment(ReadAlignment::kNotAligned, "", 0);

  aligner_.FastAlignReadsToHaplotype(haplotype, &haplotype_score, &read_scores);
  EXPECT_EQ(expected_hap_score, haplotype_score);
  EXPECT_THAT(read_scores,
              testing::UnorderedElementsAreArray(expected_read_scores));
}

// In this test haplotype has coverage zero starting at position 20.
// Make sure that haplotype is discarded.
TEST_F(FastPassAlignerTest, HaplotypeHasZeroCoverageInsideInterval_Test) {
  // Read 1: exactly matches haplotype
  // Read 2: has 'AAA' insertion relative to the haplotype
  // Read 3: matches haplotype with 2 mismatches
  // Read 4: matches haplotype with 1 mismatch
  // Read 5: matches haplotype exactly
  // Read 6: matches haplotype with 4 bases del
  //
  // 0         1         2         3         4         5         6         7
  // 01234567890123456789012---3456789012345678901234567890123456789012345678901
  //
  // ATCAAGGGAAAAAGTGCCCAGGG---CCAAATATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAGCT
  // ATCAAGGGAAAAAGTGCCCA
  //                     GGGAAACCAAATATGTTTTG
  //                            ATATGTTATGGGTTATGCAGGA
  //                                GTTTTGGGTTTTGCAGGTCA
  //                                              AGGACAAAGTATGGTT
  //                                                  CAAAGTATGGTTG---TGAGCT
  aligner_.set_reads({
      "ATCAAGGGAAAAAGTGCCCA",    // "20=", 80
      "GGGAAACCAAATATGTTTTG",    // "3=3I14=", 68
      "ATATGTTATGGGTTATGCAGGA",  // "7=1X6=1X7=", 68
      "GTTTTGGGTTTTGCAGGTCA",    // "17=1X2=", 51
      "AGGACAAAGTATGGTT",        // "16=", 64
      "CAAAGTATGGTTGTGAGCT"      // "13=4D6=", 65
  });
  AlignerOptions aligner_options;
  aligner_options.set_kmer_size(3);
  aligner_.set_ref_prefix_len(11);
  aligner_.set_ref_suffix_len(11);
  aligner_options.set_max_num_of_mismatches(2);
  aligner_.set_options(aligner_options);
  aligner_.BuildIndex();

  // reference, with ref_prefix_len = 11, and ref_suffix_len = 11
  string haplotype =
      "ATCAAGGGAAAAAGTGCCCAGGGCCAAATATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAGC"
      "T";
  int haplotype_score = 0;
  std::vector<ReadAlignment> read_scores(aligner_.get_reads().size());

  // Expected values
  int expected_hap_score = 0;

  std::vector<ReadAlignment> expected_read_scores(aligner_.get_reads().size());
  expected_read_scores[0] = ReadAlignment(0, "20=", 80);
  expected_read_scores[1] = ReadAlignment(ReadAlignment::kNotAligned, "", 0);
  expected_read_scores[2] = ReadAlignment(ReadAlignment::kNotAligned, "", 0);
  expected_read_scores[3] = ReadAlignment(ReadAlignment::kNotAligned, "", 0);
  expected_read_scores[4] = ReadAlignment(ReadAlignment::kNotAligned, "", 0);
  expected_read_scores[5] = ReadAlignment(ReadAlignment::kNotAligned, "", 0);

  aligner_.FastAlignReadsToHaplotype(haplotype, &haplotype_score, &read_scores);
  EXPECT_EQ(expected_hap_score, haplotype_score);
  EXPECT_THAT(read_scores,
              testing::UnorderedElementsAreArray(expected_read_scores));
}

TEST_F(FastPassAlignerTest, IsAlignmentNormalized_Match) {
  aligner_.set_reference("ATGCTGCACTCTCTCTCAGCTGTCACC");
  EXPECT_TRUE(aligner_.IsAlignmentNormalized(
      /*cigar=*/{
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 14)
      },
      /*ref_offset=*/7,
      /*read_sequence=*/"ACTCTCTCTCAGCT"
      ));
}

TEST_F(FastPassAlignerTest, IsAlignmentNormalized_Del_Not_Norm) {
  aligner_.set_reference("ATGCTGCACTCTCTCTCAGCTGTCACC");
  EXPECT_FALSE(aligner_.IsAlignmentNormalized(
      /*cigar=*/{
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 4),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 2),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 10),
      },
      /*ref_offset=*/7,
      /*read_sequence=*/"ACTCTCTCAGCTGT"
      ));
}

TEST_F(FastPassAlignerTest, IsAlignmentNormalized_Del_Norm) {
  aligner_.set_reference("ATGCTGCACTCTCTCTCAGCTGTCACC");
  EXPECT_TRUE(aligner_.IsAlignmentNormalized(
      /*cigar=*/{
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::DELETE, 2),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 13),
      },
      /*ref_offset=*/7,
      /*read_sequence=*/"ACTCTCTCAGCTGT"
      ));
}

TEST_F(FastPassAlignerTest, IsAlignmentNormalized_Ins_Not_Norm) {
  aligner_.set_reference("ATGCTGCACTCTCTCTCAGCTGTCACC");
  EXPECT_FALSE(aligner_.IsAlignmentNormalized(
      /*cigar=*/{
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 4),
        CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 8),
      },
      /*ref_offset=*/7,
      /*read_sequence=*/"ACTCTCTCTCTCAGCTGT"
      ));
}

TEST_F(FastPassAlignerTest, IsAlignmentNormalized_Ins_Norm) {
  aligner_.set_reference("ATGCTGCACTCTCTCTCAGCTGTCACC");
  EXPECT_TRUE(aligner_.IsAlignmentNormalized(
      /*cigar=*/{
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 1),
        CigarOp(nucleus::genomics::v1::CigarUnit::INSERT, 2),
        CigarOp(nucleus::genomics::v1::CigarUnit::ALIGNMENT_MATCH, 10),
      },
      /*ref_offset=*/7,
      /*read_sequence=*/"ACTCTCTCTCTCAGCTGT"
      ));
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
