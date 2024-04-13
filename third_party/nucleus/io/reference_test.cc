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
 *
 */

#include "third_party/nucleus/io/reference.h"

#include <memory>
#include <string>
#include <utility>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "absl/log/check.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "third_party/nucleus/io/reader_base.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/fasta.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "tensorflow/core/lib/core/status.h"

using absl::StrCat;
using std::make_pair;

namespace nucleus {

using ::testing::IsEmpty;
using ::testing::Not;
using ::testing::UnorderedElementsAre;

string TestFastaPath() { return GetTestData("test.fasta"); }

typedef std::unique_ptr<GenomeReference> CreateGenomeReferenceFunc(
    const string& fasta_path, int cache_size);

// Tests are parameterized by: reader factory, cache size.
class GenomeReferenceTest : public ::testing::TestWithParam<
    std::pair<CreateGenomeReferenceFunc*, int>> {
 protected:
  void SetUp() override {
    ref_ = (*GetParam().first)(TestFastaPath(), GetParam().second);
  }
  const GenomeReference& Ref() const { return *ref_; }

  std::unique_ptr<const GenomeReference> ref_;
};

TEST_P(GenomeReferenceTest, TestBasic) {
  EXPECT_THAT(Ref().ContigNames(),
              UnorderedElementsAre("chrM", "chr1", "chr2"));
  EXPECT_THAT(Ref().Contigs().size(), 3);

  const auto& chrm = *Ref().Contig("chrM").ValueOrDie();
  EXPECT_EQ(100, chrm.n_bases());
  EXPECT_EQ("chrM", chrm.name());
  EXPECT_EQ(0, chrm.pos_in_fasta());

  const auto& chr1 = *Ref().Contig("chr1").ValueOrDie();
  EXPECT_EQ(76, chr1.n_bases());
  EXPECT_EQ("chr1", chr1.name());
  EXPECT_EQ(1, chr1.pos_in_fasta());

  const auto& chr2 = *Ref().Contig("chr2").ValueOrDie();
  EXPECT_EQ(121, chr2.n_bases());
  EXPECT_EQ("chr2", chr2.name());
  EXPECT_EQ(2, chr2.pos_in_fasta());
}


TEST_P(GenomeReferenceTest, TestIsValidInterval) {
  // Checks that we can check that an unknown chromosome isn't valid.
  EXPECT_FALSE(Ref().IsValidInterval(MakeRange("unknown_chr", 0, 1)));

  for (const auto& chr : Ref().ContigNames()) {
    const auto n_bases = Ref().Contig(chr).ValueOrDie()->n_bases();

    EXPECT_TRUE(Ref().IsValidInterval(MakeRange(chr, 0, n_bases)));
    for (int i = 0; i < n_bases; ++i) {
      EXPECT_TRUE(Ref().IsValidInterval(MakeRange(chr, 0, i+1)));
      EXPECT_TRUE(Ref().IsValidInterval(MakeRange(chr, i, i+1)));
    }
    EXPECT_FALSE(Ref().IsValidInterval(MakeRange(chr, -10, 0)));
    EXPECT_FALSE(Ref().IsValidInterval(MakeRange(chr, -1, 0)));
    EXPECT_FALSE(Ref().IsValidInterval(MakeRange(chr, 10, 9)));
    EXPECT_FALSE(Ref().IsValidInterval(MakeRange(chr, 0, n_bases + 1)));
    EXPECT_FALSE(Ref().IsValidInterval(MakeRange(chr, 0, n_bases + 100)));
    EXPECT_FALSE(Ref().IsValidInterval(MakeRange(chr, n_bases, n_bases)));
    EXPECT_FALSE(
        Ref().IsValidInterval(MakeRange(chr, n_bases + 100, n_bases + 100)));
  }
}

TEST_P(GenomeReferenceTest, NotOKIfContigCalledWithBadName) {
  EXPECT_THAT(Ref().Contig("missing"),
              IsNotOKWithMessage("Unknown contig missing"));
}

TEST_P(GenomeReferenceTest, NotOKIfIntervalIsInvalid) {
  // Asking for bad chromosome values produces death.
  StatusOr<string> result = Ref().GetBases(MakeRange("missing", 0, 1));
  EXPECT_THAT(result,
              IsNotOKWithCodeAndMessage(static_cast<absl::StatusCode>(
                                            absl::StatusCode::kInvalidArgument),
                                        "Invalid interval"));

  // Starting before 0 is detected.
  EXPECT_THAT(Ref().GetBases(MakeRange("chrM", -1, 1)),
              IsNotOKWithMessage("Invalid interval"));

  // chr1 exists, but this range's start is beyond the chr.
  EXPECT_THAT(Ref().GetBases(MakeRange("chr1", 1000, 1010)),
              IsNotOKWithMessage("Invalid interval"));

  // chr1 exists, but this range's end is beyond the chr.
  EXPECT_THAT(Ref().GetBases(MakeRange("chr1", 0, 1010)),
              IsNotOKWithMessage("Invalid interval"));
}

TEST_P(GenomeReferenceTest, TestHasContig) {
  EXPECT_TRUE(Ref().HasContig("chrM"));
  EXPECT_TRUE(Ref().HasContig("chr1"));
  EXPECT_TRUE(Ref().HasContig("chr2"));
  EXPECT_FALSE(Ref().HasContig("chr3"));
  EXPECT_FALSE(Ref().HasContig("chr"));
  EXPECT_FALSE(Ref().HasContig(""));
}

// Checks that GetBases work in all its forms for the given arguments.
void CheckGetBases(const GenomeReference& ref,
                   const string& chrom, const int64 start, const int64 end,
                   const string& expected_bases) {
  StatusOr<string> query = ref.GetBases(MakeRange(chrom, start, end));
  ASSERT_THAT(query, IsOK());
  EXPECT_THAT(query.ValueOrDie(), expected_bases);
}


TEST_P(GenomeReferenceTest, TestReferenceBases) {
  CheckGetBases(Ref(), "chrM", 0, 100,
                "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTC"
                "GTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTG");

  CheckGetBases(Ref(), "chr1", 0, 76,
                "ACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCTACTCTCCTAAATCCCTTCT"
                "CGTCCCCATGGATGA");

  CheckGetBases(Ref(), "chr2", 0, 121,
                "CGCTNCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGAC"
                "ATCTGGTTCCTACTTCAGGGCCATAAAGCCTAAATAGCCCACACGTTCCC"
                "CTTAAATAAGACATCACGATG");
}


TEST_P(GenomeReferenceTest, TestGetBasesParts) {
  CheckGetBases(Ref(), "chrM", 0, 10, "GATCACAGGT");
  CheckGetBases(Ref(), "chrM", 0, 9, "GATCACAGG");
  CheckGetBases(Ref(), "chrM", 1, 9, "ATCACAGG");
  CheckGetBases(Ref(), "chrM", 3, 7, "CACA");
  CheckGetBases(Ref(), "chrM", 90, 100, "CGAGACGCTG");
  CheckGetBases(Ref(), "chrM", 90, 99, "CGAGACGCT");
  CheckGetBases(Ref(), "chrM", 91, 100, "GAGACGCTG");
  CheckGetBases(Ref(), "chrM", 92, 100, "AGACGCTG");
  CheckGetBases(Ref(), "chrM", 92, 99, "AGACGCT");
  CheckGetBases(Ref(), "chrM", 92, 98, "AGACGC");

  CheckGetBases(Ref(), "chrM", 0, 1, "G");
  CheckGetBases(Ref(), "chrM", 1, 2, "A");
  CheckGetBases(Ref(), "chrM", 2, 3, "T");
  CheckGetBases(Ref(), "chrM", 3, 4, "C");
  CheckGetBases(Ref(), "chrM", 4, 5, "A");
  CheckGetBases(Ref(), "chrM", 5, 6, "C");

  // crosses the boundary of the index when max_bin_size is 5
  CheckGetBases(Ref(), "chrM", 4, 6, "AC");

  // 0-bp interval requests should return the empty string.
  CheckGetBases(Ref(), "chrM", 0, 0, "");
  CheckGetBases(Ref(), "chrM", 10, 10, "");
}


static std::unique_ptr<GenomeReference> LoadWithCaseOption(
    const string& fasta, bool keep_true_case,
    int cache_size = 64 * 1024) {
  nucleus::genomics::v1::FastaReaderOptions options =
      nucleus::genomics::v1::FastaReaderOptions();
  options.set_keep_true_case(keep_true_case);
  StatusOr<std::unique_ptr<IndexedFastaReader>> fai_status =
      IndexedFastaReader::FromFile(fasta, StrCat(fasta, ".fai"),
                                   options,
                                   cache_size);
  NUCLEUS_CHECK_OK(fai_status.status());
  return std::move(fai_status.ValueOrDie());
}

static std::unique_ptr<GenomeReference> JustLoadFai(const string& fasta,
                                                    int cache_size = 64 *
                                                                     1024) {
  return LoadWithCaseOption(fasta, false, cache_size);
}

// Test with cache disabled.
INSTANTIATE_TEST_CASE_P(GRT1, GenomeReferenceTest,
                        ::testing::Values(make_pair(&JustLoadFai, 0)));

// Test with a large cache.
INSTANTIATE_TEST_CASE_P(GRT3, GenomeReferenceTest,
                        ::testing::Values(make_pair(&JustLoadFai, 64 * 1024)));

TEST(StatusOrLoadFromFile, ReturnsBadStatusIfFaiIsMissing) {
  StatusOr<std::unique_ptr<IndexedFastaReader>> result =
      IndexedFastaReader::FromFile(GetTestData("unindexed.fasta"),
                                   GetTestData("unindexed.fasta.fai"),
                                   nucleus::genomics::v1::FastaReaderOptions());
  EXPECT_THAT(result, IsNotOKWithCodeAndMessage(
                          absl::StatusCode::kNotFound,
                          "could not load fasta and/or fai for fasta"));
}

TEST(IndexedFastaReaderTest, WriteAfterCloseIsntOK) {
  auto reader = JustLoadFai(TestFastaPath());
  ASSERT_THAT(reader->Close(), IsOK());
  EXPECT_THAT(reader->GetBases(MakeRange("chrM", 0, 100)),
              IsNotOKWithCodeAndMessage(
                  absl::StatusCode::kFailedPrecondition,
                  "can't read from closed IndexedFastaReader object"));
}

TEST(IndexedFastaReaderTest, TestTrueCase) {
  auto reader = LoadWithCaseOption(TestFastaPath(), true);
  auto iterator = reader->Iterate().ValueOrDie();
  GenomeReferenceRecord r;
  StatusOr<bool> status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("chrM", r.first);
  EXPECT_EQ(
      "GATCACAGGTCTATCACCCTATTaaCCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGT"
      "GTGCACGCGATAGCATTGCGAGACGCTG",
      r.second);
}

TEST(IndexedFastaReaderTest, TestIterate) {
  auto reader = JustLoadFai(TestFastaPath());
  auto iterator = reader->Iterate().ValueOrDie();
  GenomeReferenceRecord r;
  StatusOr<bool> status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("chrM", r.first);
  EXPECT_EQ(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGT"
      "GTGCACGCGATAGCATTGCGAGACGCTG",
      r.second);
  status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("chr1", r.first);
  EXPECT_EQ(
      "ACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCTACTCTCCTAAATCCCTTCTCGTCCCCATGG"
      "ATGA",
      r.second);
  status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("chr2", r.first);
  EXPECT_EQ(
      "CGCTNCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCC"
      "ATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG",
      r.second);

  // Reading beyond the file fails.
  status = iterator->Next(&r);
  EXPECT_FALSE(status.ValueOrDie());
}

TEST(UnindexedFastaReaderTest, ReturnsBadStatusIfFileIsMissing) {
  StatusOr<std::unique_ptr<UnindexedFastaReader>> result =
      UnindexedFastaReader::FromFile(GetTestData("nonexistent.fasta"));
  EXPECT_THAT(result, IsNotOKWithCodeAndMessage(absl::StatusCode::kNotFound,
                                                "Could not open"));
}

TEST(UnindexedFastaReaderTest, IterateAfterCloseIsntOK) {
  StatusOr<std::unique_ptr<UnindexedFastaReader>> result =
      UnindexedFastaReader::FromFile(GetTestData("unindexed.fasta"));
  auto reader = std::move(result.ValueOrDie());
  auto iterator = reader->Iterate().ValueOrDie();
  ASSERT_THAT(reader->Close(), IsOK());
  GenomeReferenceRecord r;
  StatusOr<bool> status = iterator->Next(&r);
  EXPECT_THAT(iterator->Next(&r),
              IsNotOKWithCodeAndMessage(
                  absl::StatusCode::kFailedPrecondition,
                  "Cannot iterate a closed UnindexedFastaReader"));
}

TEST(UnindexedFastaReaderTest, TestMalformed) {
  StatusOr<std::unique_ptr<UnindexedFastaReader>> result =
      UnindexedFastaReader::FromFile(GetTestData("malformed.fasta"));
  auto reader = std::move(result.ValueOrDie());
  auto iterator = reader->Iterate().ValueOrDie();
  GenomeReferenceRecord r;
  EXPECT_THAT(iterator->Next(&r),
              IsNotOKWithCodeAndMessage(absl::StatusCode::kDataLoss,
                                        "Name not found in FASTA"));
}

class UnindexedFastaReaderFileTest : public ::testing::TestWithParam<string> {};

// Test a couple of files that are formatted differently but should have the
// same contents.
INSTANTIATE_TEST_CASE_P(All, UnindexedFastaReaderFileTest,
                        ::testing::Values("unindexed.fasta", "test.fasta.gz",
                                          "unindexed_emptylines.fasta"));

TEST_P(UnindexedFastaReaderFileTest, TestIterate) {
  LOG(INFO) << "testing file " << GetParam();
  StatusOr<std::unique_ptr<UnindexedFastaReader>> result =
      UnindexedFastaReader::FromFile(GetTestData(GetParam()));
  auto reader = std::move(result.ValueOrDie());
  auto iterator = reader->Iterate().ValueOrDie();
  GenomeReferenceRecord r1;
  StatusOr<bool> status = iterator->Next(&r1);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("chrM", r1.first);
  EXPECT_EQ(
      "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGT"
      "GTGCACGCGATAGCATTGCGAGACGCTG",
      r1.second);

  GenomeReferenceRecord r2;
  status = iterator->Next(&r2);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("chr1", r2.first);
  EXPECT_EQ(
      "ACCACCATCCTCCGTGAAATCAATATCCCGCACAAGAGTGCTACTCTCCTAAATCCCTTCTCGTCCCCATGG"
      "ATGA",
      r2.second);
  GenomeReferenceRecord r3;
  status = iterator->Next(&r3);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("chr2", r3.first);
  EXPECT_EQ(
      "CGCTNCGGGCCCATAACACTTGGGGGTAGCTAAAGTGAACTGTATCCGACATCTGGTTCCTACTTCAGGGCC"
      "ATAAAGCCTAAATAGCCCACACGTTCCCCTTAAATAAGACATCACGATG",
      r3.second);

  // Reading beyond the file fails.
  GenomeReferenceRecord r4;
  status = iterator->Next(&r4);
  EXPECT_FALSE(status.ValueOrDie());
}

namespace {

// Helper method to create a test sequence.
void CreateTestSeq(std::vector<genomics::v1::ContigInfo>* contigs,
                   std::vector<genomics::v1::ReferenceSequence>* seqs,
                   const string& name, const int pos_in_fasta,
                   const int range_start, const int range_end,
                   const string& bases) {
  DCHECK(pos_in_fasta >= 0 && pos_in_fasta < contigs->size());
  genomics::v1::ContigInfo* contig = &contigs->at(pos_in_fasta);
  contig->set_name(name);
  contig->set_pos_in_fasta(pos_in_fasta);
  contig->set_n_bases(range_end - range_start);
  genomics::v1::ReferenceSequence* seq = &seqs->at(pos_in_fasta);
  seq->mutable_region()->set_reference_name(name);
  seq->mutable_region()->set_start(range_start);
  seq->mutable_region()->set_end(range_end);
  seq->set_bases(bases);
}

}  // namespace

TEST(InMemoryFastaReaderTest, TestIterate) {
  int kNum = 3;
  std::vector<genomics::v1::ContigInfo> contigs(kNum);
  std::vector<genomics::v1::ReferenceSequence> seqs(kNum);
  CreateTestSeq(&contigs, &seqs, "Chr1", 0, 0, 1, "A");
  CreateTestSeq(&contigs, &seqs, "Chr2", 1, 4, 6, "CG");
  CreateTestSeq(&contigs, &seqs, "Chr3", 2, 10, 15, "AATTC");

  std::unique_ptr<InMemoryFastaReader> reader =
      std::move(InMemoryFastaReader::Create(contigs, seqs).ValueOrDie());
  auto iterator = reader->Iterate().ValueOrDie();
  GenomeReferenceRecord r;
  StatusOr<bool> status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("Chr1", r.first);
  EXPECT_EQ("A", r.second);
  status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("Chr2", r.first);
  EXPECT_EQ("CG", r.second);
  status = iterator->Next(&r);
  EXPECT_TRUE(status.ValueOrDie());
  EXPECT_EQ("Chr3", r.first);
  EXPECT_EQ("AATTC", r.second);

  // Reading beyond the file fails.
  status = iterator->Next(&r);
  EXPECT_FALSE(status.ValueOrDie());
}

}  // namespace nucleus
