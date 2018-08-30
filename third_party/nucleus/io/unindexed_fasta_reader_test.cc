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
 *
 */

#include "third_party/nucleus/io/unindexed_fasta_reader.h"

#include <utility>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/vendor/status_matchers.h"
#include "tensorflow/core/platform/logging.h"

namespace nucleus {

TEST(UnindexedFastaReaderTest, ReturnsBadStatusIfFileIsMissing) {
  StatusOr<std::unique_ptr<UnindexedFastaReader>> result =
      UnindexedFastaReader::FromFile(GetTestData("nonexistent.fasta"));
  EXPECT_THAT(result, IsNotOKWithCodeAndMessage(tensorflow::error::NOT_FOUND,
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
                  tensorflow::error::FAILED_PRECONDITION,
                  "Cannot iterate a closed UnindexedFastaReader"));
}

TEST(UnindexedFastaReaderTest, TestMalformed) {
  StatusOr<std::unique_ptr<UnindexedFastaReader>> result =
      UnindexedFastaReader::FromFile(GetTestData("malformed.fasta"));
  auto reader = std::move(result.ValueOrDie());
  auto iterator = reader->Iterate().ValueOrDie();
  GenomeReferenceRecord r;
  EXPECT_THAT(iterator->Next(&r),
              IsNotOKWithCodeAndMessage(tensorflow::error::DATA_LOSS,
                                        "Name not found in FASTA"));
}

class UnindexedFastaReaderFileTest : public ::testing::TestWithParam<string> {};

// Test a couple of files that are formatted differently but should have the
// same contents.
INSTANTIATE_TEST_CASE_P(/* prefix */, UnindexedFastaReaderFileTest,
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

}  // namespace nucleus
