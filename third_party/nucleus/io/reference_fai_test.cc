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

#include "third_party/nucleus/io/reference_fai.h"

#include <memory>
#include <utility>
#include <vector>

#include "third_party/nucleus/io/reference_test.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/vendor/status_matchers.h"

#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/logging.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>


#include "tensorflow/core/platform/test.h"

using std::make_pair;
using tensorflow::strings::StrCat;
using testing::Eq;
using testing::StartsWith;

namespace nucleus {

static std::unique_ptr<GenomeReference>
JustLoadFai(const string& fasta, int cache_size = 64 * 1024) {
  StatusOr<std::unique_ptr<GenomeReferenceFai>> fai_status =
      GenomeReferenceFai::FromFile(fasta, StrCat(fasta, ".fai"), cache_size);
  TF_CHECK_OK(fai_status.status());
  return std::move(fai_status.ValueOrDie());
}

// Test with cache disabled.
INSTANTIATE_TEST_CASE_P(GRT1, GenomeReferenceTest,
                        ::testing::Values(make_pair(&JustLoadFai, 0)));

INSTANTIATE_TEST_CASE_P(GRT2, GenomeReferenceDeathTest,
                        ::testing::Values(make_pair(&JustLoadFai, 0)));

// Test with a large cache.
INSTANTIATE_TEST_CASE_P(GRT3, GenomeReferenceTest,
                        ::testing::Values(make_pair(&JustLoadFai, 64 * 1024)));

INSTANTIATE_TEST_CASE_P(GRT4, GenomeReferenceDeathTest,
                        ::testing::Values(make_pair(&JustLoadFai, 64 * 1024)));

TEST(StatusOrLoadFromFile, ReturnsBadStatusIfFaiIsMissing) {
  StatusOr<std::unique_ptr<GenomeReferenceFai>> result =
      GenomeReferenceFai::FromFile(GetTestData("unindexed.fasta"),
                                   GetTestData("unindexed.fasta.fai"));
  EXPECT_THAT(result, IsNotOKWithCodeAndMessage(
      tensorflow::error::NOT_FOUND,
      "could not load fasta and/or fai for fasta"));
}

TEST(ReferenceFaiTest, WriteAfterCloseIsntOK) {
  auto reader = JustLoadFai(TestFastaPath());
  ASSERT_THAT(reader->Close(), IsOK());
  EXPECT_THAT(reader->GetBases(MakeRange("chrM", 0, 100)),
              IsNotOKWithCodeAndMessage(
                  tensorflow::error::FAILED_PRECONDITION,
                  "can't read from closed GenomeReferenceFai object"));
}

}  // namespace nucleus
