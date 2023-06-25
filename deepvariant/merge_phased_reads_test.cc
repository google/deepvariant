/*
 * Copyright 2023 Google LLC.
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

#include "deepvariant/merge_phased_reads.h"

#include <vector>

#include "tensorflow/core/platform/test.h"
#include "absl/hash/hash_testing.h"
#include "absl/status/status.h"
#include "absl/status/statusor.h"
#include "absl/strings/string_view.h"

namespace learning {
namespace genomics {
namespace deepvariant {

TEST(ShardRegion, SupportsAbslHash) {
  EXPECT_TRUE(absl::VerifyTypeImplementsAbslHashCorrectly({
      ShardRegion(),
      {.shard = 1, .region = 2},
      {.shard = 2, .region = 3},
      {.shard = 0, .region = 0},
      {.shard = -1, .region = -2},
      {.shard = 10, .region = 10},
      {.shard = -10, .region = -10},
      {.shard = -10, .region = 10},
      {.shard = 10, .region = -10}
  }));
}

struct ParseFilenameExpectation {
  absl::string_view filename;
  ShardedFileSpec file_spec;
  absl::StatusCode result_code = absl::StatusCode::kOk;
};

class ShardedFileSpecTest :
    public testing::TestWithParam<ParseFilenameExpectation> {
};

TEST_P(ShardedFileSpecTest, parse_sharded_file_spec_utest) {
  const ParseFilenameExpectation& exp = GetParam();
  absl::StatusOr<ShardedFileSpec> file_spec =
      parse_sharded_file_spec(exp.filename);
  EXPECT_EQ(file_spec.status().code(), exp.result_code);
  if (exp.result_code == absl::StatusCode::kOk) {
    EXPECT_EQ(file_spec.value().basename, exp.file_spec.basename);
    EXPECT_EQ(file_spec.value().nshards, exp.file_spec.nshards);
    EXPECT_EQ(file_spec.value().suffix, exp.file_spec.suffix);
  }
}

INSTANTIATE_TEST_SUITE_P(
    ShardedFileTests,
    ShardedFileSpecTest,
    ::testing::ValuesIn(std::vector<ParseFilenameExpectation>({
    {
      .filename = "/dir/foo/bar@3",
      .file_spec = {
        .basename = "/dir/foo/bar",
        .nshards = 3,
        .suffix = "",
      },
      .result_code = absl::StatusCode::kOk,
    },
    {
      .filename = "/dir/foo/bar@3.txt",
      .file_spec = {
        .basename = "/dir/foo/bar",
        .nshards = 3,
        .suffix = "txt",
      },
      .result_code = absl::StatusCode::kOk,
    },
    {
      .filename = "/dir/foo/bar_not_sharded.txt",
      .file_spec = {
        .basename = "",
        .nshards = 0,
        .suffix = "",
      },
      .result_code = absl::StatusCode::kInvalidArgument,
    },
    {
      .filename = "/dir/foo/incorrect_sharded_spec@99999999999999999.txt",
      .file_spec = {
        .basename = "",
        .nshards = 0,
        .suffix = "",
      },
      .result_code = absl::StatusCode::kOutOfRange,
    },
 })));

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
