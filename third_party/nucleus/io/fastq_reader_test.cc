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

#include "third_party/nucleus/io/fastq_reader.h"

#include <utility>
#include <vector>

#include <gmock/gmock.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"

namespace nucleus {

using std::vector;

using ::testing::Pointwise;

constexpr char kFastqFilename[] = "test_reads.fastq";
constexpr char kGzippedFastqFilename[] = "test_reads.fastq.gz";
constexpr char kBgzippedFastqFilename[] = "test_reads.bgzip.fastq.gz";

void CreateRecord(const string& id, const string& description,
                  const string& sequence, const string& quality,
                  nucleus::genomics::v1::FastqRecord* record) {
  record->set_id(id);
  record->set_description(description);
  record->set_sequence(sequence);
  record->set_quality(quality);
}

class FastqReaderTest : public ::testing::Test {
 protected:
  void SetUp() override {
    nucleus::genomics::v1::FastqRecord first, second, third, fourth;
    CreateRecord("NODESC:header", "", "GATTACA", "BB>B@FA", &first);
    CreateRecord(
        "M01321:49:000000000-A6HWP:1:1101:17009:2216", "1:N:0:1",
        "CGTTAGCGCAGGGGGCATCTTCACACTGGTGACAGGTAACCGCCGTAGTAAAGGTTCCGCCTTTCACT",
        "AAAAABF@BBBDGGGG?FFGFGHBFBFBFABBBHGGGFHHCEFGGGGG?FGFFHEDG3EFGGGHEGHG",
        &second);
    CreateRecord("FASTQ", "contains multiple spaces in description",
                 "CGGCTGGTCAGGCTGACATCGCCGCCGGCCTGCAGCGAGCCGCTGC",
                 "FAFAF;F/9;.:/;999B/9A.DFFF;-->.AAB/FC;9-@-=;=.", &third);
    CreateRecord("FASTQ_with_trailing_space", "", "CGG", "FAD", &fourth);
    golden_ = {first, second, third, fourth};
  }

  vector<nucleus::genomics::v1::FastqRecord> golden_;
};

TEST_F(FastqReaderTest, NormalIterationWorks) {
  std::unique_ptr<FastqReader> reader = std::move(
      FastqReader::FromFile(GetTestData(kFastqFilename),
                            nucleus::genomics::v1::FastqReaderOptions())
          .ValueOrDie());

  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden_));
}

TEST_F(FastqReaderTest, GzippedIterationWorks) {
  auto opts = nucleus::genomics::v1::FastqReaderOptions();
  std::unique_ptr<FastqReader> reader =
      std::move(FastqReader::FromFile(GetTestData(kGzippedFastqFilename), opts)
                    .ValueOrDie());

  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden_));
}

TEST_F(FastqReaderTest, BgzippedIterationWorks) {
  auto opts = nucleus::genomics::v1::FastqReaderOptions();
  std::unique_ptr<FastqReader> reader =
      std::move(FastqReader::FromFile(GetTestData(kBgzippedFastqFilename), opts)
                    .ValueOrDie());

  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden_));
}
}  // namespace nucleus
