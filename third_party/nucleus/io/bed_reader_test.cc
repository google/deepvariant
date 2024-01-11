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

#include "third_party/nucleus/io/bed_reader.h"

#include <utility>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"

namespace nucleus {

using std::vector;

using ::testing::Pointwise;

constexpr char kBedFilename[] = "test_regions.bed";
constexpr char kGzippedBedFilename[] = "test_regions.bed.gz";
constexpr char kMalformedBedFilename[] = "malformed.bed";
constexpr char k5ColumnBedFileName[] = "5col.bed";

class BedReaderTest : public ::testing::Test {
 protected:
  void SetUp() override {
    nucleus::genomics::v1::BedRecord first, second;
    first.set_reference_name("chr1");
    first.set_start(10);
    first.set_end(20);
    first.set_name("first");
    first.set_score(100);
    first.set_strand(nucleus::genomics::v1::BedRecord::FORWARD_STRAND);
    first.set_thick_start(12);
    first.set_thick_end(18);
    first.set_item_rgb("255,124,1");
    first.set_block_count(3);
    first.set_block_sizes("2,6,2");
    first.set_block_starts("10,12,18");

    second.set_reference_name("chr1");
    second.set_start(100);
    second.set_end(200);
    second.set_name("second");
    second.set_score(250);
    second.set_strand(nucleus::genomics::v1::BedRecord::NO_STRAND);
    second.set_thick_start(120);
    second.set_thick_end(180);
    second.set_item_rgb("252,122,12");
    second.set_block_count(2);
    second.set_block_sizes("35,40");
    second.set_block_starts("100,160");

    golden_ = {first, second};
  }

  vector<nucleus::genomics::v1::BedRecord> golden_;
};

TEST_F(BedReaderTest, NormalIterationWorks) {
  std::unique_ptr<BedReader> reader =
      std::move(BedReader::FromFile(GetTestData(kBedFilename),
                                    nucleus::genomics::v1::BedReaderOptions())
                    .ValueOrDie());

  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden_));
}

TEST_F(BedReaderTest, GzippedIterationWorks) {
  auto opts = nucleus::genomics::v1::BedReaderOptions();
  std::unique_ptr<BedReader> reader = std::move(
      BedReader::FromFile(GetTestData(kGzippedBedFilename), opts).ValueOrDie());

  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden_));
}

TEST_F(BedReaderTest, FieldRestrictionWorks) {
  auto opts = nucleus::genomics::v1::BedReaderOptions();
  opts.set_num_fields(4);
  std::unique_ptr<BedReader> reader = std::move(
      BedReader::FromFile(GetTestData(kBedFilename), opts).ValueOrDie());

  nucleus::genomics::v1::BedRecord efirst, esecond;
  efirst.set_reference_name("chr1");
  efirst.set_start(10);
  efirst.set_end(20);
  efirst.set_name("first");

  esecond.set_reference_name("chr1");
  esecond.set_start(100);
  esecond.set_end(200);
  esecond.set_name("second");

  vector<nucleus::genomics::v1::BedRecord> expected = {efirst, esecond};

  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), expected));
}

TEST_F(BedReaderTest, MalformedBedRecord) {
  std::unique_ptr<BedReader> reader =
      std::move(BedReader::FromFile(GetTestData(kMalformedBedFilename),
                                    nucleus::genomics::v1::BedReaderOptions())
                    .ValueOrDie());

  EXPECT_DEATH(as_vector(reader->Iterate()),
               "BED record has invalid number of fields");
}

TEST_F(BedReaderTest, FiveColBedRecord) {
  auto status = BedReader::FromFile(GetTestData(k5ColumnBedFileName),
                                    nucleus::genomics::v1::BedReaderOptions());
  // See https://github.com/google/deepvariant/issues/374#issuecomment-723752207
  EXPECT_TRUE(status.ok());
}

}  // namespace nucleus
