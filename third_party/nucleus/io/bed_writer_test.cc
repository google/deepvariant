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

#include "third_party/nucleus/io/bed_writer.h"

#include <utility>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/env.h"

namespace nucleus {

using std::vector;

class BedWriterTest : public ::testing::Test {
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
    header_.set_num_fields(12);
  }

  nucleus::genomics::v1::BedHeader header_;
  vector<nucleus::genomics::v1::BedRecord> golden_;
};

TEST_F(BedWriterTest, WritingWorks) {
  const string output_filename = MakeTempFile("writes_bed.bed");
  std::unique_ptr<BedWriter> writer =
      std::move(BedWriter::ToFile(output_filename, header_,
                                  nucleus::genomics::v1::BedWriterOptions())
                    .ValueOrDie());
  for (const nucleus::genomics::v1::BedRecord& record : golden_) {
    ASSERT_THAT(writer->Write(record), IsOK());
  }
  ASSERT_THAT(writer->Close(), IsOK());
  writer.reset();

  string contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &contents));
  const string kExpectedBedContent =
      "chr1\t10\t20\tfirst\t100\t+\t12\t18\t255,124,1\t3\t2,6,2\t10,12,18\n"
      "chr1\t100\t200\tsecond\t250\t.\t120\t180\t252,122,12\t2\t35,40\t100,"
      "160\n";
  EXPECT_EQ(kExpectedBedContent, contents);
}

TEST_F(BedWriterTest, WritesGzippedFiles) {
  const string output_filename = MakeTempFile("writes_bed.bed.gz");
  std::unique_ptr<BedWriter> writer =
      std::move(BedWriter::ToFile(output_filename, header_,
                                  nucleus::genomics::v1::BedWriterOptions())
                    .ValueOrDie());
  ASSERT_THAT(writer->Close(), IsOK());

  string contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &contents));
  EXPECT_THAT(IsGzipped(contents),
              "BED writer should be able to writed gzipped output");
}

TEST_F(BedWriterTest, WritingTruncatedWorks) {
  const string output_filename = MakeTempFile("writes_short_bed.bed");
  nucleus::genomics::v1::BedHeader truncated_header;
  truncated_header.set_num_fields(6);
  std::unique_ptr<BedWriter> writer =
      std::move(BedWriter::ToFile(output_filename, truncated_header,
                                  nucleus::genomics::v1::BedWriterOptions())
                    .ValueOrDie());
  for (const nucleus::genomics::v1::BedRecord& record : golden_) {
    ASSERT_THAT(writer->Write(record), IsOK());
  }
  ASSERT_THAT(writer->Close(), IsOK());
  writer.reset();

  string contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &contents));
  const string kExpectedBedContent =
      "chr1\t10\t20\tfirst\t100\t+\n"
      "chr1\t100\t200\tsecond\t250\t.\n";
  EXPECT_EQ(kExpectedBedContent, contents);
}
}  // namespace nucleus
