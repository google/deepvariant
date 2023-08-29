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

#include "third_party/nucleus/io/bedgraph_writer.h"

#include <utility>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/bedgraph_reader.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/env.h"

namespace nucleus {

using genomics::v1::BedGraphRecord;

namespace {

constexpr char kBedGraphFilename[] = "test_regions.bedgraph";

// Helper to create a BedGraphRecord for testing.
BedGraphRecord MakeTestRecord(const string& name, int64 start, int64 end,
                              double data_value) {
  nucleus::genomics::v1::BedGraphRecord r;
  r.set_reference_name(name);
  r.set_start(start);
  r.set_end(end);
  r.set_data_value(data_value);
  return r;
}

}  // namespace

TEST(BedGraphWriterTest, Writes) {
  std::vector<BedGraphRecord> expected = {
      MakeTestRecord("chr1", 10, 20, 100.1),
      MakeTestRecord("chr1", 100, 200, 250.50),
      MakeTestRecord("chr1", 300, 320, 25.13)};
  const string output_filename = MakeTempFile("writes.bedgraph");
  std::unique_ptr<BedGraphWriter> writer =
      std::move(BedGraphWriter::ToFile(output_filename).ValueOrDie());
  for (const BedGraphRecord& record : expected) {
    ASSERT_THAT(writer->Write(record), IsOK());
  }
  ASSERT_THAT(writer->Close(), IsOK());
  writer.reset();

  string contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &contents));
  const string kExpectedBedContent =
      "chr1\t10\t20\t100.1\n"
      "chr1\t100\t200\t250.5\n"
      "chr1\t300\t320\t25.13\n";
  EXPECT_EQ(kExpectedBedContent, contents);
  TF_CHECK_OK(tensorflow::Env::Default()->DeleteFile(output_filename));
}

TEST(BedGraphWriterTest, WritesGzippedFiles) {
  std::vector<BedGraphRecord> expected = {
      MakeTestRecord("chr1", 10, 20, 100.1),
      MakeTestRecord("chr1", 100, 200, 250.50),
      MakeTestRecord("chr1", 300, 320, 25.13)};
  const string output_filename = MakeTempFile("writes.bedgraph.gz");
  std::unique_ptr<BedGraphWriter> writer =
      std::move(BedGraphWriter::ToFile(output_filename).ValueOrDie());
  for (const BedGraphRecord& record : expected) {
    ASSERT_THAT(writer->Write(record), IsOK());
  }
  ASSERT_THAT(writer->Close(), IsOK());

  string contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &contents));
  EXPECT_THAT(IsGzipped(contents),
              "BED writer should be able to writed gzipped output");
  TF_CHECK_OK(tensorflow::Env::Default()->DeleteFile(output_filename));
}

TEST(BedGraphWriterTest, RoundTrip) {
  std::unique_ptr<BedGraphReader> reader = std::move(
      BedGraphReader::FromFile(GetTestData(kBedGraphFilename)).ValueOrDie());
  const string output_filename = MakeTempFile("writes.bedgraph");
  std::vector<BedGraphRecord> expected = as_vector(reader->Iterate());
  std::unique_ptr<BedGraphWriter> writer =
      std::move(BedGraphWriter::ToFile(output_filename).ValueOrDie());
  for (const BedGraphRecord& record : expected) {
    ASSERT_THAT(writer->Write(record), IsOK());
  }
  ASSERT_THAT(writer->Close(), IsOK());
  writer.reset();

  std::unique_ptr<BedGraphReader> reader2 =
      std::move(BedGraphReader::FromFile(output_filename).ValueOrDie());

  std::vector<BedGraphRecord> actual = as_vector(reader2->Iterate());

  EXPECT_THAT(actual, ::testing::Pointwise(EqualsProto(), expected));
  TF_CHECK_OK(tensorflow::Env::Default()->DeleteFile(output_filename));
}

}  // namespace nucleus
