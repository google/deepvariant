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

#include "third_party/nucleus/io/bedgraph_reader.h"

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

using ::testing::Pointwise;

namespace {
constexpr char kBedGraphFilename[] = "test_regions.bedgraph";
constexpr char kGzippedBedGraphFilename[] = "test_regions.bedgraph.gz";

// Helper to create a BedGraphRecord for testing.
genomics::v1::BedGraphRecord GetTestRecord(const string& name, int64 start,
                                           int64 end, double data_value) {
  nucleus::genomics::v1::BedGraphRecord r;
  r.set_reference_name(name);
  r.set_start(start);
  r.set_end(end);
  r.set_data_value(data_value);
  return r;
}
}  // namespace

TEST(BedGraphReaderTest, NormalIterationWorks) {
  std::vector<genomics::v1::BedGraphRecord> expected = {
      GetTestRecord("chr1", 10, 20, 100), GetTestRecord("chr1", 100, 200, 250),
      GetTestRecord("chr1", 300, 400, 150.1),
      GetTestRecord("chr1", 500, 501, 20.13)};

  std::unique_ptr<BedGraphReader> reader = std::move(
      BedGraphReader::FromFile(GetTestData(kBedGraphFilename)).ValueOrDie());
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), expected));
}

TEST(BedGraphReaderTest, GzippedIterationWorks) {
  std::unique_ptr<BedGraphReader> reader1 = std::move(
      BedGraphReader::FromFile(GetTestData(kBedGraphFilename)).ValueOrDie());

  std::unique_ptr<BedGraphReader> reader2 =
      std::move(BedGraphReader::FromFile(GetTestData(kGzippedBedGraphFilename))
                    .ValueOrDie());
  EXPECT_THAT(as_vector(reader1->Iterate()),
              Pointwise(EqualsProto(), as_vector(reader2->Iterate())));
}

}  // namespace nucleus
