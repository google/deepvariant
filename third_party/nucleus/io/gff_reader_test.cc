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

// Tests for GffReader class.

#include "third_party/nucleus/io/gff_reader.h"

#include <utility>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/gff.pb.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace nucleus {
namespace {

using nucleus::genomics::v1::GffHeader;
using nucleus::genomics::v1::GffRecord;

const char kExpectedHeaderRecord[] =
    R"(gff_version: "gff-version 3.2.1"
       sequence_regions {
         reference_name: "ctg123"
         start: 0
         end: 1497228
       }
    )";

const char kExpectedGffRecord1[] =
    R"(range {
         reference_name: "ctg123"
         start: 999
         end: 9000
       }
       source: "GenBank"
       type: "gene"
       score: 2.5
       strand: FORWARD_STRAND
       phase: 0
       attributes {
         key: "ID"
         value: "gene00001"
       }
       attributes {
         key: "Name"
         value: "EDEN"
       }
     )";

const char kExpectedGffRecord2[] =
    R"(range {
         reference_name: "ctg123"
         start: 999
         end: 1012
       }
       score: -inf
       phase: -1
     )";

TEST(GffReaderTest, ReadsExampleFile) {
  string examples_fname = nucleus::GetTestData("test_features.gff");

  auto reader_or = GffReader::FromFile(examples_fname);
  EXPECT_TRUE(reader_or.ok());

  auto reader = std::move(reader_or.ValueOrDie());
  const GffHeader& header = reader->Header();
  EXPECT_THAT(header, EqualsProto(kExpectedHeaderRecord));

  // Load the records.
  std::vector<GffRecord> gff_records = as_vector(reader->Iterate());
  EXPECT_EQ(2, gff_records.size());
  // Inspect the records.
  EXPECT_THAT(gff_records[0], EqualsProto(kExpectedGffRecord1));
  EXPECT_THAT(gff_records[1], EqualsProto(kExpectedGffRecord2));
}

// TODO: Add more comprehensive tests

}  // namespace

}  // namespace nucleus
