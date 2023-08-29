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

#include "third_party/nucleus/io/gff_writer.h"

#include <utility>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/gff.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "third_party/nucleus/core/statusor.h"
#include "google/protobuf/text_format.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/env.h"

namespace nucleus {

namespace {

using nucleus::genomics::v1::GffHeader;
using nucleus::genomics::v1::GffRecord;

const char kHeaderRecord[] =
    R"(gff_version: "gff-version 3.2.1"
       sequence_regions {
         reference_name: "ctg123"
         start: 0
         end: 1497228
       }
    )";

const char kGffRecord1[] =
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

const char kGffRecord2[] =
    R"(range {
         reference_name: "ctg123"
         start: 999
         end: 1012
       }
       score: -inf
       phase: -1
     )";

const char kExpectedGffText[] =
    "##gff-version 3.2.1\n"
    "##sequence-region ctg123 1 1497228\n"
    "ctg123\tGenBank\tgene\t1000\t9000\t2.5\t+\t0\tID=gene00001;Name=EDEN\n"
    "ctg123\t.\t.\t1000\t1012\t.\t.\t.\t\n";

TEST(GffWriterTest, WritesGffRecords) {
  GffHeader header;
  google::protobuf::TextFormat::ParseFromString(kHeaderRecord, &header);

  GffRecord record1, record2;
  google::protobuf::TextFormat::ParseFromString(kGffRecord1, &record1);
  google::protobuf::TextFormat::ParseFromString(kGffRecord2, &record2);

  string out_fname = MakeTempFile("gff_writer_test_1.gff");
  StatusOr<std::unique_ptr<GffWriter>> writer_or =
      GffWriter::ToFile(out_fname, header);

  ASSERT_THAT(writer_or.status(), IsOK());
  std::unique_ptr<GffWriter> gff_writer = std::move(writer_or.ValueOrDie());

  ASSERT_THAT(gff_writer->Write(record1), IsOK());
  ASSERT_THAT(gff_writer->Write(record2), IsOK());
  ASSERT_THAT(gff_writer->Close(), IsOK());

  string contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           out_fname, &contents));
  EXPECT_EQ(kExpectedGffText, contents);
}

}  // namespace
}  // namespace nucleus
