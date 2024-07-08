/*
 * Copyright 2017 Google LLC.
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

#include "deepvariant/postprocess_variants.h"

#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

namespace {

CallVariantsOutput CreateSingleSiteCalls(absl::string_view reference_name,
                                         int start,
                                         int end) {
  CallVariantsOutput single_site_call;
  // Add one call to fulfill the assumption of variant having one call.
  single_site_call.mutable_variant()->add_calls();
  single_site_call.mutable_variant()->set_reference_name(
      string(reference_name));
  single_site_call.mutable_variant()->set_start(start);
  single_site_call.mutable_variant()->set_end(end);
  return single_site_call;
}

CallVariantsOutput CreateSingleSiteCalls(absl::string_view reference_name,
                                         int start,
                                         int end, double quality) {
  CallVariantsOutput single_site_call =
      CreateSingleSiteCalls(reference_name, start, end);
  single_site_call.mutable_variant()->set_quality(quality);
  return single_site_call;
}

}  // namespace

TEST(ProcessSingleSiteCallTfRecords, BasicCase) {
  std::vector<nucleus::genomics::v1::ContigInfo> contigs =
      nucleus::CreateContigInfos({"chr1", "chr10"}, {0, 1000});
  std::vector<CallVariantsOutput> single_site_calls;
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 1000, 1001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr1", 1, 2));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.9));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.7));
  const string& input_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.in.tfrecord");
  const string& output_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.out.tfrecord");
  nucleus::WriteProtosToTFRecord(single_site_calls, input_tfrecord_path);

  ProcessSingleSiteCallTfRecords(contigs, {input_tfrecord_path},
                                 output_tfrecord_path,
                                 std::vector<nucleus::genomics::v1::Range>());
  std::vector<CallVariantsOutput> output =
      nucleus::ReadProtosFromTFRecord<CallVariantsOutput>(output_tfrecord_path);

  EXPECT_EQ(output.size(), 5);
  EXPECT_EQ(output[0].variant().reference_name(), "chr1");
  EXPECT_EQ(output[1].variant().reference_name(), "chr10");
  EXPECT_EQ(output[2].variant().reference_name(), "chr10");
  EXPECT_EQ(output[3].variant().reference_name(), "chr10");
  EXPECT_EQ(output[4].variant().reference_name(), "chr10");
  EXPECT_EQ(output[0].variant().start(), 1);
  EXPECT_EQ(output[1].variant().start(), 1000);
  EXPECT_EQ(output[2].variant().start(), 2000);
  EXPECT_EQ(output[3].variant().start(), 2000);
  EXPECT_EQ(output[4].variant().start(), 2000);
  EXPECT_EQ(output[0].variant().end(), 2);
  EXPECT_EQ(output[1].variant().end(), 1001);
  EXPECT_EQ(output[2].variant().end(), 2001);
  EXPECT_EQ(output[3].variant().end(), 2002);
  EXPECT_EQ(output[4].variant().end(), 2002);
  // Order of calls with the same reference, start, end should be preserved.
  EXPECT_EQ(output[3].variant().quality(), 0.9);
  EXPECT_EQ(output[4].variant().quality(), 0.7);
}

TEST(ProcessSingleSiteCallTfRecords, SplitByRange) {
  std::vector<nucleus::genomics::v1::ContigInfo> contigs =
      nucleus::CreateContigInfos({"chr1", "chr10"}, {0, 1000});
  std::vector<CallVariantsOutput> single_site_calls;
  std::vector<nucleus::genomics::v1::Range> ranges = {
      nucleus::MakeRange("chr1", 0, 100)};
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 1000, 1001));
  single_site_calls.push_back(CreateSingleSiteCalls("chr1", 1, 2));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.9));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002, 0.7));
  const string& input_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.in.tfrecord");
  const string& output_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.out.tfrecord");
  nucleus::WriteProtosToTFRecord(single_site_calls, input_tfrecord_path);

  ProcessSingleSiteCallTfRecords(contigs, {input_tfrecord_path},
                                 output_tfrecord_path, ranges);
  std::vector<CallVariantsOutput> output =
      nucleus::ReadProtosFromTFRecord<CallVariantsOutput>(output_tfrecord_path);

  EXPECT_EQ(output.size(), 1);
  EXPECT_EQ(output[0].variant().reference_name(), "chr1");
  EXPECT_EQ(output[0].variant().start(), 1);
  EXPECT_EQ(output[0].variant().end(), 2);
}

TEST(ProcessSingleSiteCallTfRecords, SplitByRangesHandleNeighbors) {
  std::vector<nucleus::genomics::v1::ContigInfo> contigs =
      nucleus::CreateContigInfos({"chr1", "chr10"}, {0, 1000});
  std::vector<CallVariantsOutput> single_site_calls;
  std::vector<nucleus::genomics::v1::Range> ranges = {
      nucleus::MakeRange("chr10", 0, 100),
      nucleus::MakeRange("chr10", 100, 200),
  };
  single_site_calls.push_back(CreateSingleSiteCalls("chr1", 1, 2));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 50, 51));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 100, 110));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 199, 200));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 199, 210));
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 200, 210));
  const string& input_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.in.tfrecord");
  const string& output_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.out.tfrecord");
  nucleus::WriteProtosToTFRecord(single_site_calls, input_tfrecord_path);

  ProcessSingleSiteCallTfRecords(contigs, {input_tfrecord_path},
                                 output_tfrecord_path, ranges);
  std::vector<CallVariantsOutput> output =
      nucleus::ReadProtosFromTFRecord<CallVariantsOutput>(output_tfrecord_path);

  EXPECT_EQ(output.size(), 4);
  EXPECT_EQ(output[0].variant().reference_name(), "chr10");
  EXPECT_EQ(output[0].variant().start(), 50);
  EXPECT_EQ(output[0].variant().end(), 51);
  EXPECT_EQ(output[1].variant().reference_name(), "chr10");
  EXPECT_EQ(output[1].variant().start(), 100);
  EXPECT_EQ(output[1].variant().end(), 110);
  EXPECT_EQ(output[2].variant().reference_name(), "chr10");
  EXPECT_EQ(output[2].variant().start(), 199);
  EXPECT_EQ(output[2].variant().end(), 200);
  EXPECT_EQ(output[3].variant().reference_name(), "chr10");
  EXPECT_EQ(output[3].variant().start(), 199);
  EXPECT_EQ(output[3].variant().end(), 210);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
