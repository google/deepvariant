/*
 * Copyright 2017 Google Inc.
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

#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "tensorflow/core/lib/core/stringpiece.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using tensorflow::StringPiece;

namespace {

CallVariantsOutput CreateSingleSiteCalls(StringPiece reference_name, int start,
                                         int end) {
  CallVariantsOutput single_site_call;
  // Add one call to fulfill the assumption of variant having one call.
  single_site_call.mutable_variant()->add_calls();
  single_site_call.mutable_variant()->set_reference_name(
      reference_name.ToString());
  single_site_call.mutable_variant()->set_start(start);
  single_site_call.mutable_variant()->set_end(end);
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
  single_site_calls.push_back(CreateSingleSiteCalls("chr10", 2000, 2002));
  const string& input_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.in.tfrecord");
  const string& output_tfrecord_path = nucleus::MakeTempFile(
      "ProessSingleSiteCallTfRecordsBasicCase.out.tfrecord");
  nucleus::WriteProtosToTFRecord(single_site_calls, input_tfrecord_path);

  ProcessSingleSiteCallTfRecords(contigs, {input_tfrecord_path},
                                 output_tfrecord_path);
  std::vector<CallVariantsOutput> output =
      nucleus::ReadProtosFromTFRecord<CallVariantsOutput>(output_tfrecord_path);

  EXPECT_EQ(output.size(), 4);
  EXPECT_EQ(output[0].variant().reference_name(), "chr1");
  EXPECT_EQ(output[1].variant().reference_name(), "chr10");
  EXPECT_EQ(output[2].variant().reference_name(), "chr10");
  EXPECT_EQ(output[3].variant().reference_name(), "chr10");
  EXPECT_EQ(output[0].variant().start(), 1);
  EXPECT_EQ(output[1].variant().start(), 1000);
  EXPECT_EQ(output[2].variant().start(), 2000);
  EXPECT_EQ(output[3].variant().start(), 2000);
  EXPECT_EQ(output[0].variant().end(), 2);
  EXPECT_EQ(output[1].variant().end(), 1001);
  EXPECT_EQ(output[2].variant().end(), 2001);
  EXPECT_EQ(output[3].variant().end(), 2002);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
