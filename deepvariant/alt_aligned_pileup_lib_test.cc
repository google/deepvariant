/*
 * Copyright 2024 Google LLC.
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

#include "deepvariant/alt_aligned_pileup_lib.h"

#include <cstdint>
#include <string>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/cigar.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

namespace {

struct TrimCigarTestData {
  int64_t ref_start;
  int64_t ref_length;
  // Cannot use absl::string_view here due to dependencies to a legacy code.
  std::vector<std::string> input_cigar;
  std::vector<std::string> expected_cigar;
  int64_t trimmed_read_start;
  int64_t trimmed_read_length;
};

google::protobuf::RepeatedPtrField<nucleus::genomics::v1::CigarUnit> CreateCigar(
    const std::vector<std::string>& cigar_elements) {
  auto cigar_vector = nucleus::MakeCigar(cigar_elements);
  google::protobuf::RepeatedPtrField<nucleus::genomics::v1::CigarUnit> cigar;
  // Add with iterator is not available in external version, so wil have to add
  // elements one by one.
  for (const auto& element : cigar_vector) {
    auto new_el = cigar.Add();
    new_el->set_operation(element.operation());
    new_el->set_operation_length(element.operation_length());
  }
  return cigar;
}

class TrimCigarTest : public testing::TestWithParam<TrimCigarTestData> {};

TEST_P(TrimCigarTest, TrimCigarTestCases) {
  const TrimCigarTestData& param = GetParam();
  // Create input cigar.
  google::protobuf::RepeatedPtrField<nucleus::genomics::v1::CigarUnit> cigar =
      CreateCigar(param.input_cigar);
  // Create expected cigar.
  ::google::protobuf::RepeatedPtrField<nucleus::genomics::v1::CigarUnit> expected_cigar =
      CreateCigar(param.expected_cigar);
  // read_start and read_length.
  int64_t read_start = param.trimmed_read_start;
  int64_t read_length = param.trimmed_read_length;
  // Output cigar.
  ::google::protobuf::RepeatedPtrField<nucleus::genomics::v1::CigarUnit> new_cigar;
  TrimCigar(cigar, param.ref_start, param.ref_length, &new_cigar, read_start,
            read_length);
  EXPECT_EQ(read_start, param.trimmed_read_start);
  EXPECT_EQ(read_length, param.trimmed_read_length);
  // EqualsProto() is not available in external version so we have to use the
  // Nucleus matcher.
  EXPECT_THAT(new_cigar,
              testing::Pointwise(nucleus::EqualsProto(), expected_cigar));
}

INSTANTIATE_TEST_SUITE_P(
    TrimCigarTests, TrimCigarTest,
    testing::ValuesIn(std::vector<TrimCigarTestData>(
        {// Trim cigar with INS.
         {10, 20, {"20M", "5I", "10M"}, {"10M", "5I", "10M"}, 10, 25},
         // Trim cigar with DEL.
         {10, 20, {"20M", "5D", "10M"}, {"10M", "5D", "5M"}, 10, 15},
         // Trim cigar with INS where res_start falls into INS.
         {22, 10, {"20M", "5I", "20M"}, {"10M"}, 27, 10},
         // Trim cigar with DEL where res_start falls into DEL.
         {22, 10, {"20M", "5D", "20M"}, {"3D", "7M"}, 20, 7},
         // Trim cigar where ref_start starts beyond the end of the read.
         {50, 20, {"20M", "5I", "10M"}, {}, 35, 0},
         // Trim cigar where ref_length goes beyond the reads end.
         {10, 40, {"20M", "5I", "10M"}, {"10M", "5I", "10M"}, 10, 25}})));

}  // namespace
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
