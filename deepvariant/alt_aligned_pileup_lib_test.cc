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
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {

namespace {

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

nucleus::genomics::v1::Read MakeRead(
    int64_t ref_start, const std::string& bases,
    const std::vector<std::string>& cigar,
    const std::vector<int>& aligned_quality,
    const std::string& read_name = "test_read") {
  nucleus::genomics::v1::Read read =
      nucleus::MakeRead("chr1", ref_start, bases, cigar, read_name);
  read.mutable_aligned_quality()->Clear();
  for (const auto& quality : aligned_quality) {
    read.mutable_aligned_quality()->Add(quality);
  }
  return read;
}

struct TrimCigarTestData {
  int64_t ref_start = 0;
  int64_t ref_length = 0;
  // Cannot use absl::string_view here due to dependencies to a legacy code.
  std::vector<std::string> input_cigar;
  std::vector<std::string> expected_cigar;
  int64_t trimmed_read_start = 0;
  int64_t trimmed_read_length = 0;
};

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

struct TrimReadTestData {
  int64_t read_ref_start = 0;
  int64_t trim_ref_start = 0;
  int64_t trim_ref_length = 0;
  // Cannot use absl::string_view here due to dependencies to a legacy code.
  std::string bases;
  std::vector<std::string> cigar;
  std::vector<int> aligned_quality;
  int64_t expected_read_ref_start = 0;
  std::string expected_bases;
  std::vector<std::string> expected_cigar;
  std::vector<int> expected_aligned_quality;
};

class TrimReadTest : public testing::TestWithParam<TrimReadTestData> {};

TEST_P(TrimReadTest, TrimReadTestCases) {
  const TrimReadTestData& param = GetParam();
  nucleus::genomics::v1::Read read = MakeRead(
      param.read_ref_start, param.bases, param.cigar, param.aligned_quality);
  nucleus::genomics::v1::Read expected_read =
      MakeRead(param.expected_read_ref_start, param.expected_bases,
               param.expected_cigar, param.expected_aligned_quality);
  nucleus::genomics::v1::Range region;
  region.set_reference_name("chr1");
  region.set_start(param.trim_ref_start);
  region.set_end(param.trim_ref_start + param.trim_ref_length);
  nucleus::genomics::v1::Read trimmed_read = TrimRead(read, region);
  EXPECT_THAT(trimmed_read, nucleus::EqualsProto(expected_read));
}

INSTANTIATE_TEST_SUITE_P(TrimReadTests, TrimReadTest,
                         testing::ValuesIn(std::vector<TrimReadTestData>({
                             {10,
                              15,
                              5,
                              "ACGTACGTAAAAAAGTGTGATC",
                              {"22M"},
                              {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                               12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22},
                              15,
                              "CGTAA",
                              {"5M"},
                              {6, 7, 8, 9, 10}},
                             {10,
                              15,
                              5,
                              "ACGTACGTAAAAAAGTGTGATC",
                              {"2M", "3I", "17M"},
                              {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                               12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22},
                              15,
                              "AAAAA",
                              {"5M"},
                              {9, 10, 11, 12, 13}},
                             {10,
                              15,
                              5,
                              "ACGTACGTAAAAAAGTGTGATC",
                              {"2M", "3D", "20M"},
                              {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                               12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22},
                              15,
                              "GTACG",
                              {"5M"},
                              {3, 4, 5, 6, 7}},
                             {10,
                              8,
                              5,
                              "ACGTACGTAAAAAAGTGTGATC",
                              {"22M"},
                              {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                               12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22},
                              10,
                              "ACG",
                              {"3M"},
                              {1, 2, 3}},
                             {10,
                              10,
                              22,
                              "ACGTACGTAAAAAAGTGTGATC",
                              {"22M"},
                              {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                               12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22},
                              10,
                              "ACGTACGTAAAAAAGTGTGATC",
                              {"22M"},
                              {1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11,
                               12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22}},
                         })));

}  // namespace
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
