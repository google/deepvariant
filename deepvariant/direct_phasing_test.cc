/*
 * Copyright 2021 Google LLC.
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

#include "deepvariant/direct_phasing.h"

#include <sys/types.h>

#include <string_view>

#include "deepvariant/protos/deepvariant.pb.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace {

using nucleus::genomics::v1::Variant;

class DirectPhasingTest : public ::testing::Test {
 public:
  using AlleleSupport =
      absl::flat_hash_map<std::string, std::vector<std::string>>;

  DeepVariantCall MakeCandidate(
      int64_t start, int64_t end,
      AlleleSupport allele_support = AlleleSupport()) {
    DeepVariantCall candidate;
    Variant* variant = candidate.mutable_variant();
    variant->set_start(start);
    variant->set_end(end);
    if (!allele_support.empty()) {
      auto allele_support_field = candidate.mutable_allele_support_ext();
      for (const auto& one_allele_support : allele_support) {
        auto& support_infos = (*allele_support_field)[one_allele_support.first];
        for (const auto& read_name : one_allele_support.second) {
          auto* read_info = support_infos.add_read_infos();
          read_info->set_read_name(read_name);
          read_info->set_is_low_quality(false);
        }
      }
    }
    return candidate;
  }
};

TEST_F(DirectPhasingTest, TestAlleleTypeFromCandidateSubstitution) {
  EXPECT_EQ(AlleleType::SUBSTITUTION,
            AlleleTypeFromCandidate("CC", MakeCandidate(100, 102)));
}

TEST_F(DirectPhasingTest, TestAlleleTypeFromCandidateDeletion) {
  EXPECT_EQ(AlleleType::DELETION,
            AlleleTypeFromCandidate("C", MakeCandidate(100, 102)));
}

TEST_F(DirectPhasingTest, TestAlleleTypeFromCandidateInsertion) {
  EXPECT_EQ(AlleleType::INSERTION,
            AlleleTypeFromCandidate("CCC", MakeCandidate(100, 101)));
}

TEST_F(DirectPhasingTest, TestAlleleTypeFromCandidateOneBaseSubstitution) {
  EXPECT_EQ(AlleleType::SUBSTITUTION,
            AlleleTypeFromCandidate("A", MakeCandidate(100, 101)));
}

struct TestNumOfSubstitutionAllelesTestCase {
  DeepVariantCall candidate;
  int expected_allele_depth;
};

TEST_F(DirectPhasingTest, TestNumOfSubstitutionAllelesMultipleSubAlleles) {
  EXPECT_EQ(2,
    NumOfSubstitutionAlleles(MakeCandidate(100, 101, {
      {"A", {"read1", "read2", "read3"}},  // SUB allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST_F(DirectPhasingTest, TestNumOfSubstitutionAllelesUncalledAllelePresent) {
  EXPECT_EQ(1,
    NumOfSubstitutionAlleles(MakeCandidate(100, 101, {
      {"UNCALLED_ALLELE", {"read1", "read2", "read3"}},  // Uncalled allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST_F(DirectPhasingTest, TestNumOfIndelAlleles2Sub1Indel) {
  EXPECT_EQ(1,
    NumOfIndelAlleles(MakeCandidate(100, 101, {
      {"A", {"read1", "read2", "read3"}},  // SUB allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST_F(DirectPhasingTest, TestNumOfIndelAllelesUncalledPresent) {
  EXPECT_EQ(2,
    NumOfIndelAlleles(MakeCandidate(100, 103, {
      {"UNCALLED_ALLELE", {"read1", "read2", "read3"}},  // Uncalled allele
      {"C", {"read4", "read5"}},          // DEL allele
      {"CCCC", {"read6", "read7"}}}         // INS allele
      )));
}

TEST_F(DirectPhasingTest, TestSubstitutionAllelesDepth2SubAlleles) {
  EXPECT_EQ(5,
    SubstitutionAllelesDepth(MakeCandidate(100, 101, {
      {"A", {"read1", "read2", "read3"}},  // SUB allele
      {"C", {"read4", "read5"}},          // SUB allele
      {"CC", {"read6", "read7"}}}         // INDEL allele
      )));
}

TEST_F(DirectPhasingTest, TestSubstitutionAllelesDepth1UncalledAnd2Indels) {
  EXPECT_EQ(0,
    SubstitutionAllelesDepth(MakeCandidate(100, 103, {
      {"UNCALLED_ALLELE", {"read1", "read2", "read3"}},  // Uncalled allele
      {"C", {"read4", "read5"}},          // DEL allele
      {"CCCC", {"read6", "read7"}}}         // INS allele
      )));
}

}  // namespace
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
