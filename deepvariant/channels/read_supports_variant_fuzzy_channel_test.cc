/*
 * Copyright 2025 Google LLC.
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
 * CONTRACT, STRICT LIABILITY, OR TORT WITHINGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#include "deepvariant/channels/read_supports_variant_fuzzy_channel.h"

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"
#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace learning {
namespace genomics {
namespace deepvariant {
namespace {

using ::nucleus::genomics::v1::Read;

// Creates a DeepVariantCall proto for testing.
DeepVariantCall CreateDvCall(
    absl::string_view ref, const std::vector<std::string>& alts,
    const std::map<std::string, std::vector<std::string>>& allele_support_map,
    const std::vector<int>& alt_ps_values = {},
    const std::vector<std::string>& ref_support_reads = {}) {
  DeepVariantCall dv_call;
  dv_call.mutable_variant()->set_reference_bases(ref.data());
  for (const auto& alt : alts) {
    dv_call.mutable_variant()->add_alternate_bases(alt);
  }
  for (const auto& pair : allele_support_map) {
    DeepVariantCall_SupportingReads supporting_reads;
    for (const auto& read_name : pair.second) {
      supporting_reads.add_read_names(read_name);
    }
    (*dv_call.mutable_allele_support())[pair.first] = supporting_reads;
  }
  if (!alt_ps_values.empty()) {
    nucleus::SetInfoField("ALT_PS", alt_ps_values, dv_call.mutable_variant());
  }
  for (const auto& read_name : ref_support_reads) {
    dv_call.add_ref_support(read_name);
  }
  return dv_call;
}

Read CreateRead(absl::string_view fragment_name, int read_number,
                 std::optional<int> hp_value = std::nullopt) {
  Read read;
  read.set_fragment_name(fragment_name.data());
  read.set_read_number(read_number);
  if (hp_value.has_value()) {
    auto& info = (*read.mutable_info())["HP"];
    info.add_values()->set_int_value(hp_value.value());
  }
  return read;
}

class ReadSupportsVariantFuzzyChannelTest : public ::testing::Test {
 protected:
  PileupImageOptions options_;
  ReadSupportsVariantFuzzyChannel channel_{10, options_};
};

TEST_F(ReadSupportsVariantFuzzyChannelTest, ExactMatch) {
  DeepVariantCall dv_call =
      CreateDvCall(/*ref=*/"A", /*alts=*/{"AC", "ACC"},
                   /*allele_support_map=*/{{"AC", {"read1/0"}},
                                          {"ACC", {"read2/0"}}},
                   /*alt_ps_values=*/{0, 1, 1});
  Read read_1 = CreateRead("read1", 0, 1);
  Read read_2 = CreateRead("read2", 0, 1);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read_1, {"AC"}), 1);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read_2, {"AC"}), 10);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read_1, {"ACC"}), 10);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read_2, {"ACC"}), 1);
}

TEST_F(ReadSupportsVariantFuzzyChannelTest, FuzzyMatch1bp) {
  DeepVariantCall dv_call =
      CreateDvCall(/*ref=*/"A", /*alts=*/{"AC", "ACC"},
                   /*allele_support_map=*/{{"ACC", {"read1/0"}}},
                   /*alt_ps_values=*/{0, 1, 1});
  Read read = CreateRead("read1", 0, 1);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read, {"AC"}), 10);
}

TEST_F(ReadSupportsVariantFuzzyChannelTest, FuzzyMatch2bp) {
  DeepVariantCall dv_call =
      CreateDvCall(/*ref=*/"A", /*alts=*/{"AC", "ACCC"},
                   /*allele_support_map=*/{{"ACCC", {"read1/0"}}},
                   /*alt_ps_values=*/{0, 1, 1});
  Read read = CreateRead("read1", 0, 1);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read, {"AC"}), 9);
}

TEST_F(ReadSupportsVariantFuzzyChannelTest, FuzzyMatchPhaseMismatch) {
  DeepVariantCall dv_call =
      CreateDvCall(/*ref=*/"A", /*alts=*/{"AC", "ACC"},
                   /*allele_support_map=*/{{"ACC", {"read1/0"}}},
                   /*alt_ps_values=*/{0, 1, 1});
  Read read = CreateRead("read1", 0, 2);
  // Phase mismatch (read HP=2, variant ALT_PS=1), so no fuzzy support.
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read, {"AC"}), 0);
}

TEST_F(ReadSupportsVariantFuzzyChannelTest, FuzzyMatchPhaseZero) {
  DeepVariantCall dv_call =
      CreateDvCall(/*ref=*/"A", /*alts=*/{"AC", "ACC"},
                   /*allele_support_map=*/{{"ACC", {"read1/0"}}},
                   /*alt_ps_values=*/{0, 0, 0});
  Read read = CreateRead("read1", 0, 1);
  // Variant phase is 0, so fuzzy match should trigger.
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read, {"AC"}), 10);
}

TEST_F(ReadSupportsVariantFuzzyChannelTest, RefSupport) {
  DeepVariantCall dv_call =
      CreateDvCall(/*ref=*/"A", /*alts=*/{"ATGC", "AA"},
                   /*allele_support_map=*/{},
                   /*alt_ps_values=*/{}, /*ref_support_reads=*/{"read1/0"});
  Read read = CreateRead("read1", 0);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read, {"ATGC"}), 0);
  EXPECT_EQ(channel_.ReadSupportsAlt(dv_call, read, {"AA"}), 10);
}

}  // namespace
}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
