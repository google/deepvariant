/*
 * Copyright 2018 Google Inc.
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

#include <map>

#include "deepvariant/realigner/fast_pass_aligner.h"
#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace learning {
namespace genomics {
namespace deepvariant {

class Aligner2Test : public ::testing::Test {
 protected:
  FastPassAligner aligner_;

  void SetUp() override {
    aligner_.set_reference(
        "ATCAAGGGAAAAAGTGCCCAGGGCCAAATATGTTTTGGGTTTTGCAGGACAAAGTATGGTTGAAACTGAG"
            "CTGAAGATATG");
  }
};

TEST_F(Aligner2Test, ReadsIndexIntegrationTest) {
  aligner_.set_reads({"AAACCC", "CTCTCT", "TGAGCTGAAG"});

  KmerIndexType expected_index = {
    {"AAA", {KmerOccurrence(ReadId(0), KmerOffset(0))}},
    {"AAC", {KmerOccurrence(ReadId(0), KmerOffset(1))}},
    {"ACC", {KmerOccurrence(ReadId(0), KmerOffset(2))}},
    {"CCC", {KmerOccurrence(ReadId(0), KmerOffset(3))}},
    {"CTC",
     {KmerOccurrence(ReadId(1), KmerOffset(0)),
      KmerOccurrence(ReadId(1), KmerOffset(2))}},
    {"TCT",
     {KmerOccurrence(ReadId(1), KmerOffset(1)),
      KmerOccurrence(ReadId(1), KmerOffset(3))}},
    {"TGA",
     {KmerOccurrence(ReadId(2), KmerOffset(0)),
      KmerOccurrence(ReadId(2), KmerOffset(5))}},
    {"GAG", {KmerOccurrence(ReadId(2), KmerOffset(1))}},
    {"AGC", {KmerOccurrence(ReadId(2), KmerOffset(2))}},
    {"GCT", {KmerOccurrence(ReadId(2), KmerOffset(3))}},
    {"CTG", {KmerOccurrence(ReadId(2), KmerOffset(4))}},
    {"GAA", {KmerOccurrence(ReadId(2), KmerOffset(6))}},
    {"AAG", {KmerOccurrence(ReadId(2), KmerOffset(7))}}};

  aligner_.set_kmer_size(3);
  aligner_.BuildIndex();

  EXPECT_EQ(aligner_.GetKmerIndex(), expected_index);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
