/*
 * Copyright 2022 Google LLC.
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

#include "third_party/nucleus/io/merge_variants.h"

#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"

namespace {

using nucleus::EqualsProto;

nucleus::genomics::v1::Variant SingleVariantCallWithLiklihood(
    const std::vector<double>& liklihoods) {
  nucleus::genomics::v1::VariantCall vcall;
  for (const auto& l : liklihoods) {
    vcall.add_genotype_likelihood(l);
  }

  nucleus::genomics::v1::Variant variant;
  *variant.add_calls() = vcall;
  return variant;
}

TEST(ZeroScaleGl, GlScaled) {
  nucleus::genomics::v1::Variant v =
      SingleVariantCallWithLiklihood({1.5, 7, 10});
  nucleus::ZeroScaleGl(&v);

  EXPECT_THAT(v, EqualsProto(SingleVariantCallWithLiklihood({-8.5, -3, 0})));
}

struct GvcfTestCase {
  std::string test_name;
  std::vector<std::string> alts;
  std::vector<double> gls;
  std::vector<double> vaf;
};

class GvcfAltAddedTest : public testing::TestWithParam<GvcfTestCase> {};

nucleus::genomics::v1::Variant BuildTestVariant(const GvcfTestCase test_case) {
  nucleus::genomics::v1::Variant variant;
  for (const auto& ab : test_case.alts) {
    variant.add_alternate_bases(ab);
  }

  nucleus::genomics::v1::VariantCall vcall;
  for (const auto& gl : test_case.gls) {
    vcall.add_genotype_likelihood(gl);
  }

  vcall.mutable_info()->insert({{"VAF", {}}});
  for (const auto& v : test_case.vaf) {
    nucleus::genomics::v1::Value val;
    val.set_number_value(v);
    vcall.mutable_info()->at("VAF").mutable_values()->Add(std::move(val));
  }

  *variant.add_calls() = vcall;

  return variant;
}

TEST_P(GvcfAltAddedTest, TestTransfromToGvcf) {
  nucleus::genomics::v1::Variant test_case = BuildTestVariant(GetParam());

  nucleus::genomics::v1::Variant expected;
  expected.MergeFrom(test_case);
  expected.add_alternate_bases("<*>");
  for (int i = 0; i < GetParam().alts.size() + 2; i++) {
    expected.mutable_calls(0)->add_genotype_likelihood(-99);
  }
  nucleus::genomics::v1::Value val;
  val.set_number_value(0);
  expected.mutable_calls(0)->mutable_info()->at("VAF").mutable_values()->Add(
      std::move(val));

  nucleus::TransfromToGvcf(&test_case);

  EXPECT_THAT(test_case, EqualsProto(expected));
}

INSTANTIATE_TEST_SUITE_P(
    GvcfTests, GvcfAltAddedTest,
    testing::Values(GvcfTestCase{"OneAlt",
                                 {"C"},
                                 {-2.0457574905606752, -0.004364805402450088,
                                  -3.0},
                                 {0.75}},
                    GvcfTestCase{"MultiAlts",
                                 {"G", "C"},
                                 {-1.1368906918484387, -0.5279124552610386,
                                  -0.5923808731731073, -0.8155431286425007,
                                  -0.8415961054266092, -1.108308924501657},
                                 {0.5, 0.1}},
                    GvcfTestCase{"MultiAlts2",
                                 {"G", "C", "T"},
                                 {-0.7956722868920258, -0.663917423732382,
                                  -1.493986734511771, -0.8202531343562444,
                                  -0.9377869397242453, -1.0415699718993066,
                                  -1.4176189291054515, -1.5795151893394743,
                                  -1.8101482990393198, -0.8139951558313916},
                                 {0.5, 0.1, 0.05}}),
    [](const testing::TestParamInfo<GvcfAltAddedTest::ParamType>& info) {
      return info.param.test_name;
    });

class GvcfNoAltAddedTest : public testing::TestWithParam<GvcfTestCase> {};

TEST_P(GvcfNoAltAddedTest, TestTransfromToGvcf) {
  nucleus::genomics::v1::Variant test_case = BuildTestVariant(GetParam());

  nucleus::genomics::v1::Variant expected;
  expected.MergeFrom(test_case);

  nucleus::TransfromToGvcf(&test_case);

  EXPECT_THAT(test_case, EqualsProto(expected));
}

INSTANTIATE_TEST_SUITE_P(
    GvcfTests, GvcfNoAltAddedTest,
    testing::Values(GvcfTestCase{"OneAlt",
                                 {"<*>"},
                                 {-2.0457574905606752, -0.004364805402450088,
                                  -3.0},
                                 {0.75}},
                    GvcfTestCase{"MultiAlts",
                                 {"G", "<*>", "C"},
                                 {-1.1368906918484387, -0.5279124552610386,
                                  -0.5923808731731073, -0.8155431286425007,
                                  -0.8415961054266092, -1.108308924501657},
                                 {0.5, 0.1}},
                    GvcfTestCase{"LastAls",
                                 {"G", "C", "T", "<*>"},
                                 {-0.7956722868920258, -0.663917423732382,
                                  -1.493986734511771, -0.8202531343562444,
                                  -0.9377869397242453, -1.0415699718993066,
                                  -1.4176189291054515, -1.5795151893394743,
                                  -1.8101482990393198, -0.8139951558313916},
                                 {0.5, 0.1, 0.05}}),
    [](const testing::TestParamInfo<GvcfAltAddedTest::ParamType>& info) {
      return info.param.test_name;
    });

// Helper method to create a test sequence.
void CreateTestSeq(
    const std::string& name, const int pos_in_fasta, const int range_start,
    const int range_end, const std::string& bases,
    std::vector<nucleus::genomics::v1::ContigInfo>* contigs,
    std::vector<nucleus::genomics::v1::ReferenceSequence>* seqs) {
  CHECK(pos_in_fasta >= 0 && pos_in_fasta < contigs->size());
  nucleus::genomics::v1::ContigInfo* contig = &contigs->at(pos_in_fasta);
  contig->set_name(name);
  contig->set_pos_in_fasta(pos_in_fasta);
  contig->set_n_bases(range_end - range_start);
  nucleus::genomics::v1::ReferenceSequence* seq = &seqs->at(pos_in_fasta);
  seq->mutable_region()->set_reference_name(name);
  seq->mutable_region()->set_start(range_start);
  seq->mutable_region()->set_end(range_end);
  seq->set_bases(bases);
}

nucleus::genomics::v1::Variant VariantCallWithStartEnd(std::string ref_name,
                                                       int start, int end,
                                                       std::string refs) {
  nucleus::genomics::v1::Variant variant;
  variant.set_reference_name(ref_name);
  variant.set_start(start);
  variant.set_end(end);
  variant.set_reference_bases(refs);
  return variant;
}

TEST(CreateRecordFromTemplateTest, TemplateCreated) {
  int kNum = 1;
  std::vector<nucleus::genomics::v1::ContigInfo> contigs(kNum);
  std::vector<nucleus::genomics::v1::ReferenceSequence> seqs(kNum);

  // Creating a InMemoryFastaReader with a test sequence.
  CreateTestSeq("chr1", 0, 0, 32, "AACCGGTTACGTTCGATTTTAAAACCCCGGGG", &contigs,
                &seqs);
  std::unique_ptr<nucleus::InMemoryFastaReader> ref = std::move(
      nucleus::InMemoryFastaReader::Create(contigs, seqs).ValueOrDie());

  nucleus::genomics::v1::Variant v = VariantCallWithStartEnd("1", 10, 15, "G");

  std::unique_ptr<nucleus::genomics::v1::Variant> t =
      nucleus::CreateRecordFromTemplate(v, 10, 12, *ref);

  EXPECT_THAT(t, testing::Pointee(
                     EqualsProto(VariantCallWithStartEnd("1", 10, 12, "G"))));
}

}  // namespace
