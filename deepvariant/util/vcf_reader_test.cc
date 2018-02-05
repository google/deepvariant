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

#include "deepvariant/util/vcf_reader.h"

#include "deepvariant/util/testing/protocol-buffer-matchers.h"
#include "deepvariant/util/test_utils.h"
#include "deepvariant/util/utils.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace nucleus {

using std::vector;

using ::testing::Eq;
using ::testing::Pointwise;
using ::testing::SizeIs;

using nucleus::genomics::v1::Variant;
using nucleus::proto::IgnoringFieldPaths;

// These files are made using the Verily converter.  See:
//  g3doc/learning/genomics/g3doc/variant-encoding.md
constexpr char kVcfSamplesFilename[] = "test_samples.vcf";
constexpr char kVcfIndexSamplesFilename[] = "test_samples.vcf.gz";
constexpr char kVcfSamplesGoldenFilename[] = "test_samples.vcf.golden.tfrecord";
constexpr char kVcfSitesFilename[] = "test_sites.vcf";
constexpr char kVcfSitesGoldenFilename[] = "test_sites.vcf.golden.tfrecord";

/* These are made using our own converter, vcf_to_tfrecord, for example:
TESTDATA=./core/testdata
blaze run -c opt internal:vcf_to_tfrecord -- \
  --input=$(pwd)/$TESTDATA/test_likelihoods_input.vcf \
  --output=$(pwd)/$TESTDATA/test_likelihoods.vcf.golden.tfrecord
*/
constexpr char kVcfLikelihoodsFilename[] = "test_likelihoods_input.vcf";
constexpr char kVcfLikelihoodsGoldenFilename[] = "test_likelihoods.vcf.golden.tfrecord";  // NOLINT
constexpr char kVcfPhasesetFilename[] = "test_phaseset.vcf";
constexpr char kVcfPhasesetGoldenFilename[] = "test_phaseset.vcf.golden.tfrecord";  // NOLINT
constexpr char kVcfAlleleDepthFilename[] = "test_allele_depth.vcf";
constexpr char kVcfAlleleDepthGoldenFilename[] = "test_allele_depth.vcf.golden.tfrecord";  // NOLINT
constexpr char kVcfVariantAlleleFrequencyFilename[] = "test_vaf.vcf";
constexpr char kVcfVariantAlleleFrequencyGoldenFilename[] = "test_vaf.vcf.golden.tfrecord";  // NOLINT

constexpr int CHR1_SIZE = 248956422;
constexpr int CHR2_SIZE = 242193529;
constexpr int CHR3_SIZE = 198295559;
constexpr int CHRX_SIZE = 156040895;

TEST(VcfFileOnlySites, IterationWorks) {
  std::unique_ptr<VcfReader> reader = std::move(
      VcfReader::FromFile(GetTestData(kVcfSitesFilename), VcfReaderOptions())
          .ValueOrDie());

  vector<Variant> golden =
      ReadProtosFromTFRecord<Variant>(GetTestData(kVcfSitesGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()),
              Pointwise(IgnoringFieldPaths({"info"}, EqualsProto()),
                        golden));
}

class VcfWithSamplesReaderTest : public ::testing::Test {
 protected:
  void SetUp() override {
    options_.set_index_mode(IndexHandlingMode::INDEX_BASED_ON_FILENAME);
    indexed_vcf_ = GetTestData(kVcfIndexSamplesFilename);
    golden_ =
        ReadProtosFromTFRecord<Variant>(GetTestData(kVcfSamplesGoldenFilename));
    ignored_fields_ = {"info", "calls.info"};
    RecreateReader();
  }

  void RecreateReader(VcfReaderOptions* options_override = nullptr) {
    VcfReaderOptions* options = options_override ? options_override : &options_;
    reader_ =
        std::move(VcfReader::FromFile(indexed_vcf_, *options).ValueOrDie());
  }

  VcfReaderOptions options_;
  string indexed_vcf_;
  vector<Variant> golden_;
  vector<string> ignored_fields_;
  std::unique_ptr<VcfReader> reader_;
};

TEST_F(VcfWithSamplesReaderTest, IterationWorksWithoutIndex) {
  // Checks that iterate() produces all of the variants in our golden file
  // in order for an unindexed VCF file.
  VcfReaderOptions options;
  RecreateReader(&options);
  EXPECT_THAT(
      as_vector(reader_->Iterate()),
      Pointwise(IgnoringFieldPaths(ignored_fields_, EqualsProto()), golden_));
}

TEST_F(VcfWithSamplesReaderTest, IterationWorksWithIndex) {
  // Checks that iterate() produces all of the variants in our golden file
  // in order for an indexed VCF file.
  EXPECT_THAT(
      as_vector(reader_->Iterate()),
      Pointwise(IgnoringFieldPaths(ignored_fields_, EqualsProto()), golden_));
}

TEST_F(VcfWithSamplesReaderTest, QueryWorks) {
  // Get all of the variants on chr1 from golden.
  vector<Variant> subgolden;
  for (Variant& v : golden_) {
    if (v.reference_name() == "chr1")
        subgolden.push_back(v);
  }

  // Compare the chr1 golden variants to those we get from Query.
  EXPECT_THAT(
      as_vector(reader_->Query(MakeRange("chr1", 0, CHR1_SIZE))),
      Pointwise(IgnoringFieldPaths(ignored_fields_, EqualsProto()), subgolden));
}

TEST_F(VcfWithSamplesReaderTest, QueryRangesIsCorrect) {
  // There's a variant at chr3:14319, test that query works exactly.
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 14318, 14319))),
              SizeIs(1));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 14317, 14318))),
              SizeIs(0));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 14319, 14320))),
              SizeIs(0));
  // Starts one before, runs to one after, so we get it.
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 14317, 14319))),
              SizeIs(1));
  // Starts 100 bp before and runs to 100 bp after, so we get it.
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 14217, 14419))),
              SizeIs(1));
  // End is far enough to grab a few more variants.
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 14318, 60000))),
              SizeIs(2));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 99999, 500000))),
              SizeIs(4));
  // There aren't any variants on the valid contig "chr4".
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr4", 9999, 50000))),
              SizeIs(0));
}

TEST_F(VcfWithSamplesReaderTest, WholeChromosomeQueries) {
  // Test a bunch of misc. queries.
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr1", 0, CHR1_SIZE))),
              SizeIs(711));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr2", 0, CHR2_SIZE))),
              SizeIs(34));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr3", 0, CHR3_SIZE))),
              SizeIs(6));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chrX", 0, CHRX_SIZE))),
              SizeIs(2));
}

TEST(VcfReaderLikelihoodsTest, MatchesGolden) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfLikelihoodsFilename),
                                    VcfReaderOptions())
                    .ValueOrDie());
  vector<Variant> golden = ReadProtosFromTFRecord<Variant>(
      GetTestData(kVcfLikelihoodsGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

TEST(VcfReaderPhasesetTest, MatchesGolden) {
  // Verify that we can still read the phaseset fields correctly.
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfPhasesetFilename),
                                    VcfReaderOptions())
                    .ValueOrDie());
  vector<Variant> golden =
      ReadProtosFromTFRecord<Variant>(GetTestData(kVcfPhasesetGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

TEST(VcfReaderAlleleDepthTest, MatchesGolden) {
  // Verify that we can still read the AD and DP fields correctly.
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfAlleleDepthFilename),
                                    VcfReaderOptions())
                    .ValueOrDie());
  vector<Variant> golden =
      ReadProtosFromTFRecord<Variant>(
      GetTestData(kVcfAlleleDepthGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

TEST(VcfReaderVariantAlleleFrequencyTest, MatchesGolden) {
  // Verify that we can still read the VAF field correctly.
  std::unique_ptr<VcfReader> reader = std::move(
      VcfReader::FromFile(GetTestData(kVcfVariantAlleleFrequencyFilename),
                          VcfReaderOptions())
          .ValueOrDie());
  vector<Variant> golden = ReadProtosFromTFRecord<Variant>(
      GetTestData(kVcfVariantAlleleFrequencyGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

}  // namespace nucleus
