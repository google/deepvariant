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
 *
 */

#include "third_party/nucleus/io/vcf_reader.h"

#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace nucleus {

using std::vector;

using ::testing::Eq;
using ::testing::Not;
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

constexpr char kValidVcfHeaderFilename[] = "test_valid_vcf_header_parsing.vcf";
constexpr char kInvalidVcfHeaderFilename[] =
    "test_invalid_vcf_header_parsing.vcf";
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

void AddTestContig(nucleus::genomics::v1::VcfHeader& header, const string& name,
                   const string& description = "", const int n_bases = 0,
                   const int pos_in_fasta = 0,
                   const std::vector<string>& kvExtra = {}) {
  nucleus::genomics::v1::ContigInfo* contig = header.add_contigs();
  contig->set_name(name);
  contig->set_description(description);
  contig->set_n_bases(n_bases);
  contig->set_pos_in_fasta(pos_in_fasta);
  for (int i = 0; i < kvExtra.size(); i += 2) {
    (*contig->mutable_extra())[kvExtra[i]] = kvExtra[i + 1];
  }
}

void AddTestFilter(nucleus::genomics::v1::VcfHeader& header, const string& id,
                   const string& description = "") {
  nucleus::genomics::v1::VcfFilterInfo* filter = header.add_filters();
  filter->set_id(id);
  filter->set_description(description);
}

void AddTestInfo(nucleus::genomics::v1::VcfHeader& header, const string& id,
                 const string& number, const string& type,
                 const string& description, const string& source = "",
                 const string& version = "") {
  nucleus::genomics::v1::VcfInfo* info = header.add_infos();
  info->set_id(id);
  info->set_number(number);
  info->set_type(type);
  info->set_description(description);
  info->set_source(source);
  info->set_version(version);
}

void AddTestFormat(nucleus::genomics::v1::VcfHeader& header, const string& id,
                   const string& number, const string& type,
                   const string& description) {
  nucleus::genomics::v1::VcfFormatInfo* format = header.add_formats();
  format->set_id(id);
  format->set_number(number);
  format->set_type(type);
  format->set_description(description);
}

void AddTestStructuredExtra(nucleus::genomics::v1::VcfHeader& header,
                            const string& key,
                            const std::vector<string>& extra_pairs) {
  nucleus::genomics::v1::VcfStructuredExtra* extra =
      header.add_structured_extras();
  extra->set_key(key);
  for (int i = 0; i < extra_pairs.size(); i += 2) {
    nucleus::genomics::v1::VcfExtra& toAdd = *extra->mutable_fields()->Add();
    toAdd.set_key(extra_pairs[i]);
    toAdd.set_value(extra_pairs[i + 1]);
  }
}

void AddTestExtra(nucleus::genomics::v1::VcfHeader& header, const string& key,
                  const string& value) {
  nucleus::genomics::v1::VcfExtra* extra = header.add_extras();
  extra->set_key(key);
  extra->set_value(value);
}

TEST(ValidVcfHeaderParsing, MatchesProto) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kValidVcfHeaderFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());
  nucleus::genomics::v1::VcfHeader header_proto;
  header_proto.set_fileformat("VCFv4.2");
  AddTestFilter(header_proto, "PASS", "All filters passed");
  AddTestInfo(header_proto, "DP", "1", "Integer",
              "Read depth of all samples summed together.");
  AddTestInfo(header_proto, "DB", "0", "Flag", "dbSNP membership", "dbSNP",
              "build 129");
  AddTestFormat(header_proto, "GT", "1", "String", "Genotype");
  AddTestFormat(header_proto, "GQ", "1", "Integer", "Genotype Quality");
  AddTestFormat(header_proto, "AD", "R", "Integer",
                "Read depth of all passing filters reads for each allele.");
  AddTestFormat(header_proto, "GL", "G", "Float",
                "Genotype likelihoods, log10 encoded");
  AddTestContig(header_proto, "Chr1", "", 50);
  AddTestContig(header_proto, "chr2", "", 81195210, 1,
                {"URL", "ftp://somewhere.org/assembly.fa", "md5", "fakemd5",
                 "species", "Homo sapiens"});
  AddTestStructuredExtra(header_proto, "META",
                         {"ID", "Assay", "Type", "String", "Number", ".",
                          "Values", "[WholeGenome, Exome]"});
  AddTestStructuredExtra(
      header_proto, "PEDIGREE",
      {"Name_0", "G0-ID", "Name_1", "G1-ID", "Name_3", "GN-ID"});
  AddTestExtra(header_proto, "pedigreeDB", "http://url.to.pedigre.es/search");
  header_proto.add_sample_names("Fido");
  header_proto.add_sample_names("Spot");
  EXPECT_THAT(reader->Header(), EqualsProto(header_proto));
}

TEST(VcfFileOnlySites, IterationWorks) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfSitesFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());

  vector<Variant> golden =
      ReadProtosFromTFRecord<Variant>(GetTestData(kVcfSitesGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

class VcfWithSamplesReaderTest : public ::testing::Test {
 protected:
  void SetUp() override {
    options_.set_index_mode(
        nucleus::genomics::v1::IndexHandlingMode::INDEX_BASED_ON_FILENAME);
    indexed_vcf_ = GetTestData(kVcfIndexSamplesFilename);
    golden_ =
        ReadProtosFromTFRecord<Variant>(GetTestData(kVcfSamplesGoldenFilename));
    RecreateReader();
  }

  void RecreateReader(
      nucleus::genomics::v1::VcfReaderOptions* options_override = nullptr) {
    nucleus::genomics::v1::VcfReaderOptions* options =
        options_override ? options_override : &options_;
    reader_ =
        std::move(VcfReader::FromFile(indexed_vcf_, *options).ValueOrDie());
  }

  nucleus::genomics::v1::VcfReaderOptions options_;
  string indexed_vcf_;
  vector<Variant> golden_;
  std::unique_ptr<VcfReader> reader_;
};

TEST_F(VcfWithSamplesReaderTest, IterationWorksWithoutIndex) {
  // Checks that iterate() produces all of the variants in our golden file
  // in order for an unindexed VCF file.
  nucleus::genomics::v1::VcfReaderOptions options;
  RecreateReader(&options);
  EXPECT_THAT(
      as_vector(reader_->Iterate()), Pointwise(EqualsProto(), golden_));
}

TEST_F(VcfWithSamplesReaderTest, IterationWorksWithIndex) {
  // Checks that iterate() produces all of the variants in our golden file
  // in order for an indexed VCF file.
  EXPECT_THAT(as_vector(reader_->Iterate()), Pointwise(EqualsProto(), golden_));
}

TEST_F(VcfWithSamplesReaderTest, FilteringInfoFieldsWorks) {
  // Checks that iterate() filters FORMAT fields out as we expect.
  nucleus::genomics::v1::VcfReaderOptions options;
  RecreateReader(&options);
  std::vector<Variant> all_infos = as_vector(reader_->Iterate());

  // No genotype (a special-cased field in all records in the test data).
  options.add_excluded_info_fields("AC");
  RecreateReader(&options);
  std::vector<Variant> no_ac = as_vector(reader_->Iterate());

  // Cut the golden data down to the first 700 records that all contain AC
  // so that pointwise comparisons work.
  int numRecordsToCompare = 700;
  golden_.resize(numRecordsToCompare);
  all_infos.resize(numRecordsToCompare);
  no_ac.resize(numRecordsToCompare);

  EXPECT_THAT(all_infos, Pointwise(EqualsProto(), golden_));

  EXPECT_THAT(no_ac, Pointwise(Not(EqualsProto()), golden_));
  EXPECT_THAT(no_ac,
              Pointwise(IgnoringFieldPaths({"info"}, EqualsProto()), golden_));
}

TEST_F(VcfWithSamplesReaderTest, FilteringFormatFieldsWorks) {
  // Checks that iterate() filters FORMAT fields out as we expect.
  nucleus::genomics::v1::VcfReaderOptions options;
  RecreateReader(&options);
  std::vector<Variant> all_formats = as_vector(reader_->Iterate());

  // No genotype (a special-cased field in all records in the test data).
  options.add_excluded_format_fields("GT");
  RecreateReader(&options);
  std::vector<Variant> no_genotype = as_vector(reader_->Iterate());

  // No GQ (a normal FORMAT field set in all records in the test data).
  options.Clear();
  options.add_excluded_format_fields("GQ");
  RecreateReader(&options);
  std::vector<Variant> no_gq = as_vector(reader_->Iterate());

  // No genotype or GQ.
  options.Clear();
  options.add_excluded_format_fields("GT");
  options.add_excluded_format_fields("GQ");
  RecreateReader(&options);
  std::vector<Variant> no_gt_or_gq = as_vector(reader_->Iterate());

  EXPECT_THAT(all_formats, Pointwise(EqualsProto(), golden_));

  EXPECT_THAT(no_genotype, Pointwise(Not(EqualsProto()), golden_));
  EXPECT_THAT(no_genotype,
              Pointwise(IgnoringFieldPaths({"calls.genotype"}, EqualsProto()),
                        golden_));

  EXPECT_THAT(no_gq, Pointwise(Not(EqualsProto()), golden_));
  EXPECT_THAT(
      no_gq,
      Pointwise(IgnoringFieldPaths({"calls.info"}, EqualsProto()), golden_));

  EXPECT_THAT(no_gt_or_gq, Pointwise(Not(EqualsProto()), golden_));
  EXPECT_THAT(no_gt_or_gq,
              Pointwise(Not(IgnoringFieldPaths({"calls.info"}, EqualsProto())),
                        golden_));
  EXPECT_THAT(
      no_gt_or_gq,
      Pointwise(Not(IgnoringFieldPaths({"calls.genotype"}, EqualsProto())),
                golden_));
  EXPECT_THAT(no_gt_or_gq,
              Pointwise(IgnoringFieldPaths({"calls.genotype", "calls.info"},
                                           EqualsProto()),
                        golden_));
}

TEST_F(VcfWithSamplesReaderTest, QueryWorks) {
  // Get all of the variants on chr1 from golden.
  vector<Variant> subgolden;
  for (Variant& v : golden_) {
    if (v.reference_name() == "chr1")
        subgolden.push_back(v);
  }

  // Compare the chr1 golden variants to those we get from Query.
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr1", 0, CHR1_SIZE))),
              Pointwise(EqualsProto(), subgolden));
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
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());
  vector<Variant> golden = ReadProtosFromTFRecord<Variant>(
      GetTestData(kVcfLikelihoodsGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

TEST(VcfReaderPhasesetTest, MatchesGolden) {
  // Verify that we can still read the phaseset fields correctly.
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfPhasesetFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());
  vector<Variant> golden =
      ReadProtosFromTFRecord<Variant>(GetTestData(kVcfPhasesetGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

TEST(VcfReaderAlleleDepthTest, MatchesGolden) {
  // Verify that we can still read the AD and DP fields correctly.
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfAlleleDepthFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
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
                          nucleus::genomics::v1::VcfReaderOptions())
          .ValueOrDie());
  vector<Variant> golden = ReadProtosFromTFRecord<Variant>(
      GetTestData(kVcfVariantAlleleFrequencyGoldenFilename));
  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), golden));
}

}  // namespace nucleus
