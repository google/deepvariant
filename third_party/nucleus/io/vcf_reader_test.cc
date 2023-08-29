/*
 * Copyright 2018 Google LLC.
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

#include <stddef.h>

#include <utility>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "google/protobuf/map.h"
#include "google/protobuf/repeated_field.h"
#include "tensorflow/core/lib/core/status.h"

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
constexpr char kVcfNoHeaderFilename[] = "vcf_no_header.vcf";
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
constexpr char kRedefinedFormatFilename[] = "test_redefined_formats.vcf";

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
  for (size_t i = 0; i < kvExtra.size(); i += 2) {
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
  for (size_t i = 0; i < extra_pairs.size(); i += 2) {
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

TEST(VcfReaderTest, FromFileWithHeader) {
  std::unique_ptr<VcfReader> reader_with_header =
      std::move(VcfReader::FromFile(GetTestData(kVcfSitesFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());

  std::unique_ptr<VcfReader> reader = std::move(
      VcfReader::FromFileWithHeader(GetTestData(kVcfNoHeaderFilename),
                                    nucleus::genomics::v1::VcfReaderOptions(),
                                    reader_with_header->Header())
          .ValueOrDie());
  EXPECT_THAT(reader->Header(), EqualsProto(reader_with_header->Header()));
  vector<Variant> expected = as_vector(reader_with_header->Iterate());

  EXPECT_THAT(as_vector(reader->Iterate()), Pointwise(EqualsProto(), expected));
}

// Fails when a VcfHeader is specified and the vcf file already has a header.
TEST(VcfReaderTest, FromFileWithHeaderFailure) {
  std::unique_ptr<VcfReader> header_reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfSitesFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());

  StatusOr<std::unique_ptr<VcfReader>> status_or =
      VcfReader::FromFileWithHeader(GetTestData(kVcfSitesFilename),
                                    nucleus::genomics::v1::VcfReaderOptions(),
                                    header_reader->Header());
  EXPECT_THAT(status_or, IsNotOKWithMessage("Unexpected header in"));
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
  AddTestFormat(header_proto, "CH", "1", "Character", "Character FORMAT field");
  AddTestFormat(header_proto, "FU", ".", "String",
                "Multiple Strings FORMAT field");
  AddTestFormat(header_proto, "F2", "2", "Character",
                "2 Characters FORMAT field");
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
  vector<Variant> actual = as_vector(reader->Iterate());
  ASSERT_EQ(2u, actual.size());
  Variant expected;
  expected.set_reference_name("Chr1");
  expected.set_start(20);
  expected.set_end(21);
  expected.add_names("DogSNP1");
  expected.set_reference_bases("A");
  expected.add_alternate_bases("T");

  nucleus::genomics::v1::Value* val;
  {
    nucleus::genomics::v1::VariantCall* call1 = expected.add_calls();
    call1->set_call_set_name("Fido");
    call1->add_genotype(0);
    call1->add_genotype(1);

    nucleus::genomics::v1::ListValue lv;
    val = lv.add_values();
    val->set_string_value("a1");
    val = lv.add_values();
    val->set_string_value("b1");
    (*call1->mutable_info())["FU"] = lv;

    nucleus::genomics::v1::ListValue lv2;
    val = lv2.add_values();
    val->set_string_value("a");
    val = lv2.add_values();
    val->set_string_value("b");
    (*call1->mutable_info())["F2"] = lv2;
  }
  {
    nucleus::genomics::v1::VariantCall* call2 = expected.add_calls();
    call2->set_call_set_name("Spot");
    call2->add_genotype(0);
    call2->add_genotype(1);
    nucleus::genomics::v1::ListValue lv;
    val = lv.add_values();
    val->set_int_value(42);
    (*call2->mutable_info())["GQ"] = lv;

    nucleus::genomics::v1::ListValue lv2;
    val = lv2.add_values();
    val->set_string_value("c1");
    val = lv2.add_values();
    val->set_string_value("d1");
    (*call2->mutable_info())["FU"] = lv2;

    nucleus::genomics::v1::ListValue lv3;
    val = lv3.add_values();
    val->set_string_value("c");
    val = lv3.add_values();
    val->set_string_value("d");
    (*call2->mutable_info())["F2"] = lv3;
  }
  // Ignore record-level info fields (DP and DB).
  actual[0].clear_info();
  EXPECT_THAT(actual[0], EqualsProto(expected));
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

TEST(VcfReaderFromStringTest, MatchesGolden) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(
          GetTestData(kVcfPhasesetFilename),
          nucleus::genomics::v1::VcfReaderOptions()).ValueOrDie());
  vector<Variant> golden =
      ReadProtosFromTFRecord<Variant>(GetTestData(kVcfPhasesetGoldenFilename));
  vector<Variant> parsed(5);
  NUCLEUS_CHECK_OK(reader->FromString(
      "Chr1\t21\tDogSNP1\tA\tT\t0\t.\t.\tGT:GQ\t0/1:.\t0/1:42", &(parsed[0])));
  NUCLEUS_CHECK_OK(reader->FromString(
      "Chr1\t22\tDogSNP2\tA\tT\t0\t.\t.\tGT:PL\t0/1:.\t0|1:50,40,60",
      &(parsed[1])));
  NUCLEUS_CHECK_OK(
      reader->FromString("Chr1\t23\tDogSNP3\tA\tT\t0\t.\t.\tGT:GL:PS\t"
                         "0/1:.:.\t0/1:-5.0,-4.0,-6.0:.",
                         &(parsed[2])));
  NUCLEUS_CHECK_OK(
      reader->FromString("Chr1\t24\tDogSNP4\tA\tT\t0\t.\t.\tGT:PL:PS\t"
                         "0|1:50,40,60:24\t0|1:50,40,60:.",
                         &(parsed[3])));
  NUCLEUS_CHECK_OK(
      reader->FromString("Chr1\t25\tDogSNP5\ta\tt\t0\t.\t.\tGT:GQ:PS:PL\t"
                         "0|1:42:24:50,40,60\t1|1:42:.:50,40,60",
                         &(parsed[4])));
  EXPECT_THAT(parsed, Pointwise(EqualsProto(), golden));
}

TEST(VcfReaderMultipleNamesTest, SplitsOnSemicolon) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(
          GetTestData(kVcfPhasesetFilename),
          nucleus::genomics::v1::VcfReaderOptions()).ValueOrDie());
  Variant v;
  NUCLEUS_CHECK_OK(reader->FromString(
      "Chr1\t21\tDogSNP1;CatSNP2\tA\tT\t0\t.\t.\tGT:GQ\t0/1:.\t0/1:42", &v));
  ASSERT_EQ(2, v.names_size());
  EXPECT_EQ("DogSNP1", v.names(0));
  EXPECT_EQ("CatSNP2", v.names(1));
}

TEST(VcfReaderDifferentPloidyTest, Simple) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(
          GetTestData(kVcfPhasesetFilename),
          nucleus::genomics::v1::VcfReaderOptions()).ValueOrDie());
  Variant v;
  NUCLEUS_CHECK_OK(reader->FromString(
      "Chr1\t21\t.\tC\t<SYMBOLIC>\t49\t.\t.\tGT:PS:GQ\t0|1:1:45\t.:.:.", &v));
  ASSERT_EQ(2, v.calls_size());
  EXPECT_THAT(v.calls(0).genotype(), testing::ElementsAre(0, 1));
  EXPECT_THAT(v.calls(1).genotype(), testing::ElementsAre(-1));

  Variant v2;
  NUCLEUS_CHECK_OK(reader->FromString(
      "Chr1\t21\t.\tC\t<SYMBOLIC>\t49\t.\t.\tGT:PS:GQ\t0|1:1:45\t0:24:42",
      &v2));
  ASSERT_EQ(2, v2.calls_size());
  EXPECT_THAT(v2.calls(0).genotype(), testing::ElementsAre(0, 1));
  EXPECT_THAT(v2.calls(1).genotype(), testing::ElementsAre(0));

}

TEST(HandleRedefinedFormatFields, MatchesProto) {
  nucleus::genomics::v1::VcfReaderOptions options;
  options.set_store_gl_and_pl_in_info_map(true);
  std::unique_ptr<VcfReader> reader = std::move(
      VcfReader::FromFile(GetTestData(kRedefinedFormatFilename), options)
          .ValueOrDie());
  vector<Variant> actual = as_vector(reader->Iterate());
  EXPECT_EQ(1, actual.size());
  Variant expected;
  expected.set_reference_name("Chr1");
  expected.set_start(20);
  expected.set_end(21);
  expected.add_names("DogSNP1");
  expected.set_reference_bases("A");
  expected.add_alternate_bases("T");
  expected.set_quality(10);

  nucleus::genomics::v1::VariantCall* call1 = expected.add_calls();
  call1->set_call_set_name("Fido");
  call1->add_genotype(0);
  call1->add_genotype(1);

  nucleus::genomics::v1::VariantCall* call2 = expected.add_calls();
  call2->set_call_set_name("Spot");
  call2->add_genotype(0);
  call2->add_genotype(1);
  nucleus::genomics::v1::ListValue lv;
  nucleus::genomics::v1::Value* val = lv.add_values();
  val->set_int_value(42);
  (*call2->mutable_info())["GL"] = lv;

  EXPECT_THAT(actual[0], EqualsProto(expected));
}

// Tests that parsing succeeds even when undefined header fields (info, contig,
// filter and format) exist.
TEST(VcfReaderTest, MissingHeaderDefinitions) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfPhasesetFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());
  Variant v1;
  // AB is an undefined FORMAT tag.
  NUCLEUS_CHECK_OK(reader->FromString(
      "Chr1\t21\tDogSNP1\tA\tT\t0\t.\t.\tGT:GQ:AB\t0/1:.:abc\t0/1:42:def",
      &v1));
  EXPECT_EQ("abc", v1.calls(0).info().at("AB").values(0).string_value());
  EXPECT_EQ("def", v1.calls(1).info().at("AB").values(0).string_value());
  Variant v2;
  // q10 is an undefined FILTER tag.
  NUCLEUS_CHECK_OK(reader->FromString(
      "Chr1\t22\tDogSNP2\tA\tT\t0\tq10\t.\tGT:PL\t0/1:.\t0|1:50,40,60", &v2));
  EXPECT_EQ("q10", v2.filter(0));

  Variant v3;
  // Chr2 is an undefined contig.
  NUCLEUS_CHECK_OK(
      reader->FromString("Chr2\t23\tDogSNP3\tA\tT\t0\t.\t.\tGT:GL:PS\t"
                         "0/1:.:.\t0/1:-5.0,-4.0,-6.0:.",
                         &v3));
  EXPECT_EQ("Chr2", v3.reference_name());
}

TEST(VcfReaderTest, HeaderAccessor) {
  std::unique_ptr<VcfReader> reader =
      std::move(VcfReader::FromFile(GetTestData(kVcfSamplesFilename),
                                    nucleus::genomics::v1::VcfReaderOptions())
                    .ValueOrDie());
  // A test to make sure header's fields can be iterated in a loop.
  for (const auto& name : reader->Header().sample_names()) {
    EXPECT_FALSE(name.empty());
  }
}

}  // namespace nucleus
