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

#include "third_party/nucleus/io/vcf_writer.h"

#include <memory>
#include <utility>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/platform/types.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "google/protobuf/repeated_field.h"
#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/platform/env.h"

namespace nucleus {

using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using std::vector;

// TODO: we should factor out a testdata.h

// Note that test_likelihoods_output.vcf is different from
// test_likelihoods_input.vcf because our VCF writer doesn't
// allow per-position control over the GL/PL formats.
constexpr char kVcfLikelihoodsVcfOutput[] = "test_likelihoods_output.vcf";
constexpr char kVcfLikelihoodsGoldenFilename[] = "test_likelihoods.vcf.golden.tfrecord";  // NOLINT

// This is the expected output header of the DogVcfWriter defined below.
constexpr char kExpectedHeaderFmt[] =
    "##fileformat=VCFv4.2\n"
    "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
    "##FILTER=<ID=RefCall,Description=\"Most likely reference\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of "
    "the interval\">\n"  // NOLINT
    "##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP "
    "membership\",Source=\"dbSNP\",Version=\"build 129\">\n"  // NOLINT
    "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele "
    "counts\">\n"  // NOLINT
    "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Frequency of each "
    "ALT allele\">\n"  // NOLINT
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype "
    "Quality\">\n"  // NOLINT
    "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth of all "
    "passing filters reads.\">\n"  // NOLINT
    "##FORMAT=<ID=MIN_DP,Number=1,Type=Integer,Description=\"Minimum DP "
    "observed within the GVCF block.\">\n"  // NOLINT
    "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth of all "
    "passing filters reads for each allele.\">\n"  // NOLINT
    "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele "
    "fractions.\">\n"  // NOLINT
    "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype "
    "likelihoods, log10 encoded\">\n"  // NOLINT
    "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype "
    "likelihoods, Phred encoded\">\n"  // NOLINT
    "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n"
    "##PEDIGREE=<Name_1=\"Val_1\",Name_2=\"Val_2\">\n"
    "##pedigreeDb=http://my.pedigre.es\n"
    "##contig=<ID=Chr1,length=50,description=\"Dog chromosome 1\">\n"
    "##contig=<ID=Chr2,length=25>\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFido\tSpot\n";

// Build the skeleton of a VCF file for some pretend variants.
// This routine will populate headers but not any records.
std::unique_ptr<VcfWriter> MakeDogVcfWriter(
    const string& fname, const bool round_qual, const bool include_gl = true,
    const std::vector<string>& excluded_infos = {},
    const std::vector<string>& excluded_formats = {},
    bool exclude_header = false) {
  nucleus::genomics::v1::VcfHeader header;
  // FILTERs. Note that the PASS filter automatically gets added even though it
  // is not present here.
  auto& filter = *header.mutable_filters()->Add();
  filter.set_id("RefCall");
  filter.set_description("Most likely reference");

  // INFOs.
  auto& infoEnd = *header.mutable_infos()->Add();
  infoEnd.set_id("END");
  infoEnd.set_number("1");
  infoEnd.set_type("Integer");
  infoEnd.set_description("Stop position of the interval");
  auto& infodbsnp = *header.mutable_infos()->Add();
  infodbsnp.set_id("DB");
  infodbsnp.set_number("0");
  infodbsnp.set_type("Flag");
  infodbsnp.set_description("dbSNP membership");
  infodbsnp.set_source("dbSNP");
  infodbsnp.set_version("build 129");
  auto& infoAc = *header.mutable_infos()->Add();
  infoAc.set_id("AC");
  infoAc.set_number("A");
  infoAc.set_type("Integer");
  infoAc.set_description("Allele counts");
  auto& infoAf = *header.mutable_infos()->Add();
  infoAf.set_id("AF");
  infoAf.set_number("A");
  infoAf.set_type("Float");
  infoAf.set_description("Frequency of each ALT allele");

  // FORMATs.
  auto& format1 = *header.mutable_formats()->Add();
  format1.set_id("GT");
  format1.set_number("1");
  format1.set_type("String");
  format1.set_description("Genotype");
  auto& format2 = *header.mutable_formats()->Add();
  format2.set_id("GQ");
  format2.set_number("1");
  format2.set_type("Integer");
  format2.set_description("Genotype Quality");
  auto& format3 = *header.mutable_formats()->Add();
  format3.set_id("DP");
  format3.set_number("1");
  format3.set_type("Integer");
  format3.set_description("Read depth of all passing filters reads.");
  auto& format4 = *header.mutable_formats()->Add();
  format4.set_id("MIN_DP");
  format4.set_number("1");
  format4.set_type("Integer");
  format4.set_description("Minimum DP observed within the GVCF block.");
  auto& format5 = *header.mutable_formats()->Add();
  format5.set_id("AD");
  format5.set_number("R");
  format5.set_type("Integer");
  format5.set_description(
      "Read depth of all passing filters reads for each allele.");
  auto& format6 = *header.mutable_formats()->Add();
  format6.set_id("VAF");
  format6.set_number("A");
  format6.set_type("Float");
  format6.set_description("Variant allele fractions.");
  if (include_gl) {
    auto& format7 = *header.mutable_formats()->Add();
    format7.set_id("GL");
    format7.set_number("G");
    format7.set_type("Float");
    format7.set_description("Genotype likelihoods, log10 encoded");
  }
  auto& format8 = *header.mutable_formats()->Add();
  format8.set_id("PL");
  format8.set_number("G");
  format8.set_type("Integer");
  format8.set_description("Genotype likelihoods, Phred encoded");
  auto& format9 = *header.mutable_formats()->Add();
  format9.set_id("PS");
  format9.set_number("1");
  format9.set_type("Integer");
  format9.set_description("Phase set");

  // Structured extras.
  auto& sExtra = *header.mutable_structured_extras()->Add();
  sExtra.set_key("PEDIGREE");
  auto& f1 = *sExtra.mutable_fields()->Add();
  f1.set_key("Name_1");
  f1.set_value("Val_1");
  auto& f2 = *sExtra.mutable_fields()->Add();
  f2.set_key("Name_2");
  f2.set_value("Val_2");

  // Unstructured extras.
  auto& extra = *header.mutable_extras()->Add();
  extra.set_key("pedigreeDb");
  extra.set_value("http://my.pedigre.es");

  // Contigs.
  auto& contig1 = *header.mutable_contigs()->Add();
  contig1.set_name("Chr1");
  contig1.set_description("Dog chromosome 1");
  contig1.set_n_bases(50);
  contig1.set_pos_in_fasta(0);
  auto& contig2 = *header.mutable_contigs()->Add();
  contig2.set_name("Chr2");
  contig2.set_n_bases(25);
  contig2.set_pos_in_fasta(1);

  // Samples.
  header.mutable_sample_names()->Add("Fido");
  header.mutable_sample_names()->Add("Spot");

  nucleus::genomics::v1::VcfWriterOptions writer_options;
  if (round_qual) {
    writer_options.set_round_qual_values(true);
  }
  for (const string& info : excluded_infos) {
    writer_options.add_excluded_info_fields(info);
  }
  for (const string& fmt : excluded_formats) {
    writer_options.add_excluded_format_fields(fmt);
  }

  writer_options.set_exclude_header(exclude_header);

  return std::move(
      VcfWriter::ToFile(fname, header, writer_options).ValueOrDie());
}

Variant MakeVariant(const vector<string>& names, const string& refName,
                    int refStart, int refEnd, const string& refBases,
                    const vector<string>& altBases) {
  Variant v;
  for (string name : names) {
    v.mutable_names()->Add(std::move(name));
  }
  v.set_reference_name(refName);
  v.set_start(refStart);
  v.set_end(refEnd);
  v.set_reference_bases(refBases);
  for (string alt : altBases) {
    v.mutable_alternate_bases()->Add(std::move(alt));
  }
  return v;
}

VariantCall MakeVariantCall(const string& callSetName, vector<int> genotypes) {
  VariantCall vc;
  vc.set_call_set_name(callSetName);
  for (int gt : genotypes) {
    vc.add_genotype(gt);
  }
  return vc;
}

// Tests that VcfWriter infers the file open mode from file path correctly.
TEST(VcfWriterTest, OpenMode) {
  EXPECT_EQ("w", string(VcfWriter::GetOpenMode("HG002.vcf")));
  EXPECT_EQ("w", string(VcfWriter::GetOpenMode("HG002")));
  EXPECT_EQ("wz", string(VcfWriter::GetOpenMode("HG002.vcf.gz")));
  EXPECT_EQ("wb", string(VcfWriter::GetOpenMode("HG002.bcf.gz")));
  EXPECT_EQ("wbu", string(VcfWriter::GetOpenMode("HG002.bcf")));
}

TEST(VcfWriterTest, WritesVCF) {
  // This test verifies that VcfWriter writes the expected VCF file
  string output_filename = MakeTempFile("writes_vcf.vcf");
  auto writer = MakeDogVcfWriter(output_filename, false);

  // A named variant with no qual, phase
  Variant v1 = MakeVariant({"DogSNP1"}, "Chr1", 20, 21, "A", {"T"});
  *v1.add_calls() = MakeVariantCall("Fido", {0, 1});
  *v1.add_calls() = MakeVariantCall("Spot", {0, 0});
  ASSERT_THAT(writer->Write(v1), IsOK());

  // A second variant, with no name, but qual, filter
  Variant v2 = MakeVariant({}, "Chr2", 10, 11, "C", {"G", "T"});
  v2.mutable_filter()->Add("PASS");
  v2.set_quality(10);
  *v2.add_calls() = MakeVariantCall("Fido", {0, 0});
  *v2.add_calls() = MakeVariantCall("Spot", {0, 1});
  ASSERT_THAT(writer->Write(v2), IsOK());

  // Another variant with two names, some missing calls,
  // and phasing.
  Variant v3 = MakeVariant({"DogSNP3", "Woof10003"},
                            "Chr2", 15, 16, "C", {"", "T"});
  v3.mutable_filter()->Add("PASS");
  v3.set_quality(10.567);
  *v3.add_calls() = MakeVariantCall("Fido", {-1, -1});
  VariantCall call2 = MakeVariantCall("Spot", {-1, 1});
  call2.set_is_phased(true);
  SetInfoField("PS", 10, &call2);
  *v3.add_calls() = call2;
  ASSERT_THAT(writer->Write(v3), IsOK());

  // Another variant that has mixed-ploidy calls--copy number variation?
  Variant v4 = MakeVariant({"DogSNP4"}, "Chr2", 17, 18, "T", {"A"});
  v4.set_quality(-1);
  *v4.add_calls() = MakeVariantCall("Fido", {0, 1, 0});
  *v4.add_calls() = MakeVariantCall("Spot", {0, 1});
  ASSERT_THAT(writer->Write(v4), IsOK());

  // A deletion variant
  Variant v5 = MakeVariant({"DogSNP5"}, "Chr2", 19, 21, "TT", {""});
  *v5.add_calls() = MakeVariantCall("Fido", {0, 1});
  *v5.add_calls() = MakeVariantCall("Spot", {0, 0});
  ASSERT_THAT(writer->Write(v5), IsOK());

  // An insertion variant
  Variant v6 = MakeVariant({"DogSNP6"}, "Chr2", 22, 22, "", {"AAA"});
  *v6.add_calls() = MakeVariantCall("Fido", {0, 0});
  *v6.add_calls() = MakeVariantCall("Spot", {1, 0});
  ASSERT_THAT(writer->Write(v6), IsOK());

  // A variant with INFO fields.
  Variant v7 = MakeVariant({"DogSNP7"}, "Chr2", 23, 24, "A", {"T", "G"});
  *v7.add_calls() = MakeVariantCall("Fido", {0, 1});
  *v7.add_calls() = MakeVariantCall("Spot", {0, 0});
  SetInfoField("DB", std::vector<bool>{true}, &v7);
  SetInfoField("AC", std::vector<int>{1, 0}, &v7);
  SetInfoField("AF", std::vector<float>{0.75, 0.0}, &v7);
  ASSERT_THAT(writer->Write(v7), IsOK());


  // A variant with a false INFO flag, which should be coded as omitting the
  // flag from the INFO fields in the VCF line.
  Variant v8 = MakeVariant({"DogSNP8"}, "Chr2", 24, 25, "A", {"T", "G"});
  *v8.add_calls() = MakeVariantCall("Fido", {0, 1});
  *v8.add_calls() = MakeVariantCall("Spot", {0, 0});
  SetInfoField("DB", std::vector<bool>{false}, &v8);
  SetInfoField("AC", std::vector<int>{1, 0}, &v8);
  SetInfoField("AF", std::vector<float>{0.75, 0.0}, &v8);
  ASSERT_THAT(writer->Write(v8), IsOK());

  // Check that the written data is as expected.
  // (Close file to guarantee flushed to disk).
  writer.reset();

  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &vcf_contents));

  const string kExpectedVcfContent =
      string(kExpectedHeaderFmt) +
      "Chr1\t21\tDogSNP1\tA\tT\t0\t.\t.\tGT\t0/1\t0/0\n"
      "Chr2\t11\t.\tC\tG,T\t10\tPASS\t.\tGT\t0/0\t0/1\n"
      "Chr2\t16\tDogSNP3;Woof10003\tC\t,T\t10.567\tPASS\t.\tGT:PS\t./"
      ".:.\t.|1:10\n"
      "Chr2\t18\tDogSNP4\tT\tA\t.\t.\t.\tGT\t0/1/0\t0/1\n"
      "Chr2\t20\tDogSNP5\tTT\t\t0\t.\t.\tGT\t0/1\t0/0\n"
      "Chr2\t23\tDogSNP6\t\tAAA\t0\t.\t.\tGT\t0/0\t1/0\n"
      "Chr2\t24\tDogSNP7\tA\tT,G\t0\t.\tDB;AC=1,0;AF=0.75,0\tGT\t0/1\t0/0\n"
      "Chr2\t25\tDogSNP8\tA\tT,G\t0\t.\tAC=1,0;AF=0.75,0\tGT\t0/1\t0/0\n";
  EXPECT_EQ(kExpectedVcfContent, vcf_contents);
}

TEST(VcfWriterTest, WritesVCFWithoutHeader) {
  string output_filename = MakeTempFile("writes_vcf.vcf");
  auto writer = MakeDogVcfWriter(output_filename, false, true, {}, {},
                                 /*exclude_header=*/true);

  Variant v1 = MakeVariant({"DogSNP1"}, "Chr1", 20, 21, "A", {"T"});
  *v1.add_calls() = MakeVariantCall("Fido", {0, 1});
  *v1.add_calls() = MakeVariantCall("Spot", {0, 0});
  ASSERT_THAT(writer->Write(v1), IsOK());

  writer = nullptr;
  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &vcf_contents));

  const string kExpectedVcfContent =
      "Chr1\t21\tDogSNP1\tA\tT\t0\t.\t.\tGT\t0/1\t0/0\n";
  EXPECT_EQ(kExpectedVcfContent, vcf_contents);
}

TEST(VcfWriterTest, RoundsVCFQuals) {
  // This test verifies that VcfWriter writes the expected VCF file with rounded
  // quality values.
  string output_filename = MakeTempFile("rounds_qualities.vcf");
  auto writer = MakeDogVcfWriter(output_filename, true);

  // A named variant with no qual, phase
  Variant v1 = MakeVariant({"DogSNP1"}, "Chr1", 20, 21, "A", {"T"});
  v1.set_quality(10.44999);
  *v1.add_calls() = MakeVariantCall("Fido", {0, 1});
  *v1.add_calls() = MakeVariantCall("Spot", {0, 0});
  ASSERT_THAT(writer->Write(v1), IsOK());

  Variant v2 = MakeVariant({}, "Chr2", 10, 11, "C", {"G", "T"});
  v2.mutable_filter()->Add("PASS");
  v2.set_quality(10.4500);
  *v2.add_calls() = MakeVariantCall("Fido", {0, 0});
  *v2.add_calls() = MakeVariantCall("Spot", {0, 1});
  ASSERT_THAT(writer->Write(v2), IsOK());

  // Check that the written data is as expected.
  // (close file to guarantee flushed to disk)
  writer.reset();

  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &vcf_contents));

  const string kExpectedVcfContent =
      string(kExpectedHeaderFmt) +
      "Chr1\t21\tDogSNP1\tA\tT\t10.4\t.\t.\tGT\t0/1\t0/0\n"
      "Chr2\t11\t.\tC\tG,T\t10.5\tPASS\t.\tGT\t0/0\t0/1\n";

  EXPECT_EQ(kExpectedVcfContent, vcf_contents);
}

TEST(VcfWriterTest, ExcludesFields) {
  // This test verifies that VcfWriter writes the expected VCF file with
  // INFO and FORMAT fields excluded.
  string output_filename = MakeTempFile("excluded_fields.vcf");
  auto writer = MakeDogVcfWriter(output_filename, true, true,
                                 {"DB", "AF"}, {"GQ"});

  // A named variant with no qual, phase
  Variant v1 = MakeVariant({"DogSNP1"}, "Chr1", 20, 21, "A", {"T"});
  v1.set_quality(10.44999);
  SetInfoField("DB", std::vector<bool>{false}, &v1);
  SetInfoField("AC", std::vector<int>{1, 0}, &v1);
  SetInfoField("AF", std::vector<float>{0.75, 0.0}, &v1);

  *v1.add_calls() = MakeVariantCall("Fido", {0, 1});
  VariantCall call2 = MakeVariantCall("Spot", {0, 0});
  SetInfoField("GQ", std::vector<float>{15}, &call2);
  *v1.add_calls() = call2;
  ASSERT_THAT(writer->Write(v1), IsOK());

  Variant v2 = MakeVariant({}, "Chr2", 10, 11, "C", {"G", "T"});
  v2.mutable_filter()->Add("PASS");
  v2.set_quality(10.4500);
  *v2.add_calls() = MakeVariantCall("Fido", {0, 0});
  *v2.add_calls() = MakeVariantCall("Spot", {0, 1});
  SetInfoField("DB", std::vector<bool>{true}, &v2);
  ASSERT_THAT(writer->Write(v2), IsOK());

  // Check that the written data is as expected.
  // (close file to guarantee flushed to disk)
  writer.reset();

  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &vcf_contents));

  const string kExpectedVcfContent =
      string(kExpectedHeaderFmt) +
      "Chr1\t21\tDogSNP1\tA\tT\t10.4\t.\tAC=1,0\tGT\t0/1\t0/0\n"
      "Chr2\t11\t.\tC\tG,T\t10.5\tPASS\t.\tGT\t0/0\t0/1\n";

  EXPECT_EQ(kExpectedVcfContent, vcf_contents);
}

TEST(VcfWriterTest, WritesVCFWithLikelihoods) {
  std::vector<Variant> variants = ReadProtosFromTFRecord<Variant>(
      GetTestData(kVcfLikelihoodsGoldenFilename));
  string out_fname = MakeTempFile("likelihoods_out.vcf");
  auto writer = MakeDogVcfWriter(out_fname, false, false);
  for (const auto& variant : variants) {
    ASSERT_THAT(writer->Write(variant), IsOK());
  }
  ASSERT_THAT(writer->Close(), IsOK());

  string expected_vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(
      tensorflow::Env::Default(), GetTestData(kVcfLikelihoodsVcfOutput),
      &expected_vcf_contents));

  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           out_fname, &vcf_contents));

  EXPECT_EQ(expected_vcf_contents, vcf_contents);
}

TEST(VcfWriterTest, WritesGzippedVCF) {
  string output_filename = MakeTempFile("writes_gzipped_vcf.vcf.gz");
  auto writer = MakeDogVcfWriter(output_filename, false);
  ASSERT_THAT(writer->Close(), IsOK());

  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &vcf_contents));
  EXPECT_THAT(IsGzipped(vcf_contents),
              "VCF writer should be able to writed gzipped output");
}

TEST(VcfWriterTest, HandlesRedefinedPL) {
  string output_filename = MakeTempFile("redefined_pl.vcf");
  nucleus::genomics::v1::VcfHeader header;

  // INFOs.
  auto& infoEnd = *header.mutable_infos()->Add();
  infoEnd.set_id("END");
  infoEnd.set_number("1");
  infoEnd.set_type("Integer");
  infoEnd.set_description("Stop position of the interval");

  // FORMATs.
  auto& format1 = *header.mutable_formats()->Add();
  format1.set_id("GT");
  format1.set_number("1");
  format1.set_type("String");
  format1.set_description("Genotype");

  auto& format2 = *header.mutable_formats()->Add();
  format2.set_id("PL");
  format2.set_number("1");
  format2.set_type("Integer");
  format2.set_description("Custom PL");

  // Contigs.
  auto& contig1 = *header.mutable_contigs()->Add();
  contig1.set_name("Chr1");
  contig1.set_description("Dog chromosome 1");
  contig1.set_n_bases(50);
  contig1.set_pos_in_fasta(0);

  // Samples.
  header.mutable_sample_names()->Add("Fido");

  nucleus::genomics::v1::VcfWriterOptions writer_options;
  writer_options.set_retrieve_gl_and_pl_from_info_map(true);

  std::unique_ptr<VcfWriter> writer = std::move(
      VcfWriter::ToFile(output_filename, header, writer_options).ValueOrDie());

  Variant v = MakeVariant({}, "Chr1", 10, 11, "C", {"G", "T"});
  v.mutable_filter()->Add("PASS");
  v.set_quality(10.5);
  auto call = MakeVariantCall("Fido", {0, 1});
  nucleus::genomics::v1::ListValue lv;
  nucleus::genomics::v1::Value* val = lv.add_values();
  val->set_int_value(42);
  (*call.mutable_info())["PL"] = lv;
  *v.add_calls() = call;
  ASSERT_THAT(writer->Write(v), IsOK());

  // Check that the written data is as expected.
  // (Close file to guarantee flushed to disk).
  writer.reset();

  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &vcf_contents));

  const string kExpectedVcfContent =
      "##fileformat=VCFv4.2\n"
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
      "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of "
      "the interval\">\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "##FORMAT=<ID=PL,Number=1,Type=Integer,Description=\"Custom PL\">\n"
      "##contig=<ID=Chr1,length=50,description=\"Dog chromosome 1\">\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFido\n"
      "Chr1\t11\t.\tC\tG,T\t10.5\tPASS\t.\tGT:PL\t0/1:42\n";

  EXPECT_EQ(kExpectedVcfContent, vcf_contents);
}

}  // namespace nucleus
