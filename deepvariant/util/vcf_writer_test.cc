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

#include "deepvariant/util/vcf_writer.h"

#include <memory>
#include <vector>

#include "deepvariant/util/genomics/variants.pb.h"
#include "deepvariant/util/protos/core.pb.h"
#include "deepvariant/util/test_utils.h"
#include "deepvariant/util/utils.h"
#include "deepvariant/util/vendor/status_matchers.h"

#include "tensorflow/core/lib/core/status.h"
#include "tensorflow/core/lib/core/stringpiece.h"
#include "tensorflow/core/platform/env.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace learning {
namespace genomics {
namespace core {

using tensorflow::StringPiece;
using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;
using std::vector;

// redacted

// Note that test_likelihoods_output.vcf is different from
// test_likelihoods_input.vcf because our VCF writer writes genotype likelihoods
// to the GL format field (whereas, input given to the reader could have used
// the PL format field, instead).
constexpr char kVcfLikelihoodsVcfOutput[] = "test_likelihoods_output.vcf";
constexpr char kVcfLikelihoodsGoldenFilename[] = "test_likelihoods.vcf.golden.tfrecord";  // NOLINT


// Build the skeleton of a VCF file for some pretend variants.
// This routine will populated headers but not any records.
std::unique_ptr<VcfWriter> MakeDogVcfWriter(StringPiece fname,
                                            const bool round_qual) {
  VcfWriterOptions writer_options;
  auto& contig1 = *writer_options.mutable_contigs()->Add();
  contig1.set_name("Chr1");
  contig1.set_description("Dog chromosome 1");
  contig1.set_n_bases(50);
  contig1.set_pos_in_fasta(0);
  auto& contig2 = *writer_options.mutable_contigs()->Add();
  contig2.set_name("Chr2");
  contig2.set_description("Dog chromosome 2");
  contig2.set_n_bases(25);
  contig2.set_pos_in_fasta(1);
  writer_options.mutable_sample_names()->Add("Fido");
  writer_options.mutable_sample_names()->Add("Spot");
  if (round_qual) {
    writer_options.set_round_qual_values(true);
  }

  return std::move(
      VcfWriter::ToFile(fname.ToString(), writer_options).ValueOrDie());
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

VariantCall MakeVariantCall(StringPiece callSetName, vector<int> genotypes) {
  VariantCall vc;
  vc.set_call_set_name(callSetName.ToString());
  for (int gt : genotypes) {
    vc.add_genotype(gt);
  }
  return vc;
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
  call2.set_phaseset("*");
  *v3.add_calls() = call2;
  ASSERT_THAT(writer->Write(v3), IsOK());

  // Another variant that has mixed-ploidy calls--copy number variation?
  Variant v4 = MakeVariant({"DogSNP4"}, "Chr2", 17, 18, "T", {"A"});
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

  // Check that the written data is as expected.
  // (close file to guarantee flushed to disk)
  writer.reset();

  string vcf_contents;
  TF_CHECK_OK(tensorflow::ReadFileToString(tensorflow::Env::Default(),
                                           output_filename, &vcf_contents));

  const char kExpectedVcfContent[] =
      "##fileformat=VCFv4.2\n"
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype "
      "Quality\">\n"  // NOLINT
      "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth of all "
      "passing filters reads.\">\n"  // NOLINT
      "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth of all "
      "passing filters reads for each allele.\">\n"  // NOLINT
      "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele "
      "fractions.\">\n"  // NOLINT
      "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods, "
      "log10 encoded\">\n"  // NOLINT
      "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype "
      "likelihoods, Phred encoded\">\n"  // NOLINT
      "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the "
      "interval\">\n"  // NOLINT
      "##contig=<ID=Chr1,length=50>\n"
      "##contig=<ID=Chr2,length=25>\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFido\tSpot\n"
      "Chr1\t21\tDogSNP1\tA\tT\t0\t.\t.\tGT\t0/1\t0/0\n"
      "Chr2\t11\t.\tC\tG,T\t10\tPASS\t.\tGT\t0/0\t0/1\n"
      "Chr2\t16\tDogSNP3;Woof10003\tC\t,T\t10.567\tPASS\t.\tGT\t./.\t.|1\n"
      "Chr2\t18\tDogSNP4\tT\tA\t0\t.\t.\tGT\t0/1/0\t0/1\n"
      "Chr2\t20\tDogSNP5\tTT\t\t0\t.\t.\tGT\t0/1\t0/0\n"
      "Chr2\t23\tDogSNP6\t\tAAA\t0\t.\t.\tGT\t0/0\t1/0\n";

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

  const char kExpectedVcfContent[] =
      "##fileformat=VCFv4.2\n"
      "##FILTER=<ID=PASS,Description=\"All filters passed\">\n"
      "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
      "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype "
      "Quality\">\n"  // NOLINT
      "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth of all "
      "passing filters reads.\">\n"  // NOLINT
      "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth of all "
      "passing filters reads for each allele.\">\n"  // NOLINT
      "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele "
      "fractions.\">\n"  // NOLINT
      "##FORMAT=<ID=GL,Number=G,Type=Float,Description=\"Genotype likelihoods, "
      "log10 encoded\">\n"  // NOLINT
      "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Genotype "
      "likelihoods, Phred encoded\">\n"  // NOLINT
      "##INFO=<ID=END,Number=1,Type=Integer,Description=\"Stop position of the "
      "interval\">\n"  // NOLINT
      "##contig=<ID=Chr1,length=50>\n"
      "##contig=<ID=Chr2,length=25>\n"
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tFido\tSpot\n"
      "Chr1\t21\tDogSNP1\tA\tT\t10.4\t.\t.\tGT\t0/1\t0/0\n"
      "Chr2\t11\t.\tC\tG,T\t10.5\tPASS\t.\tGT\t0/0\t0/1\n";

  EXPECT_EQ(kExpectedVcfContent, vcf_contents);
}

TEST(VcfWriterTest, WritesVCFWithLikelihoods) {
  std::vector<Variant> variants = ReadProtosFromTFRecord<Variant>(
      GetTestData(kVcfLikelihoodsGoldenFilename));
  string out_fname = MakeTempFile("likelihoods_out.vcf");
  auto writer = MakeDogVcfWriter(out_fname, false);
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

bool IsGzipped(StringPiece input) {
  const char gzip_magic[2] = {'\x1f', '\x8b'};
  return (input.size() >= 2 &&
          input[0] == gzip_magic[0] &&
          input[1] == gzip_magic[1]);
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


}  // namespace core
}  // namespace genomics
}  // namespace learning
