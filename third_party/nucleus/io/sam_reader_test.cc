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

#include "third_party/nucleus/io/sam_reader.h"

#include <string>
#include <utility>
#include <vector>

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/io/sam_writer.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/core/status_matchers.h"
#include "tensorflow/core/lib/core/status.h"

namespace nucleus {

using nucleus::genomics::v1::LinearAlignment;
using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;
using nucleus::genomics::v1::ReadRequirements;
using nucleus::genomics::v1::SamHeader;
using nucleus::genomics::v1::SamReaderOptions;
using nucleus::proto::IgnoringFieldPaths;
using nucleus::proto::Partially;
using std::vector;
using ::testing::IsEmpty;
using ::testing::Key;
using ::testing::Pointwise;
using ::testing::SizeIs;
using ::testing::UnorderedElementsAre;

// Constants for all filenames used in this test file.
constexpr char kSamTestFilename[] = "test.sam";
constexpr char kSamOqTestFilename[] = "test_oq.sam";
constexpr char kBamTestFilename[] = "test.bam";
constexpr char kSamGoldStandardFilename[] = "test.sam.golden.tfrecord";

// Checks if the result of converting a test sam file matches the gold standard
// record io files for SamHeader and Read. This is representative of a real life
// file conversion.
TEST(ReadBamFile, MatchesGolden) {
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamTestFilename), SamReaderOptions())
          .ValueOrDie());
  vector<Read> golden =
      ReadProtosFromTFRecord<Read>(GetTestData(kSamGoldStandardFilename));
  EXPECT_THAT(
      as_vector(reader->Iterate()),
      Pointwise(IgnoringFieldPaths({"info"}, EqualsProto()), golden));
}

TEST(SamReaderTest, TestIteration) {
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamTestFilename), SamReaderOptions())
          .ValueOrDie());
  EXPECT_THAT(as_vector(reader->Iterate()), SizeIs(6));
}

// test_oq.sam is used for this test where original scores all set to 'C'
// The test checks that if use_original_base_quality_scores is set alignment
// quality scores are taken from OQ tag and all the scores properly calculated.
TEST(SamReaderTest, TestAlignedQualityOQ) {
  SamReaderOptions samReaderOptions;
  samReaderOptions.set_use_original_base_quality_scores(true);
  samReaderOptions.set_aux_field_handling(
      SamReaderOptions::PARSE_ALL_AUX_FIELDS);
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamOqTestFilename), samReaderOptions)
          .ValueOrDie());

  // Test compares aligned_quality field that is read from test_oq.sam with
  // golden set.
  auto reads = as_vector(reader->Iterate());
  std::vector<Read> golden;
  for (const auto& record : reads) {
    Read golden_read;
    auto aq = golden_read.mutable_aligned_quality();
    // Golden set is created by copying aligned_quality from test_oq.sam
    // And then overwriting all scores with 'C'-33
    aq->CopyFrom(record.aligned_quality());

    for (auto& base_quality : *aq) {
      // Score encoding described here
      // https://samtools.github.io/hts-specs/SAMv1.pdf
      base_quality = 'C' - 33;
    }
    golden.push_back(golden_read);
  }

  EXPECT_THAT(reads, Pointwise(Partially(EqualsProto()), golden));
}

// Trying to read quality scores from OQ when OQ tag is not present.
TEST(SamReaderTest, TestAlignedQualityOQWhenTagIsNotPresent) {
  SamReaderOptions samReaderOptions;
  samReaderOptions.set_use_original_base_quality_scores(true);
  samReaderOptions.set_aux_field_handling(
      SamReaderOptions::PARSE_ALL_AUX_FIELDS);
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamTestFilename), samReaderOptions)
          .ValueOrDie());

  // Test compares aligned_quality field that is read from test_oq.sam with
  // golden set.
  auto reads = as_vector(reader->Iterate());
  std::vector<Read> golden;
  for (const auto& record : reads) {
    EXPECT_EQ(record.aligned_quality().size(), 0);
  }
}

// Test that assert is raised if aux_field_handling is not set together with
// use_original_base_quality_scores.
TEST(SamReaderTest, TestFailIfParseAuxFieldsIsNotSetWithUseOriginalOqualities) {
  SamReaderOptions samReaderOptions;
  samReaderOptions.set_use_original_base_quality_scores(true);
  ASSERT_DEATH(
      SamReader::FromFile(GetTestData(kSamTestFilename), samReaderOptions),
      "aux_field_handling must be true if use_original_quality_scores is set "
      "to true");
}

TEST(SamReaderTest, TestEmptyAuxFieldsToKeepReadsEverything) {
  SamReaderOptions samReaderOptions;
  samReaderOptions.set_aux_field_handling(
      SamReaderOptions::PARSE_ALL_AUX_FIELDS);
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamTestFilename), samReaderOptions)
          .ValueOrDie());
  auto reads = as_vector(reader->Iterate());
  EXPECT_THAT(reads[0].info(),
              UnorderedElementsAre(Key("NM"), Key("MD"), Key("AS"), Key("XS"),
                                   Key("RG")));
}

TEST(SamReaderTest, TestSetAuxFieldsToKeep) {
  SamReaderOptions samReaderOptions;
  samReaderOptions.set_aux_field_handling(
      SamReaderOptions::PARSE_ALL_AUX_FIELDS);
  samReaderOptions.add_aux_fields_to_keep("NM");  // This exists in the read.
  samReaderOptions.add_aux_fields_to_keep("FOO");  // This doesn't exist.
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamTestFilename), samReaderOptions)
          .ValueOrDie());
  auto reads = as_vector(reader->Iterate());
  EXPECT_THAT(reads[0].info(), UnorderedElementsAre(Key("NM")));
}

// Test that assert is raised if aux_fields_to_keep doesn't contain OQ when
// use_original_base_quality_scores.
TEST(SamReaderTest, TestFailAuxFieldsToKeepIsNotSetWithUseOriginalOqualities) {
  SamReaderOptions samReaderOptions;
  samReaderOptions.set_use_original_base_quality_scores(true);
  samReaderOptions.set_aux_field_handling(
      SamReaderOptions::PARSE_ALL_AUX_FIELDS);
  // aux_fields_to_keep is not empty, but doesn't contain OQ.
  samReaderOptions.add_aux_fields_to_keep("HP");
  ASSERT_DEATH(
      SamReader::FromFile(GetTestData(kSamTestFilename), samReaderOptions),
      "aux_fields_to_keep must contain OQ or be empty");
}

TEST(SamReaderTest, TestIterationRespectsReadRequirements) {
  SamReaderOptions options;
  options.mutable_read_requirements()->set_keep_unaligned(false);
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamTestFilename), options)
          .ValueOrDie());
  EXPECT_THAT(as_vector(reader->Iterate()), SizeIs(5));
}

TEST(SamReaderTest, TestSamHeaderExtraction) {
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kSamTestFilename), SamReaderOptions())
          .ValueOrDie());
  const SamHeader& header = reader->Header();
  EXPECT_EQ(header.format_version(), "1.3");
  EXPECT_EQ(header.sorting_order(), SamHeader::COORDINATE);
  EXPECT_EQ(header.alignment_grouping(), SamHeader::NONE);
  EXPECT_THAT(header.contigs(), SizeIs(493));
  EXPECT_THAT(header.read_groups(), SizeIs(1));
  const nucleus::genomics::v1::ReadGroup& rg = header.read_groups(0);
  EXPECT_EQ(rg.name(), "'Illumina3D.4");
  EXPECT_EQ(rg.sequencing_center(), "GOOG");
  EXPECT_EQ(rg.description(), "description");
  EXPECT_EQ(rg.date(), "12/10/2012");
  EXPECT_EQ(rg.flow_order(), "ACMG");
  EXPECT_EQ(rg.key_sequence(), "GATTACA");
  EXPECT_THAT(rg.program_ids(), SizeIs(1));
  EXPECT_EQ(rg.program_ids(0), "bwa");
  EXPECT_EQ(rg.predicted_insert_size(), 300);
  EXPECT_EQ(rg.platform(), "illumina");
  EXPECT_EQ(rg.platform_model(), "HiSeq");
  EXPECT_EQ(rg.platform_unit(), "abcde");
  EXPECT_THAT(rg.sample_id(), IsEmpty());
  EXPECT_THAT(header.comments(), SizeIs(1));
  EXPECT_EQ(header.comments(0), "A single line comment.");
}

TEST(SamReaderTest, TestBamSampleExtraction) {
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData(kBamTestFilename), SamReaderOptions())
          .ValueOrDie());
  const SamHeader& header = reader->Header();
  EXPECT_THAT(header.format_version(), IsEmpty());
  EXPECT_EQ(header.sorting_order(), SamHeader::UNKNOWN);
  EXPECT_EQ(header.alignment_grouping(), SamHeader::NONE);
  EXPECT_THAT(header.contigs(), SizeIs(25));
  EXPECT_THAT(header.read_groups(), SizeIs(1));
  EXPECT_EQ(header.read_groups(0).sample_id(), "NA12878");
  EXPECT_THAT(header.comments(), IsEmpty());
}

TEST(SamReaderTest, TestHeaderlessSamIsNotOkay) {
  StatusOr<std::unique_ptr<SamReader>> status = SamReader::FromFile(
      GetTestData("headerless.sam"), SamReaderOptions());
  ASSERT_EQ(status.ok(), false);
}

TEST(SamReaderTest, TestMatePosition) {
  // Write a file with one record that has mapped mate reference in FLAG field,
  // but has * as RNEXT (mate reference name).
  string output_filename(MakeTempFile("sam_reader_test.sam"));
  {
    std::unique_ptr<SamReader> reader = std::move(
        SamReader::FromFile(GetTestData(kBamTestFilename), SamReaderOptions())
            .ValueOrDie());
    std::vector<Read> reads = as_vector(reader->Iterate());
    Read r = reads[0];
    EXPECT_TRUE(r.has_next_mate_position());
    r.mutable_next_mate_position()->set_position(-1);
    r.mutable_next_mate_position()->set_reference_name("*");

    std::unique_ptr<SamWriter> writer = std::move(
        SamWriter::ToFile(output_filename, reader->Header()).ValueOrDie());
    EXPECT_THAT(writer->Write(r), IsOK());
  }

  // Read the output file, and make sure the * in mate reference name is parsed
  // as an unmapped next read.
  {
    std::unique_ptr<SamReader> reader = std::move(
        SamReader::FromFile(output_filename, SamReaderOptions()).ValueOrDie());
    std::vector<Read> reads = as_vector(reader->Iterate());
    ASSERT_EQ(1u, reads.size());
    auto read = reads[0];
    EXPECT_FALSE(read.has_next_mate_position());
  }
  TF_CHECK_OK(tensorflow::Env::Default()->DeleteFile(output_filename));
}

class SamReaderQueryTest : public ::testing::Test {
 protected:
  void SetUp() override {
    indexed_bam_ = GetTestData(kBamTestFilename);
    RecreateReader();
  }

  void RecreateReader() {
    reader_ =
        std::move(SamReader::FromFile(indexed_bam_, options_).ValueOrDie());
  }

  SamReaderOptions options_;
  string indexed_bam_;
  std::unique_ptr<SamReader> reader_;
};

/*
// samtools view \
//   learning/genomics/io/testdata/NA12878_S1.chr20.10_11mb.bam \
//   chr20:10,010,000-10,011,000 | wc -l
*/
TEST_F(SamReaderQueryTest, SimpleQueriesWork) {
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr20", 9999999, 10000000))),
              SizeIs(45));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr20", 9999999, 10000100))),
              SizeIs(106));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr20", 999999, 10000000))),
              SizeIs(45));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr20", 999999, 100000000))),
              SizeIs(106));
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr20", 999999, 2000000))),
              IsEmpty());
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr10", 9999999, 10000000))),
              IsEmpty());
  EXPECT_THAT(as_vector(reader_->Query(MakeRange("chr1", 0, 100000000))),
              IsEmpty());
}


TEST_F(SamReaderQueryTest, ThatRangeIsExactlyCorrect) {
  // Tests that our range parameter gives us exactly the read we expect.
  // In IGV this reads spans chr20:9,999,912-10,000,010
  const string read_proto =
      "fragment_name: 'HSQ1004:134:C0D8DACXX:4:1304:21341:94622'";
  const string chrom = "chr20";
  const int start_inclusive0 = 9999911;
  const int end_exclusive0 = 10000010;

  // To start, our read is present in the full interval query over it.
  EXPECT_THAT(as_vector(reader_->Query(
                  MakeRange(chrom, start_inclusive0, end_exclusive0))),
              Contains(Partially(EqualsProto<Read>(read_proto))));
  // The read occurs when it fully spans the interval.
  EXPECT_THAT(as_vector(reader_->Query(
                  MakeRange(chrom, start_inclusive0 + 1, end_exclusive0 - 1))),
              Contains(Partially(EqualsProto<Read>(read_proto))));
  // The read occurs when the read starts/ends within the interval.
  EXPECT_THAT(
      as_vector(
          reader_->Query(
              MakeRange(chrom, start_inclusive0 + 5, end_exclusive0 + 5))),
      Contains(Partially(EqualsProto<Read>(read_proto))));
  EXPECT_THAT(
      as_vector(
          reader_->Query(
              MakeRange(chrom, start_inclusive0 - 5, end_exclusive0 - 5))),
      Contains(Partially(EqualsProto<Read>(read_proto))));

  // Tests that the read occurs exactly when the interval end overlaps the first
  // base of the read from the left.
  const int left_start = start_inclusive0 - 10;
  EXPECT_THAT(
      // start_inclusive0 + 1 includes the first left base of the read.
      as_vector(reader_->Query(
          MakeRange(chrom, left_start, start_inclusive0 + 1))),
      Contains(Partially(EqualsProto<Read>(read_proto))));
  EXPECT_THAT(
      // start_inclusive0 *does not* include the first left base of the read.
      as_vector(reader_->Query(MakeRange(chrom, left_start, start_inclusive0))),
      Not(Contains(Partially(EqualsProto<Read>(read_proto)))));

  // Tests that the read occurs exactly when the interval start overlaps the
  // last base of the read from the right.
  const int right_end = end_exclusive0 + 10;
  EXPECT_THAT(
      // end_exclusive0 - 1 includes the last base of the read.
      as_vector(
          reader_->Query(MakeRange(chrom, end_exclusive0 - 1, right_end))),
      Contains(Partially(EqualsProto<Read>(read_proto))));
  EXPECT_THAT(
      // end_exclusive0 *does not* include the last base of the read.
      as_vector(reader_->Query(MakeRange(chrom, end_exclusive0, right_end))),
      Not(Contains(Partially(EqualsProto<Read>(read_proto)))));
}

/*
// samtools view kBamTestFilename chr20:10,000,000-10,000,100 \
// | cut -f 5 | sort | uniq -c
//       1 0
//       1 37
//     104 60
*/
TEST_F(SamReaderQueryTest, QueriedRespectsReadRequirements) {
  Range range = MakeRange("chr20", 9999999, 10000100);

  // Without any read requirements we have 106 reads.
  EXPECT_THAT(as_vector(reader_->Query(range)), SizeIs(106));

  // Initializing the read requirements enforces things like proper placement,
  // which cuts out one read.
  options_.mutable_read_requirements();
  RecreateReader();
  EXPECT_THAT(as_vector(reader_->Query(range)), SizeIs(105));

  // There are two reads that don't have a MAPQ60. The 0 MAPQ is flagged as
  // unmapped, and there is another with MAPQ of 37. We test below that we get
  // the right number of reads for each configuration of requirements.
  options_.mutable_read_requirements()->set_keep_unaligned(true);
  RecreateReader();
  EXPECT_THAT(as_vector(reader_->Query(range)), SizeIs(106));

  options_.mutable_read_requirements()->set_keep_unaligned(false);
  RecreateReader();
  EXPECT_THAT(as_vector(reader_->Query(range)), SizeIs(105));

  options_.mutable_read_requirements()->set_min_mapping_quality(38);
  RecreateReader();
  EXPECT_THAT(as_vector(reader_->Query(range)), SizeIs(104));
}

TEST_F(SamReaderQueryTest, ReadAfterClose) {
  ASSERT_THAT(reader_->Close(), IsOK());
  EXPECT_THAT(reader_->Iterate(),
              IsNotOKWithMessage("Cannot Iterate a closed SamReader."));
  EXPECT_THAT(reader_->Query(MakeRange("chr20", 9999999, 10000000)),
              IsNotOKWithMessage("Cannot Query a closed SamReader."));
}

TEST_F(SamReaderQueryTest, NextFailsOnReleasedIterable) {
  Read read;
  std::shared_ptr<SamIterable> it = reader_->Iterate().ValueOrDie();
  ASSERT_THAT(it->Release(), IsOK());
  EXPECT_THAT(it->Next(&read), IsNotOKWithMessage("Reader is not alive"));
}

namespace sam_reader_internal {

class ReadRequirementTest : public ::testing::Test {
 protected:
  Read read_;
  ReadRequirements reqs_;

  ReadRequirementTest() {
    read_.set_fragment_name("read1");
    read_.set_aligned_sequence("ABC");
    read_.set_number_reads(2);
    read_.set_proper_placement(true);
    LinearAlignment& aln = *read_.mutable_alignment();
    aln.set_mapping_quality(90);
    *aln.mutable_position() = MakePosition("chr1", 10);
  }
};

TEST_F(ReadRequirementTest, ThatEmptyReadFailsBecauseOfNoAlignment) {
  Read read = Read();
  EXPECT_FALSE(ReadSatisfiesRequirements(read, reqs_));
  *(read.mutable_alignment()->mutable_position()) = MakePosition("chr1", 1);
  EXPECT_TRUE(ReadSatisfiesRequirements(read, reqs_));
}

TEST_F(ReadRequirementTest, TestBaseReadSatisfiesRequirements) {
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

TEST_F(ReadRequirementTest, TestDuplicateFilter) {
  read_.set_duplicate_fragment(true);
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));
  reqs_.set_keep_duplicates(true);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

TEST_F(ReadRequirementTest, TestVenderFilter) {
  read_.set_failed_vendor_quality_checks(true);
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));
  reqs_.set_keep_failed_vendor_quality_checks(true);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

TEST_F(ReadRequirementTest, TestSecondaryAlignmentFilter) {
  read_.set_secondary_alignment(true);
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));
  reqs_.set_keep_secondary_alignments(true);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

TEST_F(ReadRequirementTest, TestSupplmentaryAlignmentFilter) {
  read_.set_supplementary_alignment(true);
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));
  reqs_.set_keep_supplementary_alignments(true);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

TEST_F(ReadRequirementTest, TestProperPlacement) {
  // We don't use reads that aren't properly placed. Here the read's mate is
  // mapped to chrX but the read is mapped to chr1. This is an improper pair.
  read_.set_proper_placement(false);
  *read_.mutable_next_mate_position() = MakePosition("chrX", 25);
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));
  // Now the read's mate is mapped to chr1 so it is properly placed.
  *read_.mutable_next_mate_position() = MakePosition("chr1", 25);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
  reqs_.set_keep_improperly_placed(true);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

TEST_F(ReadRequirementTest, TestSingleEndedProperPlacement) {
  // Singled ended reads pass.
  read_.set_number_reads(1);
  read_.set_proper_placement(false);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
  reqs_.set_keep_improperly_placed(true);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

TEST_F(ReadRequirementTest, TestMappingQuality) {
  const int min_mapq = 10;

  // There's no minimum set, so even with mapq 0 this read should pass.
  read_.mutable_alignment()->set_mapping_quality(0);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));

  // We also keep reads with the default mapping_quality.
  read_.mutable_alignment()->clear_mapping_quality();
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));

  // Setting the min_mapping_quality now rejects the read.
  reqs_.set_min_mapping_quality(min_mapq);
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));

  // Check that the min_mapping_quality calculation is correct.
  read_.mutable_alignment()->set_mapping_quality(min_mapq - 1);
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));
  read_.mutable_alignment()->set_mapping_quality(min_mapq);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));

  // A read without a alignment but otherwise good will pass even without
  // satisfying our mapping quality as long as keep_unaligned is true.
  read_.clear_alignment();
  EXPECT_FALSE(ReadSatisfiesRequirements(read_, reqs_));
  reqs_.set_keep_unaligned(true);
  EXPECT_TRUE(ReadSatisfiesRequirements(read_, reqs_));
}

}  // namespace sam_reader_internal

}  // namespace nucleus
