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

#include "third_party/nucleus/io/sam_reader.h"

#include <string>

#include "third_party/nucleus/protos/index.pb.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "third_party/nucleus/vendor/status_matchers.h"

#include "tensorflow/core/lib/core/stringpiece.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"

namespace nucleus {

using nucleus::genomics::v1::Range;
using nucleus::genomics::v1::Read;
using nucleus::genomics::v1::SamHeader;
using nucleus::genomics::v1::SamReaderOptions;
using nucleus::proto::IgnoringFieldPaths;
using nucleus::proto::Partially;
using std::vector;
using tensorflow::StringPiece;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::Pointwise;
using ::testing::SizeIs;

// Constants for all filenames used in this test file.
constexpr char kSamTestFilename[] = "test.sam";
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
  EXPECT_EQ(header.comments(0), "@CO\tA single line comment.");
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
  std::unique_ptr<SamReader> reader = std::move(
      SamReader::FromFile(GetTestData("headerless.sam"), SamReaderOptions())
          .ValueOrDie());
  auto iterator = reader->Iterate().ValueOrDie();
  Read r;
  StatusOr<bool> status =  iterator->Next(&r);
  ASSERT_EQ(status.ok(), false);
}

class SamReaderQueryTest : public ::testing::Test {
 protected:
  void SetUp() override {
    options_.set_index_mode(
        nucleus::genomics::v1::IndexHandlingMode::INDEX_BASED_ON_FILENAME);
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

}  // namespace nucleus
