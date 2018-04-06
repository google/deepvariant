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

// UnitTests for allelecounter.{h,cc}.
#include "deepvariant/allelecounter.h"
#include <numeric>

#include "deepvariant/utils.h"
#include "third_party/nucleus/io/reference_fai.h"
#include "third_party/nucleus/io/reference_test.h"
#include "third_party/nucleus/protos/position.pb.h"
#include "third_party/nucleus/testing/test_utils.h"
#include "third_party/nucleus/util/utils.h"
#include "tensorflow/core/lib/strings/strcat.h"
#include "tensorflow/core/platform/logging.h"

#include <gmock/gmock-generated-matchers.h>
#include <gmock/gmock-matchers.h>
#include <gmock/gmock-more-matchers.h>

#include "tensorflow/core/platform/test.h"
#include "third_party/nucleus/testing/protocol-buffer-matchers.h"

namespace learning {
namespace genomics {
namespace deepvariant {

using nucleus::EqualsProto;
using nucleus::GenomeReference;
using nucleus::MakePosition;
using nucleus::MakeRange;
using nucleus::TestFastaPath;
using nucleus::genomics::v1::LinearAlignment;
using nucleus::genomics::v1::Read;
using tensorflow::strings::StrCat;
using ::testing::Contains;
using ::testing::Eq;
using ::testing::IsEmpty;
using ::testing::SizeIs;
using ::testing::UnorderedPointwise;

class AlleleCounterTest : public ::testing::Test {
 protected:
  typedef std::vector<Allele> CountLiteral;

  static const int start_;
  static const int end_;
  static const char seq_[];
  static const char chr_[];
  std::unique_ptr<const GenomeReference> ref_;

  AlleleCounterTest() {
    const string& test_fasta_path = TestFastaPath();
    ref_ = std::move(nucleus::GenomeReferenceFai::FromFile(
                         test_fasta_path, StrCat(test_fasta_path, ".fai"))
                         .ValueOrDie());
    read_ = MakeRead("chr1", 1, "TCCGTxx", {"5M"});
    options_.mutable_read_requirements()->set_min_base_quality(21);
  }

  // Creates a new AlleleCount on kChr from kStart to kEnd.
  std::unique_ptr<AlleleCounter> MakeCounter() {
    return MakeCounter(chr_, start_, end_);
  }

  const int min_base_quality() {
    return options_.read_requirements().min_base_quality();
  }

  // Creates a new AlleleCount on specified chr from start to end.
  std::unique_ptr<AlleleCounter> MakeCounter(const string& chr,
                                             const int64 start,
                                             const int64 end) {
    ::nucleus::genomics::v1::Range range = MakeRange(chr, start, end);
    // redacted
    // tensorflow/compiler/xla/ptr_util.h.
    return std::unique_ptr<AlleleCounter>(
        new AlleleCounter(ref_.get(), range, options_));
  }

  // Add reads to allele_count and check that the resulting AlleleCounts are
  // those expected. The expected values are vectors of Alleles, which are
  // compared to those observed in the corresponding position in each
  // AlleleCount in allele_counter. The comparison of the alleles at position i
  // and the expected alleles is done in an order-independent way.
  void AddAndCheckReads(const std::vector<Read>& reads,
                        const std::vector<CountLiteral>& expected,
                        AlleleCounter* allele_counter) {
    ASSERT_THAT(expected.size(), Eq(allele_counter->IntervalLength()));

    // Add our reads to allele_counter.
    allele_counter->Add(reads);

    // Get a vector of our read names for further testing.
    std::vector<string> read_names;
    for (const auto& read : reads) {
      read_names.push_back(allele_counter->ReadKey(read));
    }

    // The number of reads we added is NCountedReads().
    EXPECT_THAT(reads.size(), Eq(allele_counter->NCountedReads()));

    // test that all AlleleCount objects are initialized properly.
    for (int i = 0; i < allele_counter->IntervalLength(); ++i) {
      const AlleleCount& allele_count = allele_counter->Counts()[i];
      const std::vector<Allele> alleles_sum = SumAlleleCounts(allele_count);

      // All read alleles should have count of 1 in the allele_count.
      std::vector<string> read_allele_read_names;
      for (const auto& read_name_allele : allele_count.read_alleles()) {
        EXPECT_THAT(read_names, Contains(read_name_allele.first));
        const Allele& allele = read_name_allele.second;
        EXPECT_THAT(allele.bases(), Not(IsEmpty()));
        EXPECT_THAT(allele.count(), Eq(1));
      }

      // Check that the sum per allele from SumAlleleCounts is correct.
      EXPECT_THAT(alleles_sum, UnorderedPointwise(EqualsProto(), expected[i]));

      // Check that the sum of alleles-specific counts in alleles_sum is equal
      // to the total number of reads in our read_alleles field.
      const int act_total =
          std::accumulate(expected[i].cbegin(), expected[i].cend(), 0,
                          [](const int total, const Allele& allele) {
                            return total + allele.count();
                          });
      EXPECT_THAT(act_total, Eq(TotalAlleleCounts(allele_count)));
    }
  }

  void AddNReads(const int pos, const int n, const string& base,
                 AlleleCounter* counter) {
    for (int i = 0; i < n; ++i) {
      counter->Add(MakeRead("chr1", pos, base, {"1M"}));
    }
  }

  // Same as full AddAndCheckReads() but uses standard AlleleCounter produced by
  // MakeCounter().
  void AddAndCheckReads(const std::vector<Read>& reads,
                        const std::vector<CountLiteral>& expected) {
    AddAndCheckReads(reads, expected, MakeCounter().get());
  }

  // Same as full AddAndCheckReads() but uses standard AlleleCounter produced by
  // MakeCounter() and accepts a single read.
  void AddAndCheckReads(const Read& read,
                        const std::vector<CountLiteral>& expected) {
    AddAndCheckReads(std::vector<Read>{read}, expected);
  }

  // Same as full AddAndCheckReads() but accepts a single read.
  void AddAndCheckReads(const Read& read,
                        const std::vector<CountLiteral>& expected,
                        AlleleCounter* allele_counter) {
    AddAndCheckReads(std::vector<Read>{read}, expected, allele_counter);
  }

  // Creates a test Read with a unique read name.
  Read MakeRead(const string& chr, const int start, const string& bases,
                const std::vector<string>& cigar_elements) {
    Read read = nucleus::MakeRead(chr, start, bases, cigar_elements);
    // Each read gets a unique name.
    read.set_fragment_name(StrCat("read_", ++read_name_counter_));
    return read;
  }

  int read_name_counter_ = 0;
  AlleleCounterOptions options_;
  Read read_;
};

const char AlleleCounterTest::seq_[] = "TCCGT";
const char AlleleCounterTest::chr_[] = "chr1";
const int AlleleCounterTest::start_ = 10;
const int AlleleCounterTest::end_ = 15;

TEST_F(AlleleCounterTest, TestCreate) {
  // All of these tests are designed to work with 5 bp wide interval.
  auto allele_counter = MakeCounter();
  ASSERT_THAT(allele_counter->IntervalLength(), 5);

  // Test simple properties of the counter itself. EXPECT_THAT chr_
  // (const char[]) vs. reference_name (string) doesn't compile.
  EXPECT_EQ(chr_, allele_counter->Interval().reference_name());
  EXPECT_THAT(start_, allele_counter->Interval().start());
  EXPECT_THAT(end_, allele_counter->Interval().end());
  EXPECT_THAT(end_ - start_, allele_counter->IntervalLength());
  EXPECT_THAT(0, allele_counter->NCountedReads());

  // Test that all AlleleCount objects are initialized properly.
  EXPECT_THAT(allele_counter->IntervalLength(),
              allele_counter->Counts().size());
  for (int i = 0; i < allele_counter->IntervalLength(); ++i) {
    const auto& count = allele_counter->Counts()[i];
    // EXPECT_THAT chr_ (const char[]) vs. reference_name (string) doesn't
    // compile.
    EXPECT_EQ(chr_, count.position().reference_name());
    EXPECT_THAT(false, count.position().reverse_strand());
    EXPECT_THAT(allele_counter->Interval().start() + i,
                count.position().position());
    EXPECT_THAT(string{seq_}.substr(i, 1), count.ref_base());
    EXPECT_THAT(true, Eq(count.read_alleles().empty()));
  }
}

TEST_F(AlleleCounterTest, TestAddSimpleRead) {
  for (const auto& op : {"M", "X", "="}) {
    AddAndCheckReads(MakeRead(chr_, start_, "TCCGT", {StrCat(5, op)}),
                     {
                         {MakeAllele("T", AlleleType::REFERENCE, 1)},
                         {MakeAllele("C", AlleleType::REFERENCE, 1)},
                         {MakeAllele("C", AlleleType::REFERENCE, 1)},
                         {MakeAllele("G", AlleleType::REFERENCE, 1)},
                         {MakeAllele("T", AlleleType::REFERENCE, 1)},
                     });
  }
}

TEST_F(AlleleCounterTest, TestReadSpanningBeyondInterval) {
  AddAndCheckReads(MakeRead(chr_, start_ - 2, "AATCCGTAA", {"9M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestAddRead) {
  auto allele_counter = MakeCounter();
  const string seq = seq_;
  for (int start = 0; start < static_cast<int>(seq.length()); ++start) {
    for (int end = seq.length(); end > start; --end) {
      std::vector<CountLiteral> expected;
      for (int i = 0; i < allele_counter->IntervalLength(); ++i) {
        expected.push_back(CountLiteral());
        if (i >= start && i < end) {
          expected[i].push_back(
              MakeAllele(seq.substr(i, 1), AlleleType::REFERENCE, 1));
        }
      }

      const int n_bp = end - start;
      const string read_seq = seq.substr(start, n_bp);
      const string read_cigar = StrCat(n_bp, "M");

      AddAndCheckReads(MakeRead(chr_, start_ + start, read_seq, {read_cigar}),
                       expected);
    }
  }
}

TEST_F(AlleleCounterTest, TestAddSubstitutionRead) {
  for (const size_t subi : {0, 1, 2, 3, 4}) {
    string bases = seq_;
    bases[subi] = 'A';
    std::vector<std::vector<Allele>> expected;
    for (size_t i = 0; i < bases.length(); ++i) {
      auto type = subi == i ? AlleleType::SUBSTITUTION : AlleleType::REFERENCE;
      auto allele = MakeAllele(bases.substr(i, 1), type, 1);
      expected.push_back({allele});
    }
    AddAndCheckReads(MakeRead(chr_, start_, bases, {"5M"}), expected);
  }
}

TEST_F(AlleleCounterTest, TestSimpleInsertion1) {
  AddAndCheckReads(MakeRead(chr_, start_, "TCAAACGT", {"2M", "3I", "3M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("CAAA", AlleleType::INSERTION, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSimpleInsertion2) {
  AddAndCheckReads(MakeRead(chr_, start_, "TAAACCGT", {"1M", "3I", "4M"}),
                   {
                       {MakeAllele("TAAA", AlleleType::INSERTION, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSimpleInsertion3) {
  AddAndCheckReads(MakeRead(chr_, start_, "TCCGTAAA", {"5M", "3I"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("TAAA", AlleleType::INSERTION, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestDiffInsertionSizes) {
  for (int size = 1; size < 10; ++size) {
    const string bases(size, 'A');
    AddAndCheckReads(
        MakeRead(chr_, start_, StrCat("TC", bases, {"CGT"}),
                 {"2M", StrCat(size, "I"), "3M"}),
        {
            {MakeAllele("T", AlleleType::REFERENCE, 1)},
            {MakeAllele(StrCat("C", bases), AlleleType::INSERTION, 1)},
            {MakeAllele("C", AlleleType::REFERENCE, 1)},
            {MakeAllele("G", AlleleType::REFERENCE, 1)},
            {MakeAllele("T", AlleleType::REFERENCE, 1)},
        });
  }
}

TEST_F(AlleleCounterTest, TestStartInsertion) {
  // Starting insertion at the start of our interval gets dropped.
  AddAndCheckReads(MakeRead(chr_, start_, "AAATCCGT", {"3I", "5M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });

  // Starting an insertion is recorded when the read doesn't start at the start
  // of the interval.
  AddAndCheckReads(MakeRead(chr_, start_ + 1, "AAACCGT", {"3I", "4M"}),
                   {
                       {MakeAllele("TAAA", AlleleType::INSERTION, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSimpleDeletion1) {
  AddAndCheckReads(MakeRead(chr_, start_, "TCGT", {"2M", "1D", "2M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("CC", AlleleType::DELETION, 1)},
                       {},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSimpleDeletion2) {
  AddAndCheckReads(MakeRead(chr_, start_, "TCGT", {"1M", "1D", "3M"}),
                   {
                       {MakeAllele("TC", AlleleType::DELETION, 1)},
                       {},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSimpleDeletion3) {
  AddAndCheckReads(MakeRead(chr_, start_, "TCCT", {"3M", "1D", "1M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("CG", AlleleType::DELETION, 1)},
                       {},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestDeletionSize2) {
  AddAndCheckReads(MakeRead(chr_, start_, "TGT", {"1M", "2D", "2M"}),
                   {
                       {MakeAllele("TCC", AlleleType::DELETION, 1)},
                       {},
                       {},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestDeletionSize3) {
  AddAndCheckReads(MakeRead(chr_, start_, "TT", {"1M", "3D", "1M"}),
                   {
                       {MakeAllele("TCCG", AlleleType::DELETION, 1)},
                       {},
                       {},
                       {},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestDeletionSize4) {
  AddAndCheckReads(
      MakeRead(chr_, start_, "T", {"1M", "4D"}),
      {
          {MakeAllele("TCCGT", AlleleType::DELETION, 1)}, {}, {}, {}, {},
      });
}

TEST_F(AlleleCounterTest, TestStartingDeletions) {
  // A read starting with a deletion causes us to just lose the coverage over
  // the deleted base and since the deletion conceptually belongs at previous
  // base, which is right before our interval, we lose the event as well.
  AddAndCheckReads(MakeRead(chr_, start_, "CCGT", {"1D", "4M"}),
                   {
                       {},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });

  // Test that we would have recorded the event had the read started earlier.
  AddAndCheckReads(MakeRead(chr_, start_ + 1, "CGT", {"1D", "3M"}),
                   {
                       {MakeAllele("TC", AlleleType::DELETION, 1)},
                       {},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestEndingDeletions) {
  // It's no problem to have a deletion go up to the end of the interval.
  AddAndCheckReads(MakeRead(chr_, start_, "TCCG", {"4M", "1D"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("GT", AlleleType::DELETION, 1)},
                       {},
                   });

  // We can have a deletion span off the interval, and it's handled properly.
  // Here our deletion spans 2 bp beyond the end of our interval.
  AddAndCheckReads(MakeRead(chr_, start_, "TCCG", {"4M", "3D"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("GTGA", AlleleType::DELETION, 1)},
                       {},
                   });
}

TEST_F(AlleleCounterTest, TestMultipleReads) {
  // Tests that we can add up multiple reads with different starts, cigars, and
  // ends.
  AddAndCheckReads(
      {
          MakeRead(chr_, start_, "TCCGT", {"5M"}),
          MakeRead(chr_, start_, "TCGT", {"2M", "1D", "2M"}),
          MakeRead(chr_, start_ + 2, "CGT", {"3M"}),
          MakeRead(chr_, start_, "TCCAGT", {"3M", "1I", "2M"}),
          MakeRead(chr_, start_ + 2, "CG", {"2M"}),
      },
      {
          {MakeAllele("T", AlleleType::REFERENCE, 3)},
          {MakeAllele("C", AlleleType::REFERENCE, 2),
           MakeAllele("CC", AlleleType::DELETION, 1)},
          {MakeAllele("C", AlleleType::REFERENCE, 3),
           MakeAllele("CA", AlleleType::INSERTION, 1)},
          {MakeAllele("G", AlleleType::REFERENCE, 5)},
          {MakeAllele("T", AlleleType::REFERENCE, 4)},
      });
}

TEST_F(AlleleCounterTest, TestSoftClips1) {
  AddAndCheckReads(MakeRead(chr_, start_ + 2, "AACGT", {"2S", "3M"}),
                   {
                       {},
                       {MakeAllele("CAA", AlleleType::SOFT_CLIP, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSoftClips2) {
  AddAndCheckReads(MakeRead(chr_, start_ + 1, "ACCGT", {"1S", "4M"}),
                   {
                       {MakeAllele("TA", AlleleType::SOFT_CLIP, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSoftClips3) {
  // Soft clip at the start of interval is dropped
  AddAndCheckReads(MakeRead(chr_, start_, "AATCCGT", {"2S", "5M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestSoftClips4) {
  AddAndCheckReads(MakeRead(chr_, start_, "TCCGTAA", {"5M", "2S"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("TAA", AlleleType::SOFT_CLIP, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestInsertionAtChrStart) {
  const int64 chr_start = 0;
  for (const auto op : {"2S", "2I"}) {
    AddAndCheckReads(MakeRead(chr_, chr_start, "AAAC", {op, "2M"}),
                     {
                         {MakeAllele("A", AlleleType::REFERENCE, 1)},
                         {MakeAllele("C", AlleleType::REFERENCE, 1)},
                     },
                     MakeCounter(chr_, chr_start, 2).get());
  }
}

TEST_F(AlleleCounterTest, TestAtChrEnd1) {
  // Our test are built to operate when start is 2 before the last base in the
  // genome on kChr. Load the reference and set chr_start and chr_end
  // appropriately.
  const int64 chr_end = ref_->Contig(chr_).ValueOrDie()->n_bases();
  const int64 chr_start = chr_end - 2;

  for (const auto op : {"2S", "2I"}) {
    auto type =
        op == string("2I") ? AlleleType::INSERTION : AlleleType::SOFT_CLIP;
    AddAndCheckReads(MakeRead(chr_, chr_start, "GAAA", {"2M", op}),
                     {
                         {MakeAllele("G", AlleleType::REFERENCE, 1)},
                         {MakeAllele("AAA", type, 1)},
                     },
                     MakeCounter(chr_, chr_start, chr_end).get());
  }

  // We can have a read ending in a deletion going off the chromosome and it's
  // still ok.  Will generate a warning but not crash the problem.
  AddAndCheckReads(MakeRead(chr_, chr_start, "GA", {"2M", "2D"}),
                   {
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("A", AlleleType::REFERENCE, 1)},
                   },
                   MakeCounter(chr_, chr_start, chr_end).get());

  // Here's a read with matching bases going off the chromosome, which can
  // happen with aligners that don't clip the reads down to the end of the
  // chromosome. Make sure we don't blow up.
  AddAndCheckReads(MakeRead(chr_, chr_start, "GAAAAAAA", {"8M"}),
                   {
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("A", AlleleType::REFERENCE, 1)},
                   },
                   MakeCounter(chr_, chr_start, chr_end).get());
}

TEST_F(AlleleCounterTest, TestDeletionAtChrStart) {
  const int64 chr_start = 0;
  // Deletion at the start of the chrom.
  AddAndCheckReads(MakeRead(chr_, chr_start, "CA", {"2D", "2M"}),
                   {
                       {},
                       {},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("A", AlleleType::REFERENCE, 1)},
                   },
                   MakeCounter(chr_, chr_start, 4).get());
}

TEST_F(AlleleCounterTest, TestMinBaseQualSNP) {
  for (const int bad_pos : {0, 1, 2, 3, 4}) {
    auto read = MakeRead(chr_, start_, "TCCGT", {"5M"});
    std::vector<std::vector<Allele>> expected = {
        {MakeAllele("T", AlleleType::REFERENCE, 1)},
        {MakeAllele("C", AlleleType::REFERENCE, 1)},
        {MakeAllele("C", AlleleType::REFERENCE, 1)},
        {MakeAllele("G", AlleleType::REFERENCE, 1)},
        {MakeAllele("T", AlleleType::REFERENCE, 1)},
    };
    read.set_aligned_quality(bad_pos, min_base_quality() - 1);
    expected[bad_pos].clear();
    AddAndCheckReads(read, expected);
  }
}

TEST_F(AlleleCounterTest, TestMinBaseQualInsertion) {
  // A bad base in the insertion stops us from adding that allele but it
  // preserves our good base before the insertion.
  std::vector<std::vector<Allele>> expected = {
      {MakeAllele("T", AlleleType::REFERENCE, 1)},
      {MakeAllele("C", AlleleType::REFERENCE, 1)},
      {},
      {},
      {},
  };

  for (const int bad_pos : {1, 2, 3}) {
    auto read = MakeRead(chr_, start_, "TAAAC", {"1M", "3I", "1M"});
    read.set_aligned_quality(bad_pos, min_base_quality() - 1);
    AddAndCheckReads(read, expected);
  }
}

TEST_F(AlleleCounterTest, TestMinBaseQualIndelBadInitialBase) {
  // A bad base in the insertion stops us from adding that allele but it
  // preserves our good base before the insertion.
  auto read = MakeRead(chr_, start_, "TCAAACGT", {"2M", "3I", "3M"});
  // Good case -- no bad bases.
  AddAndCheckReads(read, {
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                             {MakeAllele("CAAA", AlleleType::INSERTION, 1)},
                             {MakeAllele("C", AlleleType::REFERENCE, 1)},
                             {MakeAllele("G", AlleleType::REFERENCE, 1)},
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                         });

  // The insertion has a bad base, but the C anchor is good so instead of the
  // Cxxx allele we see C as a match.
  read.set_aligned_quality(3, min_base_quality() - 1);
  AddAndCheckReads(read, {
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                             {MakeAllele("C", AlleleType::REFERENCE, 1)},
                             {MakeAllele("C", AlleleType::REFERENCE, 1)},
                             {MakeAllele("G", AlleleType::REFERENCE, 1)},
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                         });

  // This case has both the insertion and the anchor base being bad, so there's
  // no count anywhere.
  read.set_aligned_quality(1, min_base_quality() - 1);
  AddAndCheckReads(read, {
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                             {},
                             {MakeAllele("C", AlleleType::REFERENCE, 1)},
                             {MakeAllele("G", AlleleType::REFERENCE, 1)},
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                         });

  // Interestingly, if the previous base is bad but the insertion is good we
  // get the same alleles as in the all good bases case since representing the
  // insertion requires us to encode the anchor base.
  read.set_aligned_quality(3, min_base_quality() + 1);
  AddAndCheckReads(read, {
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                             {MakeAllele("CAAA", AlleleType::INSERTION, 1)},
                             {MakeAllele("C", AlleleType::REFERENCE, 1)},
                             {MakeAllele("G", AlleleType::REFERENCE, 1)},
                             {MakeAllele("T", AlleleType::REFERENCE, 1)},
                         });
}

TEST_F(AlleleCounterTest, TestSNPIndel) {
  // Test that we get the right counts when a read contains a substitution
  // immediately followed by an indel.
  AddAndCheckReads(MakeRead(chr_, start_, "TAAAACGT", {"2M", "3I", "3M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("AAAA", AlleleType::INSERTION, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestPairedReads) {
  // Tests that we count properly reads that have the same fragment_name but
  // have different read numbers (first and second of pair, for example).
  Read read1 = MakeRead(chr_, start_, "TCCAT", {"5M"});
  Read read2 = MakeRead(chr_, start_, "TCAAT", {"5M"});
  read1.set_fragment_name("fragment");
  read1.set_read_number(0);
  read2.set_fragment_name("fragment");
  read2.set_read_number(1);
  AddAndCheckReads({read1, read2},
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 2)},
                       {MakeAllele("C", AlleleType::REFERENCE, 2)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1),
                        MakeAllele("A", AlleleType::SUBSTITUTION, 1)},
                       {MakeAllele("A", AlleleType::SUBSTITUTION, 2)},
                       {MakeAllele("T", AlleleType::REFERENCE, 2)},
                   });
}

TEST_F(AlleleCounterTest, TestCanonicalBases) {
  // We ignore an N base in our read that matches.
  AddAndCheckReads(MakeRead(chr_, start_, "TCNGT", {"5M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });

  // We ignore an N base that anchors an indel, skipping the event entirely.
  AddAndCheckReads(MakeRead(chr_, start_, "TNGT", {"2M", "1D", "2M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {},
                       {},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });

  AddAndCheckReads(MakeRead(chr_, start_, "TCNAGT", {"3M", "1I", "2M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });

  // We ignore an N base that occurs in insertion, but the anchor base is good
  // so it now appears as a reference match in the counts.
  AddAndCheckReads(MakeRead(chr_, start_, "TCCNGT", {"3M", "1I", "2M"}),
                   {
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                   });
}

TEST_F(AlleleCounterTest, TestCanonicalBasesReference) {
  // Test that we handle N reference bases correctly. This involves checking
  // that we can count reads that have a canonical base at a reference N. It
  // also checks that we handle deletion alleles that span reference N bases.
  // Note that we have to use "chr2" which contains an N in the reference
  // "CGCTNCG..." is the start of chr2.
  const string chr = "chr2";
  const int start = 2;
  const int end = start + 5;

  // Our read has an A over the N base.
  AddAndCheckReads(MakeRead(chr, start, "CTACG", {"5M"}),
                   {
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {MakeAllele("A", AlleleType::SUBSTITUTION, 1)},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                   },
                   MakeCounter(chr, start, end).get());

  // We delete away this N, and therefore the deletion allele has an N base and
  // so isn't counted.
  AddAndCheckReads(MakeRead(chr, start, "CTCG", {"2M", "1D", "2M"}),
                   {
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("T", AlleleType::REFERENCE, 1)},
                       {},
                       {MakeAllele("C", AlleleType::REFERENCE, 1)},
                       {MakeAllele("G", AlleleType::REFERENCE, 1)},
                   },
                   MakeCounter(chr, start, end).get());
}

TEST_F(AlleleCounterTest, TestCountSummaries) {
  std::unique_ptr<AlleleCounter> counter = MakeCounter("chr1", 1, 4);
  AddNReads(1, 1, "C", counter.get());
  AddNReads(1, 2, "T", counter.get());
  AddNReads(2, 3, "C", counter.get());
  AddNReads(2, 4, "T", counter.get());
  AddNReads(3, 5, "A", counter.get());
  AddNReads(3, 6, "T", counter.get());

  std::vector<AlleleCountSummary> summaries = counter->SummaryCounts();
  EXPECT_THAT(summaries, SizeIs(3));

  for (int i = 0; i < 3; ++i) {
    // Checks that the reference names and positions are correct.
    EXPECT_EQ(summaries[i].reference_name(), "chr1");
    EXPECT_EQ(summaries[i].position(), i + 1);
    EXPECT_EQ(summaries[i].ref_nonconfident_read_count(), 0);
  }

  // Itemized checks that the summary counts are correct for each position.
  EXPECT_EQ(summaries[0].ref_base(), "C");
  EXPECT_EQ(summaries[0].ref_supporting_read_count(), 1);
  EXPECT_EQ(summaries[0].total_read_count(), 3);
  EXPECT_EQ(summaries[1].ref_base(), "C");
  EXPECT_EQ(summaries[1].ref_supporting_read_count(), 3);
  EXPECT_EQ(summaries[1].total_read_count(), 7);
  EXPECT_EQ(summaries[2].ref_base(), "A");
  EXPECT_EQ(summaries[2].ref_supporting_read_count(), 5);
  EXPECT_EQ(summaries[2].total_read_count(), 11);
}

}  // namespace deepvariant
}  // namespace genomics
}  // namespace learning
