// Unit tests for the sex-chromosome haploid-calling helpers shared by the
// make_examples gVCF path and postprocess (deepvariant/native/haploid_regions.h).

#include "deepvariant/native/haploid_regions.h"

#include <set>
#include <string>
#include <vector>

#include "gtest/gtest.h"

namespace deepvariant {
namespace {

TEST(HaploidRegions, ParseHaploidContigsSplitsCommaAndWhitespace) {
  EXPECT_TRUE(ParseHaploidContigs("").empty());
  EXPECT_EQ(ParseHaploidContigs("chrX,chrY"),
            (std::set<std::string>{"chrX", "chrY"}));
  // Upstream splits on ',' then on whitespace, so both separators work.
  EXPECT_EQ(ParseHaploidContigs("X Y"), (std::set<std::string>{"X", "Y"}));
  EXPECT_EQ(ParseHaploidContigs(" chrX , chrY "),
            (std::set<std::string>{"chrX", "chrY"}));
}

TEST(HaploidRegions, ParRegionsOverlapIsHalfOpen) {
  ParRegions par;
  par.by_contig["chrX"] = {{10, 20}, {100, 200}};
  EXPECT_TRUE(par.Overlaps("chrX", 10, 11));   // at start
  EXPECT_TRUE(par.Overlaps("chrX", 19, 25));   // straddles end
  EXPECT_FALSE(par.Overlaps("chrX", 20, 25));  // [20,25) abuts [10,20) — no
  EXPECT_FALSE(par.Overlaps("chrX", 0, 10));   // [0,10) abuts [10,20) — no
  EXPECT_FALSE(par.Overlaps("chrY", 10, 11));  // wrong contig
}

TEST(HaploidRegions, IsHaploidPositionRespectsContigSetAndPar) {
  const std::set<std::string> haploid = {"chrX", "chrY"};
  ParRegions par;
  par.by_contig["chrX"] = {{1000, 2000}};

  // On a haploid contig, outside PAR → haploid.
  EXPECT_TRUE(IsHaploidPosition("chrX", 500, 501, haploid, par));
  // On a haploid contig, inside PAR → diploid.
  EXPECT_FALSE(IsHaploidPosition("chrX", 1500, 1501, haploid, par));
  // Autosome → never haploid.
  EXPECT_FALSE(IsHaploidPosition("chr20", 500, 501, haploid, par));
  // Empty haploid set → never haploid.
  EXPECT_FALSE(IsHaploidPosition("chrX", 500, 501, {}, par));
}

TEST(HaploidRegions, CorrectNonautosomeBiallelicZerosHet) {
  // PL order for 1 alt: [0/0, 0/1, 1/1]; 0/1 (index 1) is het.
  std::vector<double> like = {0.2, 0.5, 0.3};
  CorrectNonautosomeProbabilities(&like, /*n_alts=*/1);
  ASSERT_EQ(like.size(), 3u);
  EXPECT_DOUBLE_EQ(like[1], 0.0);
  EXPECT_DOUBLE_EQ(like[0], 0.4);  // 0.2 / (0.2 + 0.3)
  EXPECT_DOUBLE_EQ(like[2], 0.6);  // 0.3 / (0.2 + 0.3)
}

TEST(HaploidRegions, CorrectNonautosomeTriallelicZerosAllHets) {
  // PL order for 2 alts: 0/0,0/1,1/1,0/2,1/2,2/2.
  // Het genotypes: 0/1(1), 0/2(3), 1/2(4). Hom: 0/0(0), 1/1(2), 2/2(5).
  std::vector<double> like = {0.1, 0.2, 0.1, 0.2, 0.3, 0.1};
  CorrectNonautosomeProbabilities(&like, /*n_alts=*/2);
  EXPECT_DOUBLE_EQ(like[1], 0.0);
  EXPECT_DOUBLE_EQ(like[3], 0.0);
  EXPECT_DOUBLE_EQ(like[4], 0.0);
  // Surviving hom mass {0.1, 0.1, 0.1} renormalizes to thirds.
  EXPECT_DOUBLE_EQ(like[0], 1.0 / 3.0);
  EXPECT_DOUBLE_EQ(like[2], 1.0 / 3.0);
  EXPECT_DOUBLE_EQ(like[5], 1.0 / 3.0);
}

TEST(HaploidRegions, CorrectNonautosomeAllZeroIsSafe) {
  // Degenerate input (all het, everything zeroed) must not divide by zero.
  std::vector<double> like = {0.0, 1.0, 0.0};
  CorrectNonautosomeProbabilities(&like, /*n_alts=*/1);
  for (double v : like) EXPECT_DOUBLE_EQ(v, 0.0);
}

TEST(HaploidRegions, CorrectNonautosomeSizeMismatchIsNoOp) {
  // A vector whose length doesn't match (n_alts+1)(n_alts+2)/2 is a contract
  // violation; the function must leave it untouched rather than write OOB.
  std::vector<double> like = {0.2, 0.5};  // n_alts=1 expects 3 entries.
  CorrectNonautosomeProbabilities(&like, /*n_alts=*/1);
  ASSERT_EQ(like.size(), 2u);
  EXPECT_DOUBLE_EQ(like[0], 0.2);
  EXPECT_DOUBLE_EQ(like[1], 0.5);
}

}  // namespace
}  // namespace deepvariant
