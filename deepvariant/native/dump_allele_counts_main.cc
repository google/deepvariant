// Tiny diagnostic: load a BAM + reference, run our AlleleCounter on a
// region, and dump per-position counts (ref + alt alleles) to stdout.
// Used during make_examples parity work to check why a given upstream
// candidate isn't appearing in our pipeline.
//
// Usage:
//   dump_allele_counts <ref.fa> <reads.bam> <chr:start-end>
//
// Output (one line per AlleleCount, only positions with any alt support):
//   <chrom> <pos1> <ref_base> ref=<n> alt:<allele>=<count>[<L?>] ...
//
// pos1 is 1-based for grep-friendliness.

#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>
#include <vector>

#include "absl/strings/numbers.h"
#include "absl/strings/str_split.h"
#include "deepvariant/allelecounter.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/io/sam_reader.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"

using ::learning::genomics::deepvariant::AlleleCount;
using ::learning::genomics::deepvariant::AlleleCounter;
using ::learning::genomics::deepvariant::AlleleCounterOptions;
using ::learning::genomics::deepvariant::AlleleType;

int main(int argc, char** argv) {
  if (argc != 4) {
    std::fprintf(stderr,
                 "usage: %s <ref.fa> <reads.bam> <chr:start-end>\n", argv[0]);
    return 2;
  }
  const std::string ref_path = argv[1];
  const std::string bam_path = argv[2];
  const std::string region_str = argv[3];

  // Parse "chr20:5001580-5001650" into a Range (0-based, half-open).
  nucleus::genomics::v1::Range region;
  {
    std::vector<std::string> parts = absl::StrSplit(region_str, ':');
    if (parts.size() != 2) { std::fprintf(stderr, "bad region\n"); return 2; }
    region.set_reference_name(parts[0]);
    std::vector<std::string> se = absl::StrSplit(parts[1], '-');
    if (se.size() != 2) { std::fprintf(stderr, "bad region\n"); return 2; }
    int64_t s = 0, e = 0;
    if (!absl::SimpleAtoi(se[0], &s) || !absl::SimpleAtoi(se[1], &e) ||
        s < 1 || e < s) {
      std::fprintf(stderr, "bad region\n");  // 1-based, non-empty range
      return 2;
    }
    region.set_start(s - 1);  // 1-based input → 0-based proto.
    region.set_end(e);
  }

  auto ref_or = nucleus::IndexedFastaReader::FromFile(ref_path,
                                                      ref_path + ".fai");
  if (!ref_or.ok()) { std::fprintf(stderr, "ref open failed\n"); return 1; }
  auto ref = std::move(ref_or.ValueOrDie());

  nucleus::genomics::v1::SamReaderOptions sam_opts;
  sam_opts.mutable_read_requirements()->set_min_mapping_quality(10);
  auto sam_or = nucleus::SamReader::FromFile(bam_path, sam_opts);
  if (!sam_or.ok()) { std::fprintf(stderr, "bam open failed\n"); return 1; }
  auto sam = std::move(sam_or.ValueOrDie());

  AlleleCounterOptions ac_opts;
  ac_opts.set_partition_size(1000);
  ac_opts.mutable_read_requirements()->set_min_mapping_quality(10);
  ac_opts.mutable_read_requirements()->set_min_base_quality(10);
  ac_opts.mutable_read_requirements()->set_min_base_quality_mode(
      nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
  ac_opts.set_track_ref_reads(true);

  AlleleCounter counter(ref.get(), region, /*candidates=*/{}, ac_opts);

  auto reads_or = sam->Query(region);
  if (!reads_or.ok()) { std::fprintf(stderr, "bam query failed\n"); return 1; }
  auto& reads_iter = reads_or.ValueOrDie();
  long n_reads = 0;
  nucleus::genomics::v1::Read read;
  while (true) {
    auto next = reads_iter->Next(&read);
    if (!next.ok() || !next.ValueOrDie()) break;
    counter.Add(read, "sample");
    ++n_reads;
  }
  reads_iter->Release().IgnoreError();
  std::fprintf(stderr, "%ld reads added\n", n_reads);

  for (const auto& ac : counter.Counts()) {
    int alt_total = 0;
    for (const auto& [name, allele] : ac.read_alleles()) {
      if (allele.type() != AlleleType::REFERENCE) ++alt_total;
    }
    if (alt_total == 0) continue;  // skip pure-ref positions
    const int64_t pos1 = ac.position().position() + 1;
    std::cout << ac.position().reference_name() << '\t' << pos1
              << '\t' << ac.ref_base()
              << "\tref=" << ac.ref_supporting_read_count();
    // Aggregate alts: bases→count, marking low-quality.
    std::map<std::string, std::pair<int,int>> alts;  // bases → (hq, lq)
    for (const auto& [name, allele] : ac.read_alleles()) {
      if (allele.type() == AlleleType::REFERENCE) continue;
      auto& p = alts[allele.bases()];
      if (allele.is_low_quality()) ++p.second; else ++p.first;
    }
    for (const auto& [bases, hq_lq] : alts) {
      std::cout << "\talt:" << bases << "=" << hq_lq.first;
      if (hq_lq.second) std::cout << "+" << hq_lq.second << "L";
    }
    std::cout << '\n';
  }
  return 0;
}
