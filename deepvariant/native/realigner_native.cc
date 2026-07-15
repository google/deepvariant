#include "deepvariant/native/realigner_native.h"

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "deepvariant/realigner/debruijn_graph.h"
#include "deepvariant/realigner/fast_pass_aligner.h"
#include "deepvariant/realigner/window_selector.h"
#include "absl/log/log.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/util/utils.h"

namespace deepvariant {

namespace {

using ::learning::genomics::deepvariant::AlleleCounter;
using ::learning::genomics::deepvariant::DeBruijnGraph;
using ::learning::genomics::deepvariant::FastPassAligner;
using ::learning::genomics::deepvariant::RealignerOptions;
using ::learning::genomics::deepvariant::VariantReadsWindowSelectorCandidates;

constexpr int kRefAlignMargin = 20;

// Container of a candidate window + assigned reads — mirror of
// deepvariant/realigner/realigner.py:AssemblyRegion.
struct AssemblyRegion {
  nucleus::genomics::v1::Range region;
  std::vector<std::string> haplotypes;
  std::vector<int> read_indices;  // indices into the input reads vector
  // Read span: the minimal interval covering ALL assigned reads' alignment
  // positions. Used to extend the ref window passed to FastPassAligner so
  // reads sticking out of `region` still align cleanly.
  int64_t read_span_start = 0;
  int64_t read_span_end   = 0;
};

// Returns [read_start, read_end_exclusive) in reference coords, derived
// from alignment.position + the cigar's reference span (mirrors
// nucleus.util.utils.read_end on the Python side).
std::pair<int64_t, int64_t> ReadRefSpan(
    const nucleus::genomics::v1::Read& read) {
  const int64_t start = read.alignment().position().position();
  int64_t ref_len = 0;
  for (const auto& cu : read.alignment().cigar()) {
    using ::nucleus::genomics::v1::CigarUnit;
    switch (cu.operation()) {
      case CigarUnit::ALIGNMENT_MATCH:
      case CigarUnit::SEQUENCE_MATCH:
      case CigarUnit::SEQUENCE_MISMATCH:
      case CigarUnit::DELETE:
      case CigarUnit::SKIP:
        ref_len += cu.operation_length();
        break;
      default:
        break;  // INSERT, soft/hard clip, pad don't consume reference.
    }
  }
  return {start, start + ref_len};
}

// Merge candidate positions into windows of width 2 × min_windows_distance.
// Mirror of window_selector._candidates_to_windows. Only positions with
// `count` in [min_supporting_reads, max_supporting_reads] count as
// "candidates" — without that range filter every read mismatch becomes a
// window seed and the whole region collapses into one >max_window_size
// window that gets discarded.
std::vector<nucleus::genomics::v1::Range> CandidatesToWindows(
    const std::vector<int>& candidate_counts,
    int region_start_pos, const std::string& chrom,
    int min_windows_distance, int max_window_size,
    int min_supporting_reads, int max_supporting_reads) {
  std::vector<nucleus::genomics::v1::Range> windows;
  int start_pos = -1, end_pos = -1;
  auto add_window = [&](int s, int e) {
    nucleus::genomics::v1::Range r;
    r.set_reference_name(chrom);
    r.set_start(std::max(0, s - min_windows_distance));
    r.set_end(e + min_windows_distance);
    if (r.end() - r.start() <= max_window_size) {
      windows.push_back(std::move(r));
    }
  };
  for (int i = 0; i < static_cast<int>(candidate_counts.size()); ++i) {
    const int c = candidate_counts[i];
    if (c < min_supporting_reads || c > max_supporting_reads) continue;
    const int pos = region_start_pos + i;
    if (start_pos == -1) {
      start_pos = end_pos = pos;
    } else if (pos > end_pos + 2 * min_windows_distance) {
      add_window(start_pos, end_pos);
      start_pos = end_pos = pos;
    } else {
      end_pos = pos;
    }
  }
  if (start_pos != -1) add_window(start_pos, end_pos);
  return windows;
}

}  // namespace

RealignerOptions DefaultRealignerOptions() {
  RealignerOptions opts;
  // Window selector — defaults from realigner.py.
  auto* ws = opts.mutable_ws_config();
  ws->set_min_num_supporting_reads(2);
  ws->set_max_num_supporting_reads(300);
  ws->set_min_mapq(20);
  ws->set_min_base_quality(20);
  ws->set_min_windows_distance(80);
  ws->set_max_window_size(1000);
  // Mirrors upstream's `_MIN_ALLELE_SUPPORT = 2` in realigner.py — without
  // this, AlleleFilter() in window_selector.cc accepts singleton alleles
  // (count=1), so positions with only one supporting read can still seed
  // a candidate window. Upstream rejects them.
  ws->set_min_allele_support(2);
  // 20bp on each side. Mirrors realigner.py:_WS_REGION_EXPANSION_IN_BP.
  // Used by RealignReadsForRegion when building the WindowSelector
  // AlleleCounter (so reads that overhang the region edges still
  // contribute counts at boundary positions).
  ws->set_region_expansion_in_bp(20);
  // De-Bruijn graph — defaults from realigner.py.
  auto* dbg = opts.mutable_dbg_config();
  dbg->set_min_k(10);
  dbg->set_max_k(101);
  dbg->set_step_k(1);
  dbg->set_min_mapq(14);
  dbg->set_min_base_quality(15);
  dbg->set_min_edge_weight(2);
  dbg->set_max_num_paths(256);
  // Aligner — defaults.
  auto* aln = opts.mutable_aln_config();
  aln->set_match(4);
  aln->set_mismatch(6);
  aln->set_gap_open(8);
  aln->set_gap_extend(2);
  aln->set_k(23);
  aln->set_error_rate(0.01);
  aln->set_max_num_of_mismatches(2);
  aln->set_realignment_similarity_threshold(0.16934);
  aln->set_kmer_size(32);
  return opts;
}

std::vector<nucleus::genomics::v1::Read> RealignReadsForRegion(
    const std::vector<nucleus::genomics::v1::Read>& reads,
    const nucleus::genomics::v1::Range& region,
    const AlleleCounter& counter,
    const nucleus::GenomeReference& ref_reader,
    const RealignerOptions& options) {
  if (reads.empty()) return reads;

  // ── Step 1: candidate counts via WindowSelector ──────────────────────────
  std::vector<int> counts =
      VariantReadsWindowSelectorCandidates(counter, options.ws_config());

  // ── Step 2: merge counts into windows ────────────────────────────────────
  auto windows = CandidatesToWindows(
      counts, static_cast<int>(region.start()), region.reference_name(),
      options.ws_config().min_windows_distance(),
      options.ws_config().max_window_size(),
      options.ws_config().min_num_supporting_reads(),
      options.ws_config().max_num_supporting_reads());

  LOG(INFO) << "  realigner: " << windows.size() << " candidate windows in "
            << region.reference_name() << ":" << region.start() << "-"
            << region.end()
            << " (positions with non-zero counts: "
            << std::count_if(counts.begin(), counts.end(),
                             [](int c) { return c > 0; })
            << ")";
  if (windows.empty()) return reads;

  // ── Step 3: build DeBruijn graphs → assembled regions ────────────────────
  std::vector<AssemblyRegion> assembled;
  // We need ConstProtoPtr<Read> for DeBruijnGraph::Build.
  std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
      read_ptrs;
  read_ptrs.reserve(reads.size());
  for (const auto& r : reads) {
    read_ptrs.push_back(
        nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>(&r));
  }

  // Optional per-window CSV diagnostic output, mirroring the columns
  // upstream's DiagnosticLogger writes when --realigner_diagnostics is on:
  // window, k, n_haplotypes, n_reads_in_window. Enabled by setting
  // DV_REALIGNER_DIAG_CSV to a path; a header row is written on first call.
  static std::ofstream* diag_csv = []() -> std::ofstream* {
    const char* p = std::getenv("DV_REALIGNER_DIAG_CSV");
    if (!p || !*p) return nullptr;
    auto* f = new std::ofstream(p);
    if (!f->is_open()) return nullptr;
    *f << "window,k,n_haplotypes,n_reads_in_window,hap_hash\n";
    return f;
  }();
  // Optional second log: full haplotype strings per window. Enabled by
  // DV_REALIGNER_DIAG_HAP set to a directory; one file per window with
  // one haplotype per line.
  static const char* hap_dump_dir = std::getenv("DV_REALIGNER_DIAG_HAP");

  for (const auto& window : windows) {
    auto ref_or = ref_reader.GetBases(window);
    if (!ref_or.ok()) continue;
    const std::string ref_bases = ref_or.ValueOrDie();

    // Filter reads overlapping the window (DeBruijnGraph filters internally
    // by mapq/base_quality from dbg_config).
    std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
        win_reads;
    for (size_t i = 0; i < reads.size(); ++i) {
      if (nucleus::ReadOverlapsRegion(reads[i], window)) {
        win_reads.push_back(read_ptrs[i]);
      }
    }
    if (win_reads.empty()) continue;

    // Build the graph.
    std::vector<nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
        win_reads_copy = win_reads;
    auto graph = DeBruijnGraph::Build(ref_bases, win_reads_copy,
                                       options.dbg_config());
    std::vector<std::string> haplotypes;
    int k_used = -1;
    if (graph) {
      haplotypes = graph->CandidateHaplotypes();
      k_used = graph->KmerSize();
    }
    if (diag_csv) {
      // Hash the sorted haplotype set so we can diff against upstream
      // beyond just the count. Use a simple fold so the value is stable.
      uint64_t hap_hash = 1469598103934665603ULL;  // FNV-64 offset
      for (const auto& h : haplotypes) {
        for (unsigned char c : h) {
          hap_hash ^= c;
          hap_hash *= 1099511628211ULL;
        }
        hap_hash ^= '|';
      }
      *diag_csv << window.reference_name() << ":" << (window.start() + 1)
                << "-" << window.end() << "," << k_used << ","
                << haplotypes.size() << "," << win_reads.size()
                << "," << hap_hash << "\n";
    }
    if (hap_dump_dir && !haplotypes.empty()) {
      std::string fname = std::string(hap_dump_dir) + "/" +
          window.reference_name() + ":" +
          std::to_string(window.start() + 1) + "-" +
          std::to_string(window.end()) + ".txt";
      std::ofstream hf(fname);
      if (hf) {
        for (const auto& h : haplotypes) hf << h << "\n";
      }
    }
    if (haplotypes.empty() ||
        (haplotypes.size() == 1 && haplotypes[0] == ref_bases)) {
      continue;  // Nothing to realign in this window.
    }

    AssemblyRegion ar;
    ar.region = window;
    ar.haplotypes = std::move(haplotypes);
    assembled.push_back(std::move(ar));
  }
  if (diag_csv) diag_csv->flush();

  LOG(INFO) << "  realigner: " << assembled.size()
            << " assembled regions in "
            << region.reference_name() << ":" << region.start() << "-"
            << region.end();
  if (assembled.empty()) return reads;

  // ── Step 4: assign reads to assembled regions (max-overlap wins) ────────
  // Mirrors realigner.py:assign_reads_to_assembled_regions. For each read,
  // pick the assembled region with the maximum reference overlap; first
  // index wins in case of ties.
  std::vector<bool> read_assigned(reads.size(), false);
  std::vector<bool> ar_span_init(assembled.size(), false);
  for (size_t i = 0; i < reads.size(); ++i) {
    const auto [rs, re] = ReadRefSpan(reads[i]);
    int best_ar = -1;
    int64_t best_overlap = 0;
    for (size_t a = 0; a < assembled.size(); ++a) {
      const auto& reg = assembled[a].region;
      const int64_t lo = std::max<int64_t>(rs, reg.start());
      const int64_t hi = std::min<int64_t>(re, reg.end());
      const int64_t ov = hi - lo;
      if (ov > best_overlap) {
        best_overlap = ov;
        best_ar = static_cast<int>(a);
      }
    }
    if (best_ar < 0) continue;  // read doesn't overlap any assembled region
    auto& ar = assembled[best_ar];
    ar.read_indices.push_back(static_cast<int>(i));
    read_assigned[i] = true;
    if (!ar_span_init[best_ar]) {
      ar.read_span_start = rs;
      ar.read_span_end   = re;
      ar_span_init[best_ar] = true;
    } else {
      ar.read_span_start = std::min(ar.read_span_start, rs);
      ar.read_span_end   = std::max(ar.read_span_end,   re);
    }
  }
  for (size_t a = 0; a < assembled.size(); ++a) {
    if (!ar_span_init[a]) {
      assembled[a].read_span_start = assembled[a].region.start();
      assembled[a].read_span_end   = assembled[a].region.end();
    }
  }

  // Start with the unassigned reads (they pass through unchanged).
  std::vector<nucleus::genomics::v1::Read> out;
  out.reserve(reads.size());
  for (size_t i = 0; i < reads.size(); ++i) {
    if (!read_assigned[i]) out.push_back(reads[i]);
  }

  // ── Step 5: realign each assembled region's reads ────────────────────────
  for (const auto& ar : assembled) {
    if (ar.read_indices.empty()) continue;

    const std::string& chrom = ar.region.reference_name();
    auto contig_or = ref_reader.Contig(chrom);
    if (!contig_or.ok()) continue;
    const int64_t contig_n_bases = contig_or.ValueOrDie()->n_bases();

    // Match realigner.py: extend the ref window to the broader of the
    // assembled region and the actual reads' alignment span, plus margin.
    //   ref_start = max(0, min(read_span.start, region.start) - margin)
    //   ref_end   = min(contig_n, max(read_span.end,   region.end)   + margin)
    // The interior "window" handed to FastPassAligner stays bounded by the
    // assembled region — only the prefix/suffix grow when reads overhang.
    const int64_t span_start =
        std::min<int64_t>(ar.read_span_start, ar.region.start());
    const int64_t span_end =
        std::max<int64_t>(ar.read_span_end, ar.region.end());
    const int64_t ref_start =
        std::max<int64_t>(0, span_start - kRefAlignMargin);
    const int64_t ref_end =
        std::min<int64_t>(contig_n_bases, span_end + kRefAlignMargin);
    if (ref_end <= ar.region.end()) {
      // Mirror realigner.py:call_fast_pass_aligner — if the contig is too
      // short to form a suffix, return the region's reads unchanged. The
      // prefix can be empty (region at contig start) and FastPassAligner
      // handles that fine.
      for (int idx : ar.read_indices) out.push_back(reads[idx]);
      continue;
    }

    nucleus::genomics::v1::Range pre_range, win_range, suf_range;
    pre_range.set_reference_name(chrom);
    pre_range.set_start(ref_start);
    pre_range.set_end(ar.region.start());
    win_range = ar.region;
    suf_range.set_reference_name(chrom);
    suf_range.set_start(ar.region.end());
    suf_range.set_end(ref_end);

    auto ref_pre_or = ref_reader.GetBases(pre_range);
    auto ref_win_or = ref_reader.GetBases(win_range);
    auto ref_suf_or = ref_reader.GetBases(suf_range);
    if (!ref_pre_or.ok() || !ref_win_or.ok() || !ref_suf_or.ok()) {
      for (int idx : ar.read_indices) out.push_back(reads[idx]);
      continue;
    }
    const std::string ref_pre = ref_pre_or.ValueOrDie();
    const std::string ref_win = ref_win_or.ValueOrDie();
    const std::string ref_suf = ref_suf_or.ValueOrDie();
    const std::string ref_seq = ref_pre + ref_win + ref_suf;

    // Build the per-region read vector for the aligner.
    std::vector<nucleus::genomics::v1::Read> region_reads;
    region_reads.reserve(ar.read_indices.size());
    for (int idx : ar.read_indices) region_reads.push_back(reads[idx]);

    FastPassAligner aligner;
    auto aln_cfg = options.aln_config();
    aln_cfg.set_read_size(static_cast<int>(
        region_reads[0].aligned_sequence().size()));
    aln_cfg.set_force_alignment(false);
    aligner.set_options(aln_cfg);
    // BUG FIX (Path D Site 1, 2026-05-23): mirror upstream
    // realigner.py:call_fast_pass_aligner:779 which propagates
    // RealignerOptions.normalize_reads onto the FastPassAligner. Without
    // this, fast_pass_aligner.cc:557-568 discards any realigned alignment
    // whose CIGAR is not already left-normalized — silently dropping
    // valid shifts in T-homopolymer regions and leaving the read at its
    // original POS (the +1 DP / 1-read-off WG residual at chr12:62946475).
    aligner.set_normalize_reads(options.normalize_reads());
    aligner.set_reference(ref_seq);
    aligner.set_ref_start(chrom, static_cast<uint64_t>(ref_start));
    aligner.set_ref_prefix_len(static_cast<int>(ref_pre.size()));
    aligner.set_ref_suffix_len(static_cast<int>(ref_suf.size()));
    std::vector<std::string> haplotypes_padded;
    haplotypes_padded.reserve(ar.haplotypes.size());
    for (const auto& hap : ar.haplotypes) {
      haplotypes_padded.push_back(ref_pre + hap + ref_suf);
    }
    aligner.set_haplotypes(haplotypes_padded);

    auto realigned = aligner.AlignReads(absl::MakeConstSpan(region_reads));
    if (realigned && !realigned->empty()) {
      for (auto& r : *realigned) out.push_back(std::move(r));
    } else {
      // Fallback to original reads if alignment failed.
      for (int idx : ar.read_indices) out.push_back(reads[idx]);
    }
  }

  return out;
}

}  // namespace deepvariant
