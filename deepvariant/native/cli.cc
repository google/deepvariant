#include "deepvariant/native/cli.h"

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <system_error>
#include <thread>
#include <unordered_map>
#include <vector>

#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/flags/reflection.h"
#include "absl/flags/usage.h"
#include "absl/flags/usage_config.h"
#include <sys/sysctl.h>
#include "absl/log/globals.h"
#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "absl/strings/match.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"

// `run` subcommand flags. Reuse flags declared in the subcommand files
// (--reads, --ref, --regions, --batch_size, --num_shards, etc.) to avoid
// duplicate symbols at link time.
ABSL_FLAG(bool, include_alt_contigs, false,
          "If true, process alt/random/decoy/unplaced contigs (chr*_random, "
          "chrUn_*, etc.) in addition to canonical chromosomes. Default false "
          "to match google/deepvariant Docker behavior, which emits records "
          "only on chr1..22, chrX, chrY, chrM. Without this filter, our binary "
          "emits ~138k alt-contig records that Docker doesn't, breaking "
          "FILTER parity at WG scale.");

ABSL_FLAG(std::string, model_type, "WGS",
          "Model type: WGS, WES, PACBIO, ONT (alias ONT_R104), "
          "HYBRID_PACBIO_ILLUMINA");
ABSL_FLAG(std::string, output_vcf, "", "Output VCF path (run mode).");
ABSL_FLAG(std::string, output_gvcf, "",
          "Output gVCF path (run mode, optional).");
ABSL_FLAG(std::string, intermediate_results_dir, "/tmp/dv_run",
          "Directory for intermediate TFRecord files.");
ABSL_FLAG(std::string, model, "",
          "Path to .mlpackage model (overrides --model_type lookup).");
ABSL_FLAG(std::string, small_model_path, "",
          "Path to small_model .mlpackage. Empty = no small-model first-pass; "
          "every candidate goes through the big InceptionV3.");

// Flags owned by the subcommand .cc files — declared here for use by RunAll.
ABSL_DECLARE_FLAG(std::string, reads);
ABSL_DECLARE_FLAG(std::string, ref);
ABSL_DECLARE_FLAG(std::string, regions);
ABSL_DECLARE_FLAG(std::string, select_variant_types);
ABSL_DECLARE_FLAG(std::string, haploid_contigs);     // postprocess_main.cc
ABSL_DECLARE_FLAG(std::string, par_regions_bed);      // postprocess_main.cc
ABSL_DECLARE_FLAG(int, num_shards);
ABSL_DECLARE_FLAG(int, batch_size);
ABSL_DECLARE_FLAG(std::string, inference_backend);
ABSL_DECLARE_FLAG(std::string, ane_speculate_metal_checkpoint);
ABSL_DECLARE_FLAG(std::string, ane_speculate_metal_checkpoint_child);
ABSL_DECLARE_FLAG(std::string, ane_speculate_metal_checkpoint_parent);
ABSL_DECLARE_FLAG(std::string, ane_speculate_metal_checkpoint_somatic);
ABSL_DECLARE_FLAG(std::string, ane_speculate_metal_checkpoint_pangenome);
ABSL_DECLARE_FLAG(double, ane_speculate_confidence);

// Helper: append --ane_speculate_metal_checkpoint=... + threshold to the
// argv vector being passed to a sub-process call_variants invocation.
// `metal_ckpt` is the role-specific .dvw bundle for the GPU FP32 rerun
// (empty → call_variants will error out if backend == ane_speculate).
namespace {
inline void AppendAneSpeculateArgs(std::vector<std::string>& cv_args,
                                    const std::string& inference_backend,
                                    const std::string& metal_ckpt) {
  if (inference_backend != "ane_speculate") return;
  if (!metal_ckpt.empty()) {
    cv_args.push_back(absl::StrCat(
        "--ane_speculate_metal_checkpoint=", metal_ckpt));
  }
  cv_args.push_back(absl::StrCat(
      "--ane_speculate_confidence=",
      absl::GetFlag(FLAGS_ane_speculate_confidence)));
}

// Helper: forward the user's --select_variant_types (if any) to make_examples.
// make_examples_main.cc owns the flag and does the actual candidate filtering;
// here we just pass the value through to each per-mode make_examples worker.
inline void AppendSelectVariantTypes(std::vector<std::string>& me_args) {
  const std::string svt = absl::GetFlag(FLAGS_select_variant_types);
  if (!svt.empty()) {
    me_args.push_back(absl::StrCat("--select_variant_types=", svt));
  }
}

// Helper: forward sex-chromosome haploid-calling flags to a worker's argv.
// Called for both the make_examples and postprocess args: postprocess_main.cc
// defines the flags and applies the het-zeroing correction, and
// make_examples_main.cc declares them and uses them for haploid-aware gVCF
// reference confidence. The arg vector is the per-mode worker's command line.
inline void AppendHaploidFlags(std::vector<std::string>& args) {
  const std::string hc = absl::GetFlag(FLAGS_haploid_contigs);
  if (!hc.empty()) args.push_back(absl::StrCat("--haploid_contigs=", hc));
  const std::string par = absl::GetFlag(FLAGS_par_regions_bed);
  if (!par.empty()) args.push_back(absl::StrCat("--par_regions_bed=", par));
}
}  // namespace
ABSL_DECLARE_FLAG(std::string, checkpoint);
// Phase 9 / Step 1 — alt-aligned pileup mode (PacBio/ONT). Defined in
// make_examples_main.cc; cli.cc reads it to pick a sensible per-model
// default ("diff_channels" for PACBIO/ONT, "none" for WGS/WES) before
// passing it down to make_examples.
ABSL_DECLARE_FLAG(std::string, alt_aligned_pileup);

// DeepTrio (Step 1.5) — trio mode flags. When --reads_parent1 is set,
// run mode dispatches 3× call_variants with the appropriate child/parent
// model and writes 3 separate VCFs.
ABSL_DECLARE_FLAG(std::string, reads_parent1);
ABSL_DECLARE_FLAG(std::string, reads_parent2);
ABSL_DECLARE_FLAG(std::string, sample_name_parent1);
ABSL_DECLARE_FLAG(std::string, sample_name_parent2);
ABSL_DECLARE_FLAG(std::string, examples_child);
ABSL_DECLARE_FLAG(std::string, examples_parent1);
ABSL_DECLARE_FLAG(std::string, examples_parent2);
ABSL_DECLARE_FLAG(std::string, small_model_path_child);
ABSL_DECLARE_FLAG(std::string, small_model_path_parent);
ABSL_DECLARE_FLAG(std::string, small_model_cvo_outfile_child);
ABSL_DECLARE_FLAG(std::string, small_model_cvo_outfile_parent1);
ABSL_DECLARE_FLAG(std::string, small_model_cvo_outfile_parent2);
ABSL_FLAG(std::string, checkpoint_child, "",
          "Trio mode: model checkpoint (.dvw or .mlpackage) for child.");
ABSL_FLAG(std::string, checkpoint_parent, "",
          "Trio mode: model checkpoint shared by parent1 and parent2.");
ABSL_FLAG(std::string, output_vcf_child, "",
          "Trio mode: output VCF for the child sample.");
ABSL_FLAG(std::string, output_vcf_parent1, "",
          "Trio mode: output VCF for parent1.");
ABSL_FLAG(std::string, output_vcf_parent2, "",
          "Trio mode: output VCF for parent2.");
ABSL_FLAG(std::string, output_gvcf_child, "",
          "Trio mode: output gVCF for the child sample.");
ABSL_FLAG(std::string, output_gvcf_parent1, "",
          "Trio mode: output gVCF for parent1.");
ABSL_FLAG(std::string, output_gvcf_parent2, "",
          "Trio mode: output gVCF for parent2.");

// DeepSomatic (Step 2) — somatic mode flags. When --reads_tumor is set,
// run mode dispatches 1× call_variants on the tumor model and emits a
// single tumor VCF. tumor_only mode = no --reads_normal.
ABSL_DECLARE_FLAG(std::string, reads_tumor);
ABSL_DECLARE_FLAG(std::string, reads_normal);
ABSL_DECLARE_FLAG(std::string, sample_name_tumor);
ABSL_DECLARE_FLAG(std::string, sample_name_normal);
ABSL_DECLARE_FLAG(std::string, examples_tumor);
ABSL_DECLARE_FLAG(std::string, examples_normal);
ABSL_DECLARE_FLAG(std::string, small_model_path_somatic);
ABSL_DECLARE_FLAG(std::string, population_vcfs);
ABSL_DECLARE_FLAG(std::string, pon_filtering);
ABSL_DECLARE_FLAG(double, vsc_max_fraction_snps_for_non_target_sample);
ABSL_DECLARE_FLAG(double, vsc_max_fraction_indels_for_non_target_sample);
ABSL_DECLARE_FLAG(bool,   sort_by_alt_allele_support_somatic);
ABSL_DECLARE_FLAG(bool,   small_model_use_haplotypes);
ABSL_DECLARE_FLAG(bool,   use_direct_phasing);
ABSL_DECLARE_FLAG(std::string, small_model_cvo_outfile_tumor);
ABSL_DECLARE_FLAG(int, pileup_image_height_tumor);
ABSL_DECLARE_FLAG(int, pileup_image_height_normal);

// Pangenome-aware DV (Step 3) — When --reads_pangenome is set, run mode
// dispatches a 2-sample pangenome pipeline (pangenome=0, reads=1=main).
// Single VCF output for the reads sample.
ABSL_DECLARE_FLAG(std::string, reads_pangenome);
ABSL_DECLARE_FLAG(std::string, sample_name_pangenome);
ABSL_DECLARE_FLAG(std::string, sample_name_reads);
ABSL_DECLARE_FLAG(std::string, examples_reads);
ABSL_DECLARE_FLAG(std::string, small_model_path_pangenome);
ABSL_DECLARE_FLAG(std::string, small_model_cvo_outfile_reads);
ABSL_DECLARE_FLAG(int, pileup_image_height_pangenome);
ABSL_DECLARE_FLAG(int, pileup_image_height_reads);

namespace deepvariant {

namespace {

// Build a flag-vector for a subcommand, splicing in the given extra flags.
std::vector<char*> MakeArgv(const std::string& prog,
                             const std::vector<std::string>& extras) {
  static std::vector<std::string> storage;
  storage.clear();
  storage.push_back(prog);
  for (const auto& e : extras) storage.push_back(e);
  std::vector<char*> argv;
  for (auto& s : storage) argv.push_back(const_cast<char*>(s.c_str()));
  argv.push_back(nullptr);
  return argv;
}

// Auto-detect a sensible default for num_shards/threads. Uses
// std::thread::hardware_concurrency() (returns logical cores) and
// reserves 2 for the system (so an M4 Max with 16 cores returns 14).
// If --num_shards was set explicitly to a value > 1, that wins.
int AutoNumShards() {
  int hw = static_cast<int>(std::thread::hardware_concurrency());
  if (hw <= 0) return 1;
  if (hw <= 4) return hw;        // tiny machines: use all cores
  return std::max(1, hw - 2);    // leave headroom on bigger machines
}

int EffectiveNumShards() {
  const int explicit_n = absl::GetFlag(FLAGS_num_shards);
  // 0 (default) and 1 (=no sharding) both fall back to auto-detect.
  if (explicit_n > 1) return explicit_n;
  return AutoNumShards();
}

// IsCanonicalContig — return true if the contig name matches a canonical
// chromosome: chr1..22, chrX, chrY, chrM, chrMT (or the no-prefix forms
// 1..22, X, Y, M, MT). Reject anything with `_` (alt/random/decoy/unplaced)
// or anything that's not numeric / X / Y / M[T].
//
// Docker's run_deepvariant emits records only on canonical contigs,
// even when the BAM has reads on alt-contigs (verified empirically:
// HG002 BAM has 1.5M reads on chrUn_KI270438v1 but Docker emits 0
// records there). This helper drives our default filter to match.
bool IsCanonicalContig(absl::string_view name) {
  if (name.empty()) return false;
  // Reject anything with underscore (alt/random/decoy/unplaced).
  if (name.find('_') != absl::string_view::npos) return false;
  // Strip optional `chr` prefix.
  absl::string_view bare = name;
  if (bare.size() > 3 && bare.substr(0, 3) == "chr") bare.remove_prefix(3);
  // Sex chroms / mito.
  if (bare == "X" || bare == "Y" || bare == "M" || bare == "MT") return true;
  // Numeric 1..22.
  if (bare.empty()) return false;
  for (char c : bare) {
    if (c < '0' || c > '9') return false;
  }
  int n = 0;
  if (!absl::SimpleAtoi(bare, &n)) return false;
  return n >= 1 && n <= 22;
}

// DefaultCanonicalRegions — when --regions is empty AND
// --include_alt_contigs=false, return a comma-separated list of all
// canonical contigs from the reference's .fai index. Matches Docker's
// implicit canonical-only filter.
//
// Returns empty string on any failure (missing .fai, no canonical
// contigs found, etc.); caller falls through to the no-regions path
// in that case.
std::string DefaultCanonicalRegions(const std::string& ref_path) {
  const std::string fai_path = absl::StrCat(ref_path, ".fai");
  std::ifstream fai(fai_path);
  if (!fai) return "";
  std::vector<std::string> canonical;
  std::string line;
  while (std::getline(fai, line)) {
    const auto tab = line.find('\t');
    if (tab == std::string::npos) continue;
    const std::string name = line.substr(0, tab);
    if (IsCanonicalContig(name)) canonical.push_back(name);
  }
  if (canonical.empty()) return "";
  return absl::StrJoin(canonical, ",");
}

// CanonicalizeRegions — expand bare contig names (e.g. "chr20") to the
// explicit "chr20:1-N" form using the reference .fai. Mixed input like
// "chr20 chr21:1-100" is supported: bare names get expanded, ranges
// pass through unchanged.
//
// Why we do this: empirically, passing a bare contig name vs
// "chr20:1-64444167" through the WES pipeline produces different VCF
// record counts (19,740 vs 210,619 on chr20-full) — same Range proto
// emerges from BuildCallingRegions but somewhere downstream the bare-
// name form drops records. The bug only surfaces in WES mode at the
// full-contig scale (chr20:10M-10.1M fixture matches Docker in either
// form). Rather than chase the elusive downstream divergence, we
// canonicalize at the cli.cc boundary so every make_examples
// invocation receives the explicit-range form. F1 + FILTER parity
// already verified for the explicit form on chr20-full (210,619
// records = Docker 210,390 ± record-set drift, F1 = Docker).
std::string CanonicalizeRegions(const std::string& regions,
                                 const std::string& ref_path) {
  if (regions.empty()) return regions;
  // Build a contig length map from the .fai.
  const std::string fai_path = absl::StrCat(ref_path, ".fai");
  std::ifstream fai(fai_path);
  if (!fai) return regions;  // can't expand; pass through (callers tolerate)
  std::unordered_map<std::string, int64_t> lengths;
  std::string line;
  while (std::getline(fai, line)) {
    const auto tab = line.find('\t');
    if (tab == std::string::npos) continue;
    const std::string name = line.substr(0, tab);
    const std::string rest = line.substr(tab + 1);
    const auto tab2 = rest.find('\t');
    const std::string len_str =
        tab2 == std::string::npos ? rest : rest.substr(0, tab2);
    int64_t len;
    if (absl::SimpleAtoi(len_str, &len)) lengths[name] = len;
  }
  // Split the regions string on the same delimiters as make_examples,
  // canonicalize each token, then re-join.
  std::vector<std::string> tokens = absl::StrSplit(
      regions, absl::ByAnyChar(" \t,"), absl::SkipEmpty());
  std::vector<std::string> out;
  out.reserve(tokens.size());
  for (const auto& t : tokens) {
    if (t.find(':') != std::string::npos) {
      out.push_back(t);  // already has range
      continue;
    }
    auto it = lengths.find(t);
    if (it == lengths.end()) {
      out.push_back(t);  // unknown contig; let make_examples error out
      continue;
    }
    out.push_back(absl::StrCat(t, ":1-", it->second));
  }
  return absl::StrJoin(out, " ");
}

// EffectiveRegions — resolve the regions string to use for make_examples.
//   - If --regions is non-empty: pass through (user explicitly chose),
//     after canonicalizing bare contig names.
//   - Else if --include_alt_contigs=true: pass through empty (process all).
//   - Else: build canonical list from reference .fai (matches Docker),
//     already in explicit form via DefaultCanonicalRegions.
std::string EffectiveRegions(const std::string& user_regions,
                              const std::string& ref_path) {
  if (!user_regions.empty()) {
    return CanonicalizeRegions(user_regions, ref_path);
  }
  if (absl::GetFlag(FLAGS_include_alt_contigs)) return "";
  // DefaultCanonicalRegions also returns bare contig names; canonicalize too.
  return CanonicalizeRegions(DefaultCanonicalRegions(ref_path), ref_path);
}

// Auto-detect a sensible default for --batch_size based on physical
// RAM. The MPSGraph Inception-v3 forward pass at FP32 holds peak
// activations of ~5 MB per example mid-network plus ~100 MB of
// constant weights. Larger batches amortise the per-batch dispatch
// overhead (~50 ms) but consume proportionally more unified memory.
//
// Tiered conservative table (peak GPU footprint ≤ 50 % of physical
// RAM, leaving headroom for the OS, htslib mmap, and other tools):
//
//   < 16 GB   → batch_size 128   (8 GB Macs)
//   16-32 GB  → batch_size 512   (16 GB Macs: M1/M2/M3 Pro entry)
//   32-64 GB  → batch_size 1024  (32 GB Pro/Max, 36 GB M4 Pro)
//   ≥ 64 GB   → batch_size 2048  (64 GB+ Max/Ultra/M4 Max)
//
// User can override with --batch_size=N at any time. The auto-detect
// only kicks in when the flag is at its default value.
//
// We read RAM via sysctl(hw.memsize) which is the physical RAM in
// bytes — works on every Mac since macOS 10.0, no entitlements.
int AutoBatchSize() {
  uint64_t mem_bytes = 0;
  size_t len = sizeof(mem_bytes);
  // sysctlbyname is the macOS-portable way; #include <sys/sysctl.h> at
  // the top of the file (added below).
  if (sysctlbyname("hw.memsize", &mem_bytes, &len, nullptr, 0) != 0 ||
      mem_bytes == 0) {
    return 512;  // safe fallback
  }
  const uint64_t mem_gb = mem_bytes >> 30;  // approximate GiB
  if (mem_gb < 16) return 128;
  if (mem_gb < 32) return 512;
  if (mem_gb < 64) return 1024;
  return 2048;
}

int EffectiveBatchSize() {
  // Distinguish "user passed --batch_size on cmdline" from "default
  // value from the proto" via DefaultValue / CurrentValue string
  // comparison. (`IsSpecifiedOnCommandLine` is private in this abseil
  // version.) Edge case: a user passing exactly the default value
  // (128) gets the auto-detect path. Acceptable since 128 is the
  // smallest non-trivial value and AutoBatchSize ≥ 128 by design.
  if (auto* f = absl::FindCommandLineFlag("batch_size");
      f && f->CurrentValue() != f->DefaultValue()) {
    return absl::GetFlag(FLAGS_batch_size);
  }
  return AutoBatchSize();
}

// CanonicalModelType — canonicalize user-facing model_type aliases to the
// internal short form this file compares against. Upstream's canonical
// long-read ONT string is `ONT_R104` (also `ONT_R10` / `ONT_R9`), but every
// comparison in this file is against the short `"ONT"`. Without this, a user
// passing --model_type=ONT_R104 would be silently misrouted to the WGS/WES
// default branch. Upper-cases first to match the existing convention in this
// file (all comparisons here are against upper-case forms).
std::string CanonicalModelType(std::string mt) {
  for (char& c : mt) c = static_cast<char>(std::toupper(c));
  if (mt == "ONT_R104" || mt == "ONT_R10" || mt == "ONT_R9") return "ONT";
  return mt;
}

std::string ModelPath(const std::string& model_type) {
  if (!absl::GetFlag(FLAGS_model).empty()) {
    return absl::GetFlag(FLAGS_model);
  }
  // Default install path from deepvariant-models Homebrew formula.
  const char* prefix = std::getenv("DEEPVARIANT_MODELS_DIR");
  std::string base = prefix ? prefix : "/opt/homebrew/share/deepvariant-models";
  std::string type = model_type;
  // Normalise to lowercase.
  for (char& c : type) c = static_cast<char>(std::tolower(c));
  return absl::StrCat(base, "/", type, ".mlpackage");
}

// MergeCvoFiles — concatenate small-model CVO shards (if any) followed by the
// big-model CVO into a single `merged_cvo_path`, by raw byte copy. TFRecord
// allows naive byte concatenation because each record is self-delimiting.
//
// `small_cvo_pattern` is either empty (no small model → only big_cvo is
// copied), a single file path, or a `name@N` shard spec written one-file-per
// make_examples worker thread; in the `@N` case we expand to the conventional
// `name-NNNNN-of-NNNNN` shard files and append each in order. `big_cvo_path`
// is the single file written by call_variants.
//
// This replaces an earlier unquoted `std::system("cat ... > ...")` merge in
// the somatic/pangenome paths, which broke on paths containing spaces or shell
// metacharacters. Returns false (with LOG(ERROR)) on any I/O failure.
//
// Critical subtlety preserved from the germline merge: we read each source
// into a buffer and only write when non-empty. `operator<<(streambuf*)` sets
// the output stream's failbit when the source is empty, which silently breaks
// ALL subsequent writes — and sharded small_cvo routinely has 0-record shards.
bool MergeCvoFiles(const std::string& small_cvo_pattern,
                   const std::string& big_cvo_path,
                   const std::string& merged_cvo_path) {
  std::ofstream out(merged_cvo_path, std::ios::binary | std::ios::trunc);
  if (!out) {
    LOG(ERROR) << "Cannot open merged CVO: " << merged_cvo_path;
    return false;
  }
  bool ok = true;
  auto append_path = [&](const std::string& p) {
    std::ifstream in(p, std::ios::binary);
    if (!in) return;  // missing shard (e.g. 0-record) → skip, not an error.
    std::vector<char> buf((std::istreambuf_iterator<char>(in)),
                           std::istreambuf_iterator<char>());
    if (!buf.empty()) out.write(buf.data(), buf.size());
  };
  // small_cvo first (matching the germline / `cat small big` order), then big.
  if (!small_cvo_pattern.empty()) {
    const auto at = small_cvo_pattern.find('@');
    if (at == std::string::npos) {
      append_path(small_cvo_pattern);
    } else {
      const std::string prefix = small_cvo_pattern.substr(0, at);
      int nshard = 0;
      if (!absl::SimpleAtoi(small_cvo_pattern.substr(at + 1), &nshard) ||
          nshard <= 0) {
        LOG(ERROR) << "Bad small_cvo shard spec: " << small_cvo_pattern;
        return false;
      }
      for (int i = 0; i < nshard; ++i) {
        append_path(absl::StrCat(prefix, "-",
                                  absl::Dec(i, absl::kZeroPad5),
                                  "-of-", absl::Dec(nshard, absl::kZeroPad5)));
      }
    }
  }
  append_path(big_cvo_path);
  out.close();
  if (!out) {
    LOG(ERROR) << "Failed writing merged CVO: " << merged_cvo_path;
    ok = false;
  }
  return ok;
}

}  // namespace

// Forward decls.
int RunAllTrio(int argc, char** argv);
int RunAllSomatic(int argc, char** argv);
int RunAllPangenome(int argc, char** argv);

// ExpectsSmallModel — returns true if the model bundle for the given
// model_type declares a `trained_small_model_path` in upstream Docker's
// model.example_info.json. When this is true and the user passes an empty
// --small_model_path (resp. --small_model_path_child / _parent / _somatic),
// borderline-GQ candidates that Docker fast-paths through the small MLP go
// instead through the slower Inception-v3 path, and FILTER classification
// can drift from Docker. Long-read modes (PACBIO/ONT) regress particularly
// hard — empirically observed in B1+B2 validation 2026-05-07: ONT SNP F1
// dropped from 0.776 → 0.727 (-5%) when --small_model_path was omitted.
//
// Source of truth: tools/conversion/models/<dir>/model.example_info.json.
//   has trained_small_model_path  →  germline {WGS, ONT, PACBIO}
//                                    deepsomatic {WGS, ONT, PACBIO, FFPE_WGS}
//                                    (tumor+normal only — no tumor-only bundle
//                                     ships a small_model)
//   no trained_small_model_path   →  WES, MASSEQ, RNASEQ, HYBRID, all
//                                    tumor-only somatic, all FFPE_WES.
static bool GermlineExpectsSmallModel(const std::string& mt_upper) {
  return mt_upper == "WGS" || mt_upper == "ONT" || mt_upper == "PACBIO";
}
static bool SomaticExpectsSmallModel(const std::string& mt_upper,
                                     bool has_normal) {
  if (!has_normal) return false;  // no tumor-only bundle ships a small_model
  return mt_upper == "WGS" || mt_upper == "ONT" || mt_upper == "PACBIO" ||
         mt_upper == "FFPE_WGS";
}

// WarnIfMissingSmallModel — single-line LOG(WARNING) if `path` is empty and
// the bundle declares a small_model. `flag_name` is the user-facing flag
// (e.g., "--small_model_path"); `mt_upper` is upper-case model_type for the
// message body. No-op when path is non-empty or the bundle has no small model.
static void WarnIfMissingSmallModel(const std::string& path,
                                     const std::string& flag_name,
                                     const std::string& mt_upper,
                                     bool expects) {
  if (!path.empty() || !expects) return;
  LOG(WARNING)
      << flag_name << " is empty but model_type=" << mt_upper
      << " bundles a trained small_model in upstream Docker. "
      << "Without it every candidate goes through the big Inception-v3 "
      << "(slower) and FILTER classification may drift from Docker — "
      << "long-read modes can regress SNP F1 by several %. "
      << "Pass " << flag_name
      << "=<dir-with-layer_*_kernel.npy> (typically extracted by "
      << "tools/reference/extract_all_model_weights.sh).";
}

// LooksLikeSmallModelDir — cheap fs check: dir exists AND contains
// `layer_0_kernel.npy` (the file produced by extract_small_model_weights.sh
// for every supported small-model bundle). Matches the file the BNNS-CPU
// MLP loader will mmap at runtime.
static bool LooksLikeSmallModelDir(const std::string& dir) {
  if (dir.empty()) return false;
  struct stat st{};
  const std::string probe = absl::StrCat(dir, "/layer_0_kernel.npy");
  return ::stat(probe.c_str(), &st) == 0;
}

// AutoDiscoverGermlineSmallModel — given a `.dvw` checkpoint path, return the
// conventional sibling small-model dir if it exists, else "".
// Convention from tools/reference/extract_all_model_weights.sh:
//   <base>.dvw → <base>_small_weights/   (germline: WGS, ONT, PACBIO)
// `ckpt_path` may be empty or non-`.dvw` — in both cases we return "".
static std::string AutoDiscoverGermlineSmallModel(const std::string& ckpt_path) {
  if (ckpt_path.size() < 5) return "";
  const std::string suffix = ckpt_path.substr(ckpt_path.size() - 4);
  if (suffix != ".dvw") return "";
  const std::string base =
      ckpt_path.substr(0, ckpt_path.size() - 4);  // strip ".dvw"
  const std::string candidate = absl::StrCat(base, "_small_weights");
  if (LooksLikeSmallModelDir(candidate)) return candidate;
  return "";
}

// AutoDiscoverTrioOrSomaticSmallModel — given a `<dir>/<base>.dvw` checkpoint
// where `<base>` follows the trio/somatic naming convention, return the
// conventional sibling small-model dir if it exists, else "".
// Convention:
//   <dir>/deeptrio.<mode>_<role>.dvw → <dir>/deeptrio_<mode>_<role>_small/
//   <dir>/deepsomatic.<mode>.dvw     → <dir>/deepsomatic_<mode>_small/
// Mechanism: replace the FIRST `.` in <base> with `_`, then append `_small`.
// Returns "" if `ckpt_path` is empty, doesn't end in `.dvw`, has no `.` in
// the basename, or the candidate dir doesn't contain layer_0_kernel.npy.
static std::string AutoDiscoverTrioOrSomaticSmallModel(
    const std::string& ckpt_path) {
  if (ckpt_path.size() < 5) return "";
  if (ckpt_path.substr(ckpt_path.size() - 4) != ".dvw") return "";
  // Find the basename (start after last `/`).
  const auto slash = ckpt_path.find_last_of('/');
  const std::string parent =
      slash == std::string::npos ? "" : ckpt_path.substr(0, slash + 1);
  const std::string base = ckpt_path.substr(
      slash == std::string::npos ? 0 : slash + 1);
  // base is like "deeptrio.wgs_child.dvw" or "deepsomatic.wgs.dvw".
  // Strip ".dvw".
  const std::string base_noext = base.substr(0, base.size() - 4);
  // Replace FIRST `.` with `_`. If there's no `.`, this is not a
  // trio/somatic-style bundle and we return "".
  const auto dot = base_noext.find('.');
  if (dot == std::string::npos) return "";
  std::string flat = base_noext;
  flat[dot] = '_';
  const std::string candidate =
      absl::StrCat(parent, flat, "_small");
  if (LooksLikeSmallModelDir(candidate)) return candidate;
  return "";
}

// MaybeAutoDiscoverGermlineSmallModel — wraps AutoDiscoverGermlineSmallModel
// with the policy: only kicks in when (a) the user left the flag empty,
// (b) the bundle expects a small_model, (c) we have a checkpoint path to
// pivot off. Logs INFO when it finds a dir; the caller must still invoke
// WarnIfMissingSmallModel afterwards (with the possibly-updated path) so the
// "no small model" warning fires when discovery fails.
static void MaybeAutoDiscoverGermlineSmallModel(std::string& small_model_path,
                                                 const std::string& ckpt_path,
                                                 const std::string& flag_name,
                                                 bool expects) {
  if (!small_model_path.empty() || !expects) return;
  const std::string discovered = AutoDiscoverGermlineSmallModel(ckpt_path);
  if (discovered.empty()) return;
  LOG(INFO) << "Auto-discovered " << flag_name << "=" << discovered
            << " (sibling of --checkpoint=" << ckpt_path << ")";
  small_model_path = discovered;
}

// MaybeAutoDiscoverTrioOrSomaticSmallModel — same policy as above for the
// trio/somatic naming convention.
static void MaybeAutoDiscoverTrioOrSomaticSmallModel(
    std::string& small_model_path, const std::string& ckpt_path,
    const std::string& flag_name, bool expects) {
  if (!small_model_path.empty() || !expects) return;
  const std::string discovered =
      AutoDiscoverTrioOrSomaticSmallModel(ckpt_path);
  if (discovered.empty()) return;
  LOG(INFO) << "Auto-discovered " << flag_name << "=" << discovered
            << " (sibling of --checkpoint=" << ckpt_path << ")";
  small_model_path = discovered;
}

// EnsurePathExists — early existence check for user-supplied file/dir paths.
// Returns true (with no logging) when path is empty or `stat()` succeeds;
// returns false + LOG(ERROR) when the path is non-empty but doesn't exist.
//
// Intended use: validate --reads / --ref / --checkpoint at the top of each
// Run* dispatcher so a typo like --ref=/tmp/GRCh38.fa.bak fails in 1 ms with
// a clear "file not found" instead of failing minutes later inside Nucleus
// with "could not open SAM/FASTA reader" (cause obscured by the wrapper).
//
// Empty path is treated as "user didn't set it"; existing required-flag
// checks (LOG(ERROR) << "... required") handle that case separately, so this
// helper just no-ops on empty input.
static bool EnsurePathExists(const std::string& path,
                              const std::string& flag_name) {
  if (path.empty()) return true;
  struct stat st{};
  if (::stat(path.c_str(), &st) == 0) return true;
  LOG(ERROR) << flag_name << "=" << path
             << " not found on disk (check the path for typos).";
  return false;
}

// EnsureFastaIndexed — for --ref FASTA paths, confirm that an `.fai` sibling
// exists. Nucleus's IndexedFastaReader requires it; without one, the make_
// examples worker dies several seconds in with a generic open error. This
// catches the missing-index case in <1 ms with an actionable message
// pointing the user at `samtools faidx`.
static bool EnsureFastaIndexed(const std::string& fasta_path) {
  if (fasta_path.empty()) return true;
  const std::string fai = absl::StrCat(fasta_path, ".fai");
  struct stat st{};
  if (::stat(fai.c_str(), &st) == 0) return true;
  LOG(ERROR) << "--ref=" << fasta_path
             << " has no .fai index (expected at " << fai
             << "). Generate one with: samtools faidx " << fasta_path;
  return false;
}

// EnsureBamIndexed — for --reads BAM/CRAM paths, confirm that a sibling
// index exists (`.bai` for BAM, `.crai` for CRAM, in either samtools or
// Picard naming). Nucleus's SamReader needs the index for region queries;
// without one, the worker fails on the first `query()` call with a
// confusing "no index" error from htslib.
static bool EnsureBamIndexed(const std::string& bam_path,
                              const std::string& flag_name) {
  if (bam_path.empty()) return true;
  const auto exists = [](const std::string& p) {
    struct stat st{};
    return ::stat(p.c_str(), &st) == 0;
  };
  const std::string ext = bam_path.size() >= 4
      ? bam_path.substr(bam_path.size() - 4) : "";
  if (ext == ".bam") {
    if (exists(absl::StrCat(bam_path, ".bai"))) return true;
    if (exists(bam_path.substr(0, bam_path.size() - 4) + ".bai")) return true;
    LOG(ERROR) << flag_name << "=" << bam_path
               << " has no .bai index. Generate one with: "
               << "samtools index " << bam_path;
    return false;
  }
  if (bam_path.size() >= 5 &&
      bam_path.substr(bam_path.size() - 5) == ".cram") {
    if (exists(absl::StrCat(bam_path, ".crai"))) return true;
    if (exists(bam_path.substr(0, bam_path.size() - 5) + ".crai")) return true;
    LOG(ERROR) << flag_name << "=" << bam_path
               << " has no .crai index. Generate one with: "
               << "samtools index " << bam_path;
    return false;
  }
  // Other extensions (.sam, etc.) — skip the check; we can't enforce it.
  return true;
}

// ApplyModelFlags — appends make_examples flags from model example_info.json.
// Values mirror tools/conversion/models/<name>/model.example_info.json exactly.
static void ApplyModelFlags(const std::string& model_type,
                             std::vector<std::string>& me_args) {
  std::string mt = model_type;
  for (char& c : mt) c = static_cast<char>(std::toupper(c));

  if (mt == "PACBIO") {
    me_args.push_back("--pileup_image_width=147");
    me_args.push_back("--channel_list_preset=LONG_READ_PACBIO");
    me_args.push_back("--small_model_use_haplotypes=true");  // 106-feature model
    me_args.push_back("--min_mapping_quality=1");
    // min_base_quality intentionally NOT set for PacBio: Docker's
    // pacbio/model.example_info.json does not include this flag, so the
    // default (10) applies. ONT sets 1 explicitly; PacBio does not.
    me_args.push_back("--max_reads_per_partition=1500");
    me_args.push_back("--partition_size=25000");
    me_args.push_back("--sort_by_haplotypes=true");
    me_args.push_back("--trim_reads_for_pileup=true");
    me_args.push_back("--phase_reads=true");
    me_args.push_back("--parse_sam_aux_fields=true");
    me_args.push_back("--keep_supplementary_alignments=true");
    me_args.push_back("--realigner_enabled=false");
    // Phase 5.5d/14: enable DirectPhasing so the 106-feature haplotype
    // small_model gets DP's per-read phase output (matching upstream's
    // FeatureEncoder(haplotype, read_phases) input). Without this, our
    // small_model used BAM HP tags (whatshap haplotag), which diverge
    // from DirectPhasing at phase-block boundaries.
    me_args.push_back("--use_direct_phasing=true");
    me_args.push_back("--small_model_snp_gq_threshold=19");
    me_args.push_back("--small_model_indel_gq_threshold=22");
    me_args.push_back("--small_model_vaf_context_window_size=51");
    me_args.push_back("--vsc_min_fraction_indels=0.12");
    me_args.push_back("--vsc_min_indel_fraction_for_small_indels=0.12");
    me_args.push_back("--vsc_min_indel_fraction_for_large_indels=0.05");
    me_args.push_back("--vsc_small_indel_threshold=1");
  } else if (mt == "ONT") {
    me_args.push_back("--pileup_image_width=199");
    me_args.push_back("--channel_list_preset=LONG_READ_ONT");
    me_args.push_back("--small_model_use_haplotypes=true");  // 106-feature model
    me_args.push_back("--min_mapping_quality=1");
    me_args.push_back("--min_base_quality=1");
    me_args.push_back("--max_reads_per_partition=1500");
    me_args.push_back("--partition_size=25000");
    me_args.push_back("--sort_by_haplotypes=true");
    me_args.push_back("--trim_reads_for_pileup=true");
    me_args.push_back("--phase_reads=true");
    me_args.push_back("--parse_sam_aux_fields=true");
    me_args.push_back("--realigner_enabled=false");
    // Phase 5.5d/14: same as PACBIO — DP-fed read phases for small_model.
    me_args.push_back("--use_direct_phasing=true");
    me_args.push_back("--small_model_snp_gq_threshold=9");
    me_args.push_back("--small_model_indel_gq_threshold=17");
    me_args.push_back("--small_model_vaf_context_window_size=51");
    me_args.push_back("--vsc_min_fraction_snps=0.1");
    me_args.push_back("--vsc_min_fraction_indels=0.1");
  } else if (mt == "HYBRID_PACBIO_ILLUMINA" || mt == "HYBRID") {
    me_args.push_back("--channel_list_preset=BASE_CHANNELS");
    me_args.push_back("--trim_reads_for_pileup=true");
  } else if (mt == "MASSEQ") {
    me_args.push_back("--pileup_image_width=199");
    me_args.push_back("--channel_list_preset=MASSEQ");
    me_args.push_back("--min_mapping_quality=1");
    me_args.push_back("--max_reads_per_partition=0");
    me_args.push_back("--max_reads_for_dynamic_bases_per_region=1500");
    me_args.push_back("--partition_size=25000");
    me_args.push_back("--sort_by_haplotypes=true");
    me_args.push_back("--trim_reads_for_pileup=true");
    me_args.push_back("--phase_reads=true");
    me_args.push_back("--parse_sam_aux_fields=true");
    me_args.push_back("--realigner_enabled=false");
    me_args.push_back("--vsc_min_fraction_indels=0.12");
  } else if (mt == "RNASEQ") {
    me_args.push_back("--channel_list_preset=BASE_CHANNELS");
    me_args.push_back("--split_skip_reads=true");
    me_args.push_back("--min_mapping_quality=40");
    me_args.push_back("--max_reads_per_partition=0");
    me_args.push_back("--partition_size=10000");
  } else {
    // WGS / WES defaults.
    me_args.push_back("--realigner_enabled=true");
    // vaf_context_window=51: matches WGS/WES example_info.json.
    // Required so AlleleCounter fills all 51 VAF context positions in the
    // DeepVariantCall proto. Without this, EncodeSmallModelFeatures() reads
    // 0 for 46 of 51 positions → small model gets wrong features → GQ=20
    // borderline sites mispredicted → PASS↔NoCall FM at WG scale.
    me_args.push_back("--small_model_vaf_context_window_size=51");
  }
}

// PostprocessModelFlags — returns postprocess --flag=value args for model_type.
static std::vector<std::string> PostprocessModelFlags(
    const std::string& model_type) {
  std::vector<std::string> pp;
  std::string mt = model_type;
  for (char& c : mt) c = static_cast<char>(std::toupper(c));
  if (mt == "WES") pp.push_back("--multiallelic_mode=min");
  return pp;
}

// TrioInputDims — call_variants input shape from DeepTrio example_info.json.
struct TrioDims { int child_h; int parent_h; int channels; int width; };
static TrioDims TrioInputDims(const std::string& model_type) {
  std::string mt = model_type;
  for (char& c : mt) c = static_cast<char>(std::toupper(c));
  if (mt == "PACBIO") return {140, 140, 9, 199};
  if (mt == "ONT")    return {300, 300, 9, 199};
  if (mt == "WES")    return {300, 300, 7, 221};
  return                     {140, 140, 7, 221};  // WGS default
}

// GermlineInputDims — call_variants input shape for single-sample germline.
// Source: tools/conversion/models/<model>/model.example_info.json shape field.
// Height is always 100 for germline models.
struct GermlineDims { int channels; int width; };
static GermlineDims GermlineInputDims(const std::string& model_type) {
  std::string mt = model_type;
  for (char& c : mt) c = static_cast<char>(std::toupper(c));
  if (mt == "PACBIO")                             return {10, 147}; // base8+alt2
  if (mt == "ONT")                                return {10, 199}; // base8+alt2
  if (mt == "MASSEQ")                             return { 9, 199}; // base7+alt2
  if (mt == "HYBRID_PACBIO_ILLUMINA" ||
      mt == "HYBRID" || mt == "RNASEQ")           return { 6, 221}; // BASE_CHANNELS
  return                                                 { 7, 221}; // WGS/WES default
}

// SomaticInputDims — call_variants input shape per model_type × has_normal.
// Source: deepsomatic.<model>[_tumor_only]/model.example_info.json shape field.
struct SomaticDims { int h; int channels; int width; };
static SomaticDims SomaticInputDims(const std::string& model_type,
                                    bool has_normal) {
  std::string mt = model_type;
  for (char& c : mt) c = static_cast<char>(std::toupper(c));
  if (!has_normal) {
    // Tumor-only: h=100 for all types; channels = base+1 (allele_frequency).
    // PacBio/ONT tumor-only width=99 (narrower than TN PacBio 147).
    if (mt == "PACBIO" || mt == "ONT") return {100, 10,  99};
    return                                    {100,  8, 221};
  }
  // Tumor+normal shapes from model.example_info.json.
  if (mt == "PACBIO") return {200, 9, 147};
  if (mt == "ONT")    return {200, 9,  99};
  return                     {200, 7, 221};   // WGS/WES/FFPE_WGS/FFPE_WES
}

// SomaticModelPath — default model bundle path for somatic mode.
// Returns .mlpackage (CoreML/ane_speculate) from DEEPVARIANT_MODELS_DIR.
// Metal backend callers pass --checkpoint pointing to the .dvw in the same dir.
static std::string SomaticModelPath(const std::string& model_type,
                                    bool has_normal) {
  const char* env = std::getenv("DEEPVARIANT_MODELS_DIR");
  std::string base = env ? env
                         : "/opt/homebrew/share/deepvariant-models";
  std::string mt = model_type;
  for (char& c : mt) c = static_cast<char>(std::tolower(c));
  if (!has_normal) {
    return absl::StrCat(base, "/deepsomatic.", mt, "_tumor_only.mlpackage");
  }
  return absl::StrCat(base, "/deepsomatic.", mt, ".mlpackage");
}

// ApplySomaticModelFlags — somatic make_examples flags from
// deepsomatic.<model>[_tumor_only]/model.example_info.json flags_for_calling.
// Note: sort_by_alt_allele_support and track_ref_reads are set directly in
// make_examples_main.cc (conditioned on has_normal); not passed as flags here.
static void ApplySomaticModelFlags(const std::string& model_type,
                                    bool has_normal,
                                    std::vector<std::string>& me_args) {
  std::string mt = model_type;
  for (char& c : mt) c = static_cast<char>(std::toupper(c));

  if (!has_normal) {
    // ── Tumor-only flag dispatch ──────────────────────────────────────────
    // Mirrors deepsomatic.*_tumor_only/model.example_info.json flags_for_calling.
    // No small model for any tumor-only variant (no trained_small_model_path).
    // sort_by_alt_allele_support absent from all tumor-only JSONs → stays false
    // (handled in make_examples_main.cc).
    if (mt == "PACBIO") {
      me_args.push_back("--pileup_image_width=99");   // tumor-only width=99 not 147
      me_args.push_back("--channel_list_preset=MASSEQ");
      me_args.push_back("--alt_aligned_pileup=diff_channels");
      me_args.push_back("--sort_by_haplotypes=true");
      me_args.push_back("--phase_reads=true");
      me_args.push_back("--parse_sam_aux_fields=true");
      me_args.push_back("--trim_reads_for_pileup=true");
      me_args.push_back("--realigner_enabled=false");
      me_args.push_back("--min_mapping_quality=5");
      me_args.push_back("--partition_size=25000");
      me_args.push_back("--vsc_min_fraction_snps=0.02");
      me_args.push_back("--vsc_min_fraction_indels=0.1");
      me_args.push_back("--vsc_min_count_snps=1");
    } else if (mt == "ONT") {
      me_args.push_back("--pileup_image_width=99");
      me_args.push_back("--channel_list_preset=MASSEQ");
      me_args.push_back("--alt_aligned_pileup=diff_channels");
      me_args.push_back("--sort_by_haplotypes=true");
      me_args.push_back("--phase_reads=true");
      me_args.push_back("--parse_sam_aux_fields=true");
      me_args.push_back("--trim_reads_for_pileup=true");
      me_args.push_back("--realigner_enabled=false");
      me_args.push_back("--min_mapping_quality=5");
      me_args.push_back("--partition_size=25000");
      me_args.push_back("--vsc_min_fraction_snps=0.05");
      me_args.push_back("--vsc_min_fraction_indels=0.1");
    } else {
      // WGS/WES/FFPE_WGS/FFPE_WES tumor-only: shared vsc thresholds.
      me_args.push_back("--vsc_min_fraction_snps=0.05");
      me_args.push_back("--vsc_min_fraction_indels=0.07");
      // WGS_TO and WES_TO declare vsc_max_fraction=0.5 in their JSON;
      // FFPE_WGS_TO and FFPE_WES_TO do not. In tumor-only mode the
      // "non_target" sample doesn't exist, so this is effectively a no-op,
      // but set to match Docker's example_info.json exactly.
      if (mt == "WGS" || mt == "WES") {
        me_args.push_back("--vsc_max_fraction_snps_for_non_target_sample=0.5");
        me_args.push_back("--vsc_max_fraction_indels_for_non_target_sample=0.5");
      }
    }
    return;
  }

  // ── Tumor+normal flag dispatch ────────────────────────────────────────────
  if (mt == "PACBIO") {
    me_args.push_back("--pileup_image_width=147");
    me_args.push_back("--channel_list_preset=MASSEQ");
    me_args.push_back("--alt_aligned_pileup=diff_channels");
    me_args.push_back("--sort_by_haplotypes=true");
    me_args.push_back("--phase_reads=true");
    me_args.push_back("--parse_sam_aux_fields=true");
    me_args.push_back("--trim_reads_for_pileup=true");
    me_args.push_back("--realigner_enabled=false");
    me_args.push_back("--min_mapping_quality=5");
    me_args.push_back("--partition_size=25000");
    me_args.push_back("--vsc_min_fraction_snps=0.02");
    me_args.push_back("--vsc_min_fraction_indels=0.1");
    me_args.push_back("--vsc_min_count_snps=1");
    me_args.push_back("--small_model_snp_gq_threshold=60");
    me_args.push_back("--small_model_indel_gq_threshold=57");
    me_args.push_back("--small_model_vaf_context_window_size=51");
    me_args.push_back("--vsc_max_fraction_snps_for_non_target_sample=0.5");
    me_args.push_back("--vsc_max_fraction_indels_for_non_target_sample=0.5");
  } else if (mt == "ONT") {
    me_args.push_back("--pileup_image_width=99");
    me_args.push_back("--channel_list_preset=MASSEQ");
    me_args.push_back("--alt_aligned_pileup=diff_channels");
    me_args.push_back("--sort_by_haplotypes=true");
    me_args.push_back("--phase_reads=true");
    me_args.push_back("--parse_sam_aux_fields=true");
    me_args.push_back("--trim_reads_for_pileup=true");
    me_args.push_back("--realigner_enabled=false");
    me_args.push_back("--min_mapping_quality=5");
    me_args.push_back("--partition_size=25000");
    me_args.push_back("--vsc_min_fraction_snps=0.05");
    me_args.push_back("--vsc_min_fraction_indels=0.1");
    me_args.push_back("--small_model_snp_gq_threshold=51");
    me_args.push_back("--small_model_indel_gq_threshold=56");
    me_args.push_back("--small_model_vaf_context_window_size=51");
    // ONT uses 0.6 (not 0.5 like PacBio/WGS) per deepsomatic/ont/model.example_info.json
    me_args.push_back("--vsc_max_fraction_snps_for_non_target_sample=0.6");
    me_args.push_back("--vsc_max_fraction_indels_for_non_target_sample=0.6");
  } else if (mt == "FFPE_WGS") {
    // FFPE_WGS TN: sort_by_alt_allele_support=true + small model (in JSON).
    // No vsc_max_fraction_for_non_target_sample (NOT in FFPE_WGS JSON).
    me_args.push_back("--sort_by_alt_allele_support_somatic=true");
    me_args.push_back("--vsc_min_fraction_snps=0.029");
    me_args.push_back("--vsc_min_fraction_indels=0.05");
    me_args.push_back("--small_model_snp_gq_threshold=53");
    me_args.push_back("--small_model_indel_gq_threshold=36");
    me_args.push_back("--small_model_vaf_context_window_size=51");
  } else if (mt == "FFPE_WES") {
    // FFPE_WES TN: no sort_by_alt_allele_support, no small model,
    // no vsc_max_fraction (none declared in FFPE_WES JSON).
    me_args.push_back("--vsc_min_fraction_snps=0.029");
    me_args.push_back("--vsc_min_fraction_indels=0.05");
  } else if (mt == "WES") {
    // WES tumor+normal: vsc_max_fraction=0.5 declared; no sort_by_alt_allele,
    // no small model (WES JSON has no trained_small_model_path).
    me_args.push_back("--vsc_min_fraction_snps=0.029");
    me_args.push_back("--vsc_min_fraction_indels=0.05");
    me_args.push_back("--vsc_max_fraction_snps_for_non_target_sample=0.5");
    me_args.push_back("--vsc_max_fraction_indels_for_non_target_sample=0.5");
  } else {
    // WGS default somatic tumor+normal: sort_by_alt_allele_support=true +
    // small model GQ thresholds + vsc_max_fraction=0.5 (all in WGS JSON).
    me_args.push_back("--sort_by_alt_allele_support_somatic=true");
    me_args.push_back("--vsc_min_fraction_snps=0.029");
    me_args.push_back("--vsc_min_fraction_indels=0.05");
    me_args.push_back("--small_model_snp_gq_threshold=31");
    me_args.push_back("--small_model_indel_gq_threshold=29");
    me_args.push_back("--small_model_vaf_context_window_size=51");
    me_args.push_back("--vsc_max_fraction_snps_for_non_target_sample=0.5");
    me_args.push_back("--vsc_max_fraction_indels_for_non_target_sample=0.5");
  }
}

int RunAll(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  // Trio mode: --reads_parent1 set → dispatch the 3-sample pipeline.
  if (!absl::GetFlag(FLAGS_reads_parent1).empty()) {
    return RunAllTrio(argc, argv);
  }

  // Somatic mode: --reads_tumor set → dispatch the 2-sample (tumor+
  // normal) or 1-sample (tumor_only) pipeline. Single tumor VCF output.
  if (!absl::GetFlag(FLAGS_reads_tumor).empty()) {
    return RunAllSomatic(argc, argv);
  }

  // Pangenome-aware mode: --reads_pangenome set → 2-sample pipeline
  // (pangenome=0, reads=1=main). Single VCF output for the reads sample.
  if (!absl::GetFlag(FLAGS_reads_pangenome).empty()) {
    return RunAllPangenome(argc, argv);
  }

  const std::string model_type =
      CanonicalModelType(absl::GetFlag(FLAGS_model_type));
  const std::string reads_flag = absl::GetFlag(FLAGS_reads);
  const std::string ref_flag = absl::GetFlag(FLAGS_ref);
  const std::string output_vcf_flag = absl::GetFlag(FLAGS_output_vcf);
  const std::string user_regions = absl::GetFlag(FLAGS_regions);
  const std::string regions_flag = EffectiveRegions(user_regions, ref_flag);
  const std::string tmp_dir = absl::GetFlag(FLAGS_intermediate_results_dir);
  const int num_shards = EffectiveNumShards();

  // Ensure the intermediate results dir exists before any stage writes into it
  // (the trio/somatic/pangenome dispatchers already create it; germline did
  // not, so a non-existent --intermediate_results_dir failed mid-pipeline when
  // make_examples tried to open its TFRecord output).
  {
    std::error_code ec;
    std::filesystem::create_directories(tmp_dir, ec);
    if (ec) {
      LOG(ERROR) << "Failed to create intermediate dir " << tmp_dir << ": "
                 << ec.message();
      return 1;
    }
  }

  if (reads_flag.empty() || ref_flag.empty() || output_vcf_flag.empty()) {
    LOG(ERROR) << "Usage: deepvariant run --reads=<BAM> --ref=<FASTA> "
                  "--output_vcf=<VCF> [--model_type=WGS] [--regions=chr20]";
    return 1;
  }
  // Early-fail: catch typos in --reads / --ref / --checkpoint in <1 ms
  // instead of letting them surface minutes later as a Nucleus open error.
  // --output_vcf is intentionally NOT checked: it's an OUTPUT path that
  // postprocess will create, so its absence is expected and correct.
  if (!EnsurePathExists(reads_flag, "--reads") ||
      !EnsureBamIndexed (reads_flag, "--reads") ||
      !EnsurePathExists(ref_flag,   "--ref")   ||
      !EnsureFastaIndexed(ref_flag) ||
      !EnsurePathExists(absl::GetFlag(FLAGS_checkpoint), "--checkpoint")) {
    return 1;
  }

  // Single-process pipeline: one make_examples call (using --threads=N for
  // intra-process parallelism, writing sharded `name-NNNNN-of-NNNNN`
  // files), one call_variants call (sharded examples in, single cvo out),
  // one postprocess.
  const int n_threads = std::max(1, num_shards);
  const std::string examples_base =
      absl::StrCat(tmp_dir, "/examples.tfrecord");
  const std::string examples_pattern =
      n_threads > 1 ? absl::StrCat(examples_base, "@", n_threads)
                    : examples_base;
  const std::string cvo_pattern    = absl::StrCat(tmp_dir, "/cvo.tfrecord");
  const std::string small_cvo_base = absl::StrCat(tmp_dir, "/small_cvo.tfrecord");
  const std::string small_cvo_path =
      n_threads > 1 ? absl::StrCat(small_cvo_base, "@", n_threads)
                    : small_cvo_base;
  const std::string merged_cvo_path =
      absl::StrCat(tmp_dir, "/merged_cvo.tfrecord");
  // Phase 9 / Step 3 — gVCF intermediate non-variant TFRecord, sharded
  // per make_examples worker thread. Postprocess merges it with the
  // variant CVO stream via nucleus::MergeAndWriteVariantsAndNonVariants.
  const std::string gvcf_tfrecord_base =
      absl::StrCat(tmp_dir, "/gvcf.tfrecord");
  const std::string gvcf_tfrecord_path =
      n_threads > 1 ? absl::StrCat(gvcf_tfrecord_base, "@", n_threads)
                    : gvcf_tfrecord_base;
  const std::string gvcf_outfile = absl::GetFlag(FLAGS_output_gvcf);

  // For --inference_backend=metal, the user passes a `.dvw` weight bundle
  // via --checkpoint; for coreml (Phase 2), the bundle is a `.mlpackage`
  // resolved by --model or the default ModelPath(model_type) lookup.
  const std::string inference_backend =
      absl::GetFlag(FLAGS_inference_backend);
  const std::string user_checkpoint = absl::GetFlag(FLAGS_checkpoint);
  const std::string user_model = absl::GetFlag(FLAGS_model);
  std::string model_path;
  if (!user_checkpoint.empty()) {
    model_path = user_checkpoint;
  } else if (!user_model.empty()) {
    model_path = user_model;
  } else {
    model_path = ModelPath(model_type);
  }
  std::string small_model_path = absl::GetFlag(FLAGS_small_model_path);
  {
    std::string mt = model_type;
    for (char& c : mt) c = static_cast<char>(std::toupper(c));
    const bool expects = GermlineExpectsSmallModel(mt);
    MaybeAutoDiscoverGermlineSmallModel(small_model_path, model_path,
                                        "--small_model_path", expects);
    WarnIfMissingSmallModel(small_model_path, "--small_model_path", mt,
                            expects);
  }

  // ── Stage 1: make_examples (single in-process call, --threads=N) ─────────
  // Internally fans out N worker threads (each with its own
  // SamReader / IndexedFastaReader / ExamplesGenerator / SmallModel) writing
  // sharded `name-NNNNN-of-NNNNN` files. Downstream stages read them
  // directly via TFRecordReader's `@N` shard expansion — no end-of-stage
  // concat. One process = ~N×100 % CPU under `top`, mirroring
  // salmon/samtools.
  LOG(INFO) << "Stage 1: make_examples (--threads=" << n_threads
            << ", in-process)";
  {
    std::vector<std::string> me_args = {
        absl::StrCat("--reads=", reads_flag),
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--examples=", examples_pattern),
        absl::StrCat("--threads=", n_threads),
        "--task_id=0",
        "--num_shards=1",
    };
    if (!regions_flag.empty()) {
      me_args.push_back(absl::StrCat("--regions=", regions_flag));
    }
    if (!small_model_path.empty()) {
      me_args.push_back(absl::StrCat("--small_model=", small_model_path));
      me_args.push_back(absl::StrCat("--small_model_cvo_outfile=",
                                      small_cvo_path));
    }
    if (!gvcf_outfile.empty()) {
      me_args.push_back(absl::StrCat("--gvcf=", gvcf_tfrecord_path));
    }
    // Per-model flags from example_info.json (pileup width, channels,
    // thresholds, realigner, sort_by_haplotypes, etc.).
    ApplyModelFlags(model_type, me_args);
    // Alt-aligned pileup: user override or per-model default.
    {
      const std::string user_aap = absl::GetFlag(FLAGS_alt_aligned_pileup);
      std::string aap = user_aap;
      if (aap.empty()) {
        const std::string mt_up = [&] {
          std::string s = model_type;
          for (char& c : s) c = static_cast<char>(std::toupper(c));
          return s;
        }();
        if (mt_up == "PACBIO" || mt_up == "ONT" || mt_up == "MASSEQ") {
          aap = "diff_channels";
        } else {
          aap = "none";
        }
      }
      me_args.push_back(absl::StrCat("--alt_aligned_pileup=", aap));
    }
    // Phase 9 / Step 4c — forward --use_direct_phasing to make_examples.
    // When true, big-model candidates get is_phased + PS info field;
    // default false → byte-identical baseline.
    if (absl::GetFlag(FLAGS_use_direct_phasing)) {
      me_args.push_back("--use_direct_phasing=true");
    }
    AppendSelectVariantTypes(me_args);
    AppendHaploidFlags(me_args);  // haploid gVCF reference confidence
    auto argv_me = MakeArgv("deepvariant_make_examples", me_args);
    int n = static_cast<int>(argv_me.size()) - 1;
    if (int rc = RunMakeExamples(n, argv_me.data()); rc != 0) {
      LOG(ERROR) << "make_examples failed";
      return rc;
    }
  }

  // ── Stage 2: call_variants ────────────────────────────────────────────────
  LOG(INFO) << "Stage 2: call_variants";
  {
    const GermlineDims gdims = GermlineInputDims(model_type);
    std::vector<std::string> cv_args = {
        absl::StrCat("--examples=", examples_pattern),
        absl::StrCat("--outfile=", cvo_pattern),
        absl::StrCat("--checkpoint=", model_path),
        absl::StrCat("--batch_size=", EffectiveBatchSize()),
        absl::StrCat("--inference_backend=", inference_backend),
        absl::StrCat("--input_channels=", gdims.channels),
        absl::StrCat("--input_width=", gdims.width),
    };
    AppendAneSpeculateArgs(cv_args, inference_backend,
                           absl::GetFlag(FLAGS_ane_speculate_metal_checkpoint));
    auto argv_cv = MakeArgv("deepvariant_call_variants", cv_args);
    int n = static_cast<int>(argv_cv.size()) - 1;
    if (int rc = RunCallVariants(n, argv_cv.data()); rc != 0) {
      LOG(ERROR) << "call_variants failed";
      return rc;
    }
  }

  // ── Stage 2.5: merge small_cvo + big_cvo into a single file ──────────────
  // postprocess takes one --infile so we concatenate the (already valid)
  // TFRecord files. TFRecord allows naive byte copy since each record is
  // self-delimiting.
  // small_cvo may be a `name@N` shard spec (one file per make_examples
  // worker thread); expand and concat each shard.
  std::string postprocess_input = cvo_pattern;
  if (!small_model_path.empty()) {
    LOG(INFO) << "Stage 2.5: merge small_cvo + big_cvo → " << merged_cvo_path;
    if (!MergeCvoFiles(small_cvo_path, cvo_pattern, merged_cvo_path)) {
      return 1;
    }
    postprocess_input = merged_cvo_path;
  }

  // ── Stage 3: postprocess_variants ────────────────────────────────────────
  LOG(INFO) << "Stage 3: postprocess_variants";
  {
    std::vector<std::string> pp_args = {
        absl::StrCat("--infile=", postprocess_input),
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--output_vcf_outfile=", output_vcf_flag),
    };
    if (!gvcf_outfile.empty()) {
      pp_args.push_back(absl::StrCat("--gvcf_outfile=", gvcf_outfile));
      pp_args.push_back(absl::StrCat("--nonvariant_site_tfrecord_path=",
                                      gvcf_tfrecord_path));
    }
    // Per-model postprocess flags (e.g. WES multiallelic_mode=min).
    for (const auto& f : PostprocessModelFlags(model_type)) pp_args.push_back(f);
    AppendHaploidFlags(pp_args);
    auto argv_pp = MakeArgv("deepvariant_postprocess", pp_args);
    int n = static_cast<int>(argv_pp.size()) - 1;
    if (int rc = RunPostprocessVariants(n, argv_pp.data()); rc != 0) {
      LOG(ERROR) << "postprocess_variants failed";
      return rc;
    }
  }

  LOG(INFO) << "Done. VCF: " << output_vcf_flag;
  return 0;
}

// ──────────────────────────────────────────────────────────────────────
// Trio dispatch: one make_examples (3 sample streams), 3× call_variants
// (child + parent1 + parent2 with the appropriate child/parent model),
// 3× postprocess (one VCF per sample). Mirrors the upstream
// run_deeptrio.py command sequence at deeptrio-quick-start.md.
// ──────────────────────────────────────────────────────────────────────
int RunAllTrio(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  const std::string model_type =
      CanonicalModelType(absl::GetFlag(FLAGS_model_type));
  // Ensure intermediate_results_dir exists (may not be pre-created by caller).
  { std::system(absl::StrCat("mkdir -p '",
      absl::GetFlag(FLAGS_intermediate_results_dir), "'").c_str()); }
  const std::string ref_flag = absl::GetFlag(FLAGS_ref);
  const std::string user_regions = absl::GetFlag(FLAGS_regions);
  const std::string regions_flag = EffectiveRegions(user_regions, ref_flag);
  const std::string tmp_dir = absl::GetFlag(FLAGS_intermediate_results_dir);
  const int num_shards = EffectiveNumShards();
  const int n_threads = std::max(1, num_shards);

  const std::string reads_child   = absl::GetFlag(FLAGS_reads);
  const std::string reads_parent1 = absl::GetFlag(FLAGS_reads_parent1);
  const std::string reads_parent2 = absl::GetFlag(FLAGS_reads_parent2);
  if (reads_child.empty() || reads_parent1.empty() || reads_parent2.empty()) {
    LOG(ERROR) << "Trio: requires --reads (child), --reads_parent1, --reads_parent2";
    return 1;
  }

  const std::string out_child   = absl::GetFlag(FLAGS_output_vcf_child);
  const std::string out_parent1 = absl::GetFlag(FLAGS_output_vcf_parent1);
  const std::string out_parent2 = absl::GetFlag(FLAGS_output_vcf_parent2);
  if (out_child.empty() || out_parent1.empty() || out_parent2.empty()) {
    LOG(ERROR) << "Trio: requires --output_vcf_child, --output_vcf_parent1, "
                  "--output_vcf_parent2";
    return 1;
  }

  // Resolve per-role model checkpoints. Allow either explicit
  // --checkpoint_child / --checkpoint_parent OR fall back to legacy
  // --checkpoint (used for both — useful for smoke tests with one model).
  std::string ckpt_child  = absl::GetFlag(FLAGS_checkpoint_child);
  std::string ckpt_parent = absl::GetFlag(FLAGS_checkpoint_parent);
  if (ckpt_child.empty())  ckpt_child  = absl::GetFlag(FLAGS_checkpoint);
  if (ckpt_parent.empty()) ckpt_parent = absl::GetFlag(FLAGS_checkpoint);
  if (ckpt_child.empty() || ckpt_parent.empty()) {
    LOG(ERROR) << "Trio: requires --checkpoint_child + --checkpoint_parent "
                  "(or --checkpoint as a shared fallback)";
    return 1;
  }
  // Early-fail: catch typos in 3× --reads_* + --ref + 2× checkpoints.
  if (!EnsurePathExists(reads_child,   "--reads")         ||
      !EnsureBamIndexed (reads_child,  "--reads")         ||
      !EnsurePathExists(reads_parent1, "--reads_parent1") ||
      !EnsureBamIndexed (reads_parent1,"--reads_parent1") ||
      !EnsurePathExists(reads_parent2, "--reads_parent2") ||
      !EnsureBamIndexed (reads_parent2,"--reads_parent2") ||
      !EnsurePathExists(ref_flag,      "--ref")           ||
      !EnsureFastaIndexed(ref_flag) ||
      !EnsurePathExists(ckpt_child,    "--checkpoint_child")  ||
      !EnsurePathExists(ckpt_parent,   "--checkpoint_parent")) {
    return 1;
  }
  std::string sm_child  = absl::GetFlag(FLAGS_small_model_path_child);
  std::string sm_parent = absl::GetFlag(FLAGS_small_model_path_parent);
  {
    std::string mt = model_type;
    for (char& c : mt) c = static_cast<char>(std::toupper(c));
    // Trio bundles use the same per-mode small_model presence as germline.
    const bool expects = GermlineExpectsSmallModel(mt);
    MaybeAutoDiscoverTrioOrSomaticSmallModel(
        sm_child, ckpt_child, "--small_model_path_child", expects);
    MaybeAutoDiscoverTrioOrSomaticSmallModel(
        sm_parent, ckpt_parent, "--small_model_path_parent", expects);
    WarnIfMissingSmallModel(sm_child,  "--small_model_path_child",  mt, expects);
    WarnIfMissingSmallModel(sm_parent, "--small_model_path_parent", mt, expects);
  }

  const std::string inference_backend =
      absl::GetFlag(FLAGS_inference_backend);

  // Per-sample intermediate paths.
  struct PerSamplePaths {
    std::string role;
    std::string examples_pattern;
    std::string small_cvo_pattern;
    std::string cvo_path;
    std::string merged_cvo_path;
    std::string sm_path;        // small_model weights dir (or empty)
    std::string ckpt_path;      // big-model checkpoint
    std::string output_vcf;
    std::string output_gvcf;
  };
  // Per-role .dvw rerun bundle for ane_speculate. Child uses its own
  // model; both parents share the parent .dvw.
  const std::string ane_dvw_child =
      absl::GetFlag(FLAGS_ane_speculate_metal_checkpoint_child);
  const std::string ane_dvw_parent =
      absl::GetFlag(FLAGS_ane_speculate_metal_checkpoint_parent);

  std::array<PerSamplePaths, 3> P;
  P[0].role = "child";
  P[1].role = "parent1";
  P[2].role = "parent2";
  for (auto& p : P) {
    p.examples_pattern  = absl::StrCat(tmp_dir, "/examples_", p.role,
                                        ".tfrecord");
    if (n_threads > 1) {
      p.examples_pattern = absl::StrCat(p.examples_pattern, "@", n_threads);
    }
    p.small_cvo_pattern = absl::StrCat(tmp_dir, "/small_cvo_", p.role,
                                        ".tfrecord");
    if (n_threads > 1) {
      p.small_cvo_pattern =
          absl::StrCat(p.small_cvo_pattern, "@", n_threads);
    }
    p.cvo_path        = absl::StrCat(tmp_dir, "/cvo_", p.role, ".tfrecord");
    p.merged_cvo_path =
        absl::StrCat(tmp_dir, "/merged_cvo_", p.role, ".tfrecord");
  }
  P[0].sm_path = sm_child;   P[0].ckpt_path = ckpt_child;
  P[1].sm_path = sm_parent;  P[1].ckpt_path = ckpt_parent;
  P[2].sm_path = sm_parent;  P[2].ckpt_path = ckpt_parent;
  // Per-role ane_speculate GPU rerun bundle.
  std::array<std::string, 3> ane_dvw{ane_dvw_child, ane_dvw_parent,
                                     ane_dvw_parent};
  P[0].output_vcf = out_child;   P[0].output_gvcf = absl::GetFlag(FLAGS_output_gvcf_child);
  P[1].output_vcf = out_parent1; P[1].output_gvcf = absl::GetFlag(FLAGS_output_gvcf_parent1);
  P[2].output_vcf = out_parent2; P[2].output_gvcf = absl::GetFlag(FLAGS_output_gvcf_parent2);

  // ── Stage 1: ONE make_examples invocation produces 3 example streams.
  LOG(INFO) << "Trio Stage 1: make_examples (3-sample, --threads=" << n_threads
            << ")";
  {
    std::vector<std::string> me_args = {
        absl::StrCat("--reads=", reads_child),
        absl::StrCat("--reads_parent1=", reads_parent1),
        absl::StrCat("--reads_parent2=", reads_parent2),
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--examples_child=",   P[0].examples_pattern),
        absl::StrCat("--examples_parent1=", P[1].examples_pattern),
        absl::StrCat("--examples_parent2=", P[2].examples_pattern),
        absl::StrCat("--threads=", n_threads),
        "--task_id=0",
        "--num_shards=1",
        // Step 1.3-bis: realigner runs per-sample in the trio worker
        // (mirrors upstream's realign_reads_per_sample_multisample).
        // Closes the candidate-count gap with Docker on indel-rich
        // regions where misalignment otherwise inflates AlleleCounter
        // counts with phantom alleles.
        "--realigner_enabled=true",
    };
    if (!regions_flag.empty()) {
      me_args.push_back(absl::StrCat("--regions=", regions_flag));
    }
    if (!absl::GetFlag(FLAGS_sample_name_parent1).empty()) {
      me_args.push_back(absl::StrCat("--sample_name_parent1=",
                                      absl::GetFlag(FLAGS_sample_name_parent1)));
    }
    if (!absl::GetFlag(FLAGS_sample_name_parent2).empty()) {
      me_args.push_back(absl::StrCat("--sample_name_parent2=",
                                      absl::GetFlag(FLAGS_sample_name_parent2)));
    }
    if (!sm_child.empty()) {
      me_args.push_back(absl::StrCat("--small_model_path_child=", sm_child));
      me_args.push_back(absl::StrCat("--small_model_cvo_outfile_child=",
                                      P[0].small_cvo_pattern));
    }
    if (!sm_parent.empty()) {
      me_args.push_back(absl::StrCat("--small_model_path_parent=", sm_parent));
      me_args.push_back(absl::StrCat("--small_model_cvo_outfile_parent1=",
                                      P[1].small_cvo_pattern));
      me_args.push_back(absl::StrCat("--small_model_cvo_outfile_parent2=",
                                      P[2].small_cvo_pattern));
    }
    // Per-model pileup/read flags from example_info.json.
    ApplyModelFlags(model_type, me_args);
    // Trio vaf_context_window override: run_deeptrio.py never sets
    // small_model_vaf_context_window_size, so ALL trio models use the default
    // (5). ApplyModelFlags(WGS/PacBio/ONT) sets 51 (germline default), which
    // differs from run_deeptrio.py. Restore the default here (Abseil last-wins).
    me_args.push_back("--small_model_vaf_context_window_size=5");
    // Per-model pileup heights for trio (per-sample, stacked in make_examples).
    // WGS/PacBio: child=60 parent=40 → total 140. WES/ONT: child=100 parent=100
    // → total 300. make_examples defaults to 60/40 (WGS); override for others.
    {
      std::string mt = model_type;
      for (char& c : mt) c = static_cast<char>(std::toupper(c));
      if (mt == "WES" || mt == "ONT") {
        me_args.push_back("--pileup_image_height_child=100");
        me_args.push_back("--pileup_image_height_parent=100");
      }
      // WGS / PacBio: defaults 60/40 in make_examples_main.cc are correct.
    }
    // DeepTrio PacBio/ONT channel + width overrides.
    // ApplyModelFlags(PACBIO) sets LONG_READ_PACBIO (8ch, width=147) and
    // ApplyModelFlags(ONT) sets LONG_READ_ONT (8ch, width=199).
    // But DeepTrio PacBio/ONT models use MASSEQ preset (7ch) + alt-aligned (9)
    // with width=199.  Push overrides AFTER ApplyModelFlags; Abseil last-wins.
    {
      std::string mt = model_type;
      for (char& c : mt) c = static_cast<char>(std::toupper(c));
      if (mt == "PACBIO" || mt == "ONT") {
        // DeepTrio PacBio/ONT: channel/width overrides + trio-specific flags
        // from run_deeptrio.py (different from germline ApplyModelFlags values).
        // Abseil last-wins: these override ApplyModelFlags(PACBIO/ONT) values.
        me_args.push_back("--pileup_image_width=199");
        me_args.push_back("--channel_list_preset=MASSEQ");
        me_args.push_back("--alt_aligned_pileup=diff_channels");
        // Trio uses max_reads_for_dynamic_bases_per_region=200, not 1500
        // (run_deeptrio.py:682, 705 — germline uses 1500 via MASSEQ block).
        me_args.push_back("--max_reads_for_dynamic_bases_per_region=200");
        // discard_non_dna_regions: matches run_deeptrio.py:682,705. Flag now
        // declared in make_examples_main.cc; runtime N-region filter is a
        // future enhancement but the flag must be set for parity.
        me_args.push_back("--discard_non_dna_regions=true");
      }
      if (mt == "ONT") {
        // ONT trio overrides vs germline ONT:
        //   min_mapping_quality=5 (germline uses 1)
        //   max_reads_per_partition=500 (germline uses 1500)
        //   vsc_min_fraction_indels=0.12 (germline uses 0.1)
        // Source: run_deeptrio.py:688-706
        me_args.push_back("--min_mapping_quality=5");
        me_args.push_back("--max_reads_per_partition=500");
        me_args.push_back("--vsc_min_fraction_indels=0.12");
      }
    }
    // DeepTrio threshold overrides (upstream scripts/run_deeptrio.py).
    // WGS trio uses SNP_GQ=15 / INDEL_GQ=29; long-read models use
    // the thresholds from ApplyModelFlags() already.
    {
      std::string mt = model_type;
      for (char& c : mt) c = static_cast<char>(std::toupper(c));
      if (mt == "WGS" || mt == "WES") {
        me_args.push_back("--small_model_snp_gq_threshold=15");
        me_args.push_back("--small_model_indel_gq_threshold=29");
        // WGS trio uses realigner (different from single-sample WGS default).
        me_args.push_back("--realigner_enabled=true");
      }
    }
    // Phase 9 / Step 4c — forward --use_direct_phasing for trio path.
    if (absl::GetFlag(FLAGS_use_direct_phasing)) {
      me_args.push_back("--use_direct_phasing=true");
    }
    AppendSelectVariantTypes(me_args);
    AppendHaploidFlags(me_args);  // haploid gVCF reference confidence
    auto argv_me = MakeArgv("deepvariant_make_examples", me_args);
    int n = static_cast<int>(argv_me.size()) - 1;
    if (int rc = RunMakeExamples(n, argv_me.data()); rc != 0) {
      LOG(ERROR) << "Trio: make_examples failed";
      return rc;
    }
  }

  // ── Stage 2-3 per sample: call_variants → merge → postprocess.
  for (size_t pi = 0; pi < P.size(); ++pi) {
    auto& p = P[pi];
    LOG(INFO) << "Trio Stage 2 (" << p.role << "): call_variants";
    {
      // Per-mode pileup shape from DeepTrio example_info.json.
      const TrioDims tdims = TrioInputDims(model_type);
      const int input_h = (p.role == "child") ? tdims.child_h : tdims.parent_h;
      std::vector<std::string> cv_args = {
          absl::StrCat("--examples=", p.examples_pattern),
          absl::StrCat("--outfile=", p.cvo_path),
          absl::StrCat("--checkpoint=", p.ckpt_path),
          absl::StrCat("--batch_size=", EffectiveBatchSize()),
          absl::StrCat("--inference_backend=", inference_backend),
          absl::StrCat("--input_height=", input_h),
          absl::StrCat("--input_channels=", tdims.channels),
          absl::StrCat("--input_width=", tdims.width),
      };
      AppendAneSpeculateArgs(cv_args, inference_backend, ane_dvw[pi]);
      auto argv_cv = MakeArgv("deepvariant_call_variants", cv_args);
      int n = static_cast<int>(argv_cv.size()) - 1;
      if (int rc = RunCallVariants(n, argv_cv.data()); rc != 0) {
        LOG(ERROR) << "Trio: call_variants failed for " << p.role;
        return rc;
      }
    }

    // Stage 2.5: merge small_cvo + big cvo (per sample).
    std::string postprocess_input = p.cvo_path;
    if (!p.sm_path.empty()) {
      LOG(INFO) << "Trio Stage 2.5 (" << p.role << "): merge → "
                << p.merged_cvo_path;
      std::ofstream out(p.merged_cvo_path,
                         std::ios::binary | std::ios::trunc);
      if (!out) {
        LOG(ERROR) << "Cannot open merged CVO: " << p.merged_cvo_path;
        return 1;
      }
      auto append_path = [&](const std::string& path) {
        std::ifstream in(path, std::ios::binary);
        if (!in) return;
        // Same empty-streambuf failbit gotcha — see RunAllGermline note.
        std::vector<char> buf((std::istreambuf_iterator<char>(in)),
                               std::istreambuf_iterator<char>());
        if (!buf.empty()) out.write(buf.data(), buf.size());
      };
      auto at = p.small_cvo_pattern.find('@');
      if (at == std::string::npos) {
        append_path(p.small_cvo_pattern);
      } else {
        const std::string prefix = p.small_cvo_pattern.substr(0, at);
        int nshard = 0;
        if (!absl::SimpleAtoi(p.small_cvo_pattern.substr(at + 1), &nshard) ||
            nshard <= 0) {
          LOG(ERROR) << "Bad small_cvo shard spec: " << p.small_cvo_pattern;
          return 1;
        }
        for (int i = 0; i < nshard; ++i) {
          append_path(absl::StrCat(prefix, "-",
                                    absl::Dec(i, absl::kZeroPad5),
                                    "-of-",
                                    absl::Dec(nshard, absl::kZeroPad5)));
        }
      }
      append_path(p.cvo_path);
      out.close();
      postprocess_input = p.merged_cvo_path;
    }

    LOG(INFO) << "Trio Stage 3 (" << p.role << "): postprocess_variants";
    std::vector<std::string> pp_args = {
        absl::StrCat("--infile=", postprocess_input),
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--output_vcf_outfile=", p.output_vcf),
    };
    if (!p.output_gvcf.empty()) {
      pp_args.push_back(absl::StrCat("--gvcf_outfile=", p.output_gvcf));
    }
    AppendHaploidFlags(pp_args);
    auto argv_pp = MakeArgv("deepvariant_postprocess", pp_args);
    int n = static_cast<int>(argv_pp.size()) - 1;
    if (int rc = RunPostprocessVariants(n, argv_pp.data()); rc != 0) {
      LOG(ERROR) << "Trio: postprocess failed for " << p.role;
      return rc;
    }
    LOG(INFO) << "Trio: " << p.role << " VCF: " << p.output_vcf;
  }

  LOG(INFO) << "Trio: done. 3 VCFs at " << out_child << ", " << out_parent1
            << ", " << out_parent2;
  return 0;
}

// ──────────────────────────────────────────────────────────────────────
// DeepSomatic dispatch: one make_examples (tumor + optional normal),
// 1× call_variants on the tumor model only (normal has skip_output=true),
// 1× postprocess writing a single tumor VCF. Mirrors run_deepsomatic.py
// command sequence.
// ──────────────────────────────────────────────────────────────────────
int RunAllSomatic(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  const std::string model_type =
      CanonicalModelType(absl::GetFlag(FLAGS_model_type));
  { std::system(absl::StrCat("mkdir -p '",
      absl::GetFlag(FLAGS_intermediate_results_dir), "'").c_str()); }
  const std::string ref_flag = absl::GetFlag(FLAGS_ref);
  const std::string user_regions = absl::GetFlag(FLAGS_regions);
  const std::string regions_flag = EffectiveRegions(user_regions, ref_flag);
  const std::string tmp_dir = absl::GetFlag(FLAGS_intermediate_results_dir);
  const int num_shards = EffectiveNumShards();
  const int n_threads = std::max(1, num_shards);

  const std::string reads_tumor  = absl::GetFlag(FLAGS_reads_tumor);
  const std::string reads_normal = absl::GetFlag(FLAGS_reads_normal);
  if (reads_tumor.empty()) {
    LOG(ERROR) << "Somatic: --reads_tumor required";
    return 1;
  }
  const bool has_normal = !reads_normal.empty();

  const std::string out_vcf = absl::GetFlag(FLAGS_output_vcf);
  if (out_vcf.empty()) {
    LOG(ERROR) << "Somatic: --output_vcf required";
    return 1;
  }

  std::string ckpt = absl::GetFlag(FLAGS_checkpoint);
  if (ckpt.empty()) {
    // Auto-select model bundle based on model_type × has_normal.
    // For metal backend the user should pass --checkpoint=path/to/.dvw;
    // for coreml/ane_speculate the .mlpackage path is returned here.
    ckpt = SomaticModelPath(model_type, has_normal);
    LOG(INFO) << "Somatic: auto-selected model " << ckpt;
  }
  // Early-fail: catch typos in --reads_tumor / --reads_normal / --ref / ckpt.
  if (!EnsurePathExists(reads_tumor,  "--reads_tumor")  ||
      !EnsureBamIndexed (reads_tumor, "--reads_tumor")  ||
      !EnsurePathExists(reads_normal, "--reads_normal") ||
      !EnsureBamIndexed (reads_normal,"--reads_normal") ||
      !EnsurePathExists(ref_flag,     "--ref")          ||
      !EnsureFastaIndexed(ref_flag) ||
      !EnsurePathExists(ckpt,         "--checkpoint")) {
    return 1;
  }

  std::string sm_path =
      absl::GetFlag(FLAGS_small_model_path_somatic);
  {
    std::string mt = model_type;
    for (char& c : mt) c = static_cast<char>(std::toupper(c));
    const bool expects = SomaticExpectsSmallModel(mt, has_normal);
    MaybeAutoDiscoverTrioOrSomaticSmallModel(
        sm_path, ckpt, "--small_model_path_somatic", expects);
    WarnIfMissingSmallModel(sm_path, "--small_model_path_somatic", mt, expects);
  }

  const std::string inference_backend =
      absl::GetFlag(FLAGS_inference_backend);

  // Per-stage intermediate paths.
  const std::string examples_pattern  =
      n_threads > 1
          ? absl::StrCat(tmp_dir, "/examples_tumor.tfrecord@", n_threads)
          : absl::StrCat(tmp_dir, "/examples_tumor.tfrecord");
  const std::string small_cvo_pattern =
      n_threads > 1
          ? absl::StrCat(tmp_dir, "/small_cvo_tumor.tfrecord@", n_threads)
          : absl::StrCat(tmp_dir, "/small_cvo_tumor.tfrecord");
  const std::string cvo_path        =
      absl::StrCat(tmp_dir, "/cvo_tumor.tfrecord");
  const std::string merged_cvo_path =
      absl::StrCat(tmp_dir, "/merged_cvo_tumor.tfrecord");

  // ── Stage 1: make_examples (tumor + optional normal). ────────────
  LOG(INFO) << "Somatic Stage 1: make_examples ("
            << (has_normal ? "tumor+normal" : "tumor-only")
            << ", --threads=" << n_threads << ")";
  {
    std::vector<std::string> me_args = {
        absl::StrCat("--reads_tumor=", reads_tumor),
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--examples_tumor=", examples_pattern),
        absl::StrCat("--threads=", n_threads),
        "--task_id=0",
        "--num_shards=1",
        "--realigner_enabled=true",
    };
    if (has_normal) {
      me_args.push_back(absl::StrCat("--reads_normal=", reads_normal));
    }
    if (!regions_flag.empty()) {
      me_args.push_back(absl::StrCat("--regions=", regions_flag));
    }
    if (!absl::GetFlag(FLAGS_sample_name_tumor).empty()) {
      me_args.push_back(absl::StrCat("--sample_name_tumor=",
                                      absl::GetFlag(FLAGS_sample_name_tumor)));
    }
    if (!absl::GetFlag(FLAGS_sample_name_normal).empty()) {
      me_args.push_back(absl::StrCat("--sample_name_normal=",
                                      absl::GetFlag(FLAGS_sample_name_normal)));
    }
    if (!sm_path.empty()) {
      me_args.push_back(absl::StrCat("--small_model_path_somatic=", sm_path));
      me_args.push_back(absl::StrCat("--small_model_cvo_outfile_tumor=",
                                      small_cvo_pattern));
    }
    // Per-model flags from deepsomatic.<model>[_tumor_only]/model.example_info.json.
    ApplySomaticModelFlags(model_type, has_normal, me_args);
    // Tumor-only: forward PON VCF path for allele_frequency channel encoding.
    // Priority: explicit --population_vcfs flag > auto-discovered from
    // DEEPVARIANT_MODELS_DIR. Auto-discovery picks the correct PON per model:
    //   PACBIO/ONT → AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz
    //   WGS/WES/FFPE_* → AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz
    if (!has_normal) {
      std::string pon = absl::GetFlag(FLAGS_population_vcfs);
      if (pon.empty()) {
        // Auto-discover PON from models directory.
        const char* env = std::getenv("DEEPVARIANT_MODELS_DIR");
        std::string models_dir = env ? env : "/opt/homebrew/share/deepvariant-models";
        std::string mt_up = model_type;
        for (char& c : mt_up) c = static_cast<char>(std::toupper(c));
        const bool is_long_read = (mt_up == "PACBIO" || mt_up == "ONT");
        const std::string pon_name = is_long_read
            ? "AF_pacbio_PON_CoLoRSdb.GRCh38.AF0.05.vcf.gz"
            : "AF_ilmn_PON_DeepVariant.GRCh38.AF0.05.vcf.gz";
        pon = absl::StrCat(models_dir, "/deepsomatic_pon/", pon_name);
        // Only use auto-discovered path if file exists.
        struct stat st;
        if (stat(pon.c_str(), &st) != 0) pon.clear();
      }
      if (!pon.empty()) {
        me_args.push_back(absl::StrCat("--population_vcfs=", pon));
      }
    }
    AppendSelectVariantTypes(me_args);
    AppendHaploidFlags(me_args);  // haploid gVCF reference confidence
    auto argv_me = MakeArgv("deepvariant_make_examples", me_args);
    int n = static_cast<int>(argv_me.size()) - 1;
    if (int rc = RunMakeExamples(n, argv_me.data()); rc != 0) {
      LOG(ERROR) << "Somatic: make_examples failed";
      return rc;
    }
  }

  // ── Stage 2: call_variants on the tumor model. ────────────
  LOG(INFO) << "Somatic Stage 2: call_variants";
  {
    // Per-model input shape from deepsomatic[_tumor_only] example_info.json.
    const SomaticDims sdims = SomaticInputDims(model_type, has_normal);
    std::vector<std::string> cv_args = {
        absl::StrCat("--examples=", examples_pattern),
        absl::StrCat("--outfile=", cvo_path),
        absl::StrCat("--checkpoint=", ckpt),
        absl::StrCat("--batch_size=", EffectiveBatchSize()),
        absl::StrCat("--inference_backend=", inference_backend),
        absl::StrCat("--input_height=", sdims.h),
        absl::StrCat("--input_channels=", sdims.channels),
        absl::StrCat("--input_width=", sdims.width),
    };
    AppendAneSpeculateArgs(cv_args, inference_backend,
                           absl::GetFlag(FLAGS_ane_speculate_metal_checkpoint_somatic));
    auto argv_cv = MakeArgv("deepvariant_call_variants", cv_args);
    int n = static_cast<int>(argv_cv.size()) - 1;
    if (int rc = RunCallVariants(n, argv_cv.data()); rc != 0) {
      LOG(ERROR) << "Somatic: call_variants failed";
      return rc;
    }
  }

  // ── Stage 2.5: merge small_cvo into cvo (if SM was used). ──
  LOG(INFO) << "Somatic Stage 2.5: merge → " << merged_cvo_path;
  {
    // small_cvo_pattern is a `name@N` shard spec (n_threads>1) or a single
    // file; MergeCvoFiles expands `@N` to the per-shard files. When sm_path is
    // empty no small_cvo was written, so we merge only the big cvo.
    const std::string small = sm_path.empty() ? "" : small_cvo_pattern;
    if (!MergeCvoFiles(small, cvo_path, merged_cvo_path)) {
      LOG(ERROR) << "Somatic: merge step failed";
      return 1;
    }
  }

  // ── Stage 3: postprocess. ────────────
  LOG(INFO) << "Somatic Stage 3: postprocess_variants";
  {
    std::vector<std::string> pp_args = {
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--infile=", merged_cvo_path),
        absl::StrCat("--output_vcf_outfile=", out_vcf),
        "--process_somatic=true",
    };
    // pon_filtering: forward user flag, otherwise stay empty (matches
    // upstream's --use_default_pon_filtering=False default; auto-default
    // is opt-in via --use_default_pon_filtering=true OR by setting
    // --pon_filtering explicitly).
    {
      const std::string user_pon = absl::GetFlag(FLAGS_pon_filtering);
      if (!user_pon.empty()) {
        pp_args.push_back(absl::StrCat("--pon_filtering=", user_pon));
      }
    }
    AppendHaploidFlags(pp_args);
    auto argv_pp = MakeArgv("deepvariant_postprocess", pp_args);
    int n = static_cast<int>(argv_pp.size()) - 1;
    if (int rc = RunPostprocessVariants(n, argv_pp.data()); rc != 0) {
      LOG(ERROR) << "Somatic: postprocess_variants failed";
      return rc;
    }
  }

  LOG(INFO) << "Somatic: done. VCF at " << out_vcf;
  return 0;
}

// ──────────────────────────────────────────────────────────────────────
// Pangenome-aware DV dispatch: 2-sample make_examples
// (pangenome=0, reads=1=main); 1× call_variants on the pangenome
// model (pangenome has skip_output=true); 1× postprocess writing a
// single VCF for the reads sample. Mirrors
// run_pangenome_aware_deepvariant.py command sequence.
// ──────────────────────────────────────────────────────────────────────
int RunAllPangenome(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);
  { std::system(absl::StrCat("mkdir -p '",
      absl::GetFlag(FLAGS_intermediate_results_dir), "'").c_str()); }
  const std::string ref_flag = absl::GetFlag(FLAGS_ref);
  const std::string user_regions = absl::GetFlag(FLAGS_regions);
  const std::string regions_flag = EffectiveRegions(user_regions, ref_flag);
  const std::string tmp_dir = absl::GetFlag(FLAGS_intermediate_results_dir);
  const int num_shards = EffectiveNumShards();
  const int n_threads = std::max(1, num_shards);

  const std::string reads_main      = absl::GetFlag(FLAGS_reads);
  const std::string reads_pangenome = absl::GetFlag(FLAGS_reads_pangenome);
  if (reads_main.empty()) {
    LOG(ERROR) << "Pangenome: --reads required";
    return 1;
  }
  if (reads_pangenome.empty()) {
    LOG(ERROR) << "Pangenome: --reads_pangenome required";
    return 1;
  }

  const std::string out_vcf = absl::GetFlag(FLAGS_output_vcf);
  if (out_vcf.empty()) {
    LOG(ERROR) << "Pangenome: --output_vcf required";
    return 1;
  }

  std::string ckpt = absl::GetFlag(FLAGS_checkpoint);
  if (ckpt.empty()) {
    LOG(ERROR) << "Pangenome: --checkpoint (.dvw) required";
    return 1;
  }
  // Early-fail: catch typos in --reads / --reads_pangenome / --ref / --checkpoint.
  // --reads_pangenome can be a real BAM (synthetic reads from GBZ→BAM
  // preprocessing) so the .bai check applies normally.
  // ref_flag is declared at the top of this function (~line 1580).
  if (!EnsurePathExists(reads_main,      "--reads")           ||
      !EnsureBamIndexed (reads_main,     "--reads")           ||
      !EnsurePathExists(reads_pangenome, "--reads_pangenome") ||
      !EnsureBamIndexed (reads_pangenome,"--reads_pangenome") ||
      !EnsurePathExists(ref_flag,        "--ref")             ||
      !EnsureFastaIndexed(ref_flag) ||
      !EnsurePathExists(ckpt,            "--checkpoint")) {
    return 1;
  }

  const std::string sm_path =
      absl::GetFlag(FLAGS_small_model_path_pangenome);

  const std::string inference_backend =
      absl::GetFlag(FLAGS_inference_backend);

  // Per-stage intermediate paths (named after the reads sample).
  const std::string examples_pattern  =
      n_threads > 1
          ? absl::StrCat(tmp_dir, "/examples_reads.tfrecord@", n_threads)
          : absl::StrCat(tmp_dir, "/examples_reads.tfrecord");
  const std::string small_cvo_pattern =
      n_threads > 1
          ? absl::StrCat(tmp_dir, "/small_cvo_reads.tfrecord@", n_threads)
          : absl::StrCat(tmp_dir, "/small_cvo_reads.tfrecord");
  const std::string cvo_path        =
      absl::StrCat(tmp_dir, "/cvo_reads.tfrecord");
  const std::string merged_cvo_path =
      absl::StrCat(tmp_dir, "/merged_cvo_reads.tfrecord");

  // ── Stage 1: make_examples (reads + pangenome). ────────────
  LOG(INFO) << "Pangenome Stage 1: make_examples (--threads=" << n_threads
            << ")";
  {
    std::vector<std::string> me_args = {
        absl::StrCat("--reads=", reads_main),
        absl::StrCat("--reads_pangenome=", reads_pangenome),
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--examples_reads=", examples_pattern),
        absl::StrCat("--threads=", n_threads),
        "--task_id=0",
        "--num_shards=1",
        "--realigner_enabled=true",
    };
    if (!regions_flag.empty()) {
      me_args.push_back(absl::StrCat("--regions=", regions_flag));
    }
    if (!absl::GetFlag(FLAGS_sample_name_reads).empty()) {
      me_args.push_back(absl::StrCat("--sample_name_reads=",
                                      absl::GetFlag(FLAGS_sample_name_reads)));
    }
    if (!absl::GetFlag(FLAGS_sample_name_pangenome).empty()) {
      me_args.push_back(absl::StrCat("--sample_name_pangenome=",
                                      absl::GetFlag(FLAGS_sample_name_pangenome)));
    }
    if (!sm_path.empty()) {
      me_args.push_back(absl::StrCat("--small_model_path_pangenome=", sm_path));
      me_args.push_back(absl::StrCat("--small_model_cvo_outfile_reads=",
                                      small_cvo_pattern));
    }
    // Pangenome WGS overrides per /opt/models/pangenome_aware_deepvariant/
    // wgs/model.example_info.json:flags_for_calling. Upstream's
    // make_examples_core.py:apply_flags_for_calling reads this file at
    // runtime; we hard-code the WGS values here. Note: pangenome uses
    // the GLOBAL default vsc_min_fraction_{snps,indels} (0.12 / 0.06);
    // only min_mapping_quality is overridden to 0 (vs default 5).
    me_args.push_back("--min_mapping_quality=0");
    // Realigner SSW alignment scoring (defaults are 4/6/8/2 for WGS).
    me_args.push_back("--aln_match=2");
    me_args.push_back("--aln_mismatch=5");
    me_args.push_back("--aln_gap_open=10");
    me_args.push_back("--aln_gap_extend=1");
    me_args.push_back("--dbg_disable_graph_pruning=true");
    // Pangenome uses the DEFAULT partition_size=1000 (matching upstream:
    // run_pangenome_aware_deepvariant.py does NOT pass --partition_size, and
    // forcing 25000 in Docker errors because make_examples requires
    // --partition_size and --max_reads_per_partition to be set together).
    // Earlier we hardcoded 25000 (Phase 6 Step 3-v8) believing it matched
    // Docker, but that was wrong: with 25kb partitions, the per-partition
    // reservoir sampling (max_reads_per_partition=1500) aggressively
    // downsamples reads in high-coverage windows, dropping the few alt reads
    // at low-coverage candidate clusters (e.g. chr20:10029223-10029235, a run
    // of A>G SNPs with ~10-12 supporting reads each that Docker calls PASS but
    // 25kb-partition reservoir sampling reduced to ~1, killing the candidate).
    // partition_size=1000 mirrors Docker's per-1kb reservoir granularity.
    me_args.push_back("--partition_size=1000");
    AppendSelectVariantTypes(me_args);
    AppendHaploidFlags(me_args);  // haploid gVCF reference confidence
    auto argv_me = MakeArgv("deepvariant_make_examples", me_args);
    int n = static_cast<int>(argv_me.size()) - 1;
    if (int rc = RunMakeExamples(n, argv_me.data()); rc != 0) {
      LOG(ERROR) << "Pangenome: make_examples failed";
      return rc;
    }
  }

  // ── Stage 2: call_variants on the pangenome model. ────────────
  LOG(INFO) << "Pangenome Stage 2: call_variants";
  {
    // Pangenome WGS pileup is 200×221×7 (pangenome 100 + reads 100).
    std::vector<std::string> cv_args = {
        absl::StrCat("--examples=", examples_pattern),
        absl::StrCat("--outfile=", cvo_path),
        absl::StrCat("--checkpoint=", ckpt),
        absl::StrCat("--batch_size=", EffectiveBatchSize()),
        absl::StrCat("--inference_backend=", inference_backend),
        "--input_height=200",
        "--input_channels=7",
    };
    AppendAneSpeculateArgs(cv_args, inference_backend,
                           absl::GetFlag(FLAGS_ane_speculate_metal_checkpoint_pangenome));
    auto argv_cv = MakeArgv("deepvariant_call_variants", cv_args);
    int n = static_cast<int>(argv_cv.size()) - 1;
    if (int rc = RunCallVariants(n, argv_cv.data()); rc != 0) {
      LOG(ERROR) << "Pangenome: call_variants failed";
      return rc;
    }
  }

  // ── Stage 2.5: merge small_cvo into cvo. ──
  LOG(INFO) << "Pangenome Stage 2.5: merge → " << merged_cvo_path;
  {
    // small_cvo_pattern is a `name@N` shard spec (n_threads>1) or a single
    // file; MergeCvoFiles expands `@N` to the per-shard files. When sm_path is
    // empty no small_cvo was written, so we merge only the big cvo.
    const std::string small = sm_path.empty() ? "" : small_cvo_pattern;
    if (!MergeCvoFiles(small, cvo_path, merged_cvo_path)) {
      LOG(ERROR) << "Pangenome: merge step failed";
      return 1;
    }
  }

  // ── Stage 3: postprocess. ────────────
  LOG(INFO) << "Pangenome Stage 3: postprocess_variants";
  {
    std::vector<std::string> pp_args = {
        absl::StrCat("--ref=", ref_flag),
        absl::StrCat("--infile=", merged_cvo_path),
        absl::StrCat("--output_vcf_outfile=", out_vcf),
    };
    AppendHaploidFlags(pp_args);
    auto argv_pp = MakeArgv("deepvariant_postprocess", pp_args);
    int n = static_cast<int>(argv_pp.size()) - 1;
    if (int rc = RunPostprocessVariants(n, argv_pp.data()); rc != 0) {
      LOG(ERROR) << "Pangenome: postprocess_variants failed";
      return rc;
    }
  }

  LOG(INFO) << "Pangenome: done. VCF at " << out_vcf;
  return 0;
}

}  // namespace deepvariant

// MultiCallTool — what tool name the binary was invoked as. Set in main()
// from basename(argv[0]). Drives the per-tool help text and dispatch path.
//
// Values:
//   "deepvariant"                  — canonical binary, full subcommand suite
//   "deeptrio"                     — multi-call alias → forces trio mode
//   "deepsomatic"                  — multi-call alias → forces somatic mode
//   "pangenome-aware-deepvariant"  — multi-call alias → forces pangenome mode
//
// Mirrors upstream Google's three-binary convention (`run_deepvariant`,
// `run_deeptrio`, `run_deepsomatic`, `run_pangenome_aware_deepvariant`)
// without the disk-bloat / version-skew cost of three separate executables:
// classic Unix multi-call binary (busybox-style). Homebrew formula will
// install `deepvariant` and three symlinks pointing at it.
static std::string g_multicall_tool;  // empty until set in main()

// PrintTopLevelHelp — top-level help for whichever tool the binary was
// invoked as. Goes to stdout (it's information, not an error).
static void PrintTopLevelHelp() {
  if (g_multicall_tool == "deeptrio") {
    std::printf(
        "deeptrio — DeepTrio (child + parent1 + parent2) on Apple Silicon\n"
        "\n"
        "Usage: deeptrio --reads=<child.BAM> --reads_parent1=<BAM> --reads_parent2=<BAM> \\\n"
        "                --ref=<FASTA> --output_vcf=<child.VCF> \\\n"
        "                --output_vcf_parent1=<VCF> --output_vcf_parent2=<VCF> \\\n"
        "                [--model_type=WGS|PACBIO|ONT] [--regions=chr20]\n"
        "\n"
        "Note: the child sample uses the unsuffixed --reads / --output_vcf flags;\n"
        "parent samples use --reads_parent{1,2} / --output_vcf_parent{1,2}. The\n"
        "presence of --reads_parent1 is what triggers trio dispatch.\n"
        "\n"
        "Get all flags:           deeptrio --helpfull\n"
        "Search by name/keyword:  deeptrio --help=<substr>\n"
        "\n"
        "Equivalent canonical form: deepvariant trio <flags...>\n");
    return;
  }
  if (g_multicall_tool == "deepsomatic") {
    std::printf(
        "deepsomatic — DeepSomatic (tumor + optional normal) on Apple Silicon\n"
        "\n"
        "Usage: deepsomatic --reads_tumor=<BAM> [--reads_normal=<BAM>] \\\n"
        "                   --ref=<FASTA> --output_vcf=<VCF> \\\n"
        "                   [--model_type=WGS|PACBIO|ONT|FFPE_WGS|...] [--regions=chr20]\n"
        "\n"
        "Tumor-only:  omit --reads_normal (the model dispatch is automatic).\n"
        "\n"
        "Get all flags:           deepsomatic --helpfull\n"
        "Search by name/keyword:  deepsomatic --help=<substr>\n"
        "\n"
        "Equivalent canonical form: deepvariant somatic <flags...>\n");
    return;
  }
  if (g_multicall_tool == "pangenome-aware-deepvariant") {
    std::printf(
        "pangenome-aware-deepvariant — Pangenome-aware DV on Apple Silicon\n"
        "\n"
        "Usage: pangenome-aware-deepvariant --reads=<BAM> --reads_pangenome=<BAM> \\\n"
        "                                   --ref=<FASTA> --output_vcf=<VCF> \\\n"
        "                                   [--regions=chr20]\n"
        "\n"
        "The pangenome BAM is a GBZ-derived synthetic-haplotype BAM produced by\n"
        "the upstream Docker preprocessing step. GBZ at runtime is out of scope\n"
        "for v2; convert GBZ→BAM once via the documented Docker pipeline.\n"
        "\n"
        "Get all flags:           pangenome-aware-deepvariant --helpfull\n"
        "Search by name/keyword:  pangenome-aware-deepvariant --help=<substr>\n"
        "\n"
        "Equivalent canonical form: deepvariant pangenome <flags...>\n");
    return;
  }
  // Default: canonical `deepvariant` binary — full subcommand suite.
  std::printf(
      "deepvariant — Apple Silicon native port (v2)\n"
      "\n"
      "Usage: deepvariant <subcommand> [flags]\n"
      "\n"
      "Top-level pipelines (one BAM in, one VCF out):\n"
      "  run                 single-sample germline (WGS, WES, PACBIO, ONT, ...)\n"
      "  trio                3-sample DeepTrio (child + parent1 + parent2)\n"
      "  somatic             DeepSomatic (tumor + optional normal)\n"
      "  pangenome           pangenome-aware DeepVariant (BAM + GBZ-derived BAM)\n"
      "\n"
      "Stage-level subcommands (compose your own pipeline):\n"
      "  make_examples       BAM → tfrecord pileup examples\n"
      "  call_variants       examples → CVO via Inception-v3 / small_model\n"
      "  postprocess_variants  CVO → final VCF (+ optional gVCF)\n"
      "\n"
      "Multi-call shortcuts (Homebrew-installed symlinks; same binary, no\n"
      "version skew vs the canonical form):\n"
      "  deeptrio                       → deepvariant trio\n"
      "  deepsomatic                    → deepvariant somatic\n"
      "  pangenome-aware-deepvariant    → deepvariant pangenome\n"
      "\n"
      "Get per-subcommand flag help:  deepvariant <subcommand> --help\n");
}

// PrintVersion — print version + build metadata to stdout. Goes to stdout
// (not stderr) so the user can capture it for issue reports / CI logs.
//
// Format mirrors common Unix conventions:
//   <tool> <our-tag> (DeepVariant <upstream-version>, build <sha> <date>)
//   Apple Silicon native port — <arch>, <abi>
//
// All four substitutions are compile-time constants (DV_VERSION,
// DV_UPSTREAM_VERSION, DV_GIT_SHA, DV_BUILD_DATE) baked in via
// target_compile_definitions in CMakeLists.txt.
static void PrintVersion() {
  // Compile-time arch detection. We're arm64-only at runtime (per CLAUDE.md
  // "macOS ≥ 14, arm64 only") but emit the actual built arch for honesty.
#if defined(__aarch64__) || defined(__arm64__)
  constexpr const char* kArch = "arm64";
#elif defined(__x86_64__)
  constexpr const char* kArch = "x86_64";
#else
  constexpr const char* kArch = "unknown";
#endif
  std::printf(
      "%s %s (DeepVariant %s, build %s %s)\n"
      "Apple Silicon native port — %s, macOS\n",
      g_multicall_tool.c_str(),
      DV_VERSION, DV_UPSTREAM_VERSION,
      DV_GIT_SHA, DV_BUILD_DATE,
      kArch);
}

// DetectMultiCall — return the subcommand name to inject when the binary
// is invoked under one of its multi-call basenames. Empty string for the
// canonical `deepvariant` invocation (or any unrecognized basename).
static std::string DetectMultiCall(const char* argv0) {
  // basename(): take the substring after the last '/'.
  std::string base(argv0);
  const auto slash = base.find_last_of('/');
  if (slash != std::string::npos) base = base.substr(slash + 1);
  // Tolerate a `.exe` suffix (no-op on macOS but cheap and harmless).
  constexpr absl::string_view kExe = ".exe";
  if (base.size() > kExe.size() &&
      base.substr(base.size() - kExe.size()) == kExe) {
    base.resize(base.size() - kExe.size());
  }
  if (base == "deeptrio") return "trio";
  if (base == "deepsomatic") return "somatic";
  if (base == "pangenome-aware-deepvariant") return "pangenome";
  return "";  // canonical or unknown → no rewrite
}

int main(int argc, char** argv) {
  absl::InitializeLog();
  // Default log level: send INFO to stderr.
  absl::SetStderrThreshold(absl::LogSeverity::kInfo);

  // Multi-call dispatch — busybox-style. If the binary was invoked as
  // `deeptrio`, `deepsomatic`, or `pangenome-aware-deepvariant` (via a
  // Homebrew-installed symlink), inject the corresponding subcommand and
  // record the tool name globally so help text and SetProgramUsageMessage
  // address the right tool.
  const std::string injected = DetectMultiCall(argv[0]);
  if (!injected.empty()) {
    if      (injected == "trio")      g_multicall_tool = "deeptrio";
    else if (injected == "somatic")   g_multicall_tool = "deepsomatic";
    else if (injected == "pangenome") g_multicall_tool = "pangenome-aware-deepvariant";
  } else {
    g_multicall_tool = "deepvariant";
  }

  // Make `--help` (no args) print our flags. Abseil's default `--help`
  // matches flags whose source-file path contains the program basename,
  // which here is `deepvariant` — but our flags live in files like
  // `deepvariant/native/cli.cc`, `make_examples_main.cc`, etc., none of
  // which match "deepvariant" as a whole filename. The default behavior
  // is therefore "No flags matched". Override:
  //   - contains_helpshort_flags / contains_help_flags: match any file
  //     under `deepvariant/native/` so plain `--help` shows our ~80
  //     flags grouped by source file (and excludes absl/grpc internals).
  //   - SetProgramUsageMessage: silences the noisy
  //     "Warning: SetProgramUsageMessage() never called" emitted by
  //     ParseCommandLine, and tells the user how to drill in further.
  absl::FlagsUsageConfig usage_config;
  usage_config.contains_helpshort_flags = [](absl::string_view path) {
    return absl::StrContains(path, "deepvariant/native/");
  };
  usage_config.contains_help_flags = [](absl::string_view path) {
    return absl::StrContains(path, "deepvariant/native/");
  };
  absl::SetFlagsUsageConfig(usage_config);
  // Per-tool usage message — printed by absl::ParseCommandLine when the
  // user passes --help. We hold the storage in a function-local static so
  // it outlives the call (Abseil stores a string_view internally).
  static const std::string usage_msg = absl::StrCat(
      g_multicall_tool, " — Apple Silicon native port (v2)\n"
      "\n"
      "Get the full flag list:    ", g_multicall_tool, " --helpfull\n",
      "Get the main flags:        ", g_multicall_tool, " --help\n",
      "Search by name/keyword:    ", g_multicall_tool, " --help=<substr>");
  absl::SetProgramUsageMessage(usage_msg);

  // Multi-call binary path: dispatch directly to the injected mode.
  // The runner's own absl::ParseCommandLine handles --help / --helpfull /
  // --help=<substr>. We still intercept top-level help / version words
  // here so the user sees a per-tool synopsis (PrintTopLevelHelp branches
  // on g_multicall_tool) and version banner before any flag parsing.
  if (!injected.empty()) {
    if (argc >= 2) {
      const std::string a1(argv[1]);
      if (a1 == "-h" || a1 == "--help" || a1 == "help") {
        PrintTopLevelHelp();
        return 0;
      }
      if (a1 == "-v" || a1 == "--version" || a1 == "version") {
        PrintVersion();
        return 0;
      }
    }
    if (injected == "trio")      return deepvariant::RunAllTrio(argc, argv);
    if (injected == "somatic")   return deepvariant::RunAllSomatic(argc, argv);
    if (injected == "pangenome") return deepvariant::RunAllPangenome(argc, argv);
    // (Defensive — DetectMultiCall only returns the three names above.)
    LOG(ERROR) << "Unknown multi-call alias: " << injected;
    return 1;
  }

  // Canonical `deepvariant` invocation — subcommand-style dispatch.
  if (argc < 2) {
    PrintTopLevelHelp();
    return 0;  // No-arg invocation is informational, not an error.
  }

  const std::string sub(argv[1]);

  // Top-level help: -h, --help, help → print help to stdout, return 0.
  // (Per-subcommand --help is still handled by absl::ParseCommandLine
  // inside each Run* dispatcher.)
  if (sub == "-h" || sub == "--help" || sub == "help") {
    PrintTopLevelHelp();
    return 0;
  }

  // Top-level version: -v, --version, version → print version + build
  // metadata to stdout, return 0. Same convention as `git --version`,
  // `samtools --version`, `bcftools --version`, etc.
  if (sub == "-v" || sub == "--version" || sub == "version") {
    PrintVersion();
    return 0;
  }

  // Shift argv so subcommand sees its own flags.
  argv[1] = argv[0];
  int new_argc = argc - 1;
  char** new_argv = argv + 1;

  if (sub == "make_examples") {
    return deepvariant::RunMakeExamples(new_argc, new_argv);
  } else if (sub == "call_variants") {
    return deepvariant::RunCallVariants(new_argc, new_argv);
  } else if (sub == "postprocess_variants") {
    return deepvariant::RunPostprocessVariants(new_argc, new_argv);
  } else if (sub == "run") {
    return deepvariant::RunAll(new_argc, new_argv);
  } else if (sub == "trio") {
    return deepvariant::RunAllTrio(new_argc, new_argv);
  } else if (sub == "somatic") {
    return deepvariant::RunAllSomatic(new_argc, new_argv);
  } else if (sub == "pangenome") {
    return deepvariant::RunAllPangenome(new_argc, new_argv);
  } else {
    LOG(ERROR) << "Unknown subcommand: " << sub;
    PrintTopLevelHelp();
    return 1;
  }
}
