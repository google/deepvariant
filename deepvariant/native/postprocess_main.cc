// Native postprocess_variants — calling mode.
//
// Reads CallVariantsOutput TFRecords, groups by genomic site (multi-allelic
// merge), assigns the most-likely diploid genotype, and writes VCF with
// FORMAT fields GT:GQ:DP:AD:VAF:PL.
//
// Multi-allelic merge: upstream make_examples emits one example per
// alt-allele combination at multi-allelic sites (multi_allelic_mode =
// ADD_HET_ALT_IMAGES). Each resulting CVO carries:
//   - the same Variant (with the full alt list)
//   - cvo.alt_allele_indices.indices: which alt(s) the example tested
//   - cvo.genotype_probabilities: 3-vector
//     - if indices == [i]:   [P(0/0), P(0/(i+1)), P((i+1)/(i+1))]
//     - if indices == [i,j]: [P(other), P((i+1)/(j+1)), P(other)]
// We collect these into a likelihood table over all diploid genotypes,
// pick argmax, and emit one VCF line per site.

#include "deepvariant/native/postprocess_main.h"

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <fstream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "deepvariant/native/haploid_regions.h"
#include "deepvariant/native/haplotypes.h"
#include "deepvariant/native/tfrecord.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/log/check.h"
#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "absl/strings/numbers.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_join.h"
#include "absl/strings/str_split.h"
#include "third_party/nucleus/io/merge_variants.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/io/variant_reader.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/io/vcf_writer.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"

ABSL_FLAG(std::string, infile, "", "Input CVO TFRecord path (may be sharded).");
ABSL_DECLARE_FLAG(std::string, ref);
ABSL_DECLARE_FLAG(std::string, sample_name);
ABSL_FLAG(std::string, output_vcf_outfile, "", "Output VCF path.");
ABSL_FLAG(std::string, gvcf_outfile, "", "gVCF output path (optional).");
ABSL_FLAG(std::string, nonvariant_site_tfrecord_path, "",
          "Phase 9 / Step 3 — input non-variant TFRecord(s) produced by "
          "make_examples --gvcf=... (sharded `name@N` spec). Required when "
          "--gvcf_outfile is set; merged with the variant CVO stream via "
          "nucleus::MergeAndWriteVariantsAndNonVariants.");
ABSL_FLAG(bool, enable_temp_scaling, false,
          "Phase 8 / Tier 4 — apply post-CombineLikelihoods temperature "
          "scaling to softmax probabilities before argmax/QUAL/GQ/PL "
          "computation. Implements Guo et al. ICML 2017 calibration. "
          "Off by default to preserve baseline FILTER parity. When on, "
          "use --temp_scaling_T to set the temperature.");
ABSL_FLAG(double, temp_scaling_T, 1.0,
          "Temperature parameter for --enable_temp_scaling. T=1.0 is "
          "identity (no change). T>1 smooths probabilities (less "
          "confident, fewer PASS); T<1 sharpens (more confident, more "
          "PASS). Optimal T fit on a held-out chr21 set; ship value "
          "is determined empirically.");
ABSL_FLAG(double, qual_filter, 1.0,
          "Variants with QUAL below this become RefCall instead of PASS.");
// Default 20.0 matches upstream postprocess_variants.py default. When a
// CNN RefCall has GQ < this, upstream rewrites it to "./.": NoCall (no
// determination, low confidence). We mirror that exactly.
ABSL_FLAG(double, cnn_homref_call_min_gq, 20.0,
          "All CNN RefCalls whose GQ is less than this become ./. NoCall "
          "instead of 0/0 RefCall (matches upstream default 20.0).");
ABSL_FLAG(std::string, pon_filtering, "",
          "Optional. Only used if --process_somatic=true. Path to a Panel-of-"
          "Normals VCF. Variants whose (CHROM,POS,REF,ALT) matches the PON have "
          "PASS removed and FILTER set to PON. Mirrors upstream "
          "postprocess_variants.py:--pon_filtering. Auto-discovered by cli.cc "
          "for tumor-only modes when DEEPVARIANT_MODELS_DIR is set.");
ABSL_FLAG(bool, process_somatic, false,
          "Enable DeepSomatic-style postprocess: heterozygous (0/1) calls "
          "are reclassified as GERMLINE 0/0 (mirrors third_party/nucleus/"
          "io/vcf_writer.cc::WriteSomatic logic).");
// multiallelic_mode: CVO probability fusion for sites with >1 ALT.
// Mirrors upstream postprocess_variants.py FLAGS.multiallelic_mode from
// model example_info.json flags_for_postprocessing.
//   "product" (default/WGS): multiply probabilities across CVOs.
//   "min" (WES):             take minimum probability across kept CVOs.
ABSL_FLAG(std::string, multiallelic_mode, "product",
          "Multi-allelic CVO fusion: product (WGS default) or min (WES).");
// Sex-chromosome haploid calling (mirror of upstream postprocess_variants.py
// --haploid_contigs / --par_regions_bed). On a haploid contig (e.g. chrX/chrY
// in an XY sample) outside the pseudo-autosomal regions, heterozygous
// genotypes are disallowed: their probabilities are zeroed and the vector
// renormalized, forcing a haploid (homozygous) call.
ABSL_FLAG(std::string, haploid_contigs, "",
          "Comma/space-separated contigs to call as haploid (e.g. "
          "\"chrX,chrY\" for GRCh38, \"X,Y\" for GRCh37). Empty = all diploid.");
ABSL_FLAG(std::string, par_regions_bed, "",
          "BED of pseudo-autosomal regions exempted from haploid calling on "
          "the --haploid_contigs (they stay diploid). Empty = no exemptions.");

namespace deepvariant {

using learning::genomics::deepvariant::CallVariantsOutput;
using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;

namespace {

constexpr int kMaxPhred = 99;

// ---------------------------------------------------------------------------
// Sex-chromosome haploid calling (--haploid_contigs / --par_regions_bed).
// ParRegions / ParseHaploidContigs / LoadParRegions / IsHaploidPosition are
// shared with the make_examples gVCF path via haploid_regions.h.
// ---------------------------------------------------------------------------

// True when the variant should be called haploid: it sits on a --haploid_contig
// and does not overlap a PAR region. Mirror of postprocess_variants.py's
// `is_non_autosome(v) and not is_in_regions(v, par_regions)`.
bool IsHaploidVariant(const Variant& variant,
                      const std::set<std::string>& haploid_contigs,
                      const ParRegions& par_regions) {
  // Upstream is_in_regions -> RangeSet.variant_overlaps tests only the single
  // 0-based point variant.start (not the variant span), so a multi-base
  // variant whose start is outside every PAR interval is still corrected even
  // if its tail reaches into one. Match that exactly: probe [start, start+1).
  const int64_t start = variant.start();
  return IsHaploidPosition(variant.reference_name(), start, start + 1,
                           haploid_contigs, par_regions);
}

std::vector<std::string> ExpandShards(const std::string& spec) {
  auto at = spec.find('@');
  if (at == std::string::npos) return {spec};
  const std::string prefix = spec.substr(0, at);
  int n;
  if (!absl::SimpleAtoi(spec.substr(at + 1), &n) || n <= 0) return {spec};
  std::vector<std::string> paths;
  for (int i = 0; i < n; ++i) {
    paths.push_back(absl::StrCat(prefix, "-", absl::Dec(i, absl::kZeroPad5),
                                  "-of-", absl::Dec(n, absl::kZeroPad5)));
  }
  return paths;
}

// Mirror of nucleus/util/variant_utils.py:simplify_alleles — strips
// the longest common POSTFIX shared by all (ref, alts), leaving at
// least 1 base on every allele. Then updates variant.reference_bases,
// alternate_bases, and end. Required to match upstream's
// `merge_predictions:simplify_variant_alleles(canonical_variant)` call —
// without it we get site-extending substitutions where upstream emits
// clean SNPs (e.g. chr20:63221577 T>C encoded as a 36-bp tandem-repeat
// substitution → false overlap with neighbouring variants → spurious
// haplotype-resolution flips).
void SimplifyVariantAlleles(Variant* variant) {
  if (!variant || variant->reference_bases().empty() ||
      variant->alternate_bases_size() == 0) return;

  size_t shortest = variant->reference_bases().size();
  for (const auto& a : variant->alternate_bases()) {
    shortest = std::min(shortest, a.size());
  }
  // Find longest common postfix length, capped at shortest-1 (each allele
  // must keep at least 1 base).
  size_t common_postfix = 0;
  for (size_t i = 1; i < shortest; ++i) {
    char ref_c = variant->reference_bases()[
        variant->reference_bases().size() - i];
    bool all_same = true;
    for (const auto& a : variant->alternate_bases()) {
      if (a[a.size() - i] != ref_c) { all_same = false; break; }
    }
    if (!all_same) break;
    common_postfix = i;
  }
  if (common_postfix == 0) return;

  std::string new_ref = variant->reference_bases().substr(
      0, variant->reference_bases().size() - common_postfix);
  variant->set_reference_bases(new_ref);
  for (auto& a : *variant->mutable_alternate_bases()) {
    a = a.substr(0, a.size() - common_postfix);
  }
  variant->set_end(variant->start() + new_ref.size());
}

// Convert probability p (in [0,1]) to a phred score, capped at 99.
// Truncates toward zero (matching upstream's vcf_conversion.cc, which
// converts the double-valued Log10PErrorToPhred() into a std::vector<int>
// via implicit narrowing rather than std::round).
int ProbToPhred(double p) {
  if (p >= 1.0) return 0;
  if (p <= 0.0) return kMaxPhred;
  int phred = static_cast<int>(-10.0 * std::log10(p));
  return std::min(std::max(phred, 0), kMaxPhred);
}

// Number of diploid genotypes for a variant with `n_alts` alternates:
// 0/0, 0/1, 1/1, 0/2, 1/2, 2/2, ... = (n_alleles)*(n_alleles+1)/2.
int NumDiploidGenotypes(int n_alts) {
  const int n_alleles = n_alts + 1;
  return n_alleles * (n_alleles + 1) / 2;
}

// Return the two-allele genotype (a, b) with a <= b for the given VCF PL
// index. PL ordering: F(j/k) = k*(k+1)/2 + j  (j <= k).
std::pair<int, int> GenotypeFromPLIndex(int pl_index, int n_alts) {
  for (int k = 0; k <= n_alts; ++k) {
    for (int j = 0; j <= k; ++j) {
      const int idx = k * (k + 1) / 2 + j;
      if (idx == pl_index) return {j, k};
    }
  }
  return {0, 0};  // fallback
}

// QUAL of an alt allele. Mirrors upstream
// postprocess_variants.py:compute_quals(predictions, prediction_index=0)
// EXACTLY, including the `_QUAL_PRECISION=7` rounding step:
//
//   qual = ptrue_to_bounded_phred(min(sum(predictions[1:]), 1.0))
//        = -10 * log10(1 - sum_alt)
//   rounded_qual = round(qual, 7)
//
// The rounding is load-bearing for the AltsToRemove tie-break at
// saturated multi-allelic homref sites: there `sum_alt` is sub-ULP-
// different across alts (FP-drift between our scalar BNNS-CPU softmax
// and Docker's vectorised TF/Keras Eigen softmax), and without
// rounding the relative qual ordering can flip vs Docker. Rounding to
// 7 decimals collapses values < 5e-8 to 0 (so they tie and the first-
// iterated alt wins, matching Docker), while values ≥ 5e-8 survive at
// 1e-7 granularity (preserving Docker's pick when one alt is
// genuinely ahead). Closes the chr20 14/14 site-set diff.
double AltAlleleQual(const CallVariantsOutput& cvo) {
  if (cvo.genotype_probabilities_size() < 3) return 0.0;
  double sum_alt = 0.0;
  for (int i = 1; i < cvo.genotype_probabilities_size(); ++i) {
    sum_alt += cvo.genotype_probabilities(i);
  }
  if (sum_alt <= 0.0) return 0.0;
  if (sum_alt >= 1.0 - 1.25e-10) return kMaxPhred;
  double qual = -10.0 * std::log10(1.0 - sum_alt);
  if (qual > kMaxPhred) qual = kMaxPhred;
  // Round to 7 decimals (upstream's _QUAL_PRECISION).
  return std::round(qual * 1e7) / 1e7;
}

// Returns the set of alt-allele strings to remove from the variant.
// Mirror of postprocess_variants.py:get_alt_alleles_to_remove. An alt is
// flagged for removal when its QUAL (= phred(p_ref)) is below qual_filter.
// If every alt would be removed, the one with the highest QUAL is kept.
std::set<std::string> AltsToRemove(
    const std::vector<const CallVariantsOutput*>& cvos,
    double qual_filter) {
  std::set<std::string> to_remove;
  if (qual_filter <= 0.0 || cvos.empty()) return to_remove;
  const auto& canonical = cvos.front()->variant();
  std::string max_qual_allele;
  double max_qual = -1.0;
  for (const auto* cvo : cvos) {
    const auto& indices = cvo->alt_allele_indices().indices();
    if (indices.size() != 1) continue;
    const int idx = indices[0];
    if (idx < 0 || idx >= canonical.alternate_bases_size()) continue;
    const std::string& alt = canonical.alternate_bases(idx);
    const double qual = AltAlleleQual(*cvo);
    if (qual > max_qual) {
      max_qual = qual;
      max_qual_allele = alt;
    }
    if (qual < qual_filter) to_remove.insert(alt);
  }
  if (!max_qual_allele.empty() &&
      static_cast<int>(to_remove.size()) ==
          canonical.alternate_bases_size()) {
    to_remove.erase(max_qual_allele);  // keep the strongest one
  }
  return to_remove;
}

// Combine all CVOs for one site into a per-genotype likelihood vector.
// Mirror of postprocess_variants.py:merge_predictions "product" mode.
//
// CVOs whose alt-set intersects `alts_to_remove` are SKIPPED (they're
// "for pruned alleles"; upstream merge_predictions ignores them at line
// 1247-1248: `if is_for_pruned_allele: continue`). After pruning the
// last G alt, only the C-alt CVO contributes → predictions are exactly
// that CVO's softmax, not a multi-CVO product. This is what upstream
// does and matters for GQ at multi-allelic sites where ALL but one alt
// gets pruned.
//
// For each diploid genotype (allele1, allele2), each CVO contributes
// cvo.probs[overlap] where overlap = #{alleles in cvo's alt set}, computed
// per allele1, allele2 ∈ {ref, alt1, alt2, …}. Per-CVO contributions are
// fused by product, then normalised across all genotypes.
//
// PL ordering (VCF "G" Number): F(j/k) = k*(k+1)/2 + j  (j ≤ k).
std::vector<double> CombineLikelihoods(
    const std::vector<const CallVariantsOutput*>& cvos, int n_alts,
    const std::set<std::string>& alts_to_remove) {
  const int n_gt = NumDiploidGenotypes(n_alts);
  std::vector<double> like(n_gt, 1.0);  // multiplicative identity

  if (cvos.empty()) return like;
  // All CVOs of a site share the same `variant` (ADD_HET_ALT_IMAGES); take
  // the alt list from the first.
  const auto& alts = cvos.front()->variant().alternate_bases();

  auto pl_idx = [](int j, int k) {
    if (j > k) std::swap(j, k);
    return k * (k + 1) / 2 + j;
  };

  // Genotype 0 = REF, alleles 1..n_alts = alternate_bases[0..n_alts-1].
  // For the "in this CVO's alt set" check we need each cvo's set of alt
  // strings (from alt_allele_indices). Filter out CVOs that touch a
  // pruned alt.
  std::vector<std::set<std::string>> per_cvo_alts;
  std::vector<bool> per_cvo_kept;
  per_cvo_alts.reserve(cvos.size());
  per_cvo_kept.reserve(cvos.size());
  size_t n_kept = 0;
  for (const auto* cvo : cvos) {
    std::set<std::string> s;
    bool touches_pruned = false;
    for (int idx : cvo->alt_allele_indices().indices()) {
      if (idx >= 0 && idx < alts.size()) {
        s.insert(alts[idx]);
        if (alts_to_remove.count(alts[idx])) touches_pruned = true;
      }
    }
    per_cvo_alts.push_back(std::move(s));
    per_cvo_kept.push_back(!touches_pruned);
    if (!touches_pruned) ++n_kept;
  }

  // For every diploid genotype, fuse probabilities across kept CVOs.
  const bool use_min_mode = (absl::GetFlag(FLAGS_multiallelic_mode) == "min");
  for (int k = 0; k <= n_alts; ++k) {
    for (int j = 0; j <= k; ++j) {
      const std::string a1 = (j == 0) ? "" : alts[j - 1];  // "" = REF
      const std::string a2 = (k == 0) ? "" : alts[k - 1];
      // Collect per-CVO probability for this genotype.
      std::vector<double> cvo_probs;
      for (size_t ci = 0; ci < cvos.size(); ++ci) {
        if (!per_cvo_kept[ci]) continue;
        const auto& probs = cvos[ci]->genotype_probabilities();
        if (probs.size() < 3) continue;
        const int overlap = (a1.empty() ? 0 : per_cvo_alts[ci].count(a1)) +
                            (a2.empty() ? 0 : per_cvo_alts[ci].count(a2));
        cvo_probs.push_back(probs[overlap]);
      }
      double fused = 1.0;
      if (!cvo_probs.empty()) {
        if (use_min_mode) {
          // WES: min-probability fusion (upstream multiallelic_mode='min').
          // For each genotype, take the minimum probability across kept CVOs
          // (mirrors postprocess_variants.py::min_alt_filter).
          fused = *std::min_element(cvo_probs.begin(), cvo_probs.end());
        } else {
          // WGS default: product fusion.
          fused = 1.0;
          for (double p : cvo_probs) fused *= p;
        }
      }
      like[pl_idx(j, k)] = fused;
    }
  }

  // Normalise — only when product fusion crossed multiple kept CVOs.
  // Upstream's merge_predictions returns the raw predictions for single-
  // CVO sites and only renormalises after product fusion. For single-CVO
  // sites the FP32 softmax may saturate to exactly 1.0; renormalising by
  // the full-precision sum (=1.0+ε) sneaks the called probability
  // slightly below 1.0, which pushes ptrue_to_bounded_phred away from
  // the 99-cap and gives off-by-many GQ values.
  if (n_kept > 1) {
    double s = 0;
    for (double v : like) s += v;
    if (s <= 0.0) {
      std::fill(like.begin(), like.end(), 1.0 / n_gt);
    } else {
      for (double& v : like) v /= s;
    }
  }
  return like;
}

// Build a VcfHeader from reference contigs.
nucleus::genomics::v1::VcfHeader MakeVcfHeader(
    const std::vector<nucleus::genomics::v1::ContigInfo>& contigs,
    const std::string& sample_name) {
  nucleus::genomics::v1::VcfHeader hdr;
  hdr.set_fileformat("VCFv4.2");

  struct Filt { const char* id; const char* desc; };
  static constexpr Filt kFilters[] = {
      {"PASS",    "All filters passed"},
      {"RefCall", "Most likely homozygous reference"},
      {"LowQual", "Confidence in this variant being real is below threshold"},
      {"NoCall",
       "Site has no call due to low quality (GQ < cnn_homref_call_min_gq)"},
  };
  for (const auto& fi : kFilters) {
    auto* f = hdr.add_filters();
    f->set_id(fi.id);
    f->set_description(fi.desc);
  }
  // Somatic-only filter: GERMLINE for non-somatic variants. Mirrors
  // upstream postprocess_variants.py:2303-2308 + dv_vcf_constants.
  if (absl::GetFlag(FLAGS_process_somatic)) {
    {
      auto* f = hdr.add_filters();
      f->set_id("GERMLINE");
      f->set_description("Non somatic variants");
    }
    // PON filter: variants present in the panel of normals.
    // Only declared if --pon_filtering is set (mirrors upstream behavior:
    // header field appears only when PON filtering is active).
    if (!absl::GetFlag(FLAGS_pon_filtering).empty()) {
      auto* f = hdr.add_filters();
      f->set_id("PON");
      f->set_description("Variant present in panel of normals");
    }
  }

  // INFO fields.
  {
    auto* f = hdr.add_infos();
    f->set_id("END");
    f->set_number("1");
    f->set_type("Integer");
    f->set_description("End position (for symbolic alleles)");
  }

  // FORMAT fields. Order determines per-record column order — keep it
  // matched to upstream's gVCF (GT, GQ, [DP|MIN_DP], AD, VAF, MID, PL).
  // MIN_DP / MED_DP slot in just after GQ since gVCF reference rows
  // emit them in place of DP.
  struct Fmt {
    const char* id;
    const char* num;
    const char* type;
    const char* desc;
  };
  static constexpr Fmt fmts[] = {
      {"GT",     "1", "String",  "Genotype"},
      {"GQ",     "1", "Integer", "Conditional genotype quality"},
      {"MIN_DP", "1", "Integer", "Minimum DP observed within the gVCF block"},
      {"MED_DP", "1", "Integer", "Median DP observed within the gVCF block"},
      {"DP",     "1", "Integer", "Read depth"},
      {"AD",     "R", "Integer", "Allelic depths for ref and alt alleles"},
      {"VAF",    "A", "Float",   "Variant allele fractions"},
      {"MID",    "1", "String",  "Model identifier (small_model | deepvariant)"},
      {"PL",     "G", "Integer", "Phred-scaled genotype likelihoods"},
      // Phase 9 / Step 4c — emitted only when --use_direct_phasing=true;
      // declared unconditionally for consistent header schema.
      {"PS",     "1", "Integer", "Phase set ID (1-based position of block start)"},
  };
  for (const auto& f : fmts) {
    auto* fi = hdr.add_formats();
    fi->set_id(f.id);
    fi->set_number(f.num);
    fi->set_type(f.type);
    fi->set_description(f.desc);
  }

  // Contigs.
  for (const auto& c : contigs) {
    *hdr.add_contigs() = c;
  }
  hdr.add_sample_names(sample_name);
  return hdr;
}

}  // namespace

int RunPostprocessVariants(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  const std::string infile = absl::GetFlag(FLAGS_infile);
  const std::string outfile = absl::GetFlag(FLAGS_output_vcf_outfile);
  const std::string ref_path = absl::GetFlag(FLAGS_ref);

  if (infile.empty() || outfile.empty() || ref_path.empty()) {
    LOG(ERROR) << "Required: --infile, --output_vcf_outfile, --ref";
    return 1;
  }

  // Phase 9 / Step 3 — gVCF output. When --gvcf_outfile is set the
  // make_examples stage must have produced a non-variant Variant
  // TFRecord (one homref row per genomic position) at the path passed
  // via --nonvariant_site_tfrecord_path. After the standard variant
  // post-processing finishes (haplotype resolution + somatic GERMLINE
  // reclassification), `nucleus::MergeAndWriteVariantsAndNonVariants`
  // walks the variant + non-variant streams in lockstep, writes the
  // VCF stream, and writes the gVCF stream with each variant
  // converted to its `<NON_REF>`-extended form via TransfromToGvcf.
  const std::string gvcf_outfile = absl::GetFlag(FLAGS_gvcf_outfile);
  const std::string nonvariant_path =
      absl::GetFlag(FLAGS_nonvariant_site_tfrecord_path);
  if (!gvcf_outfile.empty() && nonvariant_path.empty()) {
    LOG(ERROR) << "--gvcf_outfile=" << gvcf_outfile
               << " requires --nonvariant_site_tfrecord_path to be set.";
    return 1;
  }

  // ── Open reference for contig order ───────────────────────────────────────
  auto ref_or = nucleus::IndexedFastaReader::FromFile(
      ref_path, absl::StrCat(ref_path, ".fai"));
  CHECK(ref_or.ok()) << "Failed to open reference: " << ref_path;
  auto ref_reader = std::move(ref_or.ValueOrDie());
  const auto& contigs = ref_reader->Contigs();

  std::map<std::string, int> contig_to_pos;
  for (int i = 0; i < static_cast<int>(contigs.size()); ++i) {
    contig_to_pos[contigs[i].name()] = i;
  }

  // ── Read all CallVariantsOutput protos ────────────────────────────────────
  const std::vector<std::string> shard_paths = ExpandShards(infile);
  std::vector<CallVariantsOutput> cvo_list;
  for (const auto& path : shard_paths) {
    auto reader = TFRecordReader::New(path);
    if (!reader) {
      LOG(WARNING) << "Cannot open shard: " << path;
      continue;
    }
    while (reader->GetNext()) {
      CallVariantsOutput cvo;
      if (!cvo.ParseFromString(reader->record())) {
        LOG(WARNING) << "Failed to parse CVO proto in " << path;
        continue;
      }
      cvo_list.push_back(std::move(cvo));
    }
    reader->Close();
  }
  LOG(INFO) << "Read " << cvo_list.size() << " CallVariantsOutput protos.";

  // ── Group CVOs by site key (chrom, pos, ref, alts) ────────────────────────
  // The variant proto is identical for all CVOs of the same site under
  // ADD_HET_ALT_IMAGES; only the alt_allele_indices differ.
  using SiteKey = std::tuple<std::string, int64_t, std::string, std::string>;
  std::map<SiteKey, std::vector<const CallVariantsOutput*>> groups;
  for (const auto& cvo : cvo_list) {
    if (!cvo.has_variant()) continue;
    const auto& v = cvo.variant();
    SiteKey k{v.reference_name(), v.start(), v.reference_bases(),
              absl::StrJoin(v.alternate_bases(), ",")};
    groups[k].push_back(&cvo);
  }
  LOG(INFO) << "Grouped into " << groups.size() << " unique sites.";

  // ── Sort sites by genomic coordinate ──────────────────────────────────────
  std::vector<SiteKey> ordered_keys;
  ordered_keys.reserve(groups.size());
  for (const auto& [k, _] : groups) ordered_keys.push_back(k);
  std::sort(ordered_keys.begin(), ordered_keys.end(),
            [&contig_to_pos](const SiteKey& a, const SiteKey& b) {
              const int pa = contig_to_pos.count(std::get<0>(a))
                                 ? contig_to_pos.at(std::get<0>(a))
                                 : INT_MAX;
              const int pb = contig_to_pos.count(std::get<0>(b))
                                 ? contig_to_pos.at(std::get<0>(b))
                                 : INT_MAX;
              if (pa != pb) return pa < pb;
              if (std::get<1>(a) != std::get<1>(b))
                return std::get<1>(a) < std::get<1>(b);
              if (std::get<2>(a) != std::get<2>(b))
                return std::get<2>(a) < std::get<2>(b);
              return std::get<3>(a) < std::get<3>(b);
            });

  // ── Open VCF writer ───────────────────────────────────────────────────────
  std::string sample_name = absl::GetFlag(FLAGS_sample_name);
  if (sample_name.empty()) sample_name = "SAMPLE";
  auto hdr = MakeVcfHeader(contigs, sample_name);
  nucleus::genomics::v1::VcfWriterOptions wr_opts;
  // Tell the writer to read PL from VariantCall.info instead of from the
  // (Float-typed) genotype_likelihood field, which lets us write Integer PL.
  wr_opts.set_retrieve_gl_and_pl_from_info_map(true);
  // Mirror upstream: print QUAL to 1 decimal (e.g. 39.4, not 39.3745).
  wr_opts.set_round_qual_values(true);
  auto writer_or = nucleus::VcfWriter::ToFile(outfile, hdr, wr_opts);
  CHECK(writer_or.ok()) << "Failed to open VCF output: " << outfile;
  auto vcf_writer = std::move(writer_or.ValueOrDie());

  const double qual_filter = absl::GetFlag(FLAGS_qual_filter);
  const double homref_min_gq = absl::GetFlag(FLAGS_cnn_homref_call_min_gq);

  // Sex-chromosome haploid calling config.
  const std::set<std::string> haploid_contigs =
      ParseHaploidContigs(absl::GetFlag(FLAGS_haploid_contigs));
  ParRegions par_regions;
  if (const std::string par_bed = absl::GetFlag(FLAGS_par_regions_bed);
      !par_bed.empty()) {
    std::string err;
    if (!LoadParRegions(par_bed, &par_regions, &err)) {
      LOG(ERROR) << err;
      return 1;
    }
  }

  int written = 0;
  int refcall = 0;
  int nocall = 0;
  // Phase 5.5d/4 — buffer variants for haplotype resolution.
  std::vector<Variant> variants_buffer;
  variants_buffer.reserve(ordered_keys.size());

  for (const auto& key : ordered_keys) {
    const auto& cvos = groups[key];
    Variant variant = cvos.front()->variant();
    const int orig_n_alts = variant.alternate_bases_size();
    const int orig_n_gt = NumDiploidGenotypes(orig_n_alts);

    // Compute alt-pruning set on the ORIGINAL alts (CVOs still reference
    // them by index). We do the actual pruning AFTER picking the
    // best genotype.
    const auto alts_to_remove = AltsToRemove(cvos, qual_filter);

    // Combine likelihoods over the ORIGINAL alt list, skipping CVOs
    // whose alt-set touches a pruned allele (mirrors upstream
    // postprocess_variants.py:merge_predictions step "is_for_pruned_allele:
    // continue"). After the call, `like[g]` for any genotype that
    // includes a pruned alt is still 1.0 (multiplicative identity since
    // every kept CVO sees `overlap=0` for pruned-alt-only genotypes — but
    // those genotypes are masked out below in any case).
    auto like = CombineLikelihoods(cvos, orig_n_alts, alts_to_remove);

    // Phase 8 / Tier 4 — temperature scaling (Guo et al. ICML 2017).
    // Applies before argmax + QUAL/GQ/PL computation. Off by default
    // (T=1.0 trivially preserves the baseline). When opt-in via
    // --enable_temp_scaling and a non-unit T, recalibrates the
    // softmax probabilities to improve expected calibration error
    // (~5-10× ECE reduction in the original CV literature). Effect
    // on F1: typically +0.02-0.10 % when T is fit on a held-out set;
    // depends on whether the baseline model is over- or under-confident
    // at borderline GQ=20 / QUAL=1 thresholds.
    static const bool kEnableTempScaling = absl::GetFlag(FLAGS_enable_temp_scaling);
    static const double kTempScalingT = absl::GetFlag(FLAGS_temp_scaling_T);
    if (kEnableTempScaling && kTempScalingT > 0.0 && kTempScalingT != 1.0) {
      const double inv_T = 1.0 / kTempScalingT;
      double sum = 0.0;
      for (size_t i = 0; i < like.size(); ++i) {
        // Pow on probabilities — avoid log(0) by clipping at the
        // same floor used for PL (1.25e-10).
        const double p = std::max(like[i], 1.25e-10);
        like[i] = std::pow(p, inv_T);
        sum += like[i];
      }
      if (sum > 0.0) {
        for (double& v : like) v /= sum;
      }
    }

    // Mask out genotypes whose alleles are in alts_to_remove. Setting
    // their likelihood to 0 makes them not selectable as argmax.
    if (!alts_to_remove.empty()) {
      for (int k = 0; k <= orig_n_alts; ++k) {
        for (int j = 0; j <= k; ++j) {
          const std::string a1 =
              (j == 0) ? "" : variant.alternate_bases(j - 1);
          const std::string a2 =
              (k == 0) ? "" : variant.alternate_bases(k - 1);
          if ((!a1.empty() && alts_to_remove.count(a1)) ||
              (!a2.empty() && alts_to_remove.count(a2))) {
            like[k * (k + 1) / 2 + j] = 0.0;
          }
        }
      }
      // Renormalise.
      double s = 0;
      for (double v : like) s += v;
      if (s > 0.0) for (double& v : like) v /= s;
    }

    // Now physically prune the variant (renumbering alts). Preserve all
    // other fields — VariantCall.info contains DP/AD/VAF set in
    // make_examples; we must NOT throw them away by replacing the proto.
    if (!alts_to_remove.empty()) {
      // Compute which original alt indices survive — index ranges from 0
      // (first alt) to n_alts-1.
      std::vector<bool> keep_alt(orig_n_alts, false);
      {
        const auto& orig_alts = variant.alternate_bases();
        for (int i = 0; i < orig_alts.size(); ++i) {
          keep_alt[i] = !alts_to_remove.count(orig_alts.Get(i));
        }
      }
      google::protobuf::RepeatedPtrField<std::string> kept_alts;
      for (const auto& a : variant.alternate_bases()) {
        if (!alts_to_remove.count(a)) *kept_alts.Add() = a;
      }
      *variant.mutable_alternate_bases() = std::move(kept_alts);

      // Mirror upstream's AlleleRemapper.reindex_allele_indexed_fields for
      // _ALT_ALLELE_INDEXED_FORMAT_FIELDS = {("AD", true), ("VAF", false),
      // ("MF", true), ("MD", true)}. AD/MF/MD have a ref entry at index 0
      // (ref_is_zero=true) so keep [0] + the kept alt slots. VAF has no ref
      // entry (ref_is_zero=false) so it just gets the kept alt slots.
      for (auto& call : *variant.mutable_calls()) {
        auto* info = call.mutable_info();
        for (const auto& field_info :
             {std::make_pair(std::string("AD"), true),
              std::make_pair(std::string("VAF"), false),
              std::make_pair(std::string("MF"), true),
              std::make_pair(std::string("MD"), true)}) {
          auto it = info->find(field_info.first);
          if (it == info->end()) continue;
          ::nucleus::genomics::v1::ListValue kept;
          const bool ref_is_zero = field_info.second;
          const auto& vals = it->second.values();
          for (int i = 0; i < vals.size(); ++i) {
            bool keep;
            if (ref_is_zero && i == 0) {
              keep = true;  // always keep the ref entry
            } else {
              const int orig_alt = ref_is_zero ? (i - 1) : i;
              keep = (orig_alt < orig_n_alts) ? keep_alt[orig_alt] : false;
            }
            if (keep) *kept.add_values() = vals.Get(i);
          }
          *it->second.mutable_values() = std::move(*kept.mutable_values());
        }
      }
    }

    const int n_alts = variant.alternate_bases_size();
    if (n_alts == 0) continue;
    const int n_gt = NumDiploidGenotypes(n_alts);

    // After pruning, remap the original-index likelihood vector down to
    // the new alt indexing. (Genotype (j, k) on pruned alts maps back to
    // (j', k') on the original alts where j', k' are the original
    // positions of the j-th and k-th non-pruned alts.)
    std::vector<int> new_to_orig(n_alts + 1);
    new_to_orig[0] = 0;
    {
      int new_pos = 1;
      for (int orig = 0; orig < orig_n_alts; ++orig) {
        if (!alts_to_remove.count(variant.alternate_bases().Get(
                std::min(new_pos - 1, n_alts - 1)))) {
          // Find the original index of variant.alternate_bases(new_pos - 1)
          // in the source CVO's alt list.
          // Since `variant` post-prune lists alts in original order, the
          // mapping for new index i is the i-th surviving original index.
        }
      }
      // Simpler reconstruction: walk pruned alts and find each in the
      // first CVO's alt list.
      const auto& orig_alts = cvos.front()->variant().alternate_bases();
      int n = 1;
      for (int i = 0; i < n_alts; ++i) {
        for (int oi = 0; oi < orig_alts.size(); ++oi) {
          if (orig_alts.Get(oi) == variant.alternate_bases(i)) {
            new_to_orig[n++] = oi + 1;
            break;
          }
        }
      }
    }
    std::vector<double> like_pruned(n_gt, 0.0);
    for (int k = 0; k <= n_alts; ++k) {
      for (int j = 0; j <= k; ++j) {
        const int oj = new_to_orig[j];
        const int ok = new_to_orig[k];
        const int new_idx = k * (k + 1) / 2 + j;
        const int orig_idx =
            std::max(oj, ok) * (std::max(oj, ok) + 1) / 2 + std::min(oj, ok);
        if (orig_idx < orig_n_gt) {
          like_pruned[new_idx] = like[orig_idx];
        }
      }
    }
    // Renormalise — but only when alts were actually pruned (the masked
    // genotypes leave the vector summing to <1). For non-pruned single-CVO
    // sites the FP32 saturation in the small_model output already means
    // predictions[0] == 1.0 exactly; renormalising by sum=1.0+ε would push
    // it below 1, which then makes ptrue_to_bounded_phred miss the 99-cap
    // and emit GQ=78 instead of 99 for very-confident homref calls.
    if (!alts_to_remove.empty()) {
      double sp = 0;
      for (double v : like_pruned) sp += v;
      if (sp > 0.0) for (double& v : like_pruned) v /= sp;
    }
    like = std::move(like_pruned);

    // Haploid correction: on a haploid contig outside the PAR, disallow
    // heterozygous genotypes before the genotype/QUAL/GQ/GL are derived from
    // `like` (mirrors merge_predictions applying
    // correct_nonautosome_probabilities to the merged probabilities).
    if (IsHaploidVariant(variant, haploid_contigs, par_regions)) {
      CorrectNonautosomeProbabilities(&like, n_alts);
    }

    // argmax genotype.
    int best = 0;
    for (int i = 1; i < n_gt; ++i) {
      if (like[i] > like[best]) best = i;
    }
    auto [j, k] = GenotypeFromPLIndex(best, n_alts);

    // QUAL = phred-scale of P(non-ref).
    //
    // Upstream's formula:
    //   qual = ptrue_to_bounded_phred(min(sum(predictions[1:]), 1.0))
    //        = phred(1 - sum(predictions[1:]))
    // *not* phred(predictions[0]) — these only agree when the prediction
    // vector sums to exactly 1.0, which it doesn't quite under FP32. Using
    // predictions[0] directly drifts QUAL by up to ~0.1 (e.g. 54.1 vs 54).
    double sum_alt = 0.0;
    for (int i = 1; i < n_gt; ++i) sum_alt += like[i];
    if (sum_alt > 1.0) sum_alt = 1.0;
    const double err_for_qual = std::max(1.0 - sum_alt, 0.0);
    double qual = (err_for_qual >= 1.0) ? 0.0
                                        : std::min(-10.0 * std::log10(err_for_qual),
                                                   static_cast<double>(kMaxPhred));
    // Mirror upstream's compute_quals: rounded_qual = round(qual, 7)
    // (postprocess_variants.py:645, _QUAL_PRECISION=7). The VCF writer
    // then rounds to 1 decimal via set_round_qual_values; this 7-decimal
    // pre-round normalises sub-ULP drift between us and Docker so the
    // 1-decimal write boundary doesn't flip QUAL by 0.1 on borderline
    // values.
    qual = std::round(qual * 1e7) / 1e7;

    // Set up the VariantCall.
    if (variant.calls_size() == 0) variant.add_calls();
    auto* call = variant.mutable_calls(0);
    call->set_call_set_name(sample_name);
    call->clear_genotype();
    call->add_genotype(j);
    call->add_genotype(k);

    // Propagate MID from any of the source CVOs. If at least one CVO in
    // this site's group was tagged as a small_model hit, use that;
    // otherwise fall back to deepvariant. (Both tags are set upstream of
    // postprocess: small_model in make_examples_main.cc, deepvariant in
    // call_variants_main.cc.)
    std::string mid;
    for (const auto* cvo : cvos) {
      for (const auto& src_call : cvo->variant().calls()) {
        auto it = src_call.info().find("MID");
        if (it != src_call.info().end() && it->second.values_size() > 0) {
          const std::string& v = it->second.values(0).string_value();
          if (v == "small_model") { mid = v; break; }
          if (mid.empty()) mid = v;
        }
      }
      if (mid == "small_model") break;
    }
    if (!mid.empty()) {
      nucleus::SetInfoField("MID", mid, call);
    }

    // GQ — mirror of postprocess_variants.py:compute_quals's
    //   gq = round(ptrue_to_bounded_phred(predictions[prediction_index]))
    // i.e. phred(1 - P_called), bounded. Different from "second-best
    // probability"; matters at the cnn_homref_call_min_gq=20 boundary.
    const double p_called = like[best];
    int gq;
    if (p_called >= 1.0) {
      gq = kMaxPhred;
    } else {
      // Mirror upstream's ptrue_to_bounded_phred: floor at 1.25e-10 (so
      // max phred is -10*log10(1.25e-10) = 99.0309) and round-to-even
      // (np.around) — std::round would split half-integer ties the wrong
      // way (35.5 → 36 instead of 36).
      const double err = std::max(1.0 - p_called, 1.25e-10);
      gq = static_cast<int>(std::nearbyint(-10.0 * std::log10(err)));
      gq = std::min(std::max(gq, 0), kMaxPhred);
    }
    nucleus::SetInfoField("GQ", gq, call);

    // GL = log10 likelihood per genotype, capped at log10(1.25e-10)
    // (mirrors upstream's perror_to_bounded_log10_perror in
    // genomics_math.py:106). Used both as the PL source and by the
    // haplotype resolver (`MaybeResolveConflictingVariants`).
    std::vector<double> gls(n_gt);
    double max_gl = -std::numeric_limits<double>::infinity();
    for (int i = 0; i < n_gt; ++i) {
      gls[i] = std::log10(std::max(like[i], 1.25e-10));
      if (gls[i] > max_gl) max_gl = gls[i];
    }
    call->clear_genotype_likelihood();
    for (double gl : gls) call->add_genotype_likelihood(gl);

    // PL = phred-scaled likelihoods. Mirrors upstream's exact flow in
    // vcf_conversion.cc:1215-1232:
    //   1. ZeroShiftLikelihoods: subtract max log10 (zero-shift).
    //   2. std::transform(..., Log10PErrorToPhred) into vector<int> →
    //      double→int via implicit narrowing = TRUNCATION (NOT
    //      std::round; the writer uses `Log10PErrorToPhred` which
    //      returns double, then `std::transform` to vector<int>).
    // Operating in LOG-space (subtract max log10 before phred) is
    // structurally different from the older PHRED-space approach
    // (compute phred[i], subtract min phred): for non-saturated
    // probabilities like=[0.6, 0.4] log-space gives PL=[0,1] (correct,
    // matches Docker), phred-space gave [0,1] too here but in general
    // diverges by 1 unit at rounding boundaries.
    std::vector<int> pl(n_gt);
    for (int i = 0; i < n_gt; ++i) {
      const double phred = -10.0 * (gls[i] - max_gl);
      int p = static_cast<int>(phred);  // truncation (matches writer).
      pl[i] = std::min(std::max(p, 0), kMaxPhred);
    }
    nucleus::SetInfoField("PL", pl, call);

    variant.set_quality(qual);

    // QUAL filter: low-confidence variants become RefCall.
    if (best == 0 || qual < qual_filter) {
      variant.add_filter("RefCall");
      ++refcall;
    } else {
      variant.add_filter("PASS");
    }

    // Mirror postprocess_variants.py:uncall_homref_gt_if_lowqual.
    // CNN RefCalls with GQ < cnn_homref_call_min_gq become "./.": NoCall.
    if (variant.filter_size() == 1 && variant.filter(0) == "RefCall" &&
        gq < homref_min_gq) {
      variant.clear_filter();
      variant.add_filter("NoCall");
      call->clear_genotype();
      call->add_genotype(-1);
      call->add_genotype(-1);
      ++nocall;
    }

    // Buffer for haplotype resolution (Phase 5.5d/4). The pre-resolution
    // GT/FILTER values (incl. uncall_homref_gt_if_lowqual above) are
    // applied first, then `MaybeResolveConflictingVariants` may rewrite
    // overlapping calls and recompute FILTER — matching upstream's
    // postprocess_variants.run_postprocess_variants_on_region order
    // (per-variant add_call_to_variant → maybe_resolve_conflicting_variants).
    // Simplify alleles (strip common postfix) — upstream
    // merge_predictions:simplify_variant_alleles. Without this our
    // tandem-repeat substitutions retain a long shared suffix which
    // makes them spuriously overlap with neighbouring variants in
    // haplotype resolution.
    SimplifyVariantAlleles(&variant);
    variants_buffer.push_back(std::move(variant));
  }

  LOG(INFO) << "Applying haplotype resolution to "
            << variants_buffer.size() << " variants ...";
  ::deepvariant::MaybeResolveConflictingVariants(&variants_buffer, qual_filter);

  const bool process_somatic = absl::GetFlag(FLAGS_process_somatic);
  // Apply the somatic GERMLINE-reclassification mutation in-place
  // (Phase 9 / Step 3 — needed because the gVCF merge path consumes
  // variants_buffer through a TFRecord round-trip rather than the
  // direct VcfWriter::WriteSomatic path).
  if (process_somatic) {
    for (auto& v : variants_buffer) {
      if (v.calls_size() == 0) continue;
      // Mirror nucleus/io/vcf_writer.cc::WriteSomatic: any non-{0/0,
      // 1/1, ./.} GT (i.e. heterozygous) gets reclassified as
      // GERMLINE 0/0. The biological assumption: a het call in a
      // tumor+normal pair is most likely a germline variant the
      // patient inherited (hom-alt would suggest LOH = somatic event).
      auto* call = v.mutable_calls(0);
      const auto& g = call->genotype();
      const bool is_homref = (g.size() == 2 && g.Get(0) == 0 && g.Get(1) == 0);
      const bool is_homalt = (g.size() == 2 && g.Get(0) == 1 && g.Get(1) == 1);
      const bool is_nocall = (g.size() == 2 && g.Get(0) < 0 && g.Get(1) < 0);
      if (!is_homref && !is_homalt && !is_nocall) {
        call->clear_genotype();
        call->add_genotype(0);
        call->add_genotype(0);
        if (v.filter_size() > 0) {
          v.clear_filter();
          v.add_filter("GERMLINE");
        }
      }
    }
  }

  // PON filtering pass — Phase 9 step (--pon_filtering, somatic only).
  // Mirrors upstream postprocess_variants.py:filter_pon. For each PASS
  // variant, look up (CHROM,POS,REF,ALT) in the PON VCF; if present,
  // remove PASS and add FILTER=PON.
  const std::string pon_path = absl::GetFlag(FLAGS_pon_filtering);
  if (process_somatic && !pon_path.empty()) {
    nucleus::genomics::v1::VcfReaderOptions pon_opts;
    auto pon_or = nucleus::VcfReader::FromFile(pon_path, pon_opts);
    CHECK(pon_or.ok()) << "PON open failed: " << pon_path;
    auto pon_reader = std::move(pon_or.ValueOrDie());

    int pon_hits = 0;
    for (auto& v : variants_buffer) {
      if (v.filter_size() == 0) continue;
      // Only check PASS variants (untouched by GERMLINE pass).
      bool has_pass = false;
      for (const auto& f : v.filter()) {
        if (f == "PASS") { has_pass = true; break; }
      }
      if (!has_pass) continue;

      // Build query range covering this site (1bp at v.start()).
      nucleus::genomics::v1::Range range;
      range.set_reference_name(v.reference_name());
      range.set_start(v.start());
      range.set_end(v.start() + 1);

      auto iter_or = pon_reader->Query(range);
      if (!iter_or.ok()) continue;
      auto iter = std::move(iter_or.ValueOrDie());

      bool match = false;
      nucleus::genomics::v1::Variant pv;
      while (true) {
        auto next_or = iter->Next(&pv);
        if (!next_or.ok() || !next_or.ValueOrDie()) break;
        if (pv.start() != v.start()) continue;
        if (pv.reference_bases() != v.reference_bases()) continue;
        // Match if any of OUR alts equal any PON alt (allow multi-allelic).
        for (const auto& our_alt : v.alternate_bases()) {
          for (const auto& pon_alt : pv.alternate_bases()) {
            if (our_alt == pon_alt) { match = true; break; }
          }
          if (match) break;
        }
        if (match) break;
      }
      if (match) {
        v.clear_filter();
        v.add_filter("PON");
        ++pon_hits;
      }
    }
    LOG(INFO) << "PON filter: " << pon_hits << " variants tagged PON.";
  }

  if (gvcf_outfile.empty()) {
    // ── Direct VCF write (no gVCF). ────────────────────────────────────
    for (const auto& v : variants_buffer) {
      auto status = vcf_writer->Write(v);
      if (!status.ok()) {
        LOG(WARNING) << "Failed to write variant at "
                     << v.reference_name() << ":" << v.start() << " — " << status;
      } else {
        ++written;
      }
    }
  } else {
    // ── gVCF merge path (Phase 9 / Step 3). ────────────────────────────
    // Round-trip variants_buffer through a temp TFRecord so we can hand
    // it to nucleus::MergeAndWriteVariantsAndNonVariants alongside the
    // sharded non-variant TFRecord written by make_examples.
    const std::string tmp_var_tfrecord =
        absl::StrCat(outfile, ".variants.tmp.tfrecord");
    {
      auto w = TFRecordWriter::New(tmp_var_tfrecord);
      CHECK(w) << "Cannot open temp variant TFRecord: " << tmp_var_tfrecord;
      for (const auto& v : variants_buffer) {
        std::string serialized;
        v.SerializeToString(&serialized);
        if (!w->WriteRecord(serialized)) {
          LOG(ERROR) << "Failed to write temp variant TFRecord: "
                     << tmp_var_tfrecord;
          std::remove(tmp_var_tfrecord.c_str());  // don't leave a partial file
          return 1;
        }
      }
      if (!w->Close()) {
        LOG(ERROR) << "Failed to flush temp variant TFRecord: "
                   << tmp_var_tfrecord;
        std::remove(tmp_var_tfrecord.c_str());  // don't leave a partial file
        return 1;
      }
    }

    // Open the variant + non-variant readers and a second VcfWriter
    // for the gVCF stream (same options + header as the main VCF).
    absl::flat_hash_map<std::string, uint32_t> contig_index_map;
    for (uint32_t i = 0; i < contigs.size(); ++i) {
      contig_index_map[contigs[i].name()] = i;
    }
    auto var_reader = nucleus::VariantReader::Open(
        tmp_var_tfrecord, /*compression=*/"", contig_index_map);
    CHECK(var_reader) << "Cannot open temp variant TFRecord for read: "
                      << tmp_var_tfrecord;

    const std::vector<std::string> nv_shards = ExpandShards(nonvariant_path);
    auto nv_reader =
        nucleus::ShardedVariantReader::Open(nv_shards, contig_index_map);
    CHECK(nv_reader) << "Cannot open non-variant TFRecord shards: "
                     << nonvariant_path;

    auto gvcf_writer_or = nucleus::VcfWriter::ToFile(gvcf_outfile, hdr, wr_opts);
    CHECK(gvcf_writer_or.ok())
        << "Failed to open gVCF output: " << gvcf_outfile;
    auto gvcf_writer = std::move(gvcf_writer_or.ValueOrDie());

    // Empty `ranges` = whole-genome (nucleus::RangesContainVariant is only
    // applied when ranges is non-empty; the make_examples region filter
    // already restricted the per-position rows to the user's --regions).
    std::vector<nucleus::genomics::v1::Range> ranges;
    nucleus::MergeAndWriteVariantsAndNonVariants(
        /*only_keep_pass=*/false, var_reader.get(), nv_reader.get(),
        vcf_writer.get(), gvcf_writer.get(), *ref_reader, ranges,
        /*process_somatic=*/process_somatic);

    written = static_cast<int>(variants_buffer.size());
    std::remove(tmp_var_tfrecord.c_str());
    LOG(INFO) << "gVCF written to " << gvcf_outfile;
  }

  // Recount filter classes after haplotype resolution (the per-variant
  // counts above may be stale where resolution rewrote GT to 0/0).
  refcall = 0; nocall = 0;
  int pass = 0;
  for (const auto& v : variants_buffer) {
    if (v.filter_size() == 0) continue;
    const std::string& f = v.filter(0);
    if (f == "RefCall") ++refcall;
    else if (f == "NoCall") ++nocall;
    else if (f == "PASS") ++pass;
  }
  LOG(INFO) << "postprocess_variants done: " << written << " VCF lines"
            << " (" << refcall << " RefCall, "
            << nocall << " NoCall, "
            << pass << " PASS).";
  return 0;
}

}  // namespace deepvariant

            