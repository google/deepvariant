// Phase 5.5d/4 — haplotype-resolution port; see haplotypes.h.

#include "deepvariant/native/haplotypes.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <string>
#include <utility>
#include <vector>

#include "absl/log/log.h"
#include "third_party/nucleus/protos/variants.pb.h"
#include "third_party/nucleus/util/utils.h"

namespace deepvariant {

namespace {

using nucleus::genomics::v1::Variant;
using nucleus::genomics::v1::VariantCall;

constexpr int kPloidy = 2;
constexpr int kMaxOverlappingVariantsToResolve = 12;

// VCF "G" PL ordering: F(j, k) = k*(k+1)/2 + j   (j ≤ k).
inline int GenotypeLikelihoodIndex(int a, int b) {
  if (a > b) std::swap(a, b);
  return b * (b + 1) / 2 + a;
}

// Number of non-ref alleles in the called genotype (0/0 → 0, 0/1 → 1,
// 1/1 or 1/2 → 2). Negative-genotype slots (./.) count as 0.
int NonrefGenotypeCount(const Variant& v) {
  if (v.calls_size() == 0) return 0;
  int n = 0;
  for (int g : v.calls(0).genotype()) if (g > 0) ++n;
  return n;
}

// True if the actual genotype calls are compatible — i.e., no reference
// position has more than `ploidy` non-ref alleles claimed across the
// covering variants. Mirrors
// `_VariantCompatibilityCalculator.all_variants_compatible`.
bool AllVariantsCompatible(const std::vector<const Variant*>& variants,
                            const std::vector<int>& nonref_counts) {
  if (variants.empty() || nonref_counts.size() != variants.size()) return true;
  // Find the union start..end of the group.
  int64_t group_start = std::numeric_limits<int64_t>::max();
  int64_t group_end = 0;
  for (const auto* v : variants) {
    group_start = std::min(group_start, (int64_t)v->start());
    group_end = std::max(group_end, (int64_t)v->end());
  }
  if (group_end <= group_start) return true;
  std::vector<int> alts_in_span(group_end - group_start, 0);
  for (size_t i = 0; i < variants.size(); ++i) {
    const Variant* v = variants[i];
    const int cnt = nonref_counts[i];
    for (int64_t pos = v->start(); pos < v->end(); ++pos) {
      alts_in_span[pos - group_start] += cnt;
    }
  }
  for (int v : alts_in_span) if (v > kPloidy) return false;
  return true;
}

// `allele_indices_with_num_alts(variant, num_alts, ploidy=2)` from
// nucleus/util/variant_utils.py.
std::vector<std::pair<int, int>> AlleleIndicesWithNumAlts(
    const Variant& v, int num_alts) {
  const int max_alt = v.alternate_bases_size();
  std::vector<std::pair<int, int>> out;
  if (num_alts == 0) {
    out.emplace_back(0, 0);
  } else if (num_alts == 1) {
    for (int i = 1; i <= max_alt; ++i) out.emplace_back(0, i);
  } else {  // num_alts == 2
    for (int i = 1; i <= max_alt; ++i) {
      for (int j = i; j <= max_alt; ++j) out.emplace_back(i, j);
    }
  }
  return out;
}

// Cartesian product of per-variant allele-indices configurations. The
// result is a list-of-lists of (a, b) pairs, one tuple per variant.
std::vector<std::vector<std::pair<int, int>>>
GetAllAlleleIndicesConfigurations(
    const std::vector<const Variant*>& variants,
    const std::vector<int>& nonref_count_config) {
  std::vector<std::vector<std::pair<int, int>>> per_variant;
  per_variant.reserve(variants.size());
  for (size_t i = 0; i < variants.size(); ++i) {
    per_variant.push_back(
        AlleleIndicesWithNumAlts(*variants[i], nonref_count_config[i]));
  }
  // Iterative Cartesian product.
  std::vector<std::vector<std::pair<int, int>>> out;
  out.push_back({});
  for (const auto& opts : per_variant) {
    std::vector<std::vector<std::pair<int, int>>> next;
    next.reserve(out.size() * opts.size());
    for (const auto& cfg : out) {
      for (const auto& o : opts) {
        auto cfg2 = cfg;
        cfg2.push_back(o);
        next.push_back(std::move(cfg2));
      }
    }
    out = std::move(next);
  }
  return out;
}

// Reads call.genotype_likelihood at the given (a, b) genotype slot.
double GenotypeLikelihood(const VariantCall& call, std::pair<int, int> ab) {
  const int idx = GenotypeLikelihoodIndex(ab.first, ab.second);
  if (idx < 0 || idx >= call.genotype_likelihood_size()) {
    return -std::numeric_limits<double>::infinity();
  }
  return call.genotype_likelihood(idx);
}

// Joint log10-likelihood = sum of per-variant GLs at the given alleles.
double AlleleIndicesConfigurationLikelihood(
    const std::vector<const Variant*>& variants,
    const std::vector<std::pair<int, int>>& cfg) {
  double total = 0.0;
  for (size_t i = 0; i < variants.size(); ++i) {
    if (variants[i]->calls_size() == 0) continue;
    total += GenotypeLikelihood(variants[i]->calls(0), cfg[i]);
  }
  return total;
}

// log10(sum(10^x)) computed in a numerically-stable way (mirror of
// genomics_math.log10sumexp).
double Log10SumExp(const std::vector<double>& xs) {
  if (xs.empty()) return -std::numeric_limits<double>::infinity();
  double m = -std::numeric_limits<double>::infinity();
  for (double x : xs) m = std::max(m, x);
  double s = 0.0;
  for (double x : xs) s += std::pow(10.0, x - m);
  return m + std::log10(s);
}

// Subtract log10sumexp so 10^x sums to 1 — mirror of
// genomics_math.normalize_log10_probs.
std::vector<double> NormalizeLog10Probs(std::vector<double> v) {
  if (v.empty()) return v;
  const double lse = Log10SumExp(v);
  for (double& x : v) x = std::min(x - lse, 0.0);
  return v;
}

// Per-variant aggregator: stores the joint LLs that touched each
// genotype slot, then `Scaled()` returns log10 marginals (subtract-max
// approximation) and `MostLikelyAllele()` returns argmax allele indices.
struct LikelihoodAggregator {
  // genotype_likelihood_index → list of LLs.
  std::vector<std::vector<double>> bucket;
  int n_alts = 0;

  static int NumLikelihoodSlots(int num_alts) {
    return GenotypeLikelihoodIndex(num_alts, num_alts) + 1;
  }

  explicit LikelihoodAggregator(int num_alts_) : n_alts(num_alts_) {
    bucket.assign(NumLikelihoodSlots(num_alts_), {});
  }

  void Add(std::pair<int, int> ab, double ll) {
    int idx = GenotypeLikelihoodIndex(ab.first, ab.second);
    if (idx >= 0 && idx < (int)bucket.size()) bucket[idx].push_back(ll);
  }

  std::vector<double> Scaled() const {
    std::vector<double> out;
    out.reserve(bucket.size());
    for (const auto& v : bucket) {
      out.push_back(v.empty() ? -std::numeric_limits<double>::infinity()
                              : Log10SumExp(v));
    }
    return NormalizeLog10Probs(std::move(out));
  }

  std::pair<int, int> MostLikelyAllele() const {
    auto s = Scaled();
    int argmax = 0;
    double m = s.empty() ? 0.0 : s[0];
    for (int i = 1; i < (int)s.size(); ++i) {
      if (s[i] > m) { m = s[i]; argmax = i; }
    }
    // Inverse of GenotypeLikelihoodIndex: walk the (a, b) pairs.
    for (int b = 0; b <= n_alts; ++b) {
      for (int a = 0; a <= b; ++a) {
        if (GenotypeLikelihoodIndex(a, b) == argmax) return {a, b};
      }
    }
    return {0, 0};
  }
};

// Recompute the FILTER field after a genotype change. Mirror of
// dv_vcf_constants.compute_filter_fields:
//   no_call → "NoCall"
//   hom_ref → "RefCall"
//   else if QUAL < min_quality → "LowQual"
//   else → "PASS"
std::string FilterFor(const Variant& v, double min_quality) {
  if (v.calls_size() == 0) return "NoCall";
  const auto& gt = v.calls(0).genotype();
  bool any = gt.size() > 0;
  bool all_neg = true, all_zero = true;
  for (int g : gt) {
    if (g != -1) all_neg = false;
    if (g != 0) all_zero = false;
  }
  if (!any || all_neg) return "NoCall";
  if (all_zero) return "RefCall";
  if (v.quality() < min_quality) return "LowQual";
  return "PASS";
}

// Group `vs` (sorted by start) into contiguous blocks of overlapping
// variants. Returns indices into `vs`.
std::vector<std::vector<size_t>> GroupOverlapping(
    const std::vector<Variant*>& vs) {
  std::vector<std::vector<size_t>> groups;
  if (vs.empty()) return groups;
  std::vector<size_t> cur{0};
  std::string prev_chrom = vs[0]->reference_name();
  int64_t prev_max_end = vs[0]->end();
  for (size_t i = 1; i < vs.size(); ++i) {
    const Variant* v = vs[i];
    if (v->reference_name() != prev_chrom ||
        (int64_t)v->start() >= prev_max_end) {
      groups.push_back(std::move(cur));
      cur = {i};
      prev_chrom = v->reference_name();
      prev_max_end = v->end();
    } else {
      cur.push_back(i);
      prev_max_end = std::max(prev_max_end, (int64_t)v->end());
    }
  }
  groups.push_back(std::move(cur));
  return groups;
}

// Apply genotype + GL update to a variant (and recompute PL info field
// + FILTER). Our VCF writer reads PL from `info["PL"]` (not from
// `genotype_likelihood`), so we must keep PL in sync after rewriting GL.
constexpr int kMaxPhred = 99;
void ApplyAlleleIndicesAndGL(Variant* v, std::pair<int, int> ab,
                              const std::vector<double>& gls,
                              double qual_filter) {
  if (v->calls_size() == 0) return;
  auto* call = v->mutable_calls(0);
  call->clear_genotype();
  call->add_genotype(ab.first);
  call->add_genotype(ab.second);
  call->clear_genotype_likelihood();
  for (double g : gls) call->add_genotype_likelihood(g);

  // Re-derive PL = round(-10 * gl_shifted_to_min0). Subtract max-gl so
  // best genotype gets PL=0.
  if (!gls.empty()) {
    double max_gl = -std::numeric_limits<double>::infinity();
    for (double g : gls) max_gl = std::max(max_gl, g);
    std::vector<int> pl(gls.size());
    for (size_t i = 0; i < gls.size(); ++i) {
      double phred = -10.0 * (gls[i] - max_gl);
      int p = (int)std::nearbyint(phred);
      pl[i] = std::min(std::max(p, 0), kMaxPhred);
    }
    auto* info_map = call->mutable_info();
    auto& pl_field = (*info_map)["PL"];
    pl_field.clear_values();
    for (int p : pl) pl_field.add_values()->set_int_value(p);
  }

  v->clear_filter();
  v->add_filter(FilterFor(*v, qual_filter));
}

// `_resolve_overlapping_variants` from haplotypes.py — takes a list of
// CONTIGUOUS overlapping variants (with non-ref calls) and rewrites
// their genotype + GL where the joint argmax agrees with the marginal
// argmax. If the algorithm punts (>12 variants, or marginals disagree
// with joint), the variants are left unchanged.
void ResolveOverlappingGroup(std::vector<Variant*>& group,
                              double qual_filter) {
  if (group.size() <= 1) return;

  std::vector<const Variant*> consts(group.begin(), group.end());
  std::vector<int> actual_counts;
  actual_counts.reserve(group.size());
  for (const Variant* v : group) actual_counts.push_back(NonrefGenotypeCount(*v));
  if (AllVariantsCompatible(consts, actual_counts)) return;

  if (group.size() > kMaxOverlappingVariantsToResolve) {
    LOG(WARNING) << "haplotypes: punting on " << group.size()
                 << " overlapping variants (> "
                 << kMaxOverlappingVariantsToResolve << ")";
    return;
  }

  // Enumerate compatible nonref-count configurations.
  // 3^N options (each variant independently 0, 1, or 2 non-ref).
  std::vector<std::vector<int>> compatible_count_configs;
  std::vector<int> cfg(group.size(), 0);
  while (true) {
    if (AllVariantsCompatible(consts, cfg)) {
      compatible_count_configs.push_back(cfg);
    }
    int i = (int)group.size() - 1;
    while (i >= 0 && cfg[i] == 2) { cfg[i] = 0; --i; }
    if (i < 0) break;
    ++cfg[i];
  }

  // For each compatible nonref-count config, enumerate allele-index
  // configurations and track joint argmax + per-variant marginals.
  std::vector<LikelihoodAggregator> aggs;
  aggs.reserve(group.size());
  for (const Variant* v : group) {
    aggs.emplace_back(v->alternate_bases_size());
  }
  std::vector<std::pair<int, int>> joint_argmax_cfg;
  double joint_argmax_ll = -std::numeric_limits<double>::infinity();
  for (const auto& nc : compatible_count_configs) {
    for (const auto& ai_cfg : GetAllAlleleIndicesConfigurations(consts, nc)) {
      double ll = AlleleIndicesConfigurationLikelihood(consts, ai_cfg);
      if (ll > joint_argmax_ll) {
        joint_argmax_ll = ll;
        joint_argmax_cfg = ai_cfg;
      }
      for (size_t i = 0; i < group.size(); ++i) aggs[i].Add(ai_cfg[i], ll);
    }
  }
  if (joint_argmax_cfg.empty()) return;  // no compatible config (should not happen)

  // Marginal argmax per variant.
  std::vector<std::pair<int, int>> marginal_cfg;
  marginal_cfg.reserve(group.size());
  for (auto& a : aggs) marginal_cfg.push_back(a.MostLikelyAllele());

  if (marginal_cfg != joint_argmax_cfg) {
    LOG(INFO) << "haplotypes: marginal vs joint disagree at "
              << group[0]->reference_name() << ":" << group[0]->start()
              << " — punting";
    return;
  }

  // Apply: rewrite genotype + GL + recompute filter.
  for (size_t i = 0; i < group.size(); ++i) {
    auto scaled = aggs[i].Scaled();
    ApplyAlleleIndicesAndGL(group[i], joint_argmax_cfg[i], scaled,
                              qual_filter);
  }
}

}  // namespace

void MaybeResolveConflictingVariants(std::vector<Variant>* variants,
                                      double qual_filter) {
  if (!variants || variants->size() < 2) return;

  std::vector<Variant*> ptrs;
  ptrs.reserve(variants->size());
  for (auto& v : *variants) ptrs.push_back(&v);

  int n_groups_total = 0, n_groups_multi = 0, n_groups_resolved = 0;
  // Group all overlapping variants.
  for (auto& group_idx : GroupOverlapping(ptrs)) {
    ++n_groups_total;
    if (group_idx.size() <= 1) continue;
    ++n_groups_multi;
    // Split each group into ref-calls (nonref count == 0) and var-calls.
    // Run resolution only on the var-calls sub-groups (mirror of upstream's
    // _maybe_resolve_mixed_calls).
    std::vector<Variant*> var_calls;
    for (size_t i : group_idx) {
      if (NonrefGenotypeCount(*ptrs[i]) > 0) var_calls.push_back(ptrs[i]);
    }
    if (var_calls.size() <= 1) continue;
    // Re-group the var-calls (some may not actually overlap each other).
    for (auto& sub_idx : GroupOverlapping(var_calls)) {
      if (sub_idx.size() <= 1) continue;
      std::vector<Variant*> sub;
      for (size_t i : sub_idx) sub.push_back(var_calls[i]);
      LOG(INFO) << "haplotypes: resolving " << sub.size()
                << " overlapping variant-calls starting at "
                << sub[0]->reference_name() << ":" << sub[0]->start();
      ResolveOverlappingGroup(sub, qual_filter);
      ++n_groups_resolved;
    }
  }
  LOG(INFO) << "haplotypes: " << n_groups_total << " total groups, "
             << n_groups_multi << " with > 1 variant, "
             << n_groups_resolved << " variant-call sub-groups resolved.";
}

}  // namespace deepvariant
