// Implementation of the 70-feature extractor for the small_model.
// Mirrors deepvariant/small_model/make_small_model_examples.py:FeatureEncoder.

#include "deepvariant/native/small_model_features.h"

#include <algorithm>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"

namespace deepvariant {

using learning::genomics::deepvariant::DeepVariantCall;
using learning::genomics::deepvariant::DeepVariantCall_ReadSupport;
using nucleus::genomics::v1::Variant;

namespace {

// Count alternate_bases entries that are not symbolic/placeholder alleles.
// Mirrors variant_utils.is_multiallelic, which default-excludes "<*>",
// "<NON_REF>", and "." before counting.
int NumNonSymbolicAlts(const Variant& v) {
  int n = 0;
  for (const auto& alt : v.alternate_bases()) {
    if (alt != "<*>" && alt != "<NON_REF>" && alt != ".") ++n;
  }
  return n;
}

// Pull the reads supporting the chosen alt allele indices into a single
// flat vector of ReadSupport pointers. If `sample_filter` is non-empty,
// only reads with matching `sample_name` are retained (mirrors upstream's
// `_filter_by_sample(read_infos, sample_name)`).
std::vector<const DeepVariantCall_ReadSupport*> GetAltReadInfos(
    const DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices,
    const std::string& sample_filter = "") {
  std::vector<const DeepVariantCall_ReadSupport*> out;
  for (int idx : alt_allele_indices) {
    if (idx < 0 || idx >= candidate.variant().alternate_bases_size()) continue;
    const auto& alt_bases = candidate.variant().alternate_bases(idx);
    auto it = candidate.allele_support_ext().find(alt_bases);
    if (it == candidate.allele_support_ext().end()) continue;
    for (const auto& r : it->second.read_infos()) {
      if (!sample_filter.empty() && r.sample_name() != sample_filter) continue;
      out.push_back(&r);
    }
  }
  return out;
}

std::vector<const DeepVariantCall_ReadSupport*> GetRefReadInfos(
    const DeepVariantCall& candidate,
    const std::string& sample_filter = "") {
  std::vector<const DeepVariantCall_ReadSupport*> out;
  for (const auto& r : candidate.ref_support_ext().read_infos()) {
    if (!sample_filter.empty() && r.sample_name() != sample_filter) continue;
    out.push_back(&r);
  }
  return out;
}

int MeanInt(const std::vector<const DeepVariantCall_ReadSupport*>& reads,
            int (*getter)(const DeepVariantCall_ReadSupport&)) {
  if (reads.empty()) return 0;
  int64_t sum = 0;
  for (const auto* r : reads) sum += getter(*r);
  return static_cast<int>(sum / static_cast<int64_t>(reads.size()));
}

int GetMQ(const DeepVariantCall_ReadSupport& r) { return r.mapping_quality(); }
int GetBQ(const DeepVariantCall_ReadSupport& r) {
  return r.average_base_quality();
}
int GetReverseStrand100(const DeepVariantCall_ReadSupport& r) {
  return r.is_reverse_strand() ? 100 : 0;
}

// SNP detection (roughly variant_utils.is_snp): every alt is a single base
// and ref is a single base.
bool IsSnp(const Variant& v, const std::set<std::string>& exclude) {
  if (v.reference_bases().size() != 1) return false;
  bool any_alt = false;
  for (const auto& a : v.alternate_bases()) {
    if (exclude.count(a)) continue;
    if (a.size() != 1) return false;
    any_alt = true;
  }
  return any_alt;
}

bool IsInsertion(const Variant& v, const std::set<std::string>& exclude) {
  bool any = false;
  for (const auto& a : v.alternate_bases()) {
    if (exclude.count(a)) continue;
    if (a.size() <= v.reference_bases().size()) return false;
    any = true;
  }
  return any;
}

bool IsDeletion(const Variant& v, const std::set<std::string>& exclude) {
  bool any = false;
  for (const auto& a : v.alternate_bases()) {
    if (exclude.count(a)) continue;
    if (a.size() >= v.reference_bases().size()) return false;
    any = true;
  }
  return any;
}

// Append 12 BaseFeatures for a particular sample-filter slice of the
// candidate's reads to `features`. Mirrors upstream's
// `FeatureEncoder.encode_base_feature` invoked over the BaseFeature
// enum in declaration order.
//
// IMPORTANT — upstream's _get_total_depth (make_small_model_examples.py:
// 292-296) is ALWAYS unfiltered: it returns
// `len(candidate.ref_support_ext.read_infos) + sum(len(r.read_infos)
// for r in candidate.allele_support_ext.values())` regardless of the
// FeatureEncoder's `sample` arg. Only `ref_read_infos_count` and
// `alt_read_infos_count` (and derivatives) are sample-filtered. So:
//   - total_depth: unfiltered
//   - alt_indices_depth: ref_count_filtered + alt_count_filtered
//   - variant_allele_frequency: 100 * alt_count_filtered / total_depth_unfiltered
//   - alt_indices_variant_allele_frequency: 100 * alt_count_filtered / alt_indices_depth
void AppendBaseFeatures(
    const DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices,
    const std::string& sample_filter,
    std::vector<float>* features) {
  auto ref_reads = GetRefReadInfos(candidate, sample_filter);
  auto alt_reads = GetAltReadInfos(candidate, alt_allele_indices, sample_filter);

  // Upstream invariant: total_depth is ALWAYS unfiltered (across all
  // samples + all alleles), even when computing per-sample features.
  int total_depth = candidate.ref_support_ext().read_infos_size();
  for (const auto& [_, support] : candidate.allele_support_ext()) {
    total_depth += support.read_infos_size();
  }

  const int n_ref = static_cast<int>(ref_reads.size());
  const int n_alt = static_cast<int>(alt_reads.size());
  const int alt_indices_depth = n_ref + n_alt;
  features->push_back(n_ref);                  // num_reads_supports_ref
  features->push_back(n_alt);                  // num_reads_supports_alt
  features->push_back(alt_indices_depth);      // alt_indices_depth
  features->push_back(total_depth);            // total_depth (unfiltered!)
  // VAF: numerator sample-filtered, denominator unfiltered total_depth.
  features->push_back(total_depth > 0 ? (100 * n_alt / total_depth) : 0);
  // alt_indices_VAF: both numerator and denominator sample-filtered.
  features->push_back(alt_indices_depth > 0
                           ? (100 * n_alt / alt_indices_depth)
                           : 0);
  features->push_back(MeanInt(ref_reads, GetMQ));
  features->push_back(MeanInt(alt_reads, GetMQ));
  features->push_back(MeanInt(ref_reads, GetBQ));
  features->push_back(MeanInt(alt_reads, GetBQ));
  features->push_back(MeanInt(ref_reads, GetReverseStrand100));
  features->push_back(MeanInt(alt_reads, GetReverseStrand100));
}

}  // namespace

std::vector<float> EncodeSmallModelFeatures(
    const DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices) {
  std::vector<float> features;
  features.reserve(kSmallModelNumFeatures);

  // Excluded alternates: those NOT in alt_allele_indices.
  std::set<std::string> exclude;
  std::set<int> indices_set(alt_allele_indices.begin(),
                             alt_allele_indices.end());
  for (int i = 0; i < candidate.variant().alternate_bases_size(); ++i) {
    if (!indices_set.count(i)) {
      exclude.insert(candidate.variant().alternate_bases(i));
    }
  }

  // ── BaseFeatures (12) — single sample, no filter ──────────────────────────
  AppendBaseFeatures(candidate, alt_allele_indices, /*sample_filter=*/"",
                      &features);

  // ── VariantFeatures (7) ───────────────────────────────────────────────────
  const auto& v = candidate.variant();
  features.push_back(IsSnp(v, exclude) ? 1 : 0);                 // is_snp
  features.push_back(IsInsertion(v, exclude) ? 1 : 0);           // is_insertion
  features.push_back(IsDeletion(v, exclude) ? 1 : 0);            // is_deletion
  // insertion_length: max(0, max over indices of (alt_len - ref_len))
  int ins_len = 0;
  for (int idx : alt_allele_indices) {
    if (idx < 0 || idx >= v.alternate_bases_size()) continue;
    int d = static_cast<int>(v.alternate_bases(idx).size()) -
            static_cast<int>(v.reference_bases().size());
    ins_len = std::max(ins_len, d);
  }
  features.push_back(std::max(0, ins_len));                      // insertion_length
  // deletion_length: max(0, max over indices of (ref_len - alt_len))
  int del_len = 0;
  for (int idx : alt_allele_indices) {
    if (idx < 0 || idx >= v.alternate_bases_size()) continue;
    int d = static_cast<int>(v.reference_bases().size()) -
            static_cast<int>(v.alternate_bases(idx).size());
    del_len = std::max(del_len, d);
  }
  features.push_back(std::max(0, del_len));                      // deletion_length
  features.push_back(NumNonSymbolicAlts(v) > 1 ? 1 : 0);         // is_multiallelic
  features.push_back(alt_allele_indices.size() > 1 ? 1 : 0);     // is_multiple_alt_alleles

  // ── VAF context (51 features, offsets -25..+25 inclusive) ─────────────────
  const auto& vaf_at_pos = candidate.allele_frequency_at_position();
  const int half = kSmallModelVafContextWindow / 2;  // 25
  for (int o = -half; o <= half; ++o) {
    const int64_t pos = v.start() + o;
    auto it = vaf_at_pos.find(static_cast<int>(pos));
    features.push_back(it != vaf_at_pos.end() ? it->second : 0);
  }

  return features;
}

// Multi-sample (trio / somatic) feature encoder. Mirrors upstream's
// `FeatureEncoder._encode_candidate_feature_dict` insertion order:
//   1. 12 BaseFeatures (combined / target-only — sample_filter="")
//   2. 12 BaseFeatures × N samples, in `sample_order` over `sample_names`
//   3. 7 VariantFeatures
//   4. 51 VAF context features
std::vector<float> EncodeSmallModelFeaturesMultiSample(
    const DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices,
    const std::vector<std::string>& sample_names,
    const std::vector<int>& sample_order) {
  const int n_samples = static_cast<int>(sample_names.size());
  const int total_features = kSmallModelNumFeatures +
      kSmallModelBaseFeaturesPerSample * n_samples;
  std::vector<float> features;
  features.reserve(total_features);

  // Excluded alternates: those NOT in alt_allele_indices.
  std::set<std::string> exclude;
  std::set<int> indices_set(alt_allele_indices.begin(),
                             alt_allele_indices.end());
  for (int i = 0; i < candidate.variant().alternate_bases_size(); ++i) {
    if (!indices_set.count(i)) {
      exclude.insert(candidate.variant().alternate_bases(i));
    }
  }

  // ── BaseFeatures (12) — combined / no sample filter ──────────────────────
  AppendBaseFeatures(candidate, alt_allele_indices, /*sample_filter=*/"",
                      &features);

  // ── Per-sample BaseFeatures (12 × N), in sample_order order ──────────────
  for (int idx : sample_order) {
    if (idx < 0 || idx >= n_samples) continue;
    AppendBaseFeatures(candidate, alt_allele_indices,
                        sample_names[idx], &features);
  }

  // ── VariantFeatures (7) ──────────────────────────────────────────────────
  const auto& v = candidate.variant();
  features.push_back(IsSnp(v, exclude) ? 1 : 0);
  features.push_back(IsInsertion(v, exclude) ? 1 : 0);
  features.push_back(IsDeletion(v, exclude) ? 1 : 0);
  int ins_len = 0;
  for (int idx : alt_allele_indices) {
    if (idx < 0 || idx >= v.alternate_bases_size()) continue;
    int d = static_cast<int>(v.alternate_bases(idx).size()) -
            static_cast<int>(v.reference_bases().size());
    ins_len = std::max(ins_len, d);
  }
  features.push_back(std::max(0, ins_len));
  int del_len = 0;
  for (int idx : alt_allele_indices) {
    if (idx < 0 || idx >= v.alternate_bases_size()) continue;
    int d = static_cast<int>(v.reference_bases().size()) -
            static_cast<int>(v.alternate_bases(idx).size());
    del_len = std::max(del_len, d);
  }
  features.push_back(std::max(0, del_len));
  features.push_back(NumNonSymbolicAlts(v) > 1 ? 1 : 0);
  features.push_back(alt_allele_indices.size() > 1 ? 1 : 0);

  // ── VAF context (51) ─────────────────────────────────────────────────────
  const auto& vaf_at_pos = candidate.allele_frequency_at_position();
  const int half = kSmallModelVafContextWindow / 2;  // 25
  for (int o = -half; o <= half; ++o) {
    const int64_t pos = v.start() + o;
    auto it = vaf_at_pos.find(static_cast<int>(pos));
    features.push_back(it != vaf_at_pos.end() ? it->second : 0);
  }

  return features;
}

// ── Haplotype-expanded feature encoder (PacBio/ONT germline, 106 features) ──
//
// Mirrors Python's SmallModelExamplesEncoder with expand_by_haplotype=True.
// After the standard 70 features, appends 36 more: for each HP in {0,1,2},
// compute 12 BaseFeatures filtering reads to those whose read_name appears in
// `read_hp_tags` with that HP value.  Reads absent from `read_hp_tags` count
// as HP=0 (unphased).
//
// `read_hp_tags` maps fragment_name+"/"+read_number → HP tag (0, 1, 2).
// Build it from the BAM reads using read.info()["HP"].

namespace {

// Compute 12 BaseFeatures for reads filtered to `hp_value`.
// `read_hp_tags` maps read_name → hp (0/1/2); absent ⟹ treated as 0.
void AppendBaseFeaturesForHP(
    const DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices,
    int8_t hp_value,
    const std::unordered_map<std::string, int8_t>& read_hp_tags,
    std::vector<float>* features) {
  // Helper: get HP for a read name (default 0 if absent).
  auto hp_of = [&](const std::string& name) -> int8_t {
    auto it = read_hp_tags.find(name);
    return (it == read_hp_tags.end()) ? 0 : it->second;
  };

  // Alt reads for this HP group.
  std::vector<const DeepVariantCall_ReadSupport*> alt_reads;
  for (int idx : alt_allele_indices) {
    if (idx < 0 || idx >= candidate.variant().alternate_bases_size()) continue;
    const auto& alt_bases = candidate.variant().alternate_bases(idx);
    auto it = candidate.allele_support_ext().find(alt_bases);
    if (it == candidate.allele_support_ext().end()) continue;
    for (const auto& r : it->second.read_infos()) {
      if (hp_of(r.read_name()) == hp_value) alt_reads.push_back(&r);
    }
  }

  // Ref reads for this HP group.
  std::vector<const DeepVariantCall_ReadSupport*> ref_reads;
  for (const auto& r : candidate.ref_support_ext().read_infos()) {
    if (hp_of(r.read_name()) == hp_value) ref_reads.push_back(&r);
  }

  // Total depth: ALWAYS unfiltered (same invariant as AppendBaseFeatures).
  int total_depth = candidate.ref_support_ext().read_infos_size();
  for (const auto& [_, support] : candidate.allele_support_ext()) {
    total_depth += support.read_infos_size();
  }

  const int n_ref = static_cast<int>(ref_reads.size());
  const int n_alt = static_cast<int>(alt_reads.size());
  const int alt_depth = n_ref + n_alt;
  features->push_back(n_ref);
  features->push_back(n_alt);
  features->push_back(alt_depth);
  features->push_back(total_depth);
  features->push_back(total_depth > 0 ? (100 * n_alt / total_depth) : 0);
  features->push_back(alt_depth > 0 ? (100 * n_alt / alt_depth) : 0);
  features->push_back(MeanInt(ref_reads, GetMQ));
  features->push_back(MeanInt(alt_reads, GetMQ));
  features->push_back(MeanInt(ref_reads, GetBQ));
  features->push_back(MeanInt(alt_reads, GetBQ));
  features->push_back(MeanInt(ref_reads, GetReverseStrand100));
  features->push_back(MeanInt(alt_reads, GetReverseStrand100));
}

}  // namespace (anonymous)

std::vector<float> EncodeSmallModelFeaturesHaplotype(
    const DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices,
    const std::unordered_map<std::string, int8_t>& read_hp_tags) {
  std::vector<float> features;
  features.reserve(kSmallModelNumFeaturesHaplotype);

  // ── Standard 70 features (same as EncodeSmallModelFeatures) ──────────────
  std::set<std::string> exclude;
  {
    std::set<int> idx_set(alt_allele_indices.begin(), alt_allele_indices.end());
    for (int i = 0; i < candidate.variant().alternate_bases_size(); ++i) {
      if (!idx_set.count(i))
        exclude.insert(candidate.variant().alternate_bases(i));
    }
  }
  AppendBaseFeatures(candidate, alt_allele_indices, /*sample_filter=*/"",
                      &features);
  const auto& v = candidate.variant();
  features.push_back(IsSnp(v, exclude) ? 1 : 0);
  features.push_back(IsInsertion(v, exclude) ? 1 : 0);
  features.push_back(IsDeletion(v, exclude) ? 1 : 0);
  int ins_len = 0, del_len = 0;
  for (int idx : alt_allele_indices) {
    if (idx < 0 || idx >= v.alternate_bases_size()) continue;
    ins_len = std::max(ins_len, static_cast<int>(v.alternate_bases(idx).size()) -
                                    static_cast<int>(v.reference_bases().size()));
    del_len = std::max(del_len, static_cast<int>(v.reference_bases().size()) -
                                    static_cast<int>(v.alternate_bases(idx).size()));
  }
  features.push_back(std::max(0, ins_len));
  features.push_back(std::max(0, del_len));
  features.push_back(NumNonSymbolicAlts(v) > 1 ? 1 : 0);
  features.push_back(static_cast<int>(alt_allele_indices.size()) > 1 ? 1 : 0);
  {
    const auto& vaf_at_pos = candidate.allele_frequency_at_position();
    const int half = kSmallModelVafContextWindow / 2;
    for (int o = -half; o <= half; ++o) {
      const int64_t pos = v.start() + o;
      auto it = vaf_at_pos.find(static_cast<int>(pos));
      features.push_back(it != vaf_at_pos.end() ? it->second : 0);
    }
  }

  // ── Haplotype-expanded block: 12 × 3 = 36 extra features ─────────────────
  // Mirrors expand_by_haplotype=True in upstream FeatureEncoder:
  //   for sample in [only_sample]:
  //     for hp in [HP_0, HP_1, HP_2]:
  //       encode 12 BaseFeatures filtered to reads with that HP tag
  for (int8_t hp = 0; hp <= 2; ++hp) {
    AppendBaseFeaturesForHP(candidate, alt_allele_indices, hp,
                             read_hp_tags, &features);
  }

  return features;
}

}  // namespace deepvariant
