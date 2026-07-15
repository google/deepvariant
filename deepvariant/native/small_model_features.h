// Compute the features that upstream's small_model takes as input.
// Feature order (must match upstream make_small_model_examples.py
// _encode_candidate_feature_dict output order so the same model
// weights produce the same predictions):
//
// Single-sample (WGS / WES / etc.) — 70 features:
//   0..11   : 12 BaseFeatures (over the target sample's reads)
//   12..18  : 7  VariantFeatures
//   19..69  : 51 VAF-context features (offset -25..+25 inclusive)
//
// Multi-sample (DeepTrio with 3 samples, DeepSomatic with 2) —
// 70 + 12 × N features. Upstream's _encode_candidate_feature_dict
// inserts dict keys in this order:
//   1. 12 BaseFeatures (combined / target-only)
//   2. For each sample in `sample_order`:
//        12 BaseFeatures filtered to that sample's reads
//   3. 7 VariantFeatures
//   4. 51 VAF context features
// → 12 + 12 × N + 7 + 51 features. For trio (N=3): 106. For somatic
// (N=2): 94.
//
// Tested against extracted upstream small_model bundles:
//   /opt/smallmodels/wgs/model.keras                     (input dim 70)
//   /opt/smallmodels/deeptrio/wgs/{child,parent}/model.keras (input dim 106)

#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "deepvariant/protos/deepvariant.pb.h"

namespace deepvariant {

constexpr int kSmallModelNumFeatures = 70;
constexpr int kSmallModelVafContextWindow = 51;
// Number of base features per sample (used for multi-sample feature dim).
constexpr int kSmallModelBaseFeaturesPerSample = 12;
// Haplotype-expanded models (PacBio/ONT germline): 3 HP groups × 12 = 36 extra.
// Total features with haplotypes = 70 + 36 = 106.
constexpr int kSmallModelNumHaplotypeFeatures = 36;
constexpr int kSmallModelNumFeaturesHaplotype =
    kSmallModelNumFeatures + kSmallModelNumHaplotypeFeatures;

// Build the 70-feature vector for a candidate, against a chosen subset of
// alt_allele_indices. Single-sample interface (WGS path).
std::vector<float> EncodeSmallModelFeatures(
    const learning::genomics::deepvariant::DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices);

// Build the 106-feature vector for haplotype-expanded models (PacBio/ONT).
// Appends 36 extra features after the standard 70: for each of HP 0, 1, 2,
// compute 12 BaseFeatures filtering allele_support reads by their HP tag.
// `read_hp_tags` maps (fragment_name + "/" + read_number) → HP value (0/1/2).
std::vector<float> EncodeSmallModelFeaturesHaplotype(
    const learning::genomics::deepvariant::DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices,
    const std::unordered_map<std::string, int8_t>& read_hp_tags);

// Build the (70 + 12*N)-feature vector for trio / somatic, where N is
// `sample_names.size()`. `sample_order` is a permutation of indices
// into `sample_names` matching the per-target rendering order
// upstream uses (for DeepTrio: child target → [0,1,2]; parent2 target
// → [2,1,0]). Reads supporting each per-sample feature group are
// filtered by `read_info.sample_name`.
std::vector<float> EncodeSmallModelFeaturesMultiSample(
    const learning::genomics::deepvariant::DeepVariantCall& candidate,
    const std::vector<int>& alt_allele_indices,
    const std::vector<std::string>& sample_names,
    const std::vector<int>& sample_order);

}  // namespace deepvariant
