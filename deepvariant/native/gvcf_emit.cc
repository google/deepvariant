// Phase 9 / Step 3 — gVCF reference-row generator implementation.
//
// Ports `make_gvcfs` from upstream's variant_caller.py:256-410 to C++.

#include "deepvariant/native/gvcf_emit.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "absl/log/log.h"

namespace deepvariant {

namespace {

constexpr double kImpossiblePLog10 = -1000.0;
constexpr double kLog10 = 2.302585092994046;  // ln(10)

// Normalise log10 probabilities so that sum(10^log10_p) == 1.
void NormalizeLog10Probs(double* a, double* b, double* c) {
  // Find max for numerical stability.
  const double m = std::max({*a, *b, *c});
  const double sum_lin =
      std::pow(10.0, *a - m) + std::pow(10.0, *b - m) + std::pow(10.0, *c - m);
  const double log_sum = std::log10(sum_lin) + m;
  *a -= log_sum;
  *b -= log_sum;
  *c -= log_sum;
}

// Phred quality of the not-best-genotype probability mass: -10 *
// log10(1 - p_ref) given log10(p_ref). Bounded at max_gq.
int Log10PtrueToPhred(double log10_p_ref, int max_gq) {
  // p_ref = 10^log10_p_ref.
  // 1 - p_ref = sum of all other genotype probs.
  // GQ = -10 * log10(1 - p_ref).
  const double p_ref = std::pow(10.0, log10_p_ref);
  if (p_ref >= 1.0) return max_gq;
  const double q = 1.0 - p_ref;
  if (q <= 0.0) return max_gq;
  const int gq = static_cast<int>(std::floor(-10.0 * std::log10(q)));
  return std::min(gq, max_gq);
}

// log10 of an effectively-impossible probability, used to drive the
// heterozygous likelihood to zero on haploid contigs (matches upstream
// variant_caller.py:IMPOSSIBLE_PROBABILITY_LOG10 = 999.0).
constexpr double kImpossibleProbabilityLog10 = 999.0;

// Compute reference-confidence likelihoods + GQ for one site, given ref/total
// read counts and per-base error rate. Returns log10 probs in [ref, het, alt].
// When `is_haploid` (a haploid contig outside the PAR), the heterozygous
// likelihood is forced to zero, mirroring variant_caller.py's is_haploid
// branch in _calc_reference_confidence.
void ReferenceConfidence(int n_ref, int n_total, double p_error,
                         bool is_haploid, double* log10_p_ref,
                         double* log10_p_het, double* log10_p_alt) {
  if (n_total <= 0) {
    // No coverage: uniform over the possible genotypes. Haploid drops het.
    *log10_p_ref = -1.0;
    *log10_p_het = is_haploid ? -kImpossibleProbabilityLog10 : -1.0;
    *log10_p_alt = -1.0;
  } else {
    const int n_alts = n_total - n_ref;
    const double logp = std::log(p_error) / kLog10;
    const double log1p = std::log1p(-p_error) / kLog10;
    *log10_p_ref = n_ref * log1p + n_alts * logp;
    *log10_p_het = is_haploid ? -kImpossibleProbabilityLog10
                              : -n_total * std::log10(2.0);
    *log10_p_alt = n_ref * logp + n_alts * log1p;
  }
  NormalizeLog10Probs(log10_p_ref, log10_p_het, log10_p_alt);
}

// Mirror upstream variant_caller.py:_quantize_gq exactly. For binsize=5:
//   raw_gq=48 → bin (48-1)//5=9 → 46
//   raw_gq=50 → bin (50-1)//5=9 → 46
// Different from a naive floor(raw/bs)*bs which would split 48 and 50 into
// separate bins (45 and 50) and emit twice as many gVCF blocks.
int QuantizeGq(int raw_gq, int binsize) {
  if (raw_gq < 1) return 0;
  if (binsize <= 1) return raw_gq;
  const int bin_number = (raw_gq - 1) / binsize;
  return bin_number * binsize + 1;
}

// Per-site computed values, used for grouping.
struct SiteEntry {
  int position;
  std::string ref_base;
  std::string ref_name;
  int n_total;
  int quantized_gq;
  int raw_gq;
  double log10_probs[3];
  bool gl_is_valid;  // true if max(log10_probs) == log10_probs[0]
};

bool IsCanonicalDnaBase(const std::string& s) {
  return s == "A" || s == "C" || s == "G" || s == "T";
}

}  // namespace

std::vector<nucleus::genomics::v1::Variant> MakeGvcfRows(
    const std::vector<learning::genomics::deepvariant::AlleleCountSummary>&
        summaries,
    const std::string& sample_name,
    double p_error, int gq_resolution, int max_gq, bool include_med_dp,
    const std::set<std::string>* haploid_contigs,
    const ParRegions* par_regions) {
  std::vector<nucleus::genomics::v1::Variant> out;
  if (summaries.empty()) return out;

  static const ParRegions kNoParRegions;
  const ParRegions& par = par_regions ? *par_regions : kNoParRegions;

  // 1. Compute per-site GQ + likelihoods.
  std::vector<SiteEntry> entries;
  entries.reserve(summaries.size());
  for (const auto& s : summaries) {
    SiteEntry e;
    e.position = s.position();
    e.ref_base = s.ref_base();
    e.ref_name = s.reference_name();
    e.n_total = s.total_read_count();
    if (!IsCanonicalDnaBase(e.ref_base)) {
      // Skip non-canonical (N, IUPAC) — upstream does the same.
      continue;
    }
    const bool is_haploid =
        haploid_contigs &&
        IsHaploidPosition(e.ref_name, e.position, e.position + 1,
                          *haploid_contigs, par);
    ReferenceConfidence(s.ref_supporting_read_count(), s.total_read_count(),
                        p_error, is_haploid, &e.log10_probs[0],
                        &e.log10_probs[1], &e.log10_probs[2]);
    e.raw_gq = Log10PtrueToPhred(e.log10_probs[0], max_gq);
    e.quantized_gq = QuantizeGq(e.raw_gq, gq_resolution);
    e.gl_is_valid =
        e.log10_probs[0] >= e.log10_probs[1] &&
        e.log10_probs[0] >= e.log10_probs[2];
    entries.push_back(std::move(e));
  }

  // 2. Group consecutive entries with same (quantized_gq, gl_is_valid). Emit
  // one merged Variant row per group when gl_is_valid; emit one Variant per
  // site when not (uncalled `./.` rows).
  //
  // is_haploid is deliberately NOT part of the grouping key: upstream
  // variant_caller.make_gvcfs groups only on (quantized_gq, has_valid_gl), so
  // keying on ploidy here would emit more blocks than upstream and break gVCF
  // parity. In practice a haploid site drops its het mass and renormalizes to a
  // higher ref probability (higher GQ), so it usually lands in a different
  // quantized-GQ bin than an adjacent diploid site and won't merge anyway.
  size_t i = 0;
  while (i < entries.size()) {
    const SiteEntry& first = entries[i];
    size_t j = i + 1;
    while (j < entries.size() &&
           entries[j].quantized_gq == first.quantized_gq &&
           entries[j].gl_is_valid == first.gl_is_valid &&
           entries[j].position == entries[j - 1].position + 1 &&
           entries[j].ref_name == first.ref_name) {
      ++j;
    }
    const SiteEntry& last = entries[j - 1];

    // Compute min_gq, min_dp, med_dp over [i, j).
    int min_gq = first.raw_gq;
    int min_dp = first.n_total;
    int min_idx = static_cast<int>(i);
    std::vector<int> dps;
    dps.reserve(j - i);
    for (size_t k = i; k < j; ++k) {
      if (entries[k].raw_gq < min_gq) {
        min_gq = entries[k].raw_gq;
        min_idx = static_cast<int>(k);
      }
      if (entries[k].n_total < min_dp) min_dp = entries[k].n_total;
      dps.push_back(entries[k].n_total);
    }
    std::sort(dps.begin(), dps.end());
    // True median: average the two middle values on even-length input, matching
    // upstream int(statistics.median(...)).
    const size_t n = dps.size();
    const int med_dp = (n % 2 == 1)
                           ? dps[n / 2]
                           : static_cast<int>((dps[n / 2 - 1] + dps[n / 2]) / 2);

    if (first.gl_is_valid) {
      // Emit ONE merged Variant for [i, j).
      nucleus::genomics::v1::Variant v;
      v.set_reference_name(first.ref_name);
      v.set_reference_bases(first.ref_base);
      v.add_alternate_bases("<*>");
      v.set_start(first.position);
      v.set_end(last.position + 1);
      auto* call = v.add_calls();
      call->set_call_set_name(sample_name);
      call->add_genotype(0);
      call->add_genotype(0);
      const auto& min_p = entries[min_idx].log10_probs;
      call->add_genotype_likelihood(min_p[0]);
      call->add_genotype_likelihood(min_p[1]);
      call->add_genotype_likelihood(min_p[2]);
      auto* info_map = call->mutable_info();
      (*info_map)["GQ"].add_values()->set_int_value(min_gq);
      (*info_map)["MIN_DP"].add_values()->set_int_value(min_dp);
      if (include_med_dp) {
        (*info_map)["MED_DP"].add_values()->set_int_value(med_dp);
      }
      // PL = phred-scaled, zero-shifted log10 likelihoods. Mirrors
      // nucleus/io/vcf_conversion.cc:1220-1228 exactly:
      //   normalized = gl - max(gl)   (ZeroShiftLikelihoods)
      //   phred = -10 * normalized    (Log10PErrorToPhred, double)
      //   pl = static_cast<int>(phred) (implicit double→int = trunc)
      {
        const double max_gl = std::max({min_p[0], min_p[1], min_p[2]});
        auto* pl_field = &(*info_map)["PL"];
        for (int g = 0; g < 3; ++g) {
          const double phred = -10.0 * (min_p[g] - max_gl);
          pl_field->add_values()->set_int_value(static_cast<int>(phred));
        }
      }
      out.push_back(std::move(v));
    } else {
      // Uncalled GT=./. for each site individually (one Variant per site).
      // Skipping merging is what upstream does — see variant_caller.py:392-410.
      for (size_t k = i; k < j; ++k) {
        nucleus::genomics::v1::Variant v_each;
        v_each.set_reference_name(entries[k].ref_name);
        v_each.set_reference_bases(entries[k].ref_base);
        v_each.add_alternate_bases("<*>");
        v_each.set_start(entries[k].position);
        v_each.set_end(entries[k].position + 1);
        auto* c = v_each.add_calls();
        c->set_call_set_name(sample_name);
        c->add_genotype(-1);
        c->add_genotype(-1);
        for (int q = 0; q < 3; ++q) {
          c->add_genotype_likelihood(entries[k].log10_probs[q]);
        }
        // PL on uncalled rows (mirrors valid-GL path above).
        {
          auto* uc_info = c->mutable_info();
          (*uc_info)["GQ"].add_values()->set_int_value(entries[k].raw_gq);
          (*uc_info)["MIN_DP"].add_values()->set_int_value(entries[k].n_total);
          if (include_med_dp) {
            (*uc_info)["MED_DP"].add_values()->set_int_value(entries[k].n_total);
          }
          const double max_gl = std::max(
              {entries[k].log10_probs[0], entries[k].log10_probs[1],
               entries[k].log10_probs[2]});
          auto* pl_field = &(*uc_info)["PL"];
          for (int g = 0; g < 3; ++g) {
            const double phred = -10.0 * (entries[k].log10_probs[g] - max_gl);
            pl_field->add_values()->set_int_value(static_cast<int>(phred));
          }
        }
        out.push_back(std::move(v_each));
      }
    }
    i = j;
  }
  return out;
}

}  // namespace deepvariant
