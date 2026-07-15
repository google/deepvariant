// Phase 9 / Step 3 — gVCF reference-row generator.
//
// Ports upstream's `make_gvcfs` algorithm from variant_caller.py to C++.
// Walks per-position AlleleCountSummary protos, computes reference
// confidence (GQ + log10 likelihoods for ref/het/homalt), groups
// consecutive sites with the same quantized GQ into single Variant
// records with the `<*>` alt allele and an END info field — exactly
// matching upstream's gVCF emission format. Output Variants are
// consumed by `nucleus::MergeAndWriteVariantsAndNonVariants` in
// postprocess to emit the final gVCF.

#pragma once

#include <set>
#include <string>
#include <vector>

#include "deepvariant/native/haploid_regions.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "third_party/nucleus/protos/variants.pb.h"

namespace deepvariant {

// Generate gVCF reference rows from per-position AlleleCountSummary.
//
// Inputs:
//   summaries: AlleleCountSummary protos in coordinate-sorted order
//              (one per genomic position in the region).
//   sample_name: emitted as VariantCall.call_set_name.
//   p_error: per-base error rate (typical: 1e-3).
//   gq_resolution: GQ binsize for grouping (typical: 1; upstream
//                  default).
//   max_gq: upper cap on GQ (typical: 50).
//   include_med_dp: emit MED_DP info field (default false).
//   haploid_contigs: contigs to call as haploid (null/empty = all diploid).
//   par_regions: PAR intervals exempted from haploid calling (null = none).
//
// Returns: coordinate-sorted Variant protos with `<*>` alt and
// `END` info field. Each row may span multiple positions if their
// quantized GQs match. On a haploid contig outside the PAR, the
// reference-confidence het likelihood is forced to zero.
std::vector<nucleus::genomics::v1::Variant> MakeGvcfRows(
    const std::vector<learning::genomics::deepvariant::AlleleCountSummary>&
        summaries,
    const std::string& sample_name,
    double p_error = 1e-3,
    int gq_resolution = 5,    // Matches upstream --gvcf_gq_binsize default.
    int max_gq = 50,
    bool include_med_dp = false,
    const std::set<std::string>* haploid_contigs = nullptr,
    const ParRegions* par_regions = nullptr);

}  // namespace deepvariant
