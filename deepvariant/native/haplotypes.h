// Phase 5.5d/4 — haplotype-resolution port from
// deepvariant/haplotypes.py.
//
// When multiple variants overlap on the reference, their genotype calls
// must be COMPATIBLE: at any reference position, the sum of non-ref
// genotype counts across covering variants must be ≤ ploidy (=2). When
// the calls would violate this (e.g. an indel called 0/1 and an inside
// SNP called 1/1 = three non-ref alleles at the SNP base), we re-search
// over the compatible (non-ref-count, allele-index) configurations,
// pick the joint argmax and the marginal argmax, and apply the result
// when they agree.
//
// This is what closes the residual 27 chr20 PASS-flips after Phase
// 5.5d/{1,2,3} — they sit on overlapping-indel + SNP groups where
// upstream forces the SNP to homref to keep ploidy=2.

#pragma once

#include <vector>

#include "third_party/nucleus/protos/variants.pb.h"

namespace deepvariant {

// Walks `variants` (sorted by chrom+start), groups overlapping ones,
// and applies haplotype resolution per group. Modifies in place.
// `qual_filter` is the threshold used when re-deriving the FILTER
// field after a genotype rewrite (mirrors qual_filter in
// `compute_filter_fields`).
void MaybeResolveConflictingVariants(
    std::vector<nucleus::genomics::v1::Variant>* variants,
    double qual_filter);

}  // namespace deepvariant
