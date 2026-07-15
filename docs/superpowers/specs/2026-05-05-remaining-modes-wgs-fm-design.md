# Design: Remaining modes + WGS FM improvement

**Date:** 2026-05-05 — **Status:** Approved

## 1. PacBio/ONT small model — 106 features (expand_by_haplotype)

**Root cause:** `kSmallModelNumFeatures=70` but PacBio model expects 106.
**Why 106:** `expand_by_haplotype=true` adds 3 HP groups × 12 base features = 36.
Standard 70 (12 base + 7 variant + 51 VAF) + 36 (HP0×12 + HP1×12 + HP2×12) = 106.

**Implementation:**
- New `AppendHaplotypeBaseFeatures(candidate, alt_indices, hp_value, features)` filters reads by HP tag and appends the same 12 base features per HP value (0, 1, 2)
- HP tags available from `candidate.allele_support` read names + a new `read_haplotypes` map passed to `EncodeSmallModelFeatures` (HP tag stored from SAM aux HP field during AlleleCounter)
- New `ABSL_FLAG(bool, small_model_use_haplotypes, false)` — auto-set for PACBIO/ONT in cli.cc
- Feature vector: 70 standard + 36 haplotype block = 106 (or 70 if flag off)
- `kSmallModelNumFeatures` stays 70; total computed at runtime

**Gate:** PacBio germline `--small_model_path=pacbio_small_weights` → 0 FM vs Docker on chr20:10M-10.1M.

## 2. WGS FM improvement — temperature scaling calibration

**Baseline:** 4,146 FM on HG002 WGS; 1,469 are PASS↔NoCall/RefCall (GQ≈20 borderline).

**Approach:** Scan T ∈ {0.6, 0.7, 0.8, 0.9} with `--enable_temp_scaling --temp_scaling_T=T` on chr20 full run. Pick T that minimises PASS↔NoCall/RefCall transitions. NoCall↔RefCall (both homref, 2,639 sites) tolerable but a lower T may also reduce those.

**Already implemented** — just needs calibration run + commit of optimal T.

## 3. DeepTrio PacBio/ONT heights

Heights already fixed (100/100 WES/ONT). PacBio trio: child=60, parent=40 = 140 (same as WGS, already correct). ONT trio: child=100, parent=100 = 300.

Test with Illumina BAM proxy; real validation needs chr1 PacBio BAM download (~5 GB).

## 4. Homebrew formula skeleton

Files: `release/homebrew/deepvariant.rb`, `release/homebrew/deepvariant-models.rb`.
Pattern mirrors existing `release/build_glnexus.sh`. Sets `DEEPVARIANT_MODELS_DIR`.
No signing/notarisation in this pass — needs Apple Developer account.
