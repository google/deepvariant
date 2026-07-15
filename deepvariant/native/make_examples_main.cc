// Native make_examples — calling mode only (no training, no labeling).
// Replaces the Python make_examples_core.py orchestration layer.
//
// Pipeline per region:
//   SamReader.Query → AlleleCounter → VariantCaller → ExamplesGenerator
//
// The heavy C++ implementations (AlleleCounter, VariantCaller, pileup image
// encoding) are fully reused from upstream; only the orchestration is new.

#include "deepvariant/native/make_examples_main.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#include "deepvariant/allelecounter.h"
#include "deepvariant/direct_phasing.h"
#include "deepvariant/make_examples_native.h"
#include "deepvariant/native/gvcf_emit.h"
#include "deepvariant/native/haploid_regions.h"
#include "deepvariant/methylation_aware_phasing.h"
#include "deepvariant/native/numpy_mt19937.h"
#include "deepvariant/native/realigner_native.h"
#include "deepvariant/native/regions.h"
#include "deepvariant/native/small_model_features.h"
#include "deepvariant/native/small_model_inference.h"
#include "deepvariant/native/dv_signpost.h"
#include "deepvariant/native/tfrecord.h"
#include "deepvariant/protos/deepvariant.pb.h"
#include "deepvariant/protos/realigner.pb.h"
#include "deepvariant/variant_calling.h"
#include "deepvariant/variant_calling_multisample.h"
#include "absl/flags/flag.h"
#include "absl/flags/parse.h"
#include "absl/log/check.h"
#include "absl/log/initialize.h"
#include "absl/log/log.h"
#include "absl/strings/str_cat.h"
#include "absl/strings/str_split.h"
#include "third_party/nucleus/io/reference.h"
#include "third_party/nucleus/io/sam_reader.h"
#include "third_party/nucleus/io/vcf_reader.h"
#include "third_party/nucleus/protos/range.pb.h"
#include "third_party/nucleus/protos/reads.pb.h"
#include "third_party/nucleus/protos/reference.pb.h"
#include "third_party/nucleus/protos/struct.pb.h"
#include "third_party/nucleus/util/proto_ptr.h"
#include "third_party/nucleus/util/utils.h"
#include <cmath>

ABSL_FLAG(std::string, gvcf, "",
          "Phase 9 / Step 3 — output non-variant TFRecord path. When "
          "non-empty, make_examples emits per-region gVCF reference "
          "rows (homref `<*>` records with GQ + MIN_DP info fields, "
          "band-coalesced) to this file alongside the regular examples "
          "output. Postprocess merges these with the variant CVOs to "
          "produce a complete gVCF. Default empty = no gVCF emission "
          "(preserves baseline). Mirrors upstream's --gvcf flag in "
          "make_examples.");
ABSL_FLAG(int32_t, gvcf_gq_binsize, 5,
          "Bin size for quantizing gVCF genotype qualities. Larger bins "
          "merge adjacent positions more aggressively, reducing the gVCF "
          "row count at the cost of GQ granularity. Mirrors upstream's "
          "--gvcf_gq_binsize default of 5.");
ABSL_FLAG(double, p_error, 1e-3,
          "Per-base sequencing error rate used by the gVCF reference "
          "confidence model. Mirrors upstream's --p_error default 0.001.");
ABSL_FLAG(bool, include_med_dp, false,
          "Emit MED_DP info field in gVCF rows (median DP per block). "
          "Mirrors upstream's --include_med_dp.");
ABSL_FLAG(bool, use_direct_phasing, false,
          "Phase 9 / Step 4 — run upstream's DirectPhasing algorithm "
          "(deepvariant/direct_phasing.{h,cc}, Boost-graph max-weight "
          "phasing) on candidates+reads per region, mark each candidate's "
          "VariantCall.is_phased and info[\"PS\"] before TFRecord emit. "
          "Default false to preserve baseline (matches our shipping "
          "default; upstream's Python default is true). Wired in both "
          "the trio (~line 1731) and solo (~line 2210) worker paths; "
          "PS info field is populated from the per-region "
          "position_to_ps map (commit fbead42f).");
ABSL_FLAG(bool, enable_methylation_calling, false,
          "Phase 9 / Step 2 — read MM/ML SAM tags for base "
          "modifications (5mC). When true, AlleleCounter computes "
          "per-allele methylation fraction (ratio of 5mC-modified "
          "to total reads supporting that allele) and the pileup "
          "image gets an extra `base_methylation` channel. Default "
          "false = no methylation-related fields emitted (matches "
          "DV WGS/WES baseline). Used for ONT/PacBio methylation "
          "calling.");
ABSL_FLAG(double, methylation_calling_threshold, 0.5,
          "Phase 9 / Step 2 — minimum methylation probability "
          "(from ML tag) for a base to be classified as 5mC. "
          "Default 0.5 matches upstream make_examples_options.py.");
ABSL_FLAG(bool, enable_methylation_aware_phasing, false,
          "Run upstream's methylation-aware phasing after DirectPhasing: "
          "methylated reference sites (alt=='.') are split out of the SNP "
          "phasing graph and used to assign a haplotype to reads left "
          "unphased by DirectPhasing, via a Wilcoxon rank-sum test on 5mC "
          "levels. Implies methylation extraction in AlleleCounter. Requires "
          "--use_direct_phasing or --small_model_use_haplotypes to have an "
          "effect. Default false = byte-identical baseline. (PacBio/ONT.)");
ABSL_FLAG(std::string, alt_aligned_pileup, "",
          "Phase 9 / Step 1 — alt-aligned pileup mode for PacBio/ONT "
          "models. One of: none, base_channels, diff_channels, rows, "
          "single_row. Default empty = inherit per-model upstream "
          "default ('diff_channels' for PACBIO/ONT, 'none' for WGS/WES). "
          "When non-empty, overrides the per-model default. Adds 2 "
          "channels for diff_channels/base_channels (7 → 9), extra "
          "rows for rows mode.");
ABSL_FLAG(int64_t, tta_seed_offset, 0,
          "Phase 8 / Tier 2 — additive offset applied to the three "
          "internal RNG seeds (make_examples opts, variant_caller, "
          "pileup_image). Default 0 = baseline (matches Docker). "
          "Non-zero: produces a different shuffle pattern in "
          "DownsampleReadIndices (when coverage > pileup height) "
          "and reservoir sampling, generating an alternative pileup "
          "view of the same region. Used by validation/run_tta.sh "
          "to orchestrate N-pass test-time augmentation.");
ABSL_FLAG(std::string, reads, "", "BAM/CRAM file with aligned reads.");
ABSL_FLAG(std::string, ref, "", "Reference FASTA (.fai index required).");
// `--examples` is the canonical pipeline filespec — defined in call_variants.
ABSL_DECLARE_FLAG(std::string, examples);
ABSL_FLAG(std::string, regions, "",
          "Whitespace-separated region strings (e.g. 'chr20 chr21:1-1000000')."
          " Empty = all contigs.");
ABSL_FLAG(std::string, exclude_regions, "",
          "Whitespace-separated regions to exclude.");
ABSL_FLAG(bool, discard_non_dna_regions, false,
          "If true, exclude reference regions containing only N bases from "
          "processing. Mirrors upstream make_examples_core.py:3382. Effective "
          "when --regions is not also set; matches Python semantics.");
ABSL_FLAG(int, task_id, 0, "0-based shard index.");
ABSL_FLAG(int, num_shards, 0,
          "Total shards. 0 or 1 means no sharding.");
ABSL_FLAG(std::string, sample_name, "",
          "Sample name (inferred from BAM header if empty).");
// Variant calling thresholds — WGS defaults.
ABSL_FLAG(int, vsc_min_count_snps, 2, "Min supporting read count for SNPs.");
ABSL_FLAG(int, vsc_min_count_indels, 2,
          "Min supporting read count for indels.");
ABSL_FLAG(double, vsc_min_fraction_snps, 0.12,
          "Min allele fraction for SNPs.");
ABSL_FLAG(double, vsc_min_fraction_indels, 0.06,
          "Min allele fraction for indels.");
ABSL_FLAG(int, partition_size, 1000,
          "AlleleCounter partition size (bp per window).");
// Default 5 mirrors upstream's make_examples_options.py
// (`--min_mapping_quality` default = 5). The candidate-emission
// AlleleCounter uses this; the WindowSelector / DBG apply their own
// stricter thresholds (20 / 14).
ABSL_FLAG(int, min_mapping_quality, 5, "Min read mapping quality.");
ABSL_FLAG(int, min_base_quality, 10, "Min base quality.");
// Mirrors make_examples_options.py's --select_variant_types: a
// whitespace-separated subset of {snps, indels, insertions, deletions,
// multi-allelics, all}. When set, only candidates whose variant matches one
// of the selectors (OR'd) are kept; empty means keep everything.
ABSL_FLAG(std::string, select_variant_types, "",
          "Whitespace-separated variant types to keep when generating "
          "examples: snps, indels, insertions, deletions, multi-allelics, "
          "all. snps/indels/insertions/deletions select bi-allelic variants "
          "of that type; multi-allelics selects any multi-allelic variant. "
          "Empty (default) keeps all candidates.");
// Sex-chromosome haploid calling. The flags are DEFINED in postprocess_main.cc
// (same multi-call binary, so defining them twice would abort at startup); we
// only read them here to drive haploid gVCF reference confidence.
ABSL_DECLARE_FLAG(std::string, haploid_contigs);
ABSL_DECLARE_FLAG(std::string, par_regions_bed);
// Small model first-pass.
ABSL_FLAG(std::string, small_model, "",
          "Path to the small_model .mlpackage. Empty = no small model "
          "(every candidate goes through the big InceptionV3 model).");
ABSL_FLAG(std::string, small_model_cvo_outfile, "",
          "TFRecord path for CVOs the small model decides directly. "
          "Read by postprocess_variants alongside the big-model CVOs.");
ABSL_FLAG(int, small_model_snp_gq_threshold, 20,
          "Min phred GQ for the small model to commit a SNP call.");
ABSL_FLAG(int, small_model_indel_gq_threshold, 28,
          "Min phred GQ for the small model to commit an indel call.");
ABSL_FLAG(bool, realigner_enabled, false,
          "Enable upstream's realigner (DeBruijnGraph + FastPassAligner) "
          "to recover candidates in indel-rich regions.");
// Realigner aligner SSW scoring params. Defaults match WGS:
// aln_match=4, aln_mismatch=6, aln_gap_open=8, aln_gap_extend=2.
// Pangenome example_info.json:flags_for_calling overrides these to
// 2/5/10/1 (more permissive matches for synthetic haplotypes vs reads).
ABSL_FLAG(int, aln_match, 4, "Realigner SSW aligner match score.");
ABSL_FLAG(int, aln_mismatch, 6, "Realigner SSW aligner mismatch penalty.");
ABSL_FLAG(int, aln_gap_open, 8, "Realigner SSW aligner gap-open penalty.");
ABSL_FLAG(int, aln_gap_extend, 2, "Realigner SSW aligner gap-extend penalty.");
// dbg_disable_graph_pruning: when true, the de-Bruijn graph pruning
// step in the realigner is skipped. Pangenome enables this to retain
// haplotype paths that would otherwise be pruned for low edge weight.
ABSL_FLAG(bool, dbg_disable_graph_pruning, false,
          "If true, skip de-Bruijn graph pruning in the realigner.");
// Per-model pileup + read-filter flags (set by cli.cc ApplyModelFlags()).
// Defaults = WGS. All values mirror upstream make_examples_options.py exactly.
ABSL_FLAG(int, pileup_image_width, 221,
          "Pileup image width. WGS/WES=221, PacBio=147, ONT/MaSeq=199.");
// Named channel preset — selects which channels are added beyond the 6 base
// channels (read_base…base_differs_from_ref):
//   WGS(default) : + insert_size(19)               → 7 ch
//   LONG_READ_PACBIO: + haplotype(7) + suppl(26)   → 8 ch (alt adds 2 → 10)
//   LONG_READ_ONT   : + haplotype(7) + fuzzy(25)   → 8 ch (alt adds 2 → 10)
//   MASSEQ          : + haplotype(7)               → 7 ch (alt adds 2 → 9)
//   BASE_CHANNELS   : no extras                    → 6 ch
ABSL_FLAG(std::string, channel_list_preset, "",
          "Channel preset: WGS, LONG_READ_PACBIO, LONG_READ_ONT, MASSEQ, "
          "BASE_CHANNELS. Empty = WGS.");
ABSL_FLAG(bool, sort_by_haplotypes, false,
          "Sort reads by HP tag in pileup (long-read models).");
ABSL_FLAG(bool, trim_reads_for_pileup, false,
          "Trim reads to pileup window before encoding.");
ABSL_FLAG(bool, phase_reads, false,
          "Phase reads using HP SAM tag.");
ABSL_FLAG(bool, parse_sam_aux_fields, false,
          "Parse auxiliary SAM fields (MM/ML for methylation, HP for phasing).");
ABSL_FLAG(bool, keep_supplementary_alignments, false,
          "Keep supplementary alignments.");
ABSL_FLAG(int, max_reads_per_partition, 1500,
          "Cap reads per partition (0 = unlimited).");
ABSL_FLAG(int, max_reads_for_dynamic_bases_per_region, -1,
          "Max reads for dynamic bases (<0 = disabled, MaSeq only).");
ABSL_FLAG(int, small_model_vaf_context_window_size, 5,
          "VAF context window for small model.");
ABSL_FLAG(double, vsc_min_indel_fraction_for_small_indels, -1.0,
          "Min allele fraction short INDELs (<0 = vsc_min_fraction_indels).");
ABSL_FLAG(double, vsc_min_indel_fraction_for_large_indels, -1.0,
          "Min allele fraction long INDELs (<0 = vsc_min_fraction_indels).");
ABSL_FLAG(int, vsc_small_indel_threshold, -1,
          "INDEL length threshold small vs large (<0 = disabled).");
ABSL_FLAG(bool, split_skip_reads, false,
          "Split reads on N CIGAR ops (RNA-seq).");
// Somatic non-target (normal) AF cap: candidates where the normal sample has
// alt VAF > threshold are skipped as clear germline het/hom.
// Default -1.0 = disabled (FFPE_WGS/FFPE_WES do not declare this in their
// model.example_info.json; WGS/WES/PacBio/ONT declare 0.5).
ABSL_FLAG(double, vsc_max_fraction_snps_for_non_target_sample, -1.0,
          "Normal AF cap for SNPs (<0 = disabled). Set 0.5 for WGS/WES/LR.");
// Sort pileup rows by alt-allele support in somatic TN mode.
// Declared by WGS + FFPE_WGS tumor+normal JSONs only; NOT by WES/FFPE_WES/
// PacBio/ONT. cli.cc sets this flag for WGS and FFPE_WGS TN only.
ABSL_FLAG(bool, sort_by_alt_allele_support_somatic, false,
          "Sort somatic pileup rows by alt support (WGS/FFPE_WGS TN only).");
ABSL_FLAG(double, vsc_max_fraction_indels_for_non_target_sample, -1.0,
          "Normal AF cap for INDELs (<0 = disabled). Set 0.5 for WGS/WES/LR.");

// Enable haplotype-expanded small model features (PacBio/ONT germline).
// When true, EncodeSmallModelFeaturesHaplotype is used instead of
// EncodeSmallModelFeatures: 70 standard + 36 HP-filtered = 106 total.
// Must be set when --small_model_path points to a 106-input model
// (pacbio_small_weights, ont_small_weights). Auto-set by cli.cc for
// PACBIO and ONT model types.
ABSL_FLAG(bool, small_model_use_haplotypes, false,
          "Use haplotype-expanded (106-feature) small model for PacBio/ONT.");

// Panel of Normals VCF for tumor-only allele_frequency pileup channel.
// Path to bgzipped+tabix-indexed VCF. When set, each tumor-only candidate's
// dv_call.allele_frequency map is populated from the PON's per-allele AF INFO
// field, enabling the 8th channel to carry population AFs as expected by
// deepsomatic.*_tumor_only models. Leave empty → default (ref=1, alts=0).
ABSL_FLAG(std::string, population_vcfs, "",
          "Panel-of-Normals VCF for tumor-only allele_frequency channel.");
ABSL_FLAG(int, threads, 1,
          "Worker threads inside this process. >1 enables true intra-process "
          "parallelism (one process showing N×100 % CPU). Each worker opens "
          "its own SamReader / IndexedFastaReader / ExamplesGenerator / "
          "SmallModel and writes to a per-thread file; results are "
          "concatenated into the final --examples / --small_model_cvo_outfile "
          "paths after all workers join.");

// ----------------------------------------------------------------------------
// DeepTrio flags (Step 1 — mirrors deeptrio/make_examples.py exactly).
// When --reads_parent1 is set, make_examples runs in trio mode: 3 samples
// (parent1 at index 0, child at index 1, parent2 at index 2; child is the
// MAIN_SAMPLE_INDEX). Each region is processed by 3 AlleleCounters keyed by
// sample_name and fed to multi_sample::VariantCaller. ExamplesGenerator
// emits 3 separate example streams (one per target sample), each rendered
// with the per-sample `order` permutation so the pileup channel-stack
// shows the target sample in slot 1.
// ----------------------------------------------------------------------------
ABSL_FLAG(std::string, reads_parent1, "",
          "Trio mode: BAM/CRAM for parent1. When set, make_examples runs "
          "as DeepTrio (3 samples: parent1, child, parent2; child = main).");
ABSL_FLAG(std::string, reads_parent2, "",
          "Trio mode: BAM/CRAM for parent2.");
ABSL_FLAG(std::string, sample_name_parent1, "",
          "Trio mode: parent1 sample name (inferred from BAM if empty).");
ABSL_FLAG(std::string, sample_name_parent2, "",
          "Trio mode: parent2 sample name (inferred from BAM if empty).");
ABSL_FLAG(int, pileup_image_height_child, 0,
          "Trio mode: pileup image height for the child sample. 0 = default "
          "(100 per upstream dt_constants.PILEUP_DEFAULT_HEIGHT_CHILD).");
ABSL_FLAG(int, pileup_image_height_parent, 0,
          "Trio mode: pileup image height for each parent sample. 0 = default "
          "(100 per upstream dt_constants.PILEUP_DEFAULT_HEIGHT_PARENT).");
ABSL_FLAG(double, downsample_fraction_child, 0.0,
          "Trio mode: downsample fraction applied to child reads (0.0 = none).");
ABSL_FLAG(double, downsample_fraction_parents, 0.0,
          "Trio mode: downsample fraction applied to both parents' reads.");
ABSL_FLAG(std::string, small_model_path_child, "",
          "Trio mode: small_model weights directory for child examples.");
ABSL_FLAG(std::string, small_model_path_parent, "",
          "Trio mode: small_model weights directory for parent examples.");
ABSL_FLAG(bool, skip_parent_calling, false,
          "Trio mode: if true, generate examples for child only "
          "(parents' SampleOptions still populated for joint candidate "
          "generation, but their example output is suppressed).");
ABSL_FLAG(std::string, examples_child, "",
          "Trio mode: examples output path for the child sample. If empty, "
          "the existing --examples flag is used as the child path.");
ABSL_FLAG(std::string, examples_parent1, "",
          "Trio mode: examples output path for the parent1 sample.");
ABSL_FLAG(std::string, examples_parent2, "",
          "Trio mode: examples output path for the parent2 sample.");
ABSL_FLAG(std::string, small_model_cvo_outfile_child, "",
          "Trio mode: small_model CVO output path for child.");
ABSL_FLAG(std::string, small_model_cvo_outfile_parent1, "",
          "Trio mode: small_model CVO output path for parent1.");
ABSL_FLAG(std::string, small_model_cvo_outfile_parent2, "",
          "Trio mode: small_model CVO output path for parent2.");

// ----------------------------------------------------------------------------
// DeepSomatic mode (Step 2):
//   tumor + normal: 2 samples — normal at index 0, tumor at index 1 (=main).
//   tumor_only:     1 sample  — tumor at index 0 (=main).
// Mirrors deepvariant/make_examples_somatic.py:tumor_normal_samples_from_flags.
// Critical somatic-specific override: vsc_min_fraction_multiplier=inf so
// candidates from the non-target (normal) sample are excluded from the
// tumor candidate set (somatic ≠ trio: we don't want normal-only variants
// in the tumor's call list).
// ----------------------------------------------------------------------------
ABSL_FLAG(std::string, reads_tumor, "",
          "Somatic mode: BAM/CRAM for the tumor sample. When set, make_examples "
          "runs as DeepSomatic (tumor + optional normal; tumor = main).");
ABSL_FLAG(std::string, reads_normal, "",
          "Somatic mode: BAM/CRAM for the normal sample. If empty, runs in "
          "tumor-only mode.");
ABSL_FLAG(std::string, sample_name_tumor, "",
          "Somatic mode: tumor sample name (inferred from BAM if empty).");
ABSL_FLAG(std::string, sample_name_normal, "",
          "Somatic mode: normal sample name (inferred from BAM if empty).");
ABSL_FLAG(int, pileup_image_height_tumor, 0,
          "Somatic mode: pileup image height for the tumor sample. 0 = default "
          "(100 per upstream dv_constants.PILEUP_DEFAULT_HEIGHT).");
ABSL_FLAG(int, pileup_image_height_normal, 0,
          "Somatic mode: pileup image height for the normal sample. 0 = default "
          "(100 per upstream dv_constants.PILEUP_DEFAULT_HEIGHT).");
ABSL_FLAG(double, downsample_fraction_tumor, 0.0,
          "Somatic mode: downsample fraction applied to tumor reads.");
ABSL_FLAG(double, downsample_fraction_normal, 0.0,
          "Somatic mode: downsample fraction applied to normal reads.");
ABSL_FLAG(std::string, small_model_path_somatic, "",
          "Somatic mode: small_model weights directory for the tumor.");
ABSL_FLAG(std::string, small_model_cvo_outfile_tumor, "",
          "Somatic mode: small_model CVO output path for the tumor sample.");
ABSL_FLAG(std::string, examples_tumor, "",
          "Somatic mode: examples output path for the tumor sample. "
          "If empty, the existing --examples flag is used.");
ABSL_FLAG(std::string, examples_normal, "",
          "Somatic mode: examples output path for the normal sample "
          "(usually unused since normal has skip_output_generation=true).");

// ----------------------------------------------------------------------------
// Pangenome-aware DV mode (Step 3):
//   2 samples — pangenome at index 0, reads at index 1 (=main).
// Mirrors deepvariant/make_examples_pangenome_aware_dv.py:
//   reads_and_pangenome_samples_from_flags. Critical pangenome-specific
//   overrides on the pangenome sample (mirrors pangenome_sample_options
//   in upstream Python at line 239):
//   - skip_output_generation=true  (only reads' examples are emitted)
//   - skip_phasing=true            (haplotype tags from reads only)
//   - skip_normalization=true      (no read normalization on synthetic haplotypes)
//   - keep_only_window_spanning_reads (drop reads not spanning the window)
//   - channels_enum_to_blank: CH_HAPLOTYPE_TAG, CH_DIFF_CHANNELS_*,
//     CH_BASE_QUALITY, CH_MAPPING_QUALITY  (5 channels blanked in pangenome rows)
//   - alt_aligned_pileup="none"    (no alt-alignment for pangenome)
// Plus pic-level flag: trim_reads_for_pileup=true (pangenome reads
// are trimmed to fit the example window).
// At runtime the --pangenome flag accepts BAM/CRAM only; GBZ input is
// out of scope for v2 (users pre-extract via Docker's
// load_gbz_into_shared_memory if needed).
// ----------------------------------------------------------------------------
ABSL_FLAG(std::string, reads_pangenome, "",
          "Pangenome mode: BAM/CRAM for the pangenome panel. When set, "
          "make_examples runs as pangenome-aware DV (pangenome + reads; "
          "reads = main). Note: GBZ input is not supported in the native "
          "binary; convert GBZ→BAM via Docker preprocessing.");
ABSL_FLAG(std::string, sample_name_pangenome, "pangenome",
          "Pangenome mode: pangenome sample name "
          "(default 'pangenome').");
ABSL_FLAG(std::string, sample_name_reads, "",
          "Pangenome mode: reads sample name (inferred from BAM if empty).");
ABSL_FLAG(int, pileup_image_height_pangenome, 0,
          "Pangenome mode: pileup image height for the pangenome sample. "
          "0 = default 100.");
ABSL_FLAG(int, pileup_image_height_reads, 0,
          "Pangenome mode: pileup image height for the reads sample. "
          "0 = default 100.");
ABSL_FLAG(double, downsample_fraction_reads, 0.0,
          "Pangenome mode: downsample fraction applied to reads.");
ABSL_FLAG(std::string, small_model_path_pangenome, "",
          "Pangenome mode: small_model weights directory for the reads "
          "sample.");
ABSL_FLAG(std::string, small_model_cvo_outfile_reads, "",
          "Pangenome mode: small_model CVO output path for the reads sample.");
ABSL_FLAG(std::string, examples_reads, "",
          "Pangenome mode: examples output path for the reads sample. "
          "If empty, the existing --examples flag is used.");
ABSL_FLAG(std::string, examples_pangenome, "",
          "Pangenome mode: examples output path for the pangenome sample "
          "(unused since pangenome has skip_output_generation=true).");

namespace deepvariant {

using namespace learning::genomics::deepvariant;  // NOLINT

namespace {

// Build the MakeExamplesOptions proto for calling mode from flags.
MakeExamplesOptions BuildOptions(const std::string& sample_name,
                                 int task_id, int num_shards) {
  MakeExamplesOptions opts;

  opts.set_reference_filename(absl::GetFlag(FLAGS_ref));
  opts.set_examples_filename(absl::GetFlag(FLAGS_examples));
  opts.set_task_id(task_id);
  opts.set_num_shards(num_shards);
  opts.set_mode(MakeExamplesOptions::CALLING);
  // TTA seed offset (default 0 = baseline, matches Docker bit-for-bit).
  // Non-zero: shifts the 3 internal RNG seeds for test-time augmentation.
  const int64_t kTtaOff = absl::GetFlag(FLAGS_tta_seed_offset);
  opts.set_random_seed(609314161 + static_cast<int>(kTtaOff));
  // Reads-per-partition cap — mirrors upstream default 1500. Long-read models
  // may set to 0 (unlimited) via --max_reads_per_partition.
  opts.set_max_reads_per_partition(absl::GetFlag(FLAGS_max_reads_per_partition));
  {
    const int mrd = absl::GetFlag(FLAGS_max_reads_for_dynamic_bases_per_region);
    if (mrd >= 0) opts.set_max_reads_for_dynamic_bases_per_region(mrd);
  }
  // Long-read behavioral flags.
  opts.set_phase_reads(absl::GetFlag(FLAGS_phase_reads));
  opts.set_parse_sam_aux_fields(absl::GetFlag(FLAGS_parse_sam_aux_fields));
  opts.set_trim_reads_for_pileup(absl::GetFlag(FLAGS_trim_reads_for_pileup));
  {
    const bool split = absl::GetFlag(FLAGS_split_skip_reads);
    if (split) opts.mutable_realigner_options()->set_split_skip_reads(true);
  }

  // Read requirements.
  nucleus::genomics::v1::ReadRequirements read_reqs;
  read_reqs.set_min_mapping_quality(absl::GetFlag(FLAGS_min_mapping_quality));
  read_reqs.set_min_base_quality(absl::GetFlag(FLAGS_min_base_quality));
  read_reqs.set_min_base_quality_mode(
      nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
  read_reqs.set_keep_supplementary_alignments(
      absl::GetFlag(FLAGS_keep_supplementary_alignments));

  // Allele counter options.
  AlleleCounterOptions ac_opts;
  ac_opts.set_partition_size(absl::GetFlag(FLAGS_partition_size));
  *ac_opts.mutable_read_requirements() = read_reqs;
  // Required so AlleleCounter actually retains REF-supporting reads in
  // each AlleleCount.read_alleles map (otherwise the small_model sees
  // num_reads_supports_ref = 0 on every candidate and is biased).
  ac_opts.set_track_ref_reads(true);
  // Phase 9 / Step 2 — methylation calling. Wires the MM/ML SAM tag
  // reader (allelecounter.cc::GetMethylationLevel + IsMethylated) to
  // populate AlleleCount.methylation_level. Default off → byte-identical
  // baseline. Per-call methylation fraction is later read from these
  // counts in postprocess to emit MF/MT/MI INFO fields.
  const bool kMethylationOn = absl::GetFlag(FLAGS_enable_methylation_calling);
  ac_opts.set_enable_methylation_calling(kMethylationOn);
  ac_opts.set_methylation_calling_threshold(
      absl::GetFlag(FLAGS_methylation_calling_threshold));
  // Methylation-aware phasing also needs AlleleCounter to extract 5mC levels
  // (allelecounter.cc gates extraction on calling OR aware-phasing), and the
  // VariantCaller emits methylated reference sites (alt='.') as candidates
  // only when the option is set on its options.
  const bool kMethylAwarePhasingOn =
      absl::GetFlag(FLAGS_enable_methylation_aware_phasing);
  ac_opts.set_enable_methylation_aware_phasing(kMethylAwarePhasingOn);
  *opts.mutable_allele_counter_options() = ac_opts;
  // Mirror the flag onto MakeExamplesOptions (used by some downstream
  // code paths, e.g. variant emission / VCF formatting).
  opts.set_enable_methylation_calling(kMethylationOn);
  opts.set_enable_methylation_aware_phasing(kMethylAwarePhasingOn);
  // Phase 9 / Step 4 — DirectPhasing options. Only used when
  // --use_direct_phasing is set; the algorithm wraps candidates +
  // reads to emit per-variant phase info (is_phased + PS).
  opts.mutable_direct_phasing_options()->set_min_alleles_to_phase(1);

  // Variant caller options.
  VariantCallerOptions vc_opts;
  // The caller emits methylated reference sites (alt='.') as candidates only
  // when this is set; they become the methylated_ref_sites for methylation-
  // aware phasing below.
  vc_opts.set_enable_methylation_aware_phasing(kMethylAwarePhasingOn);
  vc_opts.set_min_count_snps(absl::GetFlag(FLAGS_vsc_min_count_snps));
  vc_opts.set_min_count_indels(absl::GetFlag(FLAGS_vsc_min_count_indels));
  vc_opts.set_min_fraction_snps(absl::GetFlag(FLAGS_vsc_min_fraction_snps));
  vc_opts.set_min_fraction_indels(
      absl::GetFlag(FLAGS_vsc_min_fraction_indels));
  // PacBio-style size-stratified INDEL fractions (disabled by default).
  {
    const double small_f = absl::GetFlag(FLAGS_vsc_min_indel_fraction_for_small_indels);
    const double large_f = absl::GetFlag(FLAGS_vsc_min_indel_fraction_for_large_indels);
    const int    thr     = absl::GetFlag(FLAGS_vsc_small_indel_threshold);
    if (small_f >= 0.0) vc_opts.set_vsc_min_indel_fraction_for_small_indels(static_cast<float>(small_f));
    if (large_f >= 0.0) vc_opts.set_vsc_min_indel_fraction_for_large_indels(static_cast<float>(large_f));
    if (thr     >= 0)   vc_opts.set_vsc_small_indel_threshold(thr);
  }
  // VAF context window for small model — on VariantCallerOptions (not
  // SampleOptions) so variant_calling_multisample.cc uses the right window.
  {
    const int vaf_win = absl::GetFlag(FLAGS_small_model_vaf_context_window_size);
    if (vaf_win > 0) vc_opts.set_small_model_vaf_context_window_size(vaf_win);
  }
  vc_opts.set_p_error(0.001);
  vc_opts.set_max_gq(50);
  vc_opts.set_gq_resolution(1);
  vc_opts.set_ploidy(2);
  vc_opts.set_fraction_reference_sites_to_emit(0.0);
  vc_opts.set_random_seed(1260872234 + static_cast<int>(kTtaOff));
  // Required so variant_calling_multisample.cc populates ref_support_ext —
  // without it the small_model sees zero ref-supporting reads on every
  // candidate and predicts hom_ref for everything.
  vc_opts.set_track_ref_reads(true);

  // Pileup image options (WGS defaults).
  PileupImageOptions pic;
  pic.set_reference_band_height(5);
  pic.set_base_color_offset_a_and_g(40);
  pic.set_base_color_offset_t_and_c(30);
  pic.set_base_color_stride(70);
  pic.set_allele_supporting_read_alpha(1.0f);
  pic.set_allele_unsupporting_read_alpha(0.6f);
  pic.set_other_allele_supporting_read_alpha(0.6f);
  pic.set_reference_matching_read_alpha(0.2f);
  pic.set_reference_mismatching_read_alpha(1.0f);
  pic.set_indel_anchoring_base_char("*");
  pic.set_reference_alpha(0.4f);
  pic.set_reference_base_quality(60);
  pic.set_positive_strand_color(70);
  pic.set_negative_strand_color(240);
  pic.set_base_quality_cap(40);
  pic.set_mapping_quality_cap(60);
  pic.set_height(100);
  pic.set_width(absl::GetFlag(FLAGS_pileup_image_width));
  pic.set_read_overlap_buffer_bp(5);
  pic.set_multi_allelic_mode(PileupImageOptions::ADD_HET_ALT_IMAGES);
  pic.set_random_seed(2101079370 + static_cast<int>(kTtaOff));
  // Phase 9 / Step 1 — alt-aligned pileup mode. Empty flag value
  // = inherit upstream per-model default. cli.cc sets the flag from
  // model_type before invoking make_examples; here we just read it.
  {
    std::string aap = absl::GetFlag(FLAGS_alt_aligned_pileup);
    if (aap.empty()) aap = "none";
    pic.set_alt_aligned_pileup(aap);
  }
  pic.set_types_to_alt_align("indels");
  pic.set_min_non_zero_allele_frequency(0.00001f);
  *pic.mutable_read_requirements() = read_reqs;
  // Channel configuration. The 6 base channels are always present.
  // Additional channels depend on --channel_list_preset (set by cli.cc
  // ApplyModelFlags() from the model's example_info.json).
  pic.add_channels("read_base");
  pic.add_channels("base_quality");
  pic.add_channels("mapping_quality");
  pic.add_channels("strand");
  pic.add_channels("read_supports_variant");
  pic.add_channels("base_differs_from_ref");
  {
    const std::string preset = absl::GetFlag(FLAGS_channel_list_preset);
    if (preset == "LONG_READ_PACBIO") {
      // PacBio: haplotype(CH=7) + supplementary_alignment(CH=26)
      // alt_aligned_pileup=diff_channels adds 2 more → 10 total.
      pic.add_channels("haplotype");
      pic.add_channels("supplementary_alignment");
    } else if (preset == "LONG_READ_ONT") {
      // ONT: haplotype(7) + read_supports_variant_fuzzy(25)
      // alt_aligned_pileup=diff_channels adds 2 more → 10 total.
      pic.add_channels("haplotype");
      pic.add_channels("read_supports_variant_fuzzy");
    } else if (preset == "MASSEQ") {
      // MaSeq: haplotype(7); alt_aligned_pileup adds 2 more → 9 total.
      pic.add_channels("haplotype");
    } else if (preset == "BASE_CHANNELS") {
      // HYBRID / RNASeq: 6 channels only (no extras).
    } else {
      // WGS / WES (default): add insert_size → 7 channels.
      pic.add_channels("insert_size");
    }
  }
  // Methylation channel (opt-in via --enable_methylation_calling).
  if (kMethylationOn) {
    pic.add_channels("base_methylation");
  }
  // sort_by_haplotypes: long-read models sort pileup rows by HP tag.
  pic.set_sort_by_haplotypes(absl::GetFlag(FLAGS_sort_by_haplotypes));
  // Alt-aligned channels: also appear in channels() so make_examples_native.cc
  // allocates the correct buffer size (uses channels().size() as depth).
  // diff_channels and base_channels each add 2 extra channels to the pileup.
  // Must be done BEFORE set_num_channels() so the count is correct.
  {
    const std::string& aap = pic.alt_aligned_pileup();
    if (aap == "diff_channels") {
      pic.add_channels("diff_channels_alternate_allele_1");
      pic.add_channels("diff_channels_alternate_allele_2");
    } else if (aap == "base_channels") {
      pic.add_channels("base_channels_alternate_allele_1");
      pic.add_channels("base_channels_alternate_allele_2");
    }
  }
  // num_channels is derived from the channels() list above (now including
  // any alt_aligned channels); set it explicitly so downstream code (e.g.
  // MetalInception::Create) can read the declared channel count without
  // counting the repeated field.
  pic.set_num_channels(static_cast<int>(pic.channels_size()));
  *opts.mutable_pic_options() = pic;

  // Sample options. Trio mode (--reads_parent1 set) populates 3 samples
  // in upstream order [parent1, child, parent2] (mirrors deeptrio/
  // make_examples.py:trio_samples_from_flags). Single-sample mode keeps
  // the legacy single SampleOptions.
  const std::string parent1_reads = absl::GetFlag(FLAGS_reads_parent1);
  const std::string parent2_reads = absl::GetFlag(FLAGS_reads_parent2);
  const bool trio_mode = !parent1_reads.empty();

  if (trio_mode) {
    // Per-model trio defaults (mirror scripts/run_deeptrio.py:392-399):
    //   WGS:    child=60,  parent=40   → total 140 (matches Docker
    //           example_shape=[140, 221, 7])
    //   WES:    child=100, parent=100  → total 300
    //   PACBIO: child=60,  parent=40   → total 140
    //   ONT:    child=100, parent=100  → total 300
    // Users can override via --pileup_image_height_child / _parent.
    // The model_type flag is owned by cli.cc (run mode); here in
    // make_examples_main we infer it via opts.pic_options or default
    // to WGS heights. The cli.cc trio path already passes through the
    // user's --pileup_image_height_* flags so this default only
    // matters for direct `make_examples --reads_parent1=...` usage.
    int child_h  = absl::GetFlag(FLAGS_pileup_image_height_child);
    int parent_h = absl::GetFlag(FLAGS_pileup_image_height_parent);
    if (child_h  <= 0) child_h  = 60;  // DEEP_TRIO_WGS_PILEUP_HEIGHT_CHILD
    if (parent_h <= 0) parent_h = 40;  // DEEP_TRIO_WGS_PILEUP_HEIGHT_PARENT
    const double ds_child   = absl::GetFlag(FLAGS_downsample_fraction_child);
    const double ds_parents = absl::GetFlag(FLAGS_downsample_fraction_parents);
    const std::string p1_name = absl::GetFlag(FLAGS_sample_name_parent1);
    const std::string p2_name = absl::GetFlag(FLAGS_sample_name_parent2);

    auto add_sample = [&](const std::string& role, const std::string& name,
                           const std::string& reads, int height, double ds,
                           std::initializer_list<int> order,
                           bool skip_output, const std::string& small_path) {
      SampleOptions* s = opts.add_sample_options();
      s->set_role(role);
      s->set_name(name);
      if (!reads.empty()) s->add_reads_filenames(reads);
      s->set_pileup_height(height);
      *s->mutable_variant_caller_options() = vc_opts;
      // Per-sample VC options keep the same thresholds; sample_name in
      // the proto is set by upstream via make_vc_options(sample_name=…)
      // — we bake that here so multi_sample::VariantCaller can identify
      // the target sample from its own VC opts.
      s->mutable_variant_caller_options()->set_sample_name(name);
      // Trio default override (mirrors deeptrio/make_examples.py:208):
      //   FLAGS.set_default('vsc_min_fraction_multiplier', 0.67)
      // Used by multi_sample::VariantCaller::IsGoodAltAlleleWithReason
      // when re-evaluating combined-sample evidence (apply_trio_coefficient=
      // true). Lowers the joint-promotion threshold from 0.12 to 0.0804
      // so candidates supported by < 12 % in the target sample but
      // ≥ 8 % combined evidence get promoted (matches Docker's default).
      s->mutable_variant_caller_options()->set_min_fraction_multiplier(0.67f);
      for (int o : order) s->add_order(o);
      s->set_skip_output_generation(skip_output);
      if (!small_path.empty()) s->set_small_model_path(small_path);
      if (ds > 0.0) s->set_downsample_fraction(static_cast<float>(ds));
    };

    const bool skip_parents = absl::GetFlag(FLAGS_skip_parent_calling);

    // Order in `samples_in_order`: [parent1(0), child(1), parent2(2)].
    // Each sample's `order` controls the channel-stack permutation when
    // building its OWN pileup image: the target sample is placed FIRST
    // (slot 0) in its own image so the model always finds its target
    // sample at a fixed position. parent2 additionally swaps parent1↔2
    // vs child so the "other parent" is consistent at slot 2.
    add_sample("parent1", p1_name.empty() ? "parent1" : p1_name,
               parent1_reads, parent_h, ds_parents,
               {0, 1, 2}, skip_parents,
               absl::GetFlag(FLAGS_small_model_path_parent));
    add_sample("child", sample_name,
               absl::GetFlag(FLAGS_reads), child_h, ds_child,
               {0, 1, 2}, /*skip_output=*/false,
               absl::GetFlag(FLAGS_small_model_path_child));
    add_sample("parent2", p2_name.empty() ? "parent2" : p2_name,
               parent2_reads, parent_h, ds_parents,
               {2, 1, 0}, skip_parents,
               absl::GetFlag(FLAGS_small_model_path_parent));

    // MAIN_SAMPLE_INDEX = 1 (child) per deeptrio/make_examples.py:48.
    opts.set_main_sample_index(1);
    opts.set_sample_role_to_train("child");
  } else if (!absl::GetFlag(FLAGS_reads_tumor).empty()) {
    // ──────────────── DeepSomatic mode ────────────────
    // samples_in_order = [normal(0), tumor(1)] when normal provided
    //                  = [tumor(0)]            for tumor-only.
    // Mirrors deepvariant/make_examples_somatic.py:152-218.
    const std::string normal_reads = absl::GetFlag(FLAGS_reads_normal);
    const std::string tumor_reads  = absl::GetFlag(FLAGS_reads_tumor);
    const bool has_normal = !normal_reads.empty();

    // sort_by_alt_allele_support: declared by WGS + FFPE_WGS TN JSONs only.
    // WES, FFPE_WES, PacBio, ONT do NOT declare it. Tumor-only never does.
    // cli.cc passes --sort_by_alt_allele_support_somatic=true for WGS/FFPE_WGS
    // TN only (based on each model's flags_for_calling).
    if (has_normal &&
        absl::GetFlag(FLAGS_sort_by_alt_allele_support_somatic)) {
      opts.mutable_pic_options()->set_sort_by_alt_allele_support(true);
    }

    // Tumor-only: 8th channel = allele_frequency (CH_ALLELE_FREQUENCY=8).
    // Mirrors deepsomatic.*_tumor_only/model.example_info.json channels:
    //   WGS/WES/FFPE: [1,2,3,4,5,6,19,8] (base-7 WGS channels + allele_freq)
    //   PacBio/ONT:   MASSEQ 7ch + alt_aligned×2 + allele_freq = 10ch,
    //                 matching example_info shape [100, w, 10].
    // Tumor+normal models use 7 ch (WGS/WES/FFPE) or 9 ch (long-read),
    // with no allele_frequency.
    if (!has_normal) {
      opts.mutable_pic_options()->add_channels("allele_frequency");
      opts.mutable_pic_options()->set_num_channels(
          opts.pic_options().num_channels() + 1);
    }
    int tumor_h  = absl::GetFlag(FLAGS_pileup_image_height_tumor);
    int normal_h = absl::GetFlag(FLAGS_pileup_image_height_normal);
    if (tumor_h  <= 0) tumor_h  = 100;  // dv_constants.PILEUP_DEFAULT_HEIGHT
    if (normal_h <= 0) normal_h = 100;
    const double ds_tumor  = absl::GetFlag(FLAGS_downsample_fraction_tumor);
    const double ds_normal = absl::GetFlag(FLAGS_downsample_fraction_normal);
    const std::string tumor_name  = absl::GetFlag(FLAGS_sample_name_tumor);
    const std::string normal_name = absl::GetFlag(FLAGS_sample_name_normal);

    auto add_somatic_sample = [&](const std::string& role,
                                    const std::string& name,
                                    const std::string& reads, int height,
                                    double ds, std::initializer_list<int> order,
                                    bool skip_output, bool is_tumor) {
      SampleOptions* s = opts.add_sample_options();
      s->set_role(role);
      s->set_name(name);
      if (!reads.empty()) s->add_reads_filenames(reads);
      s->set_pileup_height(height);
      *s->mutable_variant_caller_options() = vc_opts;
      s->mutable_variant_caller_options()->set_sample_name(name);
      // Somatic mirrors make_examples_somatic.py:149:
      //   FLAGS.set_default('vsc_min_fraction_multiplier', float('inf'))
      // The infinity makes the joint-promotion threshold infeasible, so
      // candidates only get promoted when the TARGET sample's own
      // VAF >= min_fraction (i.e. no normal-only candidates leak into
      // the tumor list). The std::numeric_limits<float>::infinity() value
      // is preserved in the proto float field as a true IEEE infinity.
      s->mutable_variant_caller_options()->set_min_fraction_multiplier(
          std::numeric_limits<float>::infinity());
      // Somatic non-target (normal) AF cap.
      // WGS/WES/PacBio/ONT declare 0.5 in model.example_info.json →
      // cli.cc passes --vsc_max_fraction_snps/indels_for_non_target_sample=0.5.
      // FFPE_WGS/FFPE_WES do NOT declare this flag → stays at -1 (disabled).
      // Without the cap, FFPE emits germline-het candidates and GERMLINE-filters
      // them in postprocess (the correct Docker behaviour).
      {
        const double snp_cap =
            absl::GetFlag(FLAGS_vsc_max_fraction_snps_for_non_target_sample);
        const double ind_cap =
            absl::GetFlag(FLAGS_vsc_max_fraction_indels_for_non_target_sample);
        if (snp_cap >= 0.0)
          s->mutable_variant_caller_options()
              ->set_max_fraction_snps_for_non_target_sample(
                  static_cast<float>(snp_cap));
        if (ind_cap >= 0.0)
          s->mutable_variant_caller_options()
              ->set_max_fraction_indels_for_non_target_sample(
                  static_cast<float>(ind_cap));
      }
      // Adjacent VAF context window for the small_model. DeepSomatic
      // WGS uses 51; the small_model is trained with a 51-position
      // VAF context block. Used by variant_calling_multisample.cc:1160
      // → AddAdjacentAlleleFractionsAtPosition.
      s->mutable_variant_caller_options()
          ->set_small_model_vaf_context_window_size(51);
      for (int o : order) s->add_order(o);
      s->set_skip_output_generation(skip_output);
      if (is_tumor) {
        if (!absl::GetFlag(FLAGS_small_model_path_somatic).empty()) {
          s->set_small_model_path(
              absl::GetFlag(FLAGS_small_model_path_somatic));
        }
      }
      if (ds > 0.0) s->set_downsample_fraction(static_cast<float>(ds));
    };

    if (has_normal) {
      // Order in samples_in_order: [normal(0), tumor(1)]. Tumor's
      // sample.options.order = [0, 1] places normal first in the pileup
      // stack — mirrors make_examples_somatic.py:198. Normal does NOT
      // set order (upstream only assigns order on the tumor branch).
      add_somatic_sample("normal", normal_name.empty() ? "normal" : normal_name,
                          normal_reads, normal_h, ds_normal,
                          {}, /*skip_output=*/true, /*is_tumor=*/false);
      add_somatic_sample("tumor", tumor_name.empty() ? "tumor" : tumor_name,
                          tumor_reads, tumor_h, ds_tumor,
                          {0, 1}, /*skip_output=*/false, /*is_tumor=*/true);
      opts.set_main_sample_index(1);  // tumor at index 1
    } else {
      // Tumor-only: single sample at index 0, order=[0].
      add_somatic_sample("tumor", tumor_name.empty() ? "tumor" : tumor_name,
                          tumor_reads, tumor_h, ds_tumor,
                          {0}, /*skip_output=*/false, /*is_tumor=*/true);
      opts.set_main_sample_index(0);
    }
    opts.set_sample_role_to_train("tumor");
  } else if (!absl::GetFlag(FLAGS_reads_pangenome).empty()) {
    // ──────────────── Pangenome-aware DV mode ────────────────
    // samples_in_order = [pangenome(0), reads(1)]; reads = main.
    // Mirrors deepvariant/make_examples_pangenome_aware_dv.py:
    //   reads_and_pangenome_samples_from_flags (line 207-287).
    //
    // Pangenome example_info.json:flags_for_calling per
    // /opt/models/pangenome_aware_deepvariant/wgs/model.example_info.json:
    //   keep_legacy_allele_counter_behavior: true
    //   keep_only_window_spanning_haplotypes: true
    //   keep_supplementary_alignments: true
    //   min_mapping_quality: 0
    //   normalize_reads: true
    //   pileup_image_height_pangenome: 100
    //   pileup_image_height_reads: 100
    //   pileup_image_width: 221
    //   sort_by_haplotypes: true
    //   trim_reads_for_pileup: true
    //   dbg_disable_graph_pruning: true
    //   aln_match=2 / aln_mismatch=5 / aln_gap_open=10 / aln_gap_extend=1
    // Of these, the per-sample-affecting ones are applied below; pic-
    // level (sort_by_haplotypes, trim_reads_for_pileup) are applied on
    // opts.pic_options; opts-level normalize_reads is set on opts itself.
    opts.mutable_pic_options()->set_sort_by_haplotypes(true);
    opts.set_trim_reads_for_pileup(true);
    // normalize_reads is on AlleleCounterOptions, not MakeExamplesOptions.
    opts.mutable_allele_counter_options()->set_normalize_reads(true);
    // keep_legacy_allele_counter_behavior=true → AlleleCounterOptions.
    // keep_legacy_behavior=true. When true, indel bases below min_base_quality
    // cause the indel to be skipped (stricter than the new sum-of-quality
    // gate); see allelecounter.cc:215.
    opts.mutable_allele_counter_options()->set_keep_legacy_behavior(true);
    // keep_supplementary_alignments=true → ReadRequirements field. Pangenome
    // expects supplementary alignments (HPRC haplotypes can have them) to
    // be retained.
    opts.mutable_allele_counter_options()->mutable_read_requirements()
        ->set_keep_supplementary_alignments(true);

    const std::string pangenome_reads = absl::GetFlag(FLAGS_reads_pangenome);
    const std::string main_reads      = absl::GetFlag(FLAGS_reads);
    int pangenome_h = absl::GetFlag(FLAGS_pileup_image_height_pangenome);
    int reads_h     = absl::GetFlag(FLAGS_pileup_image_height_reads);
    if (pangenome_h <= 0) pangenome_h = 100;
    if (reads_h <= 0) reads_h = 100;
    const double ds_reads = absl::GetFlag(FLAGS_downsample_fraction_reads);
    const std::string reads_name      = absl::GetFlag(FLAGS_sample_name_reads);
    const std::string pangenome_name  = absl::GetFlag(FLAGS_sample_name_pangenome);

    // Reads sample (index 1, main): order=[0,1] (pangenome first, then reads).
    {
      SampleOptions* s = opts.add_sample_options();
      s->set_role("reads");
      s->set_name(reads_name.empty() ? sample_name : reads_name);
      if (!main_reads.empty()) s->add_reads_filenames(main_reads);
      s->set_pileup_height(reads_h);
      *s->mutable_variant_caller_options() = vc_opts;
      s->mutable_variant_caller_options()->set_sample_name(
          reads_name.empty() ? sample_name : reads_name);
      // Mirror upstream Python:
      //   FLAGS.set_default('vsc_min_fraction_multiplier', float('inf'))
      // — drop candidates from the non-target (pangenome) sample.
      s->mutable_variant_caller_options()->set_min_fraction_multiplier(
          std::numeric_limits<float>::infinity());
      s->add_order(0);
      s->add_order(1);
      if (ds_reads > 0.0) s->set_downsample_fraction(static_cast<float>(ds_reads));
      if (!absl::GetFlag(FLAGS_small_model_path_pangenome).empty()) {
        s->set_small_model_path(
            absl::GetFlag(FLAGS_small_model_path_pangenome));
      }
    }
    // Pangenome sample (index 0, non-target): blank channels, skip_phasing,
    // skip_normalization, keep_only_window_spanning_reads.
    {
      SampleOptions* s = opts.add_sample_options();
      s->set_role("pangenome");
      s->set_name(pangenome_name);
      s->add_reads_filenames(pangenome_reads);
      s->set_pileup_height(pangenome_h);
      *s->mutable_variant_caller_options() = vc_opts;
      s->mutable_variant_caller_options()->set_sample_name(pangenome_name);
      s->mutable_variant_caller_options()->set_min_fraction_multiplier(
          std::numeric_limits<float>::infinity());
      s->set_skip_output_generation(true);
      s->set_keep_only_window_spanning_reads(true);
      s->set_skip_phasing(true);
      s->set_skip_normalization(true);
      // Pangenome-aware DV is WGS-only (no PacBio/ONT) → alt_aligned
      // is always "none" per upstream make_examples_pangenome_aware_dv.py.
      s->set_alt_aligned_pileup("none");
      // Per upstream make_examples_pangenome_aware_dv.py:250-256,
      // pangenome rows zero out 5 channels: HAPLOTYPE_TAG (channel 8),
      // DIFF_CHANNELS_ALTERNATE_ALLELE_1 (15), _2 (16), BASE_QUALITY (2),
      // MAPPING_QUALITY (3). Enum values are from
      // deepvariant.proto:DeepVariantChannelEnum.
      s->add_channels_enum_to_blank(::learning::genomics::deepvariant::
                                       CH_HAPLOTYPE_TAG);
      s->add_channels_enum_to_blank(::learning::genomics::deepvariant::
                                       CH_DIFF_CHANNELS_ALTERNATE_ALLELE_1);
      s->add_channels_enum_to_blank(::learning::genomics::deepvariant::
                                       CH_DIFF_CHANNELS_ALTERNATE_ALLELE_2);
      s->add_channels_enum_to_blank(::learning::genomics::deepvariant::
                                       CH_BASE_QUALITY);
      s->add_channels_enum_to_blank(::learning::genomics::deepvariant::
                                       CH_MAPPING_QUALITY);
    }
    // Note: Python pushes [pangenome, reads] but we build [reads, pangenome]
    // because main_sample_index=1 must point at reads. The Python
    // samples_in_order list builds them in [pangenome, reads] order with
    // PANGENOME_SAMPLE_INDEX=0, MAIN_SAMPLE_INDEX=1; we match that
    // ordering by swapping our additions. Re-order:
    auto* mut = opts.mutable_sample_options();
    if (mut->size() == 2) std::swap(*mut->Mutable(0), *mut->Mutable(1));
    opts.set_main_sample_index(1);  // reads at index 1
    opts.set_sample_role_to_train("reads");
  } else {
    SampleOptions* sopt = opts.add_sample_options();
    sopt->set_role("sample");
    sopt->set_name(sample_name);
    sopt->add_reads_filenames(absl::GetFlag(FLAGS_reads));
    sopt->set_pileup_height(100);  // WGS default pileup height per sample.
    *sopt->mutable_variant_caller_options() = vc_opts;
    sopt->mutable_variant_caller_options()->set_sample_name(sample_name);
    opts.set_main_sample_index(0);
    opts.set_sample_role_to_train("sample");
  }

  opts.set_variant_caller(MakeExamplesOptions::VERY_SENSITIVE_CALLER);
  // realigner_enabled and phase_reads are now driven by flags.
  opts.set_realigner_enabled(absl::GetFlag(FLAGS_realigner_enabled));
  opts.set_phase_reads(absl::GetFlag(FLAGS_phase_reads));
  opts.set_stream_examples(false);

  return opts;
}

// Returns true when --reads_parent1 was set (trio mode active).
bool IsTrioMode() {
  return !absl::GetFlag(FLAGS_reads_parent1).empty();
}

// Returns true when --reads_tumor was set (somatic mode active).
bool IsSomaticMode() {
  return !absl::GetFlag(FLAGS_reads_tumor).empty();
}

// Returns true when somatic mode AND --reads_normal is also set.
bool IsSomaticTumorNormalMode() {
  return IsSomaticMode() && !absl::GetFlag(FLAGS_reads_normal).empty();
}

// Returns true when --reads_pangenome was set (pangenome-aware DV mode).
bool IsPangenomeMode() {
  return !absl::GetFlag(FLAGS_reads_pangenome).empty();
}

// DefaultRealignerOptions() with per-flag overrides. Lets pangenome
// supply its aln_match=2/aln_mismatch=5/aln_gap_open=10/aln_gap_extend=1
// + dbg_disable_graph_pruning=true via command-line flags.
::learning::genomics::deepvariant::RealignerOptions
RealignerOptionsFromFlags() {
  auto opts = DefaultRealignerOptions();
  opts.mutable_aln_config()->set_match(absl::GetFlag(FLAGS_aln_match));
  opts.mutable_aln_config()->set_mismatch(absl::GetFlag(FLAGS_aln_mismatch));
  opts.mutable_aln_config()->set_gap_open(absl::GetFlag(FLAGS_aln_gap_open));
  opts.mutable_aln_config()->set_gap_extend(absl::GetFlag(FLAGS_aln_gap_extend));
  // BUG FIX (Path D Site 1, chr12:62946475 1-read-off, 2026-05-23):
  // Mirror upstream `realigner.py:_realigner_options` (lines 420-429):
  // when --normalize_reads is true (we hardcode this true on
  // AlleleCounterOptions at line ~821), the RealignerOptions.normalize_reads
  // must also be true so FastPassAligner does NOT discard realigned
  // alignments whose CIGAR is not left-normalized. Without this, reads
  // in T-homopolymer regions (e.g. chr12:62946475 GTTTT>G in a 16-T run)
  // whose realigned CIGAR has any shiftable indel get thrown out by
  // `fast_pass_aligner.cc:557-568 IsAlignmentNormalized()` check,
  // leaving them at their original POS — losing the +1 DP contribution
  // that Docker counts (DP=27 vs ours DP=26 at this site WG-wide).
  opts.set_normalize_reads(true);
  if (absl::GetFlag(FLAGS_dbg_disable_graph_pruning)) {
    // Match upstream make_examples_core.py: dbg_disable_graph_pruning=true
    // dispatches to PruneLite() (debruijn_graph.cc:257-258), which only
    // removes orphan vertices instead of unreachable + low-weight edges.
    // Critical for pangenome at sites with adjacent insertions: keeping
    // low-weight haplotypes lets reads supporting the simple SNP
    // realign correctly instead of being absorbed by the long insertion
    // haplotype.
    opts.mutable_dbg_config()->set_disable_graph_pruning(true);
  }
  return opts;
}

// Port of upstream realigner.py:split_reads (called from realign_reads when
// --split_skip_reads is set, the RNA-seq default). Splits any read whose CIGAR
// contains a SKIP (N) operation — i.e. a spliced RNA read spanning an intron —
// into separate sub-reads, one per exonic segment, dropping the N gap. Each
// segment ≥ _MIN_SPLIT_LEN (15) aligned bases is retained, with its own start
// position and a `_p<part>` fragment-name suffix (mirrors copy_read). Without
// this, intron-spanning reads inflate the pileup with phantom reference/deletion
// evidence across the intron, degrading the pileup image so the big model emits
// ~homref (QUAL≈0.1 → NoCall) where Docker calls PASS. The native realigner set
// realigner_options.split_skip_reads=true but never acted on it; this restores
// the behavior. Constants/op-sets mirror nucleus/util/cigar.py.
static std::vector<nucleus::genomics::v1::Read> SplitReadsOnSkip(
    const std::vector<nucleus::genomics::v1::Read>& reads) {
  namespace ng = nucleus::genomics::v1;
  using ng::CigarUnit;
  constexpr int kMinSplitLen = 15;
  auto is_ref_adv = [](int op) {
    return op == CigarUnit::ALIGNMENT_MATCH || op == CigarUnit::SEQUENCE_MATCH ||
           op == CigarUnit::DELETE || op == CigarUnit::SKIP ||
           op == CigarUnit::SEQUENCE_MISMATCH;
  };
  auto is_read_adv = [](int op) {
    return op == CigarUnit::ALIGNMENT_MATCH || op == CigarUnit::SEQUENCE_MATCH ||
           op == CigarUnit::INSERT || op == CigarUnit::CLIP_SOFT ||
           op == CigarUnit::SEQUENCE_MISMATCH;
  };
  std::vector<ng::Read> out;
  out.reserve(reads.size());
  for (const auto& read : reads) {
    bool has_skip = false;
    for (const auto& c : read.alignment().cigar())
      if (c.operation() == CigarUnit::SKIP) { has_skip = true; break; }
    if (!has_skip) { out.push_back(read); continue; }

    int part = 0, read_start = 0, read_offset = 0, reference_offset = 0;
    auto make_part = [&](int p) {
      ng::Read nr;
      nr.CopyFrom(read);
      nr.clear_alignment();
      nr.clear_aligned_sequence();
      nr.clear_aligned_quality();
      auto* pos = nr.mutable_alignment()->mutable_position();
      pos->set_reference_name(read.alignment().position().reference_name());
      pos->set_reverse_strand(read.alignment().position().reverse_strand());
      nr.mutable_alignment()->set_mapping_quality(
          read.alignment().mapping_quality());
      nr.set_fragment_name(absl::StrCat(read.fragment_name(), "_p", p));
      return nr;
    };
    ng::Read new_read = make_part(part);
    const int ncig = read.alignment().cigar_size();
    for (int n = 0; n < ncig; ++n) {
      const auto& cig = read.alignment().cigar(n);
      const bool on_last = (n + 1 == ncig);
      const int op = cig.operation();
      if (is_ref_adv(op)) {
        if (new_read.alignment().position().position() == 0) {
          new_read.mutable_alignment()->mutable_position()->set_position(
              read.alignment().position().position() + reference_offset);
        }
        reference_offset += cig.operation_length();
      }
      if (is_read_adv(op)) read_offset += cig.operation_length();
      if (op != CigarUnit::SKIP) *new_read.mutable_alignment()->add_cigar() = cig;
      if (op == CigarUnit::SKIP || on_last) {
        new_read.set_aligned_sequence(
            read.aligned_sequence().substr(read_start, read_offset - read_start));
        new_read.set_aligned_quality(
            read.aligned_quality().substr(read_start, read_offset - read_start));
        if (static_cast<int>(new_read.aligned_sequence().size()) >= kMinSplitLen)
          out.push_back(new_read);
        if (!on_last) {
          read_start = read_offset;
          ++part;
          new_read = make_part(part);
        }
      }
    }
  }
  return out;
}

// Infer sample name from the first RG:SM field in the BAM header.
std::string InferSampleName(
    const nucleus::genomics::v1::SamHeader& header) {
  for (const auto& rg : header.read_groups()) {
    if (!rg.sample_id().empty()) return rg.sample_id();
  }
  return "sample";
}

// Walk a 51-bp window of AlleleCounts around the candidate and populate
// the candidate's allele_frequency_at_position map with VAF (×100, integer)
// at each position. The map is used by the small_model's VAF-context
// features (offsets −25..+25 around the variant).
void PopulateVafContext(
    DeepVariantCall* candidate,
    const std::vector<AlleleCount>& allele_counts) {
  if (allele_counts.empty()) return;
  const int64_t variant_pos = candidate->variant().start();
  const int64_t region_start = allele_counts.front().position().position();
  const int64_t local_idx = variant_pos - region_start;
  constexpr int kHalfWindow = kSmallModelVafContextWindow / 2;  // 25
  for (int o = -kHalfWindow; o <= kHalfWindow; ++o) {
    const int64_t idx = local_idx + o;
    if (idx < 0 || idx >= static_cast<int64_t>(allele_counts.size())) continue;
    const auto& ac = allele_counts[idx];
    const int depth = ac.ref_supporting_read_count() + ac.read_alleles_size();
    const int vaf = depth > 0 ? (100 * ac.read_alleles_size()) / depth : 0;
    (*candidate->mutable_allele_frequency_at_position())[
        ac.position().position()] = vaf;
  }
}

// Returns true if the (alt_idx-only) sub-variant is a SNP — used to pick
// the small_model GQ threshold (snp=20 vs indel=28).
bool IsSnpAlt(const nucleus::genomics::v1::Variant& v, int alt_idx) {
  if (alt_idx < 0 || alt_idx >= v.alternate_bases_size()) return false;
  return v.reference_bases().size() == 1 &&
         v.alternate_bases(alt_idx).size() == 1;
}

// Multi-index version: SNP iff REF is 1 base AND every alt in
// `alt_indices` is 1 base. Mirror of nucleus/util/variant_utils.is_snp(
// variant, exclude_alleles) where exclude_alleles is the complement of
// alt_indices.
bool IsSnpForIndices(const nucleus::genomics::v1::Variant& v,
                      const std::vector<int>& alt_indices) {
  if (alt_indices.empty()) return false;
  if (v.reference_bases().size() != 1) return false;
  for (int idx : alt_indices) {
    if (idx < 0 || idx >= v.alternate_bases_size()) return false;
    if (v.alternate_bases(idx).size() != 1) return false;
  }
  return true;
}

// ---------------------------------------------------------------------------
// select_variant_types candidate filtering
//
// Faithful port of make_examples_core.filter_candidates +
// nucleus/util/variant_utils. The variant-type predicates exclude the gVCF
// '<*>' allele, the '<NON_REF>' symbolic allele, and the '.' missing field
// from the alt set (variant_utils._non_excluded_alts) before classifying.
// ---------------------------------------------------------------------------

// True for alts ignored by the type predicates (variant_utils default set).
bool IsExcludedAlt(const std::string& alt) {
  return alt == "<*>" || alt == "<NON_REF>" || alt == ".";
}

// A methylated reference site: a candidate whose only alt is the missing-field
// '.', i.e. a reference position retained as a candidate because it was checked
// for methylation. Mirror of make_examples_core._is_methylated_reference_site.
bool IsMethylatedRefSite(const DeepVariantCall& c) {
  return c.variant().alternate_bases_size() == 1 &&
         c.variant().alternate_bases(0) == ".";
}

// Alt alleles that count toward variant-type classification (non-excluded).
std::vector<const std::string*> RelevantAlts(
    const nucleus::genomics::v1::Variant& v) {
  std::vector<const std::string*> alts;
  for (const auto& a : v.alternate_bases()) {
    if (!IsExcludedAlt(a)) alts.push_back(&a);
  }
  return alts;
}

// is_snp: REF is 1 bp and every non-excluded alt is 1 bp (>=1 such alt).
bool VarIsSnp(const nucleus::genomics::v1::Variant& v) {
  const auto alts = RelevantAlts(v);
  if (v.reference_bases().size() != 1 || alts.empty()) return false;
  for (const auto* a : alts) {
    if (a->size() != 1) return false;
  }
  return true;
}

// is_indel: at least one non-excluded alt, and REF>1 or some alt>1 bp.
bool VarIsIndel(const nucleus::genomics::v1::Variant& v) {
  const auto alts = RelevantAlts(v);
  if (alts.empty()) return false;
  if (v.reference_bases().size() > 1) return true;
  for (const auto* a : alts) {
    if (a->size() > 1) return true;
  }
  return false;
}

bool VarIsBiallelic(const nucleus::genomics::v1::Variant& v) {
  return RelevantAlts(v).size() == 1;
}

bool VarIsMultiallelic(const nucleus::genomics::v1::Variant& v) {
  return RelevantAlts(v).size() > 1;
}

// has_insertion/has_deletion gate on is_indel but, like variant_utils, test
// the length condition over ALL alternate_bases (not just the non-excluded).
bool VarHasInsertion(const nucleus::genomics::v1::Variant& v) {
  if (!VarIsIndel(v)) return false;
  const size_t ref_len = v.reference_bases().size();
  for (const auto& a : v.alternate_bases()) {
    if (ref_len < a.size()) return true;
  }
  return false;
}

bool VarHasDeletion(const nucleus::genomics::v1::Variant& v) {
  if (!VarIsIndel(v)) return false;
  const size_t ref_len = v.reference_bases().size();
  for (const auto& a : v.alternate_bases()) {
    if (ref_len > a.size()) return true;
  }
  return false;
}

// Parsed --select_variant_types flag. `active` is false when no selector was
// requested, in which case filtering is a no-op (keep all candidates).
struct SelectedVariantTypes {
  bool active = false;
  bool snps = false;
  bool indels = false;
  bool insertions = false;
  bool deletions = false;
  bool multiallelics = false;
  bool all = false;
};

// Parse the whitespace-separated flag. Mirrors make_examples_options.py: an
// unknown selector is a fatal command-line error (returns false + *err set).
bool ParseSelectVariantTypes(const std::string& flag,
                             SelectedVariantTypes* out, std::string* err) {
  for (absl::string_view tok :
       absl::StrSplit(flag, absl::ByAnyChar(" \t\r\n"), absl::SkipEmpty())) {
    if (tok == "snps") {
      out->snps = true;
    } else if (tok == "indels") {
      out->indels = true;
    } else if (tok == "insertions") {
      out->insertions = true;
    } else if (tok == "deletions") {
      out->deletions = true;
    } else if (tok == "multi-allelics") {
      out->multiallelics = true;
    } else if (tok == "all") {
      out->all = true;
    } else {
      *err = absl::StrCat(
          "Select variant type '", tok,
          "' not recognized. Allowed values are snps, indels, insertions, "
          "deletions, multi-allelics, all");
      return false;
    }
    out->active = true;
  }
  return true;
}

// Does the candidate's variant match any requested selector (OR'd)?
// snps/indels/insertions/deletions are bi-allelic-gated, matching the
// VARIANT_TYPE_SELECTORS table in make_examples_core.
bool CandidateSelected(const nucleus::genomics::v1::Variant& v,
                       const SelectedVariantTypes& sel) {
  if (sel.all) return true;
  if (sel.snps && VarIsSnp(v) && VarIsBiallelic(v)) return true;
  if (sel.indels && VarIsIndel(v) && VarIsBiallelic(v)) return true;
  if (sel.insertions && VarHasInsertion(v) && VarIsBiallelic(v)) return true;
  if (sel.deletions && VarHasDeletion(v) && VarIsBiallelic(v)) return true;
  if (sel.multiallelics && VarIsMultiallelic(v)) return true;
  return false;
}

// Drop candidates whose variant matches no requested selector, preserving
// order. No-op when no selector is active.
void FilterCandidatesBySelectedTypes(std::vector<DeepVariantCall>* candidates,
                                     const SelectedVariantTypes& sel) {
  if (!sel.active) return;
  candidates->erase(
      std::remove_if(candidates->begin(), candidates->end(),
                     [&](const DeepVariantCall& c) {
                       return !CandidateSelected(c.variant(), sel);
                     }),
      candidates->end());
}

// Phred = -10 * log10(p), truncated toward zero. Capped at 99.
//
// Truncation (not std::round) matches upstream's small_model
// passes_confidence_threshold(ptrue_to_bounded_phred(max_p) >= threshold)
// at the boundary: a phred of 19.5 should *fail* a threshold of 20 (which
// floor-rounds it down to 19), but std::round would push 19.5 up to 20
// and pass — flipping a candidate from big-model dispatch to a
// small_model emit.
int ProbToPhred(double p) {
  if (p <= 0.0) return 99;
  if (p >= 1.0) return 0;
  return std::min(static_cast<int>(-10.0 * std::log10(p)), 99);
}

// Build a CallVariantsOutput proto for a single (candidate, alt_idx) pair
// that the small model has resolved. We tag MID="small_model" in the
// VariantCall.info so postprocess can propagate it to the VCF.
// `alt_indices` may be a single index (single-alt CVO) or two indices
// (multi-alt combo CVO, mirrors upstream's get_set_of_allele_indices
// `multiallelic = combinations(range(N), 2)`).
CallVariantsOutput MakeSmallModelCvo(
    const DeepVariantCall& candidate, const std::vector<int>& alt_indices,
    const float* probs) {
  CallVariantsOutput cvo;
  *cvo.mutable_variant() = candidate.variant();
  for (int idx : alt_indices) cvo.mutable_alt_allele_indices()->add_indices(idx);
  // Probabilities written as double — same wire-format as the big model.
  for (int i = 0; i < 3; ++i) cvo.add_genotype_probabilities(probs[i]);
  // Tag MID in VariantCall.info["MID"]. variant_calling.cc already adds an
  // empty VariantCall, so reuse that slot rather than appending another one
  // (would trigger the VcfWriter's "calls != samples" check).
  auto* v = cvo.mutable_variant();
  if (v->calls_size() == 0) v->add_calls();
  nucleus::SetInfoField("MID", std::string("small_model"),
                         v->mutable_calls(0));
  return cvo;
}

}  // namespace

// Per-thread accumulators returned to the main thread for summing.
struct WorkerStats {
  int64_t total_candidates = 0;
  int64_t total_examples = 0;
  int64_t total_small_hits = 0;
  int64_t total_big_dispatched = 0;
};

// FillAlleleFrequencyFromPon — populate dv_call.allele_frequency map from
// a Panel-of-Normals VCF for each candidate.
// Mirrors Python's allele_frequency.add_allele_frequencies_to_candidates.
// For reads supporting an alt allele, AlleleFrequencyChannel reads the
// per-allele population AF from this map to encode the 8th pileup channel.
//
// If a candidate's position is not in the PON, sets ref=1.0, all alts=0.0
// (same as Python's fallback when population_vcf_reader is None).
static void FillAlleleFrequencyFromPon(
    std::vector<DeepVariantCall>& candidates,
    nucleus::VcfReader& pon_reader) {
  using nucleus::genomics::v1::Range;
  using nucleus::genomics::v1::Variant;
  for (auto& c : candidates) {
    const auto& v = c.variant();
    // Clear and set defaults first: ref=1.0, all ALTs=0.0.
    c.mutable_allele_frequency()->clear();
    (*c.mutable_allele_frequency())[v.reference_bases()] = 1.0f;
    for (const auto& alt : v.alternate_bases())
      (*c.mutable_allele_frequency())[alt] = 0.0f;

    Range range;
    range.set_reference_name(v.reference_name());
    range.set_start(v.start());
    range.set_end(v.end());

    auto it_or = pon_reader.Query(range);
    if (!it_or.ok()) continue;
    auto it = it_or.ValueOrDie();

    Variant pon_v;
    while (true) {
      auto next_or = it->Next(&pon_v);
      if (!next_or.ok() || !next_or.ValueOrDie()) break;
      if (pon_v.reference_bases() != v.reference_bases()) continue;

      // Find AF INFO field (per-allele, one value per ALT in PON entry).
      auto af_it = pon_v.info().find("AF");
      if (af_it == pon_v.info().end()) continue;
      const auto& af_vals = af_it->second.values();

      float sum_alt_af = 0.0f;
      for (int i = 0; i < pon_v.alternate_bases_size(); ++i) {
        const std::string& pon_alt = pon_v.alternate_bases(i);
        float af = (i < af_vals.size() && af_vals[i].has_number_value())
                       ? static_cast<float>(af_vals[i].number_value()) : 0.0f;
        // Map only PON alts that match a candidate alt.
        for (const auto& cand_alt : v.alternate_bases()) {
          if (pon_alt == cand_alt) {
            (*c.mutable_allele_frequency())[cand_alt] = af;
            sum_alt_af += af;
          }
        }
      }
      // Recompute ref AF = 1 - sum(matched alt AFs).
      (*c.mutable_allele_frequency())[v.reference_bases()] =
          std::max(0.0f, 1.0f - sum_alt_af);
      break;  // Use first matching PON entry at this position.
    }
  }
}

int RunMakeExamples(int argc, char** argv) {
  absl::ParseCommandLine(argc, argv);

  const std::string reads_path = absl::GetFlag(FLAGS_reads);
  const std::string ref_path = absl::GetFlag(FLAGS_ref);
  const std::string examples_path = absl::GetFlag(FLAGS_examples);

  if (ref_path.empty()) {
    LOG(ERROR) << "Required: --ref";
    return 1;
  }

  // Validate --select_variant_types up front so a typo fails fast rather than
  // silently keeping everything (the pre-fix behavior) or erroring mid-region.
  SelectedVariantTypes selected_types;
  {
    std::string err;
    if (!ParseSelectVariantTypes(absl::GetFlag(FLAGS_select_variant_types),
                                 &selected_types, &err)) {
      LOG(ERROR) << err;
      return 1;
    }
  }

  // Sex-chromosome haploid calling config — drives haploid gVCF reference
  // confidence on the --haploid_contigs (outside the PAR).
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
  if (reads_path.empty() && !IsSomaticMode()) {
    LOG(ERROR) << "Required: --reads (or --reads_tumor for somatic mode)";
    return 1;
  }
  // Trio mode: at least one of --examples / --examples_child must be set.
  if (IsTrioMode()) {
    const std::string ex_child = absl::GetFlag(FLAGS_examples_child);
    if (examples_path.empty() && ex_child.empty()) {
      LOG(ERROR)
          << "Trio mode requires --examples_child (or --examples as alias).";
      return 1;
    }
  } else if (IsSomaticMode()) {
    // Somatic: tumor's examples are mandatory; normal has skip_output=true.
    const std::string ex_tumor = absl::GetFlag(FLAGS_examples_tumor);
    if (examples_path.empty() && ex_tumor.empty()) {
      LOG(ERROR)
          << "Somatic mode requires --examples_tumor (or --examples as alias).";
      return 1;
    }
  } else if (IsPangenomeMode()) {
    // Pangenome: reads' examples are mandatory; pangenome has skip_output=true.
    const std::string ex_reads = absl::GetFlag(FLAGS_examples_reads);
    if (examples_path.empty() && ex_reads.empty()) {
      LOG(ERROR)
          << "Pangenome mode requires --examples_reads (or --examples).";
      return 1;
    }
  } else if (examples_path.empty()) {
    LOG(ERROR) << "Required: --examples";
    return 1;
  }

  const int task_id = absl::GetFlag(FLAGS_task_id);
  const int num_shards = std::max(1, absl::GetFlag(FLAGS_num_shards));
  const int n_threads = std::max(1, absl::GetFlag(FLAGS_threads));

  // ── Open shared reference (only used for header / contigs / sample name
  //     inference). Per-thread workers reopen their own IndexedFastaReader
  //     so AlleleCounter calls into htslib stay thread-local. ──────────────
  auto ref_or = nucleus::IndexedFastaReader::FromFile(
      ref_path, absl::StrCat(ref_path, ".fai"));
  CHECK(ref_or.ok()) << "Failed to open reference: " << ref_path;
  auto ref_reader_main = std::move(ref_or.ValueOrDie());

  // ── Infer sample name from the BAM header (cheap, single read). ──────────
  // For somatic mode, use the tumor BAM (no --reads needed); for trio we
  // also infer from the (child) --reads. Fall back to a default name if
  // the user passed neither.
  nucleus::genomics::v1::SamReaderOptions sam_opts;
  sam_opts.mutable_read_requirements()->set_min_mapping_quality(
      absl::GetFlag(FLAGS_min_mapping_quality));
  // Phase 5.5d/15 — propagate keep_supplementary_alignments to SamReader.
  // Without this, sam_reader.cc::PartialReadSatisfiesRequirements rejects
  // supplementary alignments at the BAM-read source level, regardless of
  // what we set later on `read_reqs` (which only flows into AlleleCounter
  // and PileupImage). For PACBIO/ONT, supplementary reads carry
  // significant pileup depth at chimeric-alignment regions; dropping
  // them at the source produced 8-10× DP underflow vs Docker
  // (e.g. chr20:62642 our DP=7 vs Docker DP=55) and missed candidates
  // entirely (1342 Docker-only PASS sites including homopolymer indels).
  sam_opts.mutable_read_requirements()->set_keep_supplementary_alignments(
      absl::GetFlag(FLAGS_keep_supplementary_alignments));
  // --parse_sam_aux_fields must reach the SamReader itself: ParseAuxFields is a
  // no-op unless aux_field_handling == PARSE_ALL_AUX_FIELDS, and without it the
  // reader never populates read.info() (MM/ML base modifications, HP, ...).
  // Setting it on MakeExamplesOptions alone (above) is not enough. Default off
  // keeps the byte-identical baseline (no aux fields parsed).
  if (absl::GetFlag(FLAGS_parse_sam_aux_fields)) {
    sam_opts.set_aux_field_handling(
        nucleus::genomics::v1::SamReaderOptions::PARSE_ALL_AUX_FIELDS);
  }
  {
    std::string probe_bam = reads_path;
    if (probe_bam.empty() && IsSomaticMode()) {
      probe_bam = absl::GetFlag(FLAGS_reads_tumor);
    }
    if (!probe_bam.empty()) {
      auto sam_or = nucleus::SamReader::FromFile(probe_bam, sam_opts);
      CHECK(sam_or.ok()) << "Failed to open BAM: " << probe_bam;
      auto tmp_reader = std::move(sam_or.ValueOrDie());
      std::string sn = absl::GetFlag(FLAGS_sample_name);
      if (sn.empty()) {
        sn = InferSampleName(tmp_reader->Header());
        LOG(INFO) << "Inferred sample name: " << sn;
        absl::SetFlag(&FLAGS_sample_name, sn);
      }
    }
  }
  const std::string sample_name = absl::GetFlag(FLAGS_sample_name);

  // ── Build MakeExamplesOptions ─────────────────────────────────────────────
  const MakeExamplesOptions opts = BuildOptions(sample_name, task_id, num_shards);

  // ── Build calling regions ─────────────────────────────────────────────────
  const auto& contigs = ref_reader_main->Contigs();
  std::vector<std::string> inc_regions, exc_regions;
  {
    const std::string regions_str = absl::GetFlag(FLAGS_regions);
    if (!regions_str.empty()) {
      inc_regions = absl::StrSplit(regions_str, absl::ByAnyChar(" \t,"),
                                   absl::SkipEmpty());
    }
    const std::string excl_str = absl::GetFlag(FLAGS_exclude_regions);
    if (!excl_str.empty()) {
      exc_regions = absl::StrSplit(excl_str, absl::ByAnyChar(" \t,"),
                                   absl::SkipEmpty());
    }
  }
  auto all_regions = BuildCallingRegions(contigs, inc_regions, exc_regions);
  // Partition into chunks of partition_size bp (default 1000), then shard.
  // Mirrors upstream's `regions.partition()` step. Required for realigner
  // window-set parity: each chunk runs the WindowSelector + DBG
  // independently, and adjacent chunks emit overlapping windows at the
  // chunk boundary — without partitioning we'd merge windows across chunk
  // boundaries that upstream keeps separate.
  const int64_t partition_size_bp =
      static_cast<int64_t>(absl::GetFlag(FLAGS_partition_size));
  auto partitioned = PartitionRegions(all_regions, partition_size_bp);
  auto shard_regions = ShardRegions(partitioned, task_id, num_shards);

  LOG(INFO) << "Processing " << shard_regions.size() << " regions (shard "
            << task_id << "/" << num_shards << ", threads=" << n_threads
            << ")";

  const std::string small_path = absl::GetFlag(FLAGS_small_model);
  const std::string small_cvo_path =
      absl::GetFlag(FLAGS_small_model_cvo_outfile);
  const int snp_gq_threshold = absl::GetFlag(FLAGS_small_model_snp_gq_threshold);
  const int indel_gq_threshold =
      absl::GetFlag(FLAGS_small_model_indel_gq_threshold);
  if (!small_path.empty() && small_cvo_path.empty()) {
    LOG(ERROR) << "--small_model requires --small_model_cvo_outfile";
    return 1;
  }

  // ── Atomic region cursor: workers fetch_add to claim regions. ────────────
  std::atomic<size_t> next_region{0};

  // Per-thread output paths. We use the standard `name-NNNNN-of-NNNNN`
  // shard naming so downstream stages that already understand the `@N`
  // shard spec (call_variants / postprocess via TFRecordReader) can read
  // the per-thread files directly — no end-of-stage concat needed.
  //
  // examples_path: if it already carries an `@N` suffix, we honour the
  // caller's N; otherwise we synthesise `examples_path@n_threads` and
  // shard from it. n_threads==1 collapses to the plain path.
  std::string examples_spec = examples_path;
  std::string small_cvo_spec = small_cvo_path;
  if (n_threads > 1) {
    if (examples_spec.find('@') == std::string::npos) {
      examples_spec = absl::StrCat(examples_path, "@", n_threads);
    }
    if (!small_cvo_spec.empty() &&
        small_cvo_spec.find('@') == std::string::npos) {
      small_cvo_spec = absl::StrCat(small_cvo_path, "@", n_threads);
    }
  }
  auto thread_examples_path = [&](int t) {
    return n_threads == 1 ? examples_path : ShardName(examples_spec, t);
  };
  auto thread_small_cvo_path = [&](int t) {
    return n_threads == 1 ? small_cvo_path : ShardName(small_cvo_spec, t);
  };
  // Phase 9 / Step 3 — gVCF output sharding. Same pattern as small_cvo.
  const std::string gvcf_path_top = absl::GetFlag(FLAGS_gvcf);
  std::string gvcf_spec = gvcf_path_top;
  if (n_threads > 1 && !gvcf_spec.empty() &&
      gvcf_spec.find('@') == std::string::npos) {
    gvcf_spec = absl::StrCat(gvcf_path_top, "@", n_threads);
  }
  auto thread_gvcf_path_top = [&](int t) {
    return n_threads == 1 ? gvcf_path_top : ShardName(gvcf_spec, t);
  };

  // ──────────────────────────────────────────────────────────────────
  // Trio worker — mirrors deeptrio/make_examples.py's per-region loop.
  // Opens 3 SamReaders, builds 3 AlleleCounters per region, runs
  // multi_sample::VariantCaller once per target sample, generates
  // examples with the target's `order` permutation, and writes per-
  // sample small_cvo + examples streams. The single-sample path
  // below is unchanged (preserves the WGS chr20 100% FILTER parity
  // gate already achieved at 5.5d/10).
  // ──────────────────────────────────────────────────────────────────
  // Multi-sample worker: handles BOTH trio (3 samples: parent1, child,
  // parent2) AND somatic (1-2 samples: tumor[, normal]). Same processing
  // pipeline; only role names + flag plumbing differ.
  auto run_trio_worker = [&](int tid, WorkerStats* out_stats) {
    auto t_ref_or = nucleus::IndexedFastaReader::FromFile(
        ref_path, absl::StrCat(ref_path, ".fai"));
    CHECK(t_ref_or.ok()) << "thread " << tid << ": ref reopen failed";
    auto ref_reader = std::move(t_ref_or.ValueOrDie());

    // Per-role context. Up to 3 sample slots; trio uses all 3 (parent1,
    // child, parent2), somatic uses 1-2 (tumor[, normal]).
    struct SampleCtx {
      std::string role;
      std::string name;
      std::vector<int> order;          // pileup channel-stack permutation
      int pileup_height = 100;
      bool skip_output = false;
      std::unique_ptr<nucleus::SamReader> sam_reader;
      std::unique_ptr<SmallModel> small_model;
      std::unique_ptr<TFRecordWriter> small_cvo_writer;
      std::string examples_path;
      // Per-target call_variants_outputs counters reported back as stats.
      int64_t total_candidates = 0;
      int64_t total_examples = 0;
      int64_t total_small_hits = 0;
      int64_t total_big_dispatched = 0;
    };
    const int n_samples = opts.sample_options_size();
    std::array<SampleCtx, 3> ctx;
    for (int s = 0; s < n_samples; ++s) {
      const auto& so = opts.sample_options(s);
      ctx[s].role = so.role();
      ctx[s].name = so.name();
      ctx[s].pileup_height = so.pileup_height();
      ctx[s].skip_output = so.skip_output_generation();
      for (int o : so.order()) ctx[s].order.push_back(o);
      if (so.reads_filenames_size() > 0) {
        auto sr_or = nucleus::SamReader::FromFile(so.reads_filenames(0),
                                                    sam_opts);
        CHECK(sr_or.ok()) << "thread " << tid << " " << ctx[s].role
                           << ": BAM reopen failed: " << so.reads_filenames(0);
        ctx[s].sam_reader = std::move(sr_or.ValueOrDie());
      }
    }

    // Per-role examples path lookup. Handles both trio roles
    // (parent1/child/parent2) and somatic roles (tumor/normal).
    auto multi_examples_path = [&](const std::string& role) -> std::string {
      std::string base;
      if (role == "child")
        base = absl::GetFlag(FLAGS_examples_child).empty()
                   ? examples_path
                   : absl::GetFlag(FLAGS_examples_child);
      else if (role == "parent1")
        base = absl::GetFlag(FLAGS_examples_parent1);
      else if (role == "parent2")
        base = absl::GetFlag(FLAGS_examples_parent2);
      else if (role == "tumor")
        base = absl::GetFlag(FLAGS_examples_tumor).empty()
                   ? examples_path
                   : absl::GetFlag(FLAGS_examples_tumor);
      else if (role == "normal")
        base = absl::GetFlag(FLAGS_examples_normal);
      else if (role == "reads")
        base = absl::GetFlag(FLAGS_examples_reads).empty()
                   ? examples_path
                   : absl::GetFlag(FLAGS_examples_reads);
      else if (role == "pangenome")
        base = absl::GetFlag(FLAGS_examples_pangenome);
      if (base.empty()) return "";
      return n_threads == 1 ? base : ShardName(base, tid);
    };
    auto multi_small_cvo_path = [&](const std::string& role) -> std::string {
      std::string base;
      if (role == "child")
        base = absl::GetFlag(FLAGS_small_model_cvo_outfile_child);
      else if (role == "parent1")
        base = absl::GetFlag(FLAGS_small_model_cvo_outfile_parent1);
      else if (role == "parent2")
        base = absl::GetFlag(FLAGS_small_model_cvo_outfile_parent2);
      else if (role == "tumor")
        base = absl::GetFlag(FLAGS_small_model_cvo_outfile_tumor);
      else if (role == "reads")
        base = absl::GetFlag(FLAGS_small_model_cvo_outfile_reads);
      // normal/pangenome: skip_output=true → no CVO
      if (base.empty()) return "";
      return n_threads == 1 ? base : ShardName(base, tid);
    };

    const std::string sm_child     = absl::GetFlag(FLAGS_small_model_path_child);
    const std::string sm_parent    = absl::GetFlag(FLAGS_small_model_path_parent);
    const std::string sm_somatic   = absl::GetFlag(FLAGS_small_model_path_somatic);
    const std::string sm_pangenome = absl::GetFlag(FLAGS_small_model_path_pangenome);

    auto sm_path_for_role = [&](const std::string& role) -> const std::string& {
      static const std::string empty;
      if (role == "child")   return sm_child;
      if (role == "parent1" || role == "parent2") return sm_parent;
      if (role == "tumor")   return sm_somatic;
      if (role == "reads")   return sm_pangenome;
      return empty;
    };

    std::unordered_map<std::string, std::string> example_filenames;
    for (int s = 0; s < n_samples; ++s) {
      auto& c = ctx[s];
      c.examples_path = multi_examples_path(c.role);
      if (!c.skip_output && !c.examples_path.empty()) {
        example_filenames[c.role] = c.examples_path;
      }
      const std::string& sm_path = sm_path_for_role(c.role);
      if (!sm_path.empty() && !c.skip_output) {
        c.small_model = SmallModel::Load(sm_path);
        CHECK(c.small_model) << "thread " << tid << " " << c.role
                              << ": small_model load failed: " << sm_path;
        const std::string scp = multi_small_cvo_path(c.role);
        if (!scp.empty()) {
          c.small_cvo_writer = TFRecordWriter::New(scp);
          CHECK(c.small_cvo_writer)
              << "thread " << tid << " " << c.role
              << ": small CVO writer open failed: " << scp;
        }
      }
    }

    // Per-thread PON VcfReader for tumor-only allele_frequency channel.
    // Opened once per thread (VcfReader is NOT thread-safe — each thread
    // needs its own handle). Empty path → pon_reader stays null → defaults.
    std::unique_ptr<nucleus::VcfReader> pon_reader;
    {
      const std::string pon_path = absl::GetFlag(FLAGS_population_vcfs);
      if (!pon_path.empty()) {
        nucleus::genomics::v1::VcfReaderOptions pon_opts;
        auto pon_or = nucleus::VcfReader::FromFile(pon_path, pon_opts);
        CHECK(pon_or.ok()) << "thread " << tid
                           << ": PON VCF open failed: " << pon_path;
        pon_reader = std::move(pon_or.ValueOrDie());
      }
    }

    multi_sample::VariantCaller caller(
        opts.sample_options(opts.main_sample_index()).variant_caller_options());

    ExamplesGenerator generator(opts, example_filenames);

    while (true) {
      const size_t i = next_region.fetch_add(1, std::memory_order_relaxed);
      if (i >= shard_regions.size()) break;
      const auto& region = shard_regions[i];
      LOG(INFO) << "Trio region: " << region.reference_name() << ":"
                << region.start() << "-" << region.end();

      // Per-sample: query reads, reservoir-sample, run realigner per
      // sample (mirrors upstream's realign_reads_per_sample_multisample
      // — each sample's reads are re-aligned independently against
      // assembled haplotypes; trio joint_realignment is a future
      // optimization but per-sample matches Docker's default).
      std::array<std::vector<nucleus::genomics::v1::Read>, 3> reads_per_sample_v;
      const int max_rpp = static_cast<int>(opts.max_reads_per_partition());
      const bool realigner_enabled = absl::GetFlag(FLAGS_realigner_enabled);
      for (int s = 0; s < n_samples; ++s) {
        if (!ctx[s].sam_reader) continue;
        auto reads_or = ctx[s].sam_reader->Query(region);
        if (!reads_or.ok()) {
          LOG(WARNING) << "Query failed for " << ctx[s].role << " "
                       << region.reference_name() << ":" << region.start()
                       << "-" << region.end() << " — " << reads_or.status();
          continue;
        }
        auto& reads_iter = reads_or.ValueOrDie();
        std::vector<nucleus::genomics::v1::Read> raw_reads;
        nucleus::genomics::v1::Read tmp_read;
        while (true) {
          auto next = reads_iter->Next(&tmp_read);
          if (!next.ok() || !next.ValueOrDie()) break;
          raw_reads.push_back(tmp_read);
        }
        reads_iter->Release().IgnoreError();
        if (max_rpp > 0 && raw_reads.size() > static_cast<size_t>(max_rpp)) {
          // BUG FIX (2026-05-10): the previous stable_sort by
          // (POS, fragment_name, read_number) was added in Phase 5.5d/10
          // as a "shard-count-independence guard", but it CHANGED the
          // input order to reservoir sampling vs Docker. Docker reads
          // BAM-naturally ordered (POS only, secondary by file offset),
          // and our sort by (POS, fragment_name, read_number) reorders
          // same-POS reads → reservoir picks different reads → ±1-4 read
          // DP differences at WG scale on ~79 % of FILTER-mismatch sites.
          // Removed for full Docker compatibility (user directive 2026-05-10).
          ::deepvariant::npr::NumpyMt19937 region_rng(opts.random_seed());
          auto sampled = ::deepvariant::npr::ReservoirSamplePtrs(
              raw_reads, max_rpp, region_rng);
          std::vector<nucleus::genomics::v1::Read> kept;
          kept.reserve(sampled.size());
          for (const auto* p : sampled) kept.push_back(*p);
          raw_reads = std::move(kept);
        }

        // Realign per-sample (Step 1.3-bis). Matches upstream's
        // make_examples_core.py:realign_reads_per_sample_multisample:
        // each sample's reads are reassembled against per-sample
        // de Bruijn graph haplotypes, eliminating misalignment-induced
        // phantom alleles that inflate the AlleleCounter Counts.
        //
        // Pangenome exception: upstream's `can_realign` (make_examples_
        // core.py:2208) returns False for `role == 'pangenome'` — synthetic
        // haplotypes are pre-aligned to the GBZ graph, so re-running our
        // realigner on them produces phantom alt alleles that diverge
        // from Docker's pangenome AlleleCount.
        if (realigner_enabled && ctx[s].role != "pangenome") {
          const auto realigner_opts = RealignerOptionsFromFlags();
          const int expand_bp =
              realigner_opts.ws_config().region_expansion_in_bp();
          auto contig_or = ref_reader->Contig(region.reference_name());
          const int64_t contig_n =
              contig_or.ok() ? contig_or.ValueOrDie()->n_bases()
                             : static_cast<int64_t>(region.end()) + expand_bp;
          nucleus::genomics::v1::Range ws_region;
          ws_region.set_reference_name(region.reference_name());
          ws_region.set_start(std::max<int64_t>(
              0, static_cast<int64_t>(region.start()) - expand_bp));
          ws_region.set_end(std::min<int64_t>(
              contig_n, static_cast<int64_t>(region.end()) + expand_bp));

          AlleleCounterOptions ws_ac_opts;
          ws_ac_opts.set_partition_size(
              opts.allele_counter_options().partition_size());
          ws_ac_opts.mutable_read_requirements()->set_min_mapping_quality(
              realigner_opts.ws_config().min_mapq());
          ws_ac_opts.mutable_read_requirements()->set_min_base_quality(
              realigner_opts.ws_config().min_base_quality());
          ws_ac_opts.mutable_read_requirements()->set_min_base_quality_mode(
              nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
          AlleleCounter pre(ref_reader.get(), ws_region, /*positions=*/{},
                             ws_ac_opts);
          for (const auto& r : raw_reads) pre.Add(r, ctx[s].name);
          reads_per_sample_v[s] = RealignReadsForRegion(
              raw_reads, ws_region, pre, *ref_reader, realigner_opts);
        } else {
          reads_per_sample_v[s] = std::move(raw_reads);
        }
      }

      // Build 3 AlleleCounters keyed by sample_name. Two-pass: probe
      // per sample (no candidate positions) → compute per-sample
      // candidate_positions via the multi-sample VariantCaller →
      // rebuild each AlleleCounter with its own candidate_positions.
      std::array<std::unique_ptr<AlleleCounter>, 3> counters;
      // PER-SAMPLE candidate positions (not the union, mirrors upstream
      // make_examples_core.py:2898 — `sample.variant_caller.get_candidate_
      // positions(allele_counters, sample_name)` runs per sample with a
      // single target_sample, so each sample's candidate_positions are
      // determined by THAT sample's evidence.
      //
      // Why this matters: with track_ref_reads=ON, ref reads are added to
      // AlleleCount.read_alleles ONLY at the sample's own candidate
      // positions. If a sample has no alt evidence at a position (e.g.
      // parent2 at an indel only seen in parent1+child), that position
      // is NOT a candidate for parent2 → parent2's read_alleles is empty
      // there → the candidate's ref_support_ext does not include parent2
      // reads at that position → small_model features for parent2 are 0.
      //
      // Our previous code used the UNION across all samples, which forced
      // every sample to track ref reads at every union-candidate position.
      // That inflated the small_model's combined-block total_depth (which
      // sums across all 3 samples in ref_support_ext) and produced wrong
      // SM probabilities at sites with asymmetric per-sample coverage.

      // Step 1: per-sample probe to compute per-sample candidate positions.
      std::array<std::unique_ptr<AlleleCounter>, 3> probes;
      for (int s = 0; s < n_samples; ++s) {
        if (!ctx[s].sam_reader) continue;
        probes[s] = std::make_unique<AlleleCounter>(
            ref_reader.get(), region, /*positions=*/std::vector<int>{},
            opts.allele_counter_options());
        for (const auto& r : reads_per_sample_v[s]) {
          probes[s]->Add(r, ctx[s].name);
        }
      }

      // Step 2: build per-sample candidate_positions via the multi-sample
      // VariantCaller's CallPositionsFromAlleleCounts (mirrors upstream's
      // sample.variant_caller.get_candidate_positions invocation per sample).
      std::array<std::vector<int>, 3> per_sample_cand_positions;
      {
        std::unordered_map<std::string, AlleleCounter*> probe_map;
        for (int s = 0; s < n_samples; ++s) {
          if (probes[s]) probe_map[ctx[s].name] = probes[s].get();
        }
        for (int s = 0; s < n_samples; ++s) {
          if (!probes[s]) continue;
          per_sample_cand_positions[s] =
              caller.CallPositionsFromAlleleCounts(
                  probe_map, ctx[s].name, ctx[s].role);
        }
      }

      // Step 3: rebuild each AlleleCounter with its OWN candidate_positions.
      for (int s = 0; s < n_samples; ++s) {
        if (!ctx[s].sam_reader) continue;
        counters[s] = std::make_unique<AlleleCounter>(
            ref_reader.get(), region, per_sample_cand_positions[s],
            opts.allele_counter_options());
        for (const auto& r : reads_per_sample_v[s]) {
          counters[s]->Add(r, ctx[s].name);
        }
      }
      probes = {};  // free probe memory

      // Build the unordered_map<sample_name, AlleleCounter*> map for
      // multi_sample::VariantCaller.
      std::unordered_map<std::string, AlleleCounter*> ac_map;
      for (int s = 0; s < n_samples; ++s) {
        if (counters[s]) ac_map[ctx[s].name] = counters[s].get();
      }

      // For each target sample (child, parent1, parent2): generate
      // candidates with the multi-sample API, run small_model dispatch,
      // emit examples + CVOs. Skip parents when --skip_parent_calling.
      for (int s = 0; s < n_samples; ++s) {
        SampleCtx& C = ctx[s];
        if (C.skip_output) continue;

        std::vector<DeepVariantCall> candidates =
            caller.CallsFromAlleleCounts(ac_map, C.name, C.role);
        // DirectPhasing (is_phased/PS, below) must see the FULL candidate set so
        // its SNP phasing graph matches upstream, which phases in
        // candidates_in_region *before* filter_candidates. Capture the
        // pre-prune set here (only when phasing is on) and phase that; the
        // select_variant_types prune + small_model dispatch then run on the
        // filtered set, and is_phased/PS is applied to the surviving
        // big_candidates by position. The multi-sample path still filters
        // before small_model so pruned types are not emitted via the
        // small_model CVO.
        std::vector<DeepVariantCall> phasing_candidates;
        if (absl::GetFlag(FLAGS_use_direct_phasing)) {
          phasing_candidates = candidates;  // full set, pre-prune
        }
        FilterCandidatesBySelectedTypes(&candidates, selected_types);
        if (candidates.empty()) continue;
        C.total_candidates += candidates.size();

        // Tumor-only allele_frequency channel: fill from PON VCF when present.
        // Mirrors Python's add_allele_frequencies_to_candidates called from
        // make_examples_core.py:2380 when 'allele_frequency' is in channels.
        if (pon_reader) {
          FillAlleleFrequencyFromPon(candidates, *pon_reader);
        }

        // VAF context — uses the target sample's AlleleCounts.
        if (C.small_model && counters[s]) {
          const auto& allele_counts = counters[s]->Counts();
          for (auto& c : candidates) PopulateVafContext(&c, allele_counts);
        }

        // Small-model dispatch (per alt-set), same as single-sample path.
        std::vector<DeepVariantCall> big_candidates;
        if (C.small_model) {
          for (auto& c : candidates) {
            const int n_alts = c.variant().alternate_bases_size();
            std::vector<std::vector<int>> alt_idx_sets;
            for (int i = 0; i < n_alts; ++i) alt_idx_sets.push_back({i});
            for (int i = 0; i < n_alts; ++i)
              for (int j = i + 1; j < n_alts; ++j)
                alt_idx_sets.push_back({i, j});

            // Trio small_model is multi-sample (106 features = 70
            // single-sample + 12 × 3 per-sample). Encode with the
            // target's `order` so per-sample feature blocks come in
            // the same insertion order upstream's Python uses.
            std::vector<std::string> sample_names_in_order;
            sample_names_in_order.reserve(3);
            for (int s2 = 0; s2 < 3; ++s2) {
              sample_names_in_order.push_back(ctx[s2].name);
            }

            bool any_failed = false;
            c.clear_make_examples_alt_allele_indices();
            for (const auto& idx_set : alt_idx_sets) {
              const auto features = EncodeSmallModelFeaturesMultiSample(
                  c, idx_set, sample_names_in_order, C.order);
              float probs[3] = {0, 0, 0};
              bool pred_ok =
                  C.small_model->Predict(features.data(), 1, probs);
              bool accept = false;
              if (pred_ok) {
                // Mirror upstream's _MAX_CONFIDENCE = 1 - 1e-7 clamp
                // (inference.py:46). When our BNNS-CPU saturates to exactly
                // 1.0 in FP32, ProbToPhred(1.0 - 1.0) = ProbToPhred(0) = 0
                // → GQ=0 → reject, while Docker's Eigen gives p < 1.0 →
                // GQ ≥ threshold → accept. Clamp to 1-1e-7 so saturated
                // p=1.0 maps to GQ=70 (same decision as Docker).
                const float max_p = std::min(
                    std::max({probs[0], probs[1], probs[2]}),
                    1.0f - 1e-7f);
                const int gq = ProbToPhred(1.0 - max_p);
                const int threshold = IsSnpForIndices(c.variant(), idx_set)
                                        ? snp_gq_threshold
                                        : indel_gq_threshold;
                accept = (gq >= threshold);
              }
              if (accept) {
                CallVariantsOutput cvo =
                    MakeSmallModelCvo(c, idx_set, probs);
                std::string serialized;
                cvo.SerializeToString(&serialized);
                if (C.small_cvo_writer) {
                  C.small_cvo_writer->WriteRecord(serialized);
                }
                ++C.total_small_hits;
              } else {
                auto* aai = c.add_make_examples_alt_allele_indices();
                for (int idx : idx_set) aai->add_indices(idx);
                any_failed = true;
              }
            }
            if (any_failed) {
              big_candidates.push_back(c);
              ++C.total_big_dispatched;
            }
          }
        } else {
          big_candidates = candidates;
          C.total_big_dispatched += candidates.size();
        }

        if (big_candidates.empty()) continue;

        // Phase 9 / Step 4b — DirectPhasing per-region (trio path).
        // Mirrors the single-sample wire-up at line ~1999, using only
        // the target sample's reads (reads_per_sample_v[s]). Each
        // target sample (child / parent1 / parent2) gets phased
        // independently against its own read pool — same semantic
        // as upstream's per-sample DirectPhasing invocation in
        // make_examples_core.py.
        if (absl::GetFlag(FLAGS_use_direct_phasing)) {
          std::vector<
              nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
              dp_read_ptrs;
          dp_read_ptrs.reserve(reads_per_sample_v[s].size());
          for (auto& r : reads_per_sample_v[s]) dp_read_ptrs.emplace_back(&r);
          ::learning::genomics::deepvariant::DirectPhasing dp(
              opts.direct_phasing_options());
          // Phase the full pre-prune candidate set (not just big_candidates)
          // so the phasing graph matches upstream; is_phased/PS is then applied
          // to the surviving big_candidates by position below.
          auto so = dp.PhaseReads(absl::MakeSpan(phasing_candidates),
                                   absl::MakeSpan(dp_read_ptrs));
          if (so.ok()) {
            const auto phased = dp.GetPhasedVariants();
            int64_t current_ps = -1;
            std::map<int64_t, int64_t> position_to_ps;
            for (const auto& pv : phased) {
              if (pv.is_first_in_block) current_ps = pv.position;
              if (current_ps >= 0 && pv.phase_1_bases != pv.phase_2_bases) {
                position_to_ps[pv.position] = current_ps;
              }
            }
            for (auto& c : big_candidates) {
              const int64_t pos = c.variant().start();
              auto it = position_to_ps.find(pos);
              if (it == position_to_ps.end()) continue;
              if (c.variant().calls_size() == 0) continue;
              auto* call = c.mutable_variant()->mutable_calls(0);
              call->set_is_phased(true);
              // Phase 9 / Step 4c — emit PS info field. PS = position of
              // first variant in block (1-based, VCF convention). Mirrors
              // upstream's stitch_phase_sets first-pass per-region output.
              const int ps_id = static_cast<int>(it->second + 1);
              nucleus::SetInfoField("PS", ps_id, call);
            }
          }
        }

        // ExamplesGenerator: 3 sample read vectors in upstream order
        // [parent1, child, parent2], rendered with this target's order
        // permutation. C.order tells the generator which slot of the
        // 3-sample array to put in slot 1 (target), 0, 2.
        std::vector<nucleus::ConstProtoPtr<DeepVariantCall>> cand_ptrs;
        cand_ptrs.reserve(big_candidates.size());
        for (auto& c : big_candidates) {
          cand_ptrs.push_back(
              nucleus::ConstProtoPtr<DeepVariantCall>(&c));
        }
        std::array<std::vector<nucleus::ConstProtoPtr<
            nucleus::genomics::v1::Read>>, 3> per_sample_ptrs;
        for (int q = 0; q < 3; ++q) {
          per_sample_ptrs[q].reserve(reads_per_sample_v[q].size());
          for (auto& r : reads_per_sample_v[q]) {
            per_sample_ptrs[q].push_back(
                nucleus::ConstProtoPtr<nucleus::genomics::v1::Read>(&r));
          }
        }
        std::vector<std::vector<nucleus::ConstProtoPtr<
            nucleus::genomics::v1::Read>>> reads_per_sample = {
                per_sample_ptrs[0], per_sample_ptrs[1], per_sample_ptrs[2]};
        std::vector<float> mean_coverage = {0.0f, 0.0f, 0.0f};
        std::vector<int> image_shape;

        auto stats = generator.WriteExamplesInRegion(
            absl::MakeSpan(cand_ptrs), absl::MakeSpan(reads_per_sample),
            absl::MakeSpan(C.order), C.role,
            absl::MakeSpan(mean_coverage), &image_shape);
        auto n_it = stats.find("n_examples");
        if (n_it != stats.end()) C.total_examples += n_it->second;
      }
    }  // end while next_region

    generator.SignalShardFinished();
    for (auto& c : ctx) {
      if (c.small_cvo_writer) c.small_cvo_writer->Close();
    }

    // Aggregate per-sample stats into the worker totals (sum across
    // the 3 samples — the postprocess stage will re-bucket them later).
    int64_t tot_cand = 0, tot_ex = 0, tot_small = 0, tot_big = 0;
    for (auto& c : ctx) {
      tot_cand  += c.total_candidates;
      tot_ex    += c.total_examples;
      tot_small += c.total_small_hits;
      tot_big   += c.total_big_dispatched;
    }
    out_stats->total_candidates    = tot_cand;
    out_stats->total_examples      = tot_ex;
    out_stats->total_small_hits    = tot_small;
    out_stats->total_big_dispatched = tot_big;
  };

  // Worker function: opens its own SamReader/IndexedFastaReader/
  // ExamplesGenerator/SmallModel, then loops fetching regions from
  // `next_region` until the queue is exhausted. Writes only to its own
  // per-thread files; no inter-thread mutation.
  auto run_worker = [&](int tid, WorkerStats* out_stats) {
    if (IsTrioMode() || IsSomaticMode() || IsPangenomeMode()) {
      // Multi-sample worker handles trio (3 samples), somatic (1-2),
      // and pangenome-aware (2). Dispatched on opts.sample_options_size().
      run_trio_worker(tid, out_stats);
      return;
    }
    auto t_ref_or = nucleus::IndexedFastaReader::FromFile(
        ref_path, absl::StrCat(ref_path, ".fai"));
    CHECK(t_ref_or.ok()) << "thread " << tid << ": ref reopen failed";
    auto ref_reader = std::move(t_ref_or.ValueOrDie());

    auto t_sam_or = nucleus::SamReader::FromFile(reads_path, sam_opts);
    CHECK(t_sam_or.ok()) << "thread " << tid << ": BAM reopen failed";
    auto sam_reader = std::move(t_sam_or.ValueOrDie());

    // Use the multi-sample VariantCaller (as upstream does for every sample
    // count) rather than the vcf_candidate_importer caller: only the former
    // emits methylated reference sites (alt='.'), which methylation-aware
    // phasing consumes. For a single sample the candidate set is otherwise
    // identical (validated against the WGS baseline).
    multi_sample::VariantCaller caller(
        opts.sample_options(0).variant_caller_options());

    const std::unordered_map<std::string, std::string> example_filenames = {
        {"sample", thread_examples_path(tid)}};
    ExamplesGenerator generator(opts, example_filenames);

    std::unique_ptr<SmallModel> small_model;
    std::unique_ptr<TFRecordWriter> small_cvo_writer;
    if (!small_path.empty()) {
      small_model = SmallModel::Load(small_path);
      CHECK(small_model) << "thread " << tid << ": small_model load failed";
      small_cvo_writer = TFRecordWriter::New(thread_small_cvo_path(tid));
      CHECK(small_cvo_writer)
          << "thread " << tid << ": small CVO writer open failed";
    }
    // Phase 9 / Step 3 — gVCF non-variant TFRecord writer (one per worker
    // thread, sharded). Postprocess reads via ShardedVariantReader.
    std::unique_ptr<TFRecordWriter> gvcf_writer;
    if (!gvcf_path_top.empty()) {
      gvcf_writer = TFRecordWriter::New(thread_gvcf_path_top(tid));
      CHECK(gvcf_writer)
          << "thread " << tid << ": gvcf writer open failed";
    }

    int64_t total_candidates = 0;
    int64_t total_examples = 0;
    int64_t total_small_hits = 0;
    int64_t total_big_dispatched = 0;

    // Hoist loop-invariant flag reads (and the RealignerOptions proto, which
    // is built purely from immutable flags) out of the per-region loop. absl
    // flags are immutable after parse, so these values are identical on every
    // iteration; computing them once is a pure-perf win with no value change.
    const bool kRealignerEnabled = absl::GetFlag(FLAGS_realigner_enabled);
    const bool kSplitSkipReads = absl::GetFlag(FLAGS_split_skip_reads);
    const bool kUseDirectPhasing = absl::GetFlag(FLAGS_use_direct_phasing);
    const bool kSmallModelUseHaplotypes =
        absl::GetFlag(FLAGS_small_model_use_haplotypes);
    const RealignerOptions kRealignerOpts = RealignerOptionsFromFlags();

    while (true) {
      const size_t i = next_region.fetch_add(1, std::memory_order_relaxed);
      if (i >= shard_regions.size()) break;
      const auto& region = shard_regions[i];
    // Periodic progress: log every 1000 regions (cheap modulo only). `i` is
    // the shared 0-based region index, so this fires once per 1000 regions
    // across all worker threads combined; counters are thread-local running
    // totals for this worker.
    if (i % 1000 == 0) {
      LOG(INFO) << "Progress: processed " << i << "/" << shard_regions.size()
                << " regions (thread " << tid << " examples so far="
                << total_examples << ")";
    }
    LOG(INFO) << "Region: " << region.reference_name() << ":"
              << region.start() << "-" << region.end();
    const std::string region_str =
        absl::StrCat(region.reference_name(), ":", region.start(),
                     "-", region.end());
    DV_SIGNPOST_INTERVAL_BEGIN(RegionTotal, region_str.c_str());

    // Query reads.
    DV_SIGNPOST_INTERVAL_BEGIN(BamQuery, region_str.c_str());
    auto reads_or = sam_reader->Query(region);
    if (!reads_or.ok()) {
      LOG(WARNING) << "Query failed for " << region.reference_name() << ":"
                   << region.start() << "-" << region.end()
                   << " — " << reads_or.status();
      continue;
    }
    auto& reads_iter = reads_or.ValueOrDie();

    // Collect reads into a vector (AlleleCounter needs random access).
    std::vector<nucleus::genomics::v1::Read> reads;
    nucleus::genomics::v1::Read tmp_read;
    while (true) {
      auto next = reads_iter->Next(&tmp_read);
      if (!next.ok() || !next.ValueOrDie()) break;
      reads.push_back(tmp_read);
    }
    reads_iter->Release().IgnoreError();
    DV_SIGNPOST_INTERVAL_END(BamQuery);

    // Match upstream make_examples_core.py:partition_reads_etc, which
    // applies Algorithm-R reservoir sampling to cap reads per partition
    // at `max_reads_per_partition` (default 1500). Without this cap,
    // high-coverage regions (chr20:31185000-31186000 has 5686 reads
    // post-filter) blow up the pileup-image evidence and produce a
    // different DP/AD/VAF than Docker → different small_model dispatch
    // → different deepvariant softmax → FILTER drift. The RNG is a
    // NumPy-compatible mt19937 (numpy_mt19937.h) seeded with
    // opts.random_seed (609314161, the upstream default), reset per
    // region — matches `np.random.RandomState(seed)` in
    // make_examples_core.py:2134.
    //
    // SHARD-COUNT INDEPENDENCE GUARD (2026-05-01): the upstream Python
    // pipeline runs as a single process per shard, so partition-level
    // determinism is automatic. Our native port runs N worker threads
    // in one process, all sharing the BAM via per-thread SamReader
    // instances. Empirically (chr20 trio HG002 today, num_shards=4 vs
    // num_shards=14 at Phase 5.5d/10) we observe a 0.2-0.3 % PASS-set
    // delta between the two configurations, traceable to reservoir-
    // sampling output differing across thread loads. To eliminate this
    // we stable-sort the read vector by (POS, fragment_name,
    // read_number) BEFORE feeding it to the reservoir. BAM is
    // coordinate-sorted, so reads naturally arrive in increasing POS
    // order from htslib; the secondary key (fragment_name + read_number)
    // disambiguates within-position reads deterministically. If htslib
    // already returns reads in this exact order (the BAM standard
    // guarantee), this sort is a no-op (stable sort preserves relative
    // order on equal keys); if any thread-related state introduces
    // sub-position reordering, the sort imposes the canonical order.
    // Docker's pysam.AlignmentFile.fetch returns reads in BAM-sorted
    // order, so ours-after-sort matches Docker's order.
    const int max_rpp = static_cast<int>(opts.max_reads_per_partition());
    if (max_rpp > 0 && reads.size() > static_cast<size_t>(max_rpp)) {
      const size_t orig_n = reads.size();
      // BUG FIX (2026-05-10): the previous stable_sort by
      // (POS, fragment_name, read_number) here was added in Phase 5.5d/10
      // as a "shard-count-independence guard", but it CHANGED the
      // input order to reservoir sampling vs Docker. Docker reads BAM
      // in natural order (POS only, secondary by file offset, NOT by
      // fragment_name/read_number), and our sort reordered same-POS
      // reads → reservoir picks different reads → ±1-4 read DP
      // differences at WG scale on ~79 % of FILTER-mismatch sites
      // (HG002 WG: 4,146 FM, of which ~3,200 trace to this sort).
      //
      // BAM is coordinate-sorted at the file level, and htslib's
      // SamReader::Query() iterates within a region in the BAM's
      // natural order — same as pysam.AlignmentFile.fetch which
      // Docker uses. So removing the sort makes our reservoir input
      // bit-identical to Docker's at the read level.
      //
      // The "shard-count independence" rationale doesn't apply here
      // anyway: each region is processed by a single thread that owns
      // its own SamReader, so the read-load order is deterministic
      // per region regardless of thread count.
      ::deepvariant::npr::NumpyMt19937 region_rng(opts.random_seed());
      auto sampled =
          ::deepvariant::npr::ReservoirSamplePtrs(reads, max_rpp, region_rng);
      std::vector<nucleus::genomics::v1::Read> kept;
      kept.reserve(sampled.size());
      for (const auto* p : sampled) kept.push_back(*p);
      reads = std::move(kept);
      LOG(INFO) << "  reservoir-sampled " << orig_n << " → " << reads.size()
                 << " reads (max_reads_per_partition=" << max_rpp << ")";
    }

    LOG(INFO) << "  read " << reads.size() << " reads from BAM";
    // Upstream make_examples_core.py only aborts a zero-coverage region early
    // when gVCF is disabled; with a gvcf_writer we must still fall through to
    // emit per-position reference-confidence rows (the gVCF block below tolerates
    // empty reads, and probe_candidates.empty() then `continue`s with no work).
    if (reads.empty() && gvcf_writer == nullptr) continue;

    // RNA-seq: split reads on N (SKIP) CIGAR ops into per-exon sub-reads
    // before candidate discovery / realignment / pileup. Mirrors upstream
    // realigner.py:realign_reads → split_reads (gated by --split_skip_reads,
    // the RNASEQ example_info default). Must run before the AlleleCounter so
    // intron gaps don't pollute the pileup image.
    if (kSplitSkipReads) {
      const size_t before = reads.size();
      reads = SplitReadsOnSkip(reads);
      LOG(INFO) << "  split_skip_reads: " << before << " → " << reads.size()
                << " reads (split on N CIGAR)";
    }

    // ── Optional: realign reads through assembled haplotypes ─────────────
    // Done before any AlleleCounter pass so candidate sweep + ref read
    // tracking see the realigned reads (matches upstream's flow).
    DV_SIGNPOST_INTERVAL_BEGIN(Realigner, region_str.c_str());
    std::vector<nucleus::genomics::v1::Read> working_reads;
    if (kRealignerEnabled) {
      // Pre-scan AlleleCounter for the WindowSelector. Upstream
      // (realigner/window_selector.py:_candidates_from_reads) builds
      // a *dedicated* AlleleCounter with WindowSelector-specific
      // requirements (ws_min_mapq=20, ws_min_base_quality=20) over an
      // expanded region (region_expansion_in_bp=20). The candidate-
      // emission AlleleCounter further down uses the looser
      // make_examples thresholds (10/10) — they're separate counters.
      const auto& realigner_opts = kRealignerOpts;
      const int expand_bp = realigner_opts.ws_config().region_expansion_in_bp();
      auto contig_or = ref_reader->Contig(region.reference_name());
      const int64_t contig_n =
          contig_or.ok() ? contig_or.ValueOrDie()->n_bases() :
          static_cast<int64_t>(region.end()) + expand_bp;
      nucleus::genomics::v1::Range ws_region;
      ws_region.set_reference_name(region.reference_name());
      ws_region.set_start(std::max<int64_t>(0,
          static_cast<int64_t>(region.start()) - expand_bp));
      ws_region.set_end(std::min<int64_t>(contig_n,
          static_cast<int64_t>(region.end()) + expand_bp));

      AlleleCounterOptions ws_ac_opts;
      ws_ac_opts.set_partition_size(opts.allele_counter_options().partition_size());
      ws_ac_opts.mutable_read_requirements()->set_min_mapping_quality(
          realigner_opts.ws_config().min_mapq());
      ws_ac_opts.mutable_read_requirements()->set_min_base_quality(
          realigner_opts.ws_config().min_base_quality());
      ws_ac_opts.mutable_read_requirements()->set_min_base_quality_mode(
          nucleus::genomics::v1::ReadRequirements::ENFORCED_BY_CLIENT);
      // track_ref_reads stays false — the WindowSelector doesn't use ref
      // reads (AlleleFilter() rejects REFERENCE alleles).
      AlleleCounter pre(ref_reader.get(), ws_region, /*candidates=*/{},
                        ws_ac_opts);
      for (const auto& r : reads) pre.Add(r, sample_name);
      working_reads =
          RealignReadsForRegion(reads, ws_region, pre, *ref_reader,
                                 realigner_opts);
      // Optional: dump (qname, contig, pos, mapq, cigar, seq) of
      // post-realigner reads per chunk so we can side-by-side diff
      // against upstream's --emit_realigned_reads BAM.
      static const char* dump_dir = std::getenv("DV_REALIGNED_READS_TSV");
      if (dump_dir) {
        std::string fname = std::string(dump_dir) + "/" +
            region.reference_name() + ":" +
            std::to_string(region.start()) + "-" +
            std::to_string(region.end()) + ".tsv";
        std::ofstream rf(fname);
        if (rf) {
          for (const auto& r : working_reads) {
            rf << r.fragment_name() << '/' << r.read_number() << '\t'
               << r.alignment().position().reference_name() << '\t'
               << r.alignment().position().position() << '\t'
               << r.alignment().mapping_quality() << '\t';
            for (const auto& cu : r.alignment().cigar()) {
              rf << cu.operation_length();
              switch (cu.operation()) {
                using ::nucleus::genomics::v1::CigarUnit;
                case CigarUnit::ALIGNMENT_MATCH:    rf << 'M'; break;
                case CigarUnit::INSERT:             rf << 'I'; break;
                case CigarUnit::DELETE:             rf << 'D'; break;
                case CigarUnit::SKIP:               rf << 'N'; break;
                case CigarUnit::CLIP_SOFT:          rf << 'S'; break;
                case CigarUnit::CLIP_HARD:          rf << 'H'; break;
                case CigarUnit::PAD:                rf << 'P'; break;
                case CigarUnit::SEQUENCE_MATCH:     rf << '='; break;
                case CigarUnit::SEQUENCE_MISMATCH:  rf << 'X'; break;
                default: rf << '?';
              }
            }
            rf << '\t' << r.aligned_sequence() << '\n';
          }
        }
      }
    } else {
      // No realigner: hand the read vector off by move instead of deep-
      // copying the entire vector of Read protos (the default WGS/WES path).
      // `reads` is not read again after this point — the only later
      // reference (the candidates/reads size log below) is repointed to
      // `working_reads`, whose size is identical (the realigner path also
      // preserves read count), so this is value- and output-neutral.
      working_reads = std::move(reads);
    }
    DV_SIGNPOST_INTERVAL_END(Realigner);

    // First pass: find candidate positions (no ref-read tracking yet).
    DV_SIGNPOST_INTERVAL_BEGIN(AlleleCounterProbe, region_str.c_str());
    AlleleCounter probe(ref_reader.get(), region, {},
                        opts.allele_counter_options());
    for (const auto& r : working_reads) probe.Add(r, sample_name);
    DV_SIGNPOST_INTERVAL_END(AlleleCounterProbe);

    // Phase 9 / Step 3 — gVCF non-variant TFRecord emission. Per-position
    // reference-confidence rows are written for every region (regardless
    // of whether candidates exist). Postprocess merges these with the
    // variant TFRecord via nucleus::MergeAndWriteVariantsAndNonVariants.
    if (gvcf_writer) {
      auto summaries = probe.SummaryCounts(0, 0);
      auto gvcf_rows = MakeGvcfRows(
          summaries, sample_name,
          absl::GetFlag(FLAGS_p_error),
          absl::GetFlag(FLAGS_gvcf_gq_binsize),
          /*max_gq=*/50,
          absl::GetFlag(FLAGS_include_med_dp),
          &haploid_contigs, &par_regions);
      for (const auto& v : gvcf_rows) {
        std::string serialized;
        v.SerializeToString(&serialized);
        gvcf_writer->WriteRecord(serialized);
      }
    }

    auto probe_candidates = caller.CallsFromAlleleCounts(
        {{sample_name, &probe}}, sample_name, "sample");
    if (probe_candidates.empty()) continue;

    // Second pass: rerun AlleleCounter with the candidate positions known
    // up-front. AlleleCounter only retains REF-supporting reads in its
    // read_alleles map at positions that appear in this list (when
    // track_ref_reads=true). Without this two-pass shape the small_model
    // sees num_reads_supports_ref = 0 on every candidate.
    std::vector<int> candidate_positions;
    candidate_positions.reserve(probe_candidates.size());
    for (const auto& c : probe_candidates) {
      candidate_positions.push_back(static_cast<int>(c.variant().start()));
    }
    std::sort(candidate_positions.begin(), candidate_positions.end());
    candidate_positions.erase(
        std::unique(candidate_positions.begin(), candidate_positions.end()),
        candidate_positions.end());

    DV_SIGNPOST_INTERVAL_BEGIN(AlleleCounterMain, region_str.c_str());
    AlleleCounter counter(ref_reader.get(), region, candidate_positions,
                          opts.allele_counter_options());
    for (const auto& r : working_reads) counter.Add(r, sample_name);
    DV_SIGNPOST_INTERVAL_END(AlleleCounterMain);

    std::vector<DeepVariantCall> candidates =
        caller.CallsFromAlleleCounts({{sample_name, &counter}}, sample_name,
                                     "sample");
    if (candidates.empty()) continue;

    // Phase 5.5d/14 — DirectPhasing runs BEFORE small_model dispatch so the
    // 106-feature haplotype-expanded small_model (PacBio/ONT) sees the same
    // per-read phase that upstream's FeatureEncoder does. Upstream order
    // (make_examples_core.py):
    //   1. direct_phasing.phase_reads(candidates, reads) → read_phases dict
    //   2. small_model invoked with FeatureEncoder(haplotype, read_phases)
    //   3. variant phasing via dp.GetPhasedVariants() → is_phased + PS
    // Pre-fix: small_model used BAM HP tags (whatshap haplotag from BAM
    // PG line) — these can disagree with DirectPhasing's per-region output
    // at phase-block boundaries. Sites where BAM HP=0 (unphased) but
    // DirectPhasing assigns HP=1/2 produce different 106-feature vectors,
    // flipping small_model GQ across the dispatch threshold.
    // After-fix: DP output keyed by `fragment_name + "/" + read_number`
    // (matches allelecounter.cc::ReadKey) overrides BAM HP tags. Only run
    // when --use_direct_phasing OR --small_model_use_haplotypes is set;
    // otherwise we'd waste cycles on WGS/WES paths where it has no effect.
    ::learning::genomics::deepvariant::DirectPhasing dp(
        opts.direct_phasing_options());
    bool dp_ran = false;
    std::unordered_map<std::string, int8_t> read_hp_tags;
    {
      const bool need_phasing =
          kUseDirectPhasing || kSmallModelUseHaplotypes;
      if (need_phasing) {
        // BUG FIX (chr20:23.97-23.99M PacBio FN cluster, 0e15ddb2 diagnosis):
        // upstream make_examples_core.py:2308-2317 expands the region by 20%
        // (`PHASE_READS_REGION_PADDING_PCT = 20`) before fetching reads for
        // DirectPhasing. Reads spanning region boundaries provide the
        // SNP-graph context that lets DP correctly split reads across HP=1
        // vs HP=2 in dense haplotype blocks. Without padding, DP collapses
        // all 49 alt-supporting reads at chr20:23973486 onto a single
        // haplotype, so the small_model sees "100% reads on one HP, other
        // HP empty" → predicts homref (Docker correctly splits → predicts
        // HET).
        //
        // Fix: re-fetch reads from a 20%-padded region for DP. We do NOT
        // change `working_reads` itself (still used downstream for the
        // candidate-emitting AlleleCounter, pileup encoder, etc.) — only
        // the DP input set. Upstream uses raw BAM reads (not realigned)
        // for DP; we mirror that by re-Querying the SAM reader.
        const int64_t region_len =
            static_cast<int64_t>(region.end()) -
            static_cast<int64_t>(region.start());
        const int64_t pad = std::max<int64_t>(1, region_len * 20 / 100);
        // contig_n in this scope: re-derive locally (the outer realigner
        // block's contig_n is out of scope here when realigner is disabled).
        auto contig_or_dp = ref_reader->Contig(region.reference_name());
        const int64_t contig_n_dp =
            contig_or_dp.ok()
                ? contig_or_dp.ValueOrDie()->n_bases()
                : static_cast<int64_t>(region.end()) + pad;
        nucleus::genomics::v1::Range padded_region;
        padded_region.set_reference_name(region.reference_name());
        padded_region.set_start(std::max<int64_t>(0,
            static_cast<int64_t>(region.start()) - pad));
        padded_region.set_end(std::min<int64_t>(contig_n_dp,
            static_cast<int64_t>(region.end()) + pad));

        std::vector<nucleus::genomics::v1::Read> phasing_reads;
        auto phasing_or = sam_reader->Query(padded_region);
        if (phasing_or.ok()) {
          auto& phasing_iter = phasing_or.ValueOrDie();
          nucleus::genomics::v1::Read tmp;
          while (true) {
            auto next = phasing_iter->Next(&tmp);
            if (!next.ok() || !next.ValueOrDie()) break;
            phasing_reads.push_back(std::move(tmp));
          }
          phasing_iter->Release().IgnoreError();
        }
        // Defensive fallback: if the padded query returned nothing (e.g.,
        // sam_reader transient failure), fall through to working_reads.
        const auto& dp_reads_src =
            phasing_reads.empty() ? working_reads : phasing_reads;

        std::vector<
            nucleus::ConstProtoPtr<const nucleus::genomics::v1::Read>>
            dp_read_ptrs;
        dp_read_ptrs.reserve(dp_reads_src.size());
        for (const auto& r : dp_reads_src) dp_read_ptrs.emplace_back(&r);
        // Methylation-aware phasing (PacBio/ONT): split methylated reference
        // sites (alt='.') out of the SNP phasing graph, phase the SNP
        // candidates, then use the methylated sites to phase reads DirectPhasing
        // left unphased. Mirror of make_examples_core.py's
        // enable_methylation_aware_phasing block. When the flag is off no
        // methylated-ref-site candidates exist, so this is identical to phasing
        // the full candidate set.
        const bool meth_aware =
            absl::GetFlag(FLAGS_enable_methylation_aware_phasing);
        std::vector<DeepVariantCall> methylated_ref_sites;
        std::vector<DeepVariantCall> snp_candidates;
        const std::vector<DeepVariantCall>* dp_candidates = &candidates;
        if (meth_aware) {
          for (const auto& c : candidates) {
            if (IsMethylatedRefSite(c)) methylated_ref_sites.push_back(c);
            else snp_candidates.push_back(c);
          }
          dp_candidates = &snp_candidates;
        }
        auto so = dp.PhaseReads(absl::MakeSpan(*dp_candidates),
                                  absl::MakeSpan(dp_read_ptrs));
        if (so.ok()) {
          dp_ran = true;
          std::vector<int> phases = so.ValueOrDie();
          // Refine per-read phases using 5mC levels at the methylated reference
          // sites (Wilcoxon rank-sum vote over reads DirectPhasing left at HP_0).
          if (meth_aware && !methylated_ref_sites.empty()) {
            auto meth = PerformMethylationAwarePhasing(
                absl::MakeSpan(dp_reads_src), phases, methylated_ref_sites,
                /*max_iter=*/10);
            std::vector<int>& meth_phases = std::get<0>(meth);
            // std::get<1>(meth) holds the per-site p-values, used upstream only
            // for --phasing_error_stats_output, which the native port does not
            // emit; intentionally dropped here.
            if (meth_phases.size() == phases.size()) {
              int rephased = 0;
              for (size_t i = 0; i < phases.size(); ++i) {
                if (meth_phases[i] != phases[i]) ++rephased;
              }
              LOG(INFO) << "Methylation-aware phasing in " << region_str << ": "
                        << methylated_ref_sites.size() << " methylated ref sites, "
                        << rephased << " reads rephased";
              phases = std::move(meth_phases);
            }
          }
          // phases[i] corresponds to dp_reads_src[i]: 0/1/2.
          for (size_t i = 0;
               i < dp_reads_src.size() && i < phases.size(); ++i) {
            if (phases[i] == 0) continue;  // HP_0 default; skip to save memory
            const auto& r = dp_reads_src[i];
            read_hp_tags[r.fragment_name() + "/" +
                         std::to_string(r.read_number())] =
                static_cast<int8_t>(phases[i]);
          }
        }
      }
      // Fallback: if DP didn't run or failed and we still need haplotype
      // features, use BAM HP tags (whatshap haplotag in PacBio BAMs).
      // Guards against regression for users running --small_model_use_haplotypes
      // without --use_direct_phasing on a pre-haplotagged BAM.
      if (!dp_ran && kSmallModelUseHaplotypes) {
        for (const auto& r : working_reads) {
          auto hp_it = r.info().find("HP");
          if (hp_it == r.info().end() ||
              hp_it->second.values().empty()) continue;
          const auto& hp_val = hp_it->second.values(0);
          if (!hp_val.has_number_value()) continue;
          const int8_t hp = static_cast<int8_t>(hp_val.number_value());
          if (hp == 0) continue;
          read_hp_tags[r.fragment_name() + "/" +
                       std::to_string(r.read_number())] = hp;
        }
      }
    }

    // select_variant_types pruning runs AFTER DirectPhasing so the phasing
    // graph — and the per-read HP tags it produces, which feed the haplotype
    // small_model's 106-feature vector — is built from the full candidate set.
    // This matches upstream's order: filter_candidates runs in process() only
    // after candidates_in_region has already phased. Pruning earlier would drop
    // the SNP backbone DirectPhasing relies on and silently change phasing.
    FilterCandidatesBySelectedTypes(&candidates, selected_types);
    if (candidates.empty()) continue;
    total_candidates += candidates.size();

    // Small-model first-pass dispatch. Mirror of upstream
    // `SmallModelVariantCaller.call_variants` + `make_small_model_examples.
    // get_set_of_allele_indices`:
    //   - For each candidate, enumerate the FULL set of alt-allele-indices:
    //       biallelic   = [(0,), (1,), …, (N-1,)]
    //       multiallelic = list(combinations(range(N), 2))
    //   - Run small_model on each (candidate, alt_indices) PAIR.
    //   - PER-PAIR pass/fail: if pass → emit small_model CVO; if fail →
    //     append to candidate.make_examples_alt_allele_indices so big_model
    //     generates an example for that specific alt-set only. Multiple
    //     pairs from the same candidate can split between small/big.
    DV_SIGNPOST_INTERVAL_BEGIN(SmallModel, region_str.c_str());
    std::vector<DeepVariantCall> big_candidates;
    if (small_model) {
      // Populate VAF context for every candidate (the small model's 51
      // VAF-context features need it; the big model doesn't, but it's cheap
      // and keeps both paths producing the same DeepVariantCall shape).
      const auto& allele_counts = counter.Counts();
      for (auto& c : candidates) {
        PopulateVafContext(&c, allele_counts);
      }

      // read_hp_tags is now built above (Phase 5.5d/14): DirectPhasing
      // output overrides BAM HP tags when DP runs successfully.
      const bool use_haplotypes = kSmallModelUseHaplotypes;

      for (auto& c : candidates) {
        const int n_alts = c.variant().alternate_bases_size();
        // Build the list of alt-index sets to query: single + combinations.
        std::vector<std::vector<int>> alt_idx_sets;
        for (int i = 0; i < n_alts; ++i) alt_idx_sets.push_back({i});
        for (int i = 0; i < n_alts; ++i) {
          for (int j = i + 1; j < n_alts; ++j) {
            alt_idx_sets.push_back({i, j});
          }
        }

        // Per alt-index-set: predict + decide.
        bool any_failed = false;
        c.clear_make_examples_alt_allele_indices();
        for (const auto& idx_set : alt_idx_sets) {
          const auto features = use_haplotypes
              ? EncodeSmallModelFeaturesHaplotype(c, idx_set, read_hp_tags)
              : EncodeSmallModelFeatures(c, idx_set);
          float probs[3] = {0, 0, 0};
          bool pred_ok = small_model->Predict(features.data(), 1, probs);
          bool accept = false;
          if (pred_ok) {
            // Mirror upstream _MAX_CONFIDENCE clamp (see trio path comment).
            const float max_p = std::min(
                std::max({probs[0], probs[1], probs[2]}),
                1.0f - 1e-7f);
            const int gq = ProbToPhred(1.0 - max_p);
            const int threshold = IsSnpForIndices(c.variant(), idx_set)
                                    ? snp_gq_threshold
                                    : indel_gq_threshold;
            accept = (gq >= threshold);
          }
          if (accept) {
            // Single-alt CVOs use idx_set[0]; the multi-alt (i, j) set is
            // emitted with both indices so postprocess merge can route it
            // correctly.
            CallVariantsOutput cvo = MakeSmallModelCvo(c, idx_set, probs);
            std::string serialized;
            cvo.SerializeToString(&serialized);
            small_cvo_writer->WriteRecord(serialized);
            ++total_small_hits;
          } else {
            // Failed → big model generates an example for this exact
            // alt-index-set only.
            auto* aai = c.add_make_examples_alt_allele_indices();
            for (int idx : idx_set) aai->add_indices(idx);
            any_failed = true;
          }
        }
        if (any_failed) {
          big_candidates.push_back(c);
          ++total_big_dispatched;
        }
      }
    } else {
      big_candidates = candidates;
      total_big_dispatched += candidates.size();
    }

    DV_SIGNPOST_INTERVAL_END(SmallModel);
    if (big_candidates.empty()) continue;

    // Phase 9 / Step 4b + 5.5d/14 — variant phasing now reuses the `dp`
    // object built before small_model dispatch (no second PhaseReads
    // call). Walks GetPhasedVariants() to mark each big_candidate's
    // VariantCall.is_phased = true and emit PS info field. The phase
    // set ID is per-region (= start of block); cross-region stitching
    // is a follow-up that mirrors upstream's stitch_phase_sets.
    if (dp_ran && kUseDirectPhasing) {
      const auto phased = dp.GetPhasedVariants();
      // Walk in order; track current phase set (= start of block).
      int64_t current_ps = -1;
      std::map<int64_t, int64_t> position_to_ps;
      for (const auto& pv : phased) {
        if (pv.is_first_in_block) current_ps = pv.position;
        if (current_ps >= 0 && pv.phase_1_bases != pv.phase_2_bases) {
          position_to_ps[pv.position] = current_ps;
        }
      }
      for (auto& c : big_candidates) {
        const int64_t pos = c.variant().start();
        auto it = position_to_ps.find(pos);
        if (it == position_to_ps.end()) continue;
        if (c.variant().calls_size() == 0) continue;
        auto* call = c.mutable_variant()->mutable_calls(0);
        call->set_is_phased(true);
        // Phase 9 / Step 4c — emit PS info field. PS = position of
        // first variant in block (1-based, VCF convention). Mirrors
        // upstream's stitch_phase_sets first-pass per-region output.
        const int ps_id = static_cast<int>(it->second + 1);
        nucleus::SetInfoField("PS", ps_id, call);
      }
    }

    // Wrap in ConstProtoPtr for ExamplesGenerator API.
    std::vector<nucleus::ConstProtoPtr<DeepVariantCall>> cand_ptrs;
    cand_ptrs.reserve(big_candidates.size());
    for (auto& c : big_candidates) {
      cand_ptrs.push_back(nucleus::ConstProtoPtr<DeepVariantCall>(&c));
    }

    std::vector<nucleus::ConstProtoPtr<nucleus::genomics::v1::Read>> read_ptrs;
    read_ptrs.reserve(working_reads.size());
    for (auto& r : working_reads) {
      read_ptrs.push_back(
          nucleus::ConstProtoPtr<nucleus::genomics::v1::Read>(&r));
    }
    std::vector<std::vector<
        nucleus::ConstProtoPtr<nucleus::genomics::v1::Read>>>
        reads_per_sample = {read_ptrs};

    std::vector<int> sample_order = {0};
    std::vector<float> mean_coverage = {0.0f};
    std::vector<int> image_shape;

    LOG(INFO) << "  candidates=" << candidates.size()
              << " reads=" << working_reads.size();

    DV_SIGNPOST_INTERVAL_BEGIN(PileupEncode, region_str.c_str());
    auto stats = generator.WriteExamplesInRegion(
        absl::MakeSpan(cand_ptrs), absl::MakeSpan(reads_per_sample),
        absl::MakeSpan(sample_order), "sample",
        absl::MakeSpan(mean_coverage), &image_shape);
    DV_SIGNPOST_INTERVAL_END(PileupEncode);

    auto n_it = stats.find("n_examples");
    if (n_it != stats.end()) total_examples += n_it->second;
    DV_SIGNPOST_INTERVAL_END(RegionTotal);
    }  // end while next_region

    generator.SignalShardFinished();
    if (small_cvo_writer) small_cvo_writer->Close();
    if (gvcf_writer) gvcf_writer->Close();

    out_stats->total_candidates = total_candidates;
    out_stats->total_examples = total_examples;
    out_stats->total_small_hits = total_small_hits;
    out_stats->total_big_dispatched = total_big_dispatched;
  };  // end run_worker lambda

  // ── Dispatch workers ─────────────────────────────────────────────────────
  std::vector<WorkerStats> stats_per_thread(n_threads);
  if (n_threads == 1) {
    run_worker(0, &stats_per_thread[0]);
  } else {
    std::vector<std::thread> workers;
    workers.reserve(n_threads);
    for (int t = 0; t < n_threads; ++t) {
      workers.emplace_back([&, t] { run_worker(t, &stats_per_thread[t]); });
    }
    for (auto& w : workers) w.join();
  }

  // ── Sum per-thread stats ─────────────────────────────────────────────────
  WorkerStats agg;
  for (const auto& s : stats_per_thread) {
    agg.total_candidates    += s.total_candidates;
    agg.total_examples      += s.total_examples;
    agg.total_small_hits    += s.total_small_hits;
    agg.total_big_dispatched += s.total_big_dispatched;
  }

  // No end-of-stage concat: workers wrote sharded `name-NNNNN-of-NNNNN`
  // files that downstream stages read directly via TFRecordReader's
  // `@N` shard spec expansion.

  LOG(INFO) << "make_examples done: " << agg.total_candidates << " candidates, "
            << agg.total_examples << " examples written"
            << " (small_model_hits=" << agg.total_small_hits
            << ", big_model_dispatched=" << agg.total_big_dispatched
            << ", threads=" << n_threads << ").";
  return 0;
}

}  // namespace deepvariant
