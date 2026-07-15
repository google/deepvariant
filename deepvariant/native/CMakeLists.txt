# deepvariant/native — Phase 2: call_variants (TFRecord + Core ML inference).

# ---------------------------------------------------------------------------
# dv_tfrecord — TFRecord reader/writer (C++ 17, no TF runtime)
# ---------------------------------------------------------------------------
add_library(dv_tfrecord STATIC tfrecord.cc)
target_include_directories(dv_tfrecord PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_tfrecord PUBLIC
  absl::crc32c
  proto_dv
  proto_nucleus
)

# ---------------------------------------------------------------------------
# dv_weights — .dvw mmap loader (Phase 5.5; consumed by Metal/BNNS path)
# ---------------------------------------------------------------------------
add_library(dv_weights STATIC dv_weights.cc)
target_include_directories(dv_weights PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_weights PUBLIC
  absl::log
)

# ---------------------------------------------------------------------------
# dv_bnns_finalize — Phase 5.5 deterministic CPU dense + softmax
# (sequential FP32 reduction, designed to bit-match TF CPU output)
# ---------------------------------------------------------------------------
add_library(dv_bnns_finalize STATIC bnns_finalize.mm)
set_source_files_properties(bnns_finalize.mm PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(dv_bnns_finalize PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_bnns_finalize PUBLIC
  dv_weights
  absl::log
)

# ---------------------------------------------------------------------------
# dv_metal_conv_serial — Phase 5.5c deterministic-reduction-order Conv2D
# kernel (compiled at runtime via newLibraryWithSource:). Used to
# selectively replace MPSGraph convolution2D for layers that drift
# beyond the FILTER-threshold sensitivity vs Docker.
# ---------------------------------------------------------------------------
add_library(dv_metal_conv_serial STATIC metal_conv_serial.mm)
set_source_files_properties(metal_conv_serial.mm PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(dv_metal_conv_serial PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_metal_conv_serial PUBLIC
  dv_metal_conv_kahan       # Path B: MetalConvSerial::Encode delegates to
                            # MetalConvKahan when DV_METAL_KAHAN=1 is set.
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# dv_metal_conv_kahan — Phase 5.5e/Path B Kahan-compensated Conv2D
# kernel. Same dispatch shape as conv_serial; per-thread accumulator
# uses TwoSum compensation to bound per-step error at O(ε² · |sum|)
# instead of O(ε · |sum|). Provably bit-deterministic across reduction
# orders within ~1 ULP — sufficient to match Docker FILTER classes
# regardless of their AVX-512 vs scalar reduction strategy.
# ---------------------------------------------------------------------------
add_library(dv_metal_conv_kahan STATIC metal_conv_kahan.mm)
set_source_files_properties(metal_conv_kahan.mm PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(dv_metal_conv_kahan PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
# NB: dv_metal_conv_kahan only needs ConvDesc from metal_conv_serial.h (a
# compile-time include available via target_include_directories), not any
# compiled MetalConvSerial symbol — so there is no link edge back to
# dv_metal_conv_serial. Adding one would create a static-library cycle
# (serial -> kahan is the only real edge, via MetalConvKahan::Create) that
# Apple ld reports as "ignoring duplicate libraries".
target_link_libraries(dv_metal_conv_kahan PUBLIC
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# dv_metal_det_mixed — Phase 8/Tier 6.0 deterministic Inception block
# dispatch. Wraps MetalConvSerial + MetalBnRelu + MetalAvgPool +
# MetalConcat into per-block builders + dispatchers (Mixed_5b first,
# scaled to all 11 blocks 5b–7c).
# ---------------------------------------------------------------------------
add_library(dv_metal_det_mixed STATIC metal_det_mixed.mm)
set_source_files_properties(metal_det_mixed.mm PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(dv_metal_det_mixed PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_metal_det_mixed PUBLIC
  dv_weights
  dv_metal_conv_serial
  dv_metal_det_kernels
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# dv_metal_det_kernels — Phase 5.5e deterministic AvgPool / Concat /
# GlobalAvgPool kernels needed to extend Phase 5.5c to the full
# Inception-v3 stack (blocks 5b–7c + global-avg-pool). All embed
# kernel source as strings + compile at runtime; same per-thread
# strict-serial accumulation contract as conv_serial_fp32.
# ---------------------------------------------------------------------------
add_library(dv_metal_det_kernels STATIC
  metal_avg_pool.mm
  metal_concat.mm
  metal_global_avg_pool.mm
  metal_bn_relu.mm
)
set_source_files_properties(
  metal_avg_pool.mm metal_concat.mm metal_global_avg_pool.mm
  metal_bn_relu.mm
  PROPERTIES COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(dv_metal_det_kernels PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_metal_det_kernels PUBLIC
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# dv_metal_inference — Phase 5.5 MPSGraph + BNNS Inception-v3 backend
# (Obj-C++ + Metal Performance Shaders Graph; consumes .dvw weights)
# ---------------------------------------------------------------------------
add_library(dv_metal_inference STATIC metal_inference.mm)
set_source_files_properties(metal_inference.mm PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(dv_metal_inference PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_metal_inference PUBLIC
  dv_weights
  dv_metal_conv_serial
  dv_metal_det_kernels
  dv_metal_det_mixed
  absl::log
  "-framework Metal"
  "-framework MetalPerformanceShadersGraph"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# dv_coreml — Obj-C++ Core ML inference wrapper
# (Obj-C++ requires macOS frameworks; no Python, no TF)
# ---------------------------------------------------------------------------
add_library(dv_coreml STATIC coreml_inference.mm)
set_source_files_properties(coreml_inference.mm PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(dv_coreml PUBLIC
  "${CMAKE_SOURCE_DIR}"
)
target_link_libraries(dv_coreml PUBLIC
  "-framework CoreML"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# dv_small_model — deterministic FP32 BNNS-CPU MLP runner for the
# small_model (70 → 750 → 750 → 3). Reads weights from .npy files
# (Phase 5.5d/7); replaces the Core ML path which had ~0.005-0.01 drift
# vs Docker's TF/Keras output and caused cross-MID FILTER flips at
# threshold-borderline sites.
# ---------------------------------------------------------------------------
add_library(dv_small_model STATIC
  small_model_inference.mm
  small_model_features.cc
)
target_include_directories(dv_small_model PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_small_model PUBLIC
  proto_dv
  proto_nucleus
  absl::log
)

# ---------------------------------------------------------------------------
# dv_call_variants_lib — the call_variants logic (C++, no Obj-C)
# ---------------------------------------------------------------------------
add_library(dv_call_variants_lib STATIC call_variants_main.cc)
target_include_directories(dv_call_variants_lib PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${ABSL_PREFIX}/include"
  "${CMAKE_SOURCE_DIR}/cmake/tf_stubs"
)
target_compile_options(dv_call_variants_lib PRIVATE
  "-include${CMAKE_SOURCE_DIR}/cmake/tf_stubs/tensorflow/core/platform/tf_compat.h"
)
target_link_libraries(dv_call_variants_lib PUBLIC
  dv_tfrecord
  dv_coreml
  dv_metal_inference
  dv_bnns_finalize
  proto_dv
  proto_nucleus
  absl::flags
  absl::flags_parse
  absl::log
  absl::strings
)

# ---------------------------------------------------------------------------
# Phase 2 smoke test
# (also in tests/native/ — added later after this lib compiles)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Phase 3 — TF-free replacement for nucleus::ExampleWriter
# (originally implemented against tensorflow::io::RecordWriter)
# ---------------------------------------------------------------------------
add_library(dv_example_writer STATIC
  "${CMAKE_SOURCE_DIR}/patches/example_writer_macos.cc"
)
target_include_directories(dv_example_writer PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${CMAKE_SOURCE_DIR}/cmake/tf_stubs"
  "${ABSL_PREFIX}/include"
)
target_compile_options(dv_example_writer PRIVATE
  "-include${CMAKE_SOURCE_DIR}/cmake/tf_stubs/tensorflow/core/platform/tf_compat.h"
)
target_link_libraries(dv_example_writer PUBLIC
  dv_tfrecord
  proto_nucleus
  absl::log
  absl::strings
  absl::status
)

# ---------------------------------------------------------------------------
# Phase 3 — make_examples orchestration
# ---------------------------------------------------------------------------
set(ME_INCLUDE_DIRS
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${CMAKE_SOURCE_DIR}/cmake/tf_stubs"
  "${ABSL_PREFIX}/include"
  "${BOOST_INCLUDE_DIR}"
)
set(ME_COMPILE_OPTS
  "-include${CMAKE_SOURCE_DIR}/cmake/tf_stubs/tensorflow/core/platform/tf_compat.h"
)

# ---------------------------------------------------------------------------
# dv_haploid_regions — sex-chromosome haploid-calling helpers shared by
# make_examples (gVCF) and postprocess. The header is dependency-free and
# header-only; only LoadParRegions lives here because it uses zlib to read
# plaintext-or-gzipped PAR BEDs.
# ---------------------------------------------------------------------------
add_library(dv_haploid_regions STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/haploid_regions.cc"
)
target_include_directories(dv_haploid_regions PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(dv_haploid_regions PUBLIC
  absl::strings
  z  # zlib (also pulled in transitively by htslib); for gzopen/gzgets
)

add_library(dv_make_examples_lib STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/make_examples_main.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/regions.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/realigner_native.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/gvcf_emit.cc"
)
target_include_directories(dv_make_examples_lib PUBLIC ${ME_INCLUDE_DIRS})
target_compile_options(dv_make_examples_lib PRIVATE ${ME_COMPILE_OPTS})
target_link_libraries(dv_make_examples_lib PUBLIC
  dv_make_examples_native
  dv_allelecounter
  dv_tfrecord
  dv_small_model
  dv_direct_phasing       # Phase 9 / Step 4 — per-region read phasing
  dv_haploid_regions      # haploid gVCF reference confidence (PAR BED loader)
  dv_methylation_aware_phasing  # --enable_methylation_aware_phasing
  realigner
  proto_dv
  proto_nucleus
  absl::flags
  absl::flags_parse
  absl::log
  absl::strings
  nucleus_io
)

# ---------------------------------------------------------------------------
# dv_methylation_aware_phasing — faithful C++ port of upstream's
# methylation_aware_phasing.{cc,h} (Wilcoxon rank-sum 5mC haplotype voting).
#
# It is unit-tested here (so it cannot bitrot against the native proto
# definitions) and linked into dv_make_examples_lib, which invokes
# PerformMethylationAwarePhasing from the single-sample phasing block under
# --enable_methylation_aware_phasing (mirroring make_examples_core.py's
# methylated-ref-site partitioning + post-DirectPhasing call).
# ---------------------------------------------------------------------------
add_library(dv_methylation_aware_phasing STATIC
  "${CMAKE_SOURCE_DIR}/deepvariant/methylation_aware_phasing.cc"
)
target_include_directories(dv_methylation_aware_phasing PUBLIC ${ME_INCLUDE_DIRS})
target_compile_options(dv_methylation_aware_phasing PRIVATE ${ME_COMPILE_OPTS})
target_link_libraries(dv_methylation_aware_phasing PUBLIC
  proto_dv
  proto_nucleus
  absl::log
  absl::strings
  absl::flat_hash_map   # raw_hash_set internals referenced by the .cc
  absl::flat_hash_set
)

# ---------------------------------------------------------------------------
# Phase 3 — postprocess_variants orchestration
# ---------------------------------------------------------------------------
add_library(dv_postprocess_lib STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/postprocess_main.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/haplotypes.cc"
)
target_include_directories(dv_postprocess_lib PUBLIC ${ME_INCLUDE_DIRS})
target_compile_options(dv_postprocess_lib PRIVATE ${ME_COMPILE_OPTS})
target_link_libraries(dv_postprocess_lib PUBLIC
  dv_tfrecord
  dv_haploid_regions      # haploid het-zeroing (PAR BED loader)
  proto_dv
  proto_nucleus
  absl::flags
  absl::flags_parse
  absl::log
  absl::strings
  nucleus_io
)

# ---------------------------------------------------------------------------
# debug_metal — dev-only diagnostic for Phase 5.5 MPSGraph divergence.
# Compares Metal stem_s1a output to the hand-computed reference for an
# all-zeros input.
# ---------------------------------------------------------------------------
add_executable(debug_metal "${CMAKE_CURRENT_SOURCE_DIR}/debug_metal_main.cc")
target_include_directories(debug_metal PRIVATE ${ME_INCLUDE_DIRS})
target_compile_options(debug_metal PRIVATE ${ME_COMPILE_OPTS})
target_link_libraries(debug_metal PRIVATE
  dv_metal_inference  # transitively links dv_weights (PUBLIC); do not repeat it
)

# ---------------------------------------------------------------------------
# microtest_numpy_rng — Phase 5.5d/3 verification that NumpyMt19937 +
# BoundedLemireUint32 reproduce NumPy 1.24's
# np.random.RandomState(seed).randint(...) bit-for-bit.
# ---------------------------------------------------------------------------
add_executable(microtest_numpy_rng
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_numpy_rng.cc"
)
target_include_directories(microtest_numpy_rng PRIVATE
  "${CMAKE_SOURCE_DIR}"
)

# ---------------------------------------------------------------------------
# microtest_neon_cigar_classify — A2.2 verification that the NEON
# M-block byte classifier (use_base, is_low_quality, is_ref, canonical)
# produces output byte-identical to the scalar reference across all
# (read,ref) byte pairs, all qual boundary values, both legacy and
# non-legacy CanBasesBeUsed semantics, and all length tails 0..1024.
# ---------------------------------------------------------------------------
add_executable(microtest_neon_cigar_classify
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_neon_cigar_classify.cc"
)
target_include_directories(microtest_neon_cigar_classify PRIVATE
  "${CMAKE_SOURCE_DIR}"
)

# ---------------------------------------------------------------------------
# microtest_bnns_stem — Option-2 PoC: scalar BNNS-CPU stem_s1a vs TF
# Docker AVX-512 reference. Decides whether to invest in a full
# BNNS-CPU Inception-v3 forward pass for borderline-only re-evaluation.
# ---------------------------------------------------------------------------
add_executable(microtest_bnns_stem
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_bnns_stem.cc"
)
target_include_directories(microtest_bnns_stem PRIVATE
  "${CMAKE_SOURCE_DIR}"
)
target_link_libraries(microtest_bnns_stem PRIVATE
  dv_weights
)

# ---------------------------------------------------------------------------
# microtest_conv_serial — Phase 5.5c hand-verifiable test for the
# deterministic Conv2D kernel. Compares GPU dispatch output against a
# scalar (kh,kw,c_in)-order CPU reference using std::fma — bit-exact
# match expected on healthy build.
# ---------------------------------------------------------------------------
add_executable(microtest_conv_serial
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_conv_serial.mm"
)
set_source_files_properties(
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_conv_serial.mm" PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(microtest_conv_serial PRIVATE
  ${ME_INCLUDE_DIRS}
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(microtest_conv_serial PRIVATE
  dv_metal_conv_serial
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# microtest_det_mixed5b — Phase 8/Tier 6.0 validation: det Mixed_5b
# block dispatch vs TF reference (/tmp/dv_per_layer/{stem_mp5a,5b}.npy).
# ---------------------------------------------------------------------------
add_executable(microtest_det_mixed5b
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_det_mixed5b.mm"
)
set_source_files_properties(
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_det_mixed5b.mm" PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(microtest_det_mixed5b PRIVATE
  ${ME_INCLUDE_DIRS}
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(microtest_det_mixed5b PRIVATE
  dv_metal_det_mixed
  dv_weights
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# microtest_det_inception — full chain of 11 Mixed_X det blocks vs TF
# reference (per-block max_abs/mean_abs/max_rel).
# ---------------------------------------------------------------------------
add_executable(microtest_det_inception
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_det_inception.mm"
)
set_source_files_properties(
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_det_inception.mm" PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(microtest_det_inception PRIVATE
  ${ME_INCLUDE_DIRS}
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(microtest_det_inception PRIVATE
  dv_metal_det_mixed
  dv_weights
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# microtest_conv_kahan — Phase 5.5e/Path B hand-verifiable test for the
# Kahan-compensated Conv2D kernel. Compares GPU dispatch output against
# a scalar Kahan reference and reports the precision improvement vs
# basic-FMA scalar reference.
# ---------------------------------------------------------------------------
add_executable(microtest_conv_kahan
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_conv_kahan.mm"
)
set_source_files_properties(
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_conv_kahan.mm" PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(microtest_conv_kahan PRIVATE
  ${ME_INCLUDE_DIRS}
  "${CMAKE_SOURCE_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(microtest_conv_kahan PRIVATE
  dv_metal_conv_kahan
  absl::log
  "-framework Metal"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# extract_pileup_npy — dev-only profiling tool for Phase 5.5c per-layer
# drift work. Reads N pileup images from a TFRecord (or `name@N` shard
# spec) and writes a NumPy `(N, 100, 221, 7)` FP32 NHWC array. The
# resulting `.npy` is consumed by `dump_tf_per_layer.py` (Docker) and
# `debug_metal --compare-to-reference` to profile drift on real data.
# ---------------------------------------------------------------------------
add_executable(extract_pileup_npy
  "${CMAKE_CURRENT_SOURCE_DIR}/extract_pileup_npy_main.cc"
)
target_include_directories(extract_pileup_npy PRIVATE ${ME_INCLUDE_DIRS})
target_link_libraries(extract_pileup_npy PRIVATE dv_tfrecord)

# extract_pileup_at_pos — locate a specific (chrom, start, ref, alt)
# in an examples.tfrecord and dump that single pileup as a (1,100,221,7)
# NHWC FP32 .npy. Used to byte-compare our pileup vs Docker's at a
# PASS-flip site (Phase 5.5c diagnostic).
add_executable(extract_pileup_at_pos
  "${CMAKE_CURRENT_SOURCE_DIR}/extract_pileup_at_pos_main.cc"
)
target_include_directories(extract_pileup_at_pos PRIVATE ${ME_INCLUDE_DIRS})
target_link_libraries(extract_pileup_at_pos PRIVATE dv_tfrecord)

# ---------------------------------------------------------------------------
# microtest_metal — hand-verifiable MPSGraph conv micro-tests for
# Phase 5.5a investigation. Builds tiny graphs (1×1 and 3×3) with
# inputs / weights small enough to compute the expected output by
# pencil-and-paper, and prints PASS/FAIL per test.
# ---------------------------------------------------------------------------
add_executable(microtest_metal
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_main.mm"
)
set_source_files_properties(
  "${CMAKE_CURRENT_SOURCE_DIR}/microtest_main.mm" PROPERTIES
  COMPILE_FLAGS "-fobjc-arc"
)
target_include_directories(microtest_metal PRIVATE ${ME_INCLUDE_DIRS})
target_link_libraries(microtest_metal PRIVATE
  "-framework Metal"
  "-framework MetalPerformanceShadersGraph"
  "-framework Foundation"
)

# ---------------------------------------------------------------------------
# dump_cvo — dev-only TFRecord dumper for CallVariantsOutput protos.
# Used to diff our candidate set against upstream's during parity work.
# ---------------------------------------------------------------------------
add_executable(dump_cvo
  "${CMAKE_CURRENT_SOURCE_DIR}/dump_cvo_main.cc"
)
target_include_directories(dump_cvo PRIVATE ${ME_INCLUDE_DIRS})
target_compile_options(dump_cvo PRIVATE ${ME_COMPILE_OPTS})
target_link_libraries(dump_cvo PRIVATE
  dv_tfrecord
  proto_dv
  proto_nucleus
  absl::log
  absl::strings
  absl::hash
)

# ---------------------------------------------------------------------------
# dump_allele_counts — dev-only AlleleCount dumper.
# Useful for diffing per-position ref/alt counts between our pipeline
# and upstream during candidate-set parity work.
# ---------------------------------------------------------------------------
add_executable(dump_allele_counts
  "${CMAKE_CURRENT_SOURCE_DIR}/dump_allele_counts_main.cc"
)
target_include_directories(dump_allele_counts PRIVATE ${ME_INCLUDE_DIRS})
target_compile_options(dump_allele_counts PRIVATE ${ME_COMPILE_OPTS})
target_link_libraries(dump_allele_counts PRIVATE
  dv_allelecounter
  proto_dv
  proto_nucleus
  nucleus_io
  absl::log
  absl::strings
  absl::hash
  absl::flat_hash_map
  absl::flat_hash_set
  absl::raw_hash_set
)

# ---------------------------------------------------------------------------
# deepvariant — main binary (CLI dispatcher)
# ---------------------------------------------------------------------------
add_executable(deepvariant
  "${CMAKE_CURRENT_SOURCE_DIR}/cli.cc"
)
target_include_directories(deepvariant PRIVATE ${ME_INCLUDE_DIRS})
target_compile_options(deepvariant PRIVATE ${ME_COMPILE_OPTS})

# Version metadata — captured at configure time and re-captured at every
# build via a small `git rev-parse` shim, so `deepvariant --version` always
# reflects the actual built tree (not just whatever was checked out when
# CMake was first run). Three sources:
#
#   DV_VERSION           — hand-bumped semver-ish tag for this fork
#                          (currently "v2-applesilicon"). Bumped at release.
#   DV_UPSTREAM_VERSION  — Google DeepVariant version we mirror
#                          (currently "1.10.0"). Bumped when we re-baseline.
#   DV_GIT_SHA           — short SHA of HEAD at configure time. Re-checked
#                          on every build via add_custom_target below; if
#                          the SHA drifts, cli.cc gets recompiled.
#   DV_BUILD_DATE        — ISO 8601 date of the build (UTC).
execute_process(
  COMMAND git rev-parse --short HEAD
  WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
  OUTPUT_VARIABLE DV_GIT_SHA
  OUTPUT_STRIP_TRAILING_WHITESPACE
  ERROR_QUIET)
if(NOT DV_GIT_SHA)
  set(DV_GIT_SHA "unknown")
endif()
string(TIMESTAMP DV_BUILD_DATE "%Y-%m-%d" UTC)
target_compile_definitions(deepvariant PRIVATE
  DV_VERSION="v2-applesilicon"
  DV_UPSTREAM_VERSION="1.10.0"
  DV_GIT_SHA="${DV_GIT_SHA}"
  DV_BUILD_DATE="${DV_BUILD_DATE}"
)
target_link_libraries(deepvariant PRIVATE
  dv_make_examples_lib
  dv_postprocess_lib
  dv_call_variants_lib
  dv_example_writer
  absl::flags
  absl::flags_parse
  absl::log
  absl::log_initialize
)

# Multi-call binary symlinks — busybox-style. Same binary, dispatched by
# basename(argv[0]) inside cli.cc::DetectMultiCall. Mirrors upstream's
# three-binary convention (run_deepvariant / run_deeptrio / run_deepsomatic /
# run_pangenome_aware_deepvariant) without the disk-bloat / version-skew
# cost of three separate executables.
#
# After build: build-macos/bin/{deeptrio,deepsomatic,pangenome-aware-deepvariant}
# all point at build-macos/bin/deepvariant.
#
# After `make install`: same layout under ${CMAKE_INSTALL_PREFIX}/bin/.
# Homebrew formula will use `bin.install_symlink "deepvariant" => "deeptrio"`
# (etc.) to mirror this in the final bottle.
add_custom_command(TARGET deepvariant POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          deepvariant
          "$<TARGET_FILE_DIR:deepvariant>/deeptrio"
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          deepvariant
          "$<TARGET_FILE_DIR:deepvariant>/deepsomatic"
  COMMAND ${CMAKE_COMMAND} -E create_symlink
          deepvariant
          "$<TARGET_FILE_DIR:deepvariant>/pangenome-aware-deepvariant"
  COMMENT "Creating multi-call binary symlinks (deeptrio, deepsomatic, pangenome-aware-deepvariant)"
  VERBATIM)
