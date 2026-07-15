# Phase 1 gate tests — smoke tests that nucleus_io + realigner link correctly.
# These run under ctest: cmake --build && ctest -V must pass.
# Excludes upstream tests that depend on test_utils.h (heavy TF I/O) which
# are deferred to Phase 1b once all TF stubs are complete.

set(GATE_INCLUDE_DIRS
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${CMAKE_SOURCE_DIR}/cmake/tf_stubs"
  "${ABSL_PREFIX}/include"
  "${BOOST_INCLUDE_DIR}"
  "${googletest_SOURCE_DIR}/googletest/include"
  "${googletest_SOURCE_DIR}/googlemock/include"
)
set(GATE_COMPILE_OPTS
  "-include${CMAKE_SOURCE_DIR}/cmake/tf_stubs/tensorflow/core/platform/tf_compat.h"
)

function(phase1_test name src lib)
  add_executable(${name} "${CMAKE_CURRENT_SOURCE_DIR}/${src}")
  target_include_directories(${name} PRIVATE ${GATE_INCLUDE_DIRS})
  target_compile_options(${name} PRIVATE ${GATE_COMPILE_OPTS})
  target_link_libraries(${name} PRIVATE ${lib} gtest_main gmock)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# Like phase1_test, but the source lives outside tests/native/ (it is the
# upstream test, shared verbatim). Pass an absolute path to the .cc.
function(phase1_test_abs name abs_src lib)
  add_executable(${name} "${abs_src}")
  target_include_directories(${name} PRIVATE ${GATE_INCLUDE_DIRS})
  target_compile_options(${name} PRIVATE ${GATE_COMPILE_OPTS})
  target_link_libraries(${name} PRIVATE ${lib} gtest_main gmock)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# Upstream methylation_aware_phasing unit test, run verbatim against the
# native build of the algorithm (dv_methylation_aware_phasing).
phase1_test_abs(methylation_aware_phasing_test
  "${CMAKE_SOURCE_DIR}/deepvariant/methylation_aware_phasing_test.cc"
  dv_methylation_aware_phasing)

# Header-only haploid-calling helpers (deepvariant/native/haploid_regions.h):
# needs only absl + gtest, no dv library.
add_executable(haploid_regions_test
  "${CMAKE_CURRENT_SOURCE_DIR}/haploid_regions_test.cc")
target_include_directories(haploid_regions_test PRIVATE ${GATE_INCLUDE_DIRS})
target_link_libraries(haploid_regions_test PRIVATE
  absl::strings gtest_main gmock)
add_test(NAME haploid_regions_test COMMAND haploid_regions_test)

phase1_test(nucleus_io_smoke      nucleus_io_smoke_test.cc      nucleus_io)
phase1_test(realigner_smoke       realigner_smoke_test.cc       realigner)
phase1_test(call_variants_smoke   call_variants_smoke_test.cc   dv_call_variants_lib)
phase1_test(small_model_smoke     small_model_smoke_test.cc     dv_small_model)
phase1_test(dv_weights_smoke      dv_weights_smoke_test.cc      dv_weights)
phase1_test(metal_inference_smoke metal_inference_smoke_test.cc dv_metal_inference)
phase1_test(bnns_finalize_smoke   bnns_finalize_smoke_test.cc   dv_bnns_finalize)

# Self-contained determinism / fidelity microtests (defined as executables in
# deepvariant/native/CMakeLists.txt). Each returns non-zero on mismatch, so it
# is a real ctest gate. The fixture-dependent microtests (det_mixed5b,
# det_inception, bnns_stem — they read /tmp/dv_per_layer/*.npy) are left
# unregistered since they need pre-generated per-layer dumps.
foreach(t
    microtest_numpy_rng
    microtest_neon_cigar_classify
    microtest_conv_serial
    microtest_conv_kahan)
  if(TARGET ${t})
    add_test(NAME ${t} COMMAND ${t})
  endif()
endforeach()
