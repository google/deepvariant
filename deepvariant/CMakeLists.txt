# deepvariant/ — upstream C++ libraries compiled TF-free for the native port.
# Does NOT compile training-only, Python bindings, or test files.

set(DV_INCLUDE_DIRS
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${CMAKE_SOURCE_DIR}/cmake/tf_stubs"
  "${ABSL_PREFIX}/include"
  "${BOOST_INCLUDE_DIR}"
)
set(DV_COMPILE_OPTS
  "-include${CMAKE_SOURCE_DIR}/cmake/tf_stubs/tensorflow/core/platform/tf_compat.h"
)

# ---------------------------------------------------------------------------
# Helper: compile one deepvariant static library
# ---------------------------------------------------------------------------
function(dv_library name)
  cmake_parse_arguments(ARG "" "" "SRCS;DEPS" ${ARGN})
  add_library(${name} STATIC ${ARG_SRCS})
  target_include_directories(${name} PUBLIC ${DV_INCLUDE_DIRS})
  target_compile_options(${name} PRIVATE ${DV_COMPILE_OPTS})
  if(ARG_DEPS)
    target_link_libraries(${name} PUBLIC ${ARG_DEPS})
  endif()
endfunction()

# ---------------------------------------------------------------------------
# dv_utils — deepvariant/utils.cc
# ---------------------------------------------------------------------------
dv_library(dv_utils
  SRCS utils.cc
  DEPS proto_dv absl::strings
)

# ---------------------------------------------------------------------------
# dv_channels — all pileup channel implementations
# ---------------------------------------------------------------------------
file(GLOB DV_CHANNEL_SRCS CONFIGURE_DEPENDS "channels/*.cc")
list(FILTER DV_CHANNEL_SRCS EXCLUDE REGEX "_test\\.cc$")
dv_library(dv_channels
  SRCS ${DV_CHANNEL_SRCS}
  DEPS proto_dv absl::strings absl::log absl::check nucleus_io
)

# ---------------------------------------------------------------------------
# dv_pileup — pileup_channel_lib + pileup_image_native + alt_aligned_pileup
# ---------------------------------------------------------------------------
dv_library(dv_pileup
  SRCS
    pileup_channel_lib.cc
    pileup_image_native.cc
    alt_aligned_pileup_lib.cc
  DEPS
    dv_channels
    dv_utils
    realigner  # for fast_pass_aligner used in alt-aligned pileup
    proto_dv
    proto_nucleus
    absl::algorithm
    absl::flat_hash_map
    absl::flat_hash_set
    absl::btree
    absl::log
    absl::strings
    nucleus_io
)

# ---------------------------------------------------------------------------
# dv_allelecounter — AlleleCounter + VariantCaller
# ---------------------------------------------------------------------------
dv_library(dv_allelecounter
  SRCS
    allelecounter.cc
    variant_calling.cc
    variant_calling_multisample.cc
  DEPS
    dv_utils
    proto_dv
    proto_nucleus
    absl::log
    absl::strings
    nucleus_io
)

# ---------------------------------------------------------------------------
# dv_stream_examples_stub — StreamExamples compiled with Boost IPC headers.
# stream_examples_ is always nullptr in native mode (options.stream_examples()
# returns false), so these methods are never called at runtime.
# ---------------------------------------------------------------------------
add_library(dv_stream_examples STATIC stream_examples.cc)
target_include_directories(dv_stream_examples PUBLIC
  ${DV_INCLUDE_DIRS}
)
target_compile_options(dv_stream_examples PRIVATE ${DV_COMPILE_OPTS})
target_link_libraries(dv_stream_examples PUBLIC
  proto_dv
  absl::log
  absl::strings
)

# ---------------------------------------------------------------------------
# dv_make_examples_native — ExamplesGenerator (the C++ pileup encoder).
# ---------------------------------------------------------------------------
dv_library(dv_make_examples_native
  SRCS make_examples_native.cc
  DEPS
    dv_pileup
    dv_stream_examples
    proto_dv
    proto_nucleus
    absl::flat_hash_map
    absl::flat_hash_set
    absl::log
    absl::strings
    nucleus_io
    re2::re2
)

# ---------------------------------------------------------------------------
# dv_direct_phasing — DirectPhasing for phase_reads (optional for v1.0)
# ---------------------------------------------------------------------------
dv_library(dv_direct_phasing
  SRCS direct_phasing.cc
  DEPS
    dv_allelecounter
    proto_dv
    proto_nucleus
    absl::log
    absl::strings
    nucleus_io
)
