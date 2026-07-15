# deepvariant/realigner — Smith-Waterman + de Bruijn graph realigner.
#
# Depends on: libssw, abseil, proto_nucleus, proto_dv.
# No TF, no pybind11.

add_library(realigner STATIC
  "${CMAKE_CURRENT_SOURCE_DIR}/debruijn_graph.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/fast_pass_aligner.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/ssw.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/window_selector.cc"
)

target_include_directories(realigner PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${CMAKE_SOURCE_DIR}/cmake/tf_stubs"
  "${ABSL_PREFIX}/include"
  "${BOOST_INCLUDE_DIR}"
)
target_compile_options(realigner PRIVATE
  "-include${CMAKE_SOURCE_DIR}/cmake/tf_stubs/tensorflow/core/platform/tf_compat.h"
)

target_link_libraries(realigner PUBLIC
  ssw
  nucleus_io
  proto_dv
  proto_nucleus
  absl::strings
  absl::log
  absl::check
  absl::flat_hash_map
  absl::flat_hash_set
)

# Tests deferred to tests/native/ smoke tests.

