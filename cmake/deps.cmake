# deps.cmake — All external C++ dependencies (no TensorFlow).
#
# All major deps use Homebrew (already installed) via find_package.
# Only libssw (not in Homebrew) uses FetchContent.
#
# Homebrew versions on this machine:
#   htslib  1.18        (req: 1.18)
#   abseil  20260107.1  (req: ≥ 20240722; API-compatible)
#   protobuf 34.1       (req: 21.9; API-compatible for generated code)
#
# Pangenome deps (gbwt, gbwtgraph, sdsl-lite, libdivsufsort, libhandlegraph)
# are deferred until Phase 3 (pangenome-aware DeepVariant port).

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)
set(FETCHCONTENT_UPDATES_DISCONNECTED ON)

# ---------------------------------------------------------------------------
# htslib 1.18 — Homebrew (avoids autoconf complexity on macOS)
# ---------------------------------------------------------------------------
find_program(BREW_EXECUTABLE brew REQUIRED)
execute_process(
  COMMAND ${BREW_EXECUTABLE} --prefix htslib
  OUTPUT_VARIABLE HTSLIB_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT HTSLIB_PREFIX)
  message(FATAL_ERROR "htslib not found — run: brew install htslib")
endif()

add_library(htslib::htslib STATIC IMPORTED)
find_library(HTSLIB_LIB NAMES libhts.a hts PATHS "${HTSLIB_PREFIX}/lib" REQUIRED)
set_target_properties(htslib::htslib PROPERTIES
  IMPORTED_LOCATION "${HTSLIB_LIB}"
  INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_PREFIX}/include"
)
# Resolve libdeflate via Homebrew rather than a hardcoded /opt/homebrew path,
# so the build works under a non-default Homebrew prefix or keg-only layout.
execute_process(
  COMMAND ${BREW_EXECUTABLE} --prefix libdeflate
  OUTPUT_VARIABLE LIBDEFLATE_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
find_library(LIBDEFLATE_LIB NAMES libdeflate.a deflate
  PATHS "${LIBDEFLATE_PREFIX}/lib" REQUIRED)
target_link_libraries(htslib::htslib INTERFACE
  "-framework CoreFoundation"
  "${LIBDEFLATE_LIB}"
  z bz2 lzma curl
)
message(STATUS "htslib: ${HTSLIB_LIB}")

# ---------------------------------------------------------------------------
# abseil — Homebrew (no FetchContent; avoids hash management)
# ---------------------------------------------------------------------------
execute_process(
  COMMAND ${BREW_EXECUTABLE} --prefix abseil
  OUTPUT_VARIABLE ABSL_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT ABSL_PREFIX)
  message(FATAL_ERROR "abseil not found — run: brew install abseil")
endif()
list(APPEND CMAKE_PREFIX_PATH "${ABSL_PREFIX}")
find_package(absl REQUIRED)
message(STATUS "abseil: ${ABSL_PREFIX}")

# ---------------------------------------------------------------------------
# protobuf — Homebrew (no FetchContent; avoids hash management)
# ---------------------------------------------------------------------------
execute_process(
  COMMAND ${BREW_EXECUTABLE} --prefix protobuf
  OUTPUT_VARIABLE PROTOBUF_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT PROTOBUF_PREFIX)
  message(FATAL_ERROR "protobuf not found — run: brew install protobuf")
endif()
list(APPEND CMAKE_PREFIX_PATH "${PROTOBUF_PREFIX}")
find_package(protobuf REQUIRED)
message(STATUS "protobuf: ${PROTOBUF_PREFIX}")

# Homebrew's protoc.
find_program(PROTOC protoc HINTS "${PROTOBUF_PREFIX}/bin" REQUIRED)
message(STATUS "protoc: ${PROTOC}")

# ---------------------------------------------------------------------------
# libssw 1.2.5 — Smith-Waterman aligner (realigner/)
# ---------------------------------------------------------------------------
FetchContent_Declare(
  libssw
  URL      https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/v1.2.5.tar.gz
  URL_HASH SHA256=b294c0cb6f0f3d578db11b4112a88b20583b9d4190b0a9cf04d83bb6a8704d9a
)
# libssw ships no CMakeLists, so MakeAvailable just populates ${libssw_SOURCE_DIR}
# (no add_subdirectory) — and avoids the deprecated single-arg
# FetchContent_Populate(libssw) call.
FetchContent_MakeAvailable(libssw)

# OVERLAY: replace the vendored sse2neon.h (Ratcliff/NVIDIA early version,
# 8798 lines, missing fixes) with the modern DLTcollab fork (11744 lines,
# improved fidelity for edge cases like _mm_slli_si128 byte-shifts).
# This reduces realigner SSW score drift between native arm64 (compile-time
# SSE→NEON) and Docker on Rosetta (runtime SSE→ARM translation), which
# was the source of 105/120 PASS-flips on chr20:26-31Mb pericentromere.
# See PORT_LOG 2026-05-07 "PASS-flip root-cause analysis".
if(EXISTS "${CMAKE_SOURCE_DIR}/release/vendored/sse2neon.h")
  configure_file(
    "${CMAKE_SOURCE_DIR}/release/vendored/sse2neon.h"
    "${libssw_SOURCE_DIR}/src/sse2neon.h"
    COPYONLY)
  message(STATUS "libssw: overlaid modern sse2neon.h from release/vendored/")
endif()

# libssw has no CMakeLists — define targets here.
add_library(ssw STATIC
  "${libssw_SOURCE_DIR}/src/ssw.c"
  "${libssw_SOURCE_DIR}/src/ssw.h"
  "${libssw_SOURCE_DIR}/src/ssw_cpp.cpp"
  "${libssw_SOURCE_DIR}/src/ssw_cpp.h"
)
# deepvariant/realigner/ssw.h uses #include "src/ssw_cpp.h",
# so the PARENT of src/ must be on the include path, not just src/.
target_include_directories(ssw PUBLIC "${libssw_SOURCE_DIR}")
# Apple Clang/arm64: SSW uses SSE2 intrinsics guarded by __SSE2__ —
# arm64 does not have SSE2; the fallback scalar path is used automatically.
set(DV_LIBSSW_DIR "${libssw_SOURCE_DIR}" CACHE INTERNAL "libssw source root")

# ---------------------------------------------------------------------------
# re2 — Homebrew
# ---------------------------------------------------------------------------
execute_process(
  COMMAND ${BREW_EXECUTABLE} --prefix re2
  OUTPUT_VARIABLE RE2_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT RE2_PREFIX)
  message(FATAL_ERROR "re2 not found — run: brew install re2")
endif()
add_library(re2::re2 STATIC IMPORTED)
find_library(RE2_LIB NAMES libre2.a re2 PATHS "${RE2_PREFIX}/lib" REQUIRED)
set_target_properties(re2::re2 PROPERTIES
  IMPORTED_LOCATION "${RE2_LIB}"
  INTERFACE_INCLUDE_DIRECTORIES "${RE2_PREFIX}/include"
)
message(STATUS "re2: ${RE2_LIB}")

# ---------------------------------------------------------------------------
# Boost — Homebrew (for debruijn_graph.h in realigner/)
# ---------------------------------------------------------------------------
execute_process(
  COMMAND ${BREW_EXECUTABLE} --prefix boost
  OUTPUT_VARIABLE BOOST_PREFIX
  OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NOT BOOST_PREFIX)
  message(FATAL_ERROR "boost not found — run: brew install boost")
endif()
message(STATUS "boost: ${BOOST_PREFIX}")
set(BOOST_INCLUDE_DIR "${BOOST_PREFIX}/include" CACHE INTERNAL "")

# ---------------------------------------------------------------------------
# GoogleTest — FetchContent (no standalone Homebrew package)
# ---------------------------------------------------------------------------
FetchContent_Declare(
  googletest
  URL      https://github.com/google/googletest/archive/refs/tags/v1.14.0.tar.gz
  URL_HASH SHA256=8ad598c73ad796e0d8280b082cebd82a630d73e73cd3c70057938a6501bba5d7
)
set(INSTALL_GTEST OFF)
FetchContent_MakeAvailable(googletest)
set(GTEST_PREFIX "${googletest_SOURCE_DIR}" CACHE INTERNAL "")

# ---------------------------------------------------------------------------
# zlib — guaranteed present on macOS (from Xcode SDK)
# ---------------------------------------------------------------------------
find_package(ZLIB REQUIRED)
