cmake_minimum_required(VERSION 3.27)
project(deepvariant VERSION 1.10.0 LANGUAGES CXX OBJCXX C)

# ---------------------------------------------------------------------------
# Guards: macOS arm64 only.
# ---------------------------------------------------------------------------
if(NOT APPLE OR NOT CMAKE_SYSTEM_PROCESSOR STREQUAL "arm64")
  message(FATAL_ERROR "This build targets macOS arm64 only.")
endif()
if(CMAKE_SYSTEM_VERSION VERSION_LESS "23")  # macOS 14 = Darwin 23.x
  message(FATAL_ERROR "macOS 14 (Sonoma) or newer required.")
endif()

# ---------------------------------------------------------------------------
# Language standards
# ---------------------------------------------------------------------------
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_OBJCXX_STANDARD 17)
set(CMAKE_OBJCXX_STANDARD_REQUIRED ON)

# Visibility: match TF convention — default hidden.
set(CMAKE_C_VISIBILITY_PRESET   hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

# All depedencies built as STATIC.
set(BUILD_SHARED_LIBS OFF)

# Default build type.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

# Build output goes to a single directory for easy inspection.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# ---------------------------------------------------------------------------
# Compiler flags — arm64, Clang (Apple Clang 21+)
# ---------------------------------------------------------------------------
add_compile_options(
  -arch arm64
  # Target the Apple-M microarchitecture for scheduling + the full ARMv8 ISA
  # (CRC32/dotprod/etc.). ISA/scheduling only — fast-math stays OFF globally,
  # so the deterministic conv_serial/conv_kahan/bnns reductions remain
  # bit-identical (guarded by the determinism microtests in ctest).
  -mcpu=apple-m1
  -Wall
  -Wextra
  -Wno-unused-parameter
  -Wno-missing-field-initializers
)

# ---------------------------------------------------------------------------
# Module path
# ---------------------------------------------------------------------------
list(PREPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# ---------------------------------------------------------------------------
# External dependencies (order matters: protos before nucleus)
# ---------------------------------------------------------------------------
include(deps)    # FetchContent / find_package for htslib, abseil, protobuf, ssw
include(protos)  # compile DeepVariant + nucleus + TF-example protos

# ---------------------------------------------------------------------------
# Core libraries (TF-free)
# ---------------------------------------------------------------------------
add_subdirectory(third_party/nucleus)
add_subdirectory(deepvariant/realigner)
add_subdirectory(deepvariant)          # upstream C++ libs (Phase 3)
add_subdirectory(deepvariant/native)   # runtime binary (Phase 2+)

# ---------------------------------------------------------------------------
# Tests (Phase 1 gate: ctest -V must pass)
# ---------------------------------------------------------------------------
enable_testing()
include(CTest)
add_subdirectory(tests/native)  # thin wrappers around upstream C++ test code
