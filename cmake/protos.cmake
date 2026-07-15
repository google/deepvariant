# protos.cmake — compile nucleus + deepvariant .proto files (no TF framework needed).
#
# nucleus/protos/ is self-contained:
#   example.proto  → imports feature.proto → defines tf.train.Example in namespace tensorflow
#   feature.proto  → no imports
# deepvariant/protos/ imports nucleus protos + google.protobuf builtins.
#
# All .pb.h / .pb.cc generated files land in ${CMAKE_BINARY_DIR}/proto_gen/,
# added to INTERFACE_INCLUDE_DIRECTORIES of proto_nucleus and proto_dv targets.

# PROTOC is set by deps.cmake (Homebrew protoc).
if(NOT PROTOC)
  find_program(PROTOC protoc REQUIRED HINTS "${PROTOBUF_PREFIX}/bin")
endif()

set(PROTO_GEN_DIR "${CMAKE_BINARY_DIR}/proto_gen")
file(MAKE_DIRECTORY "${PROTO_GEN_DIR}")

# dv_proto_compile(OUT_VAR PROTO_FILE PROTO_ROOT)
# PROTO_ROOT must be the directory you pass as --proto_path to protoc.
# Output files mirror the relative path under PROTO_ROOT inside PROTO_GEN_DIR.
function(dv_proto_compile OUT_SRC_VAR PROTO_FILE PROTO_ROOT)
  file(RELATIVE_PATH _rel "${PROTO_ROOT}" "${PROTO_FILE}")
  string(REGEX REPLACE "\\.proto$" ".pb.cc" _cc_rel "${_rel}")
  string(REGEX REPLACE "\\.proto$" ".pb.h"  _hh_rel "${_rel}")
  set(_cc "${PROTO_GEN_DIR}/${_cc_rel}")
  set(_hh "${PROTO_GEN_DIR}/${_hh_rel}")

  # Ensure output subdirectory exists.
  cmake_path(GET _cc PARENT_PATH _out_dir)
  file(MAKE_DIRECTORY "${_out_dir}")

  add_custom_command(
    OUTPUT  "${_cc}" "${_hh}"
    COMMAND "${PROTOC}"
            "--proto_path=${PROTO_ROOT}"
            "--cpp_out=${PROTO_GEN_DIR}"
            "${PROTO_FILE}"
    DEPENDS "${PROTO_FILE}" "${PROTOC}"
    VERBATIM
  )
  set(${OUT_SRC_VAR} "${${OUT_SRC_VAR}}" "${_cc}" PARENT_SCOPE)
endfunction()

# ---------------------------------------------------------------------------
# 1. nucleus protos (self-contained, no TF imports)
# ---------------------------------------------------------------------------
set(NUCLEUS_PROTO_ROOT "${CMAKE_SOURCE_DIR}/third_party/nucleus/protos")
file(GLOB NUCLEUS_PROTOS CONFIGURE_DEPENDS "${NUCLEUS_PROTO_ROOT}/*.proto")

set(NUCLEUS_PB_SRCS)
foreach(_p ${NUCLEUS_PROTOS})
  # proto_path = repo root so "third_party/nucleus/protos/..." resolves correctly.
  dv_proto_compile(NUCLEUS_PB_SRCS "${_p}" "${CMAKE_SOURCE_DIR}")
endforeach()

add_library(proto_nucleus STATIC ${NUCLEUS_PB_SRCS})
target_include_directories(proto_nucleus PUBLIC
  "${PROTO_GEN_DIR}"
  "${ABSL_PREFIX}/include"   # protobuf headers include absl/* transitively
)
target_link_libraries(proto_nucleus PUBLIC
  protobuf::libprotobuf
  absl::base
)

# Alias for compat with targets that link proto_tf_example separately.
# In our build the tf.train.Example type is in proto_nucleus (nucleus/protos/example.proto).
add_library(proto_tf_example ALIAS proto_nucleus)

# ---------------------------------------------------------------------------
# 2. deepvariant protos
# ---------------------------------------------------------------------------
set(DV_PROTO_ROOT "${CMAKE_SOURCE_DIR}/deepvariant/protos")
file(GLOB DV_PROTOS CONFIGURE_DEPENDS "${DV_PROTO_ROOT}/*.proto")

set(DV_PB_SRCS)
foreach(_p ${DV_PROTOS})
  # proto_path = repo root for both DV and nucleus imports.
  dv_proto_compile(DV_PB_SRCS "${_p}" "${CMAKE_SOURCE_DIR}")
endforeach()

add_library(proto_dv STATIC ${DV_PB_SRCS})
target_include_directories(proto_dv PUBLIC
  "${PROTO_GEN_DIR}"
  "${ABSL_PREFIX}/include"
)
target_link_libraries(proto_dv PUBLIC protobuf::libprotobuf proto_nucleus)
