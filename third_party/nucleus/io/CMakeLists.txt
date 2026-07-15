# nucleus_io — SAM/VCF/FASTA/BED I/O, TFRecord I/O (POSIX patch).
#
# Files excluded from this build (TF replacement pending Phase 2/3):
#   gfile.cc          → replaced by patches/gfile_macos.cc
#   tfrecord_reader.cc → replaced by patches/tfrecord_reader_macos.cc
#   tfrecord_writer.cc → replaced by patches/tfrecord_writer_macos.cc
#   example_writer.cc  → excluded (TF API; Phase 3 will provide a replacement)
#   gbz_reader.cc      → excluded (pangenome; deferred to Phase 3)
#   *_pybind.cc        → excluded (pybind11 / TF Python bindings; not needed)
#   *_test.cc          → excluded (test targets added below)

set(NUCLEUS_IO_SRCS
  "${CMAKE_CURRENT_SOURCE_DIR}/bed_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/bed_writer.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/bedgraph_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/bedgraph_writer.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/fastq_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/fastq_writer.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/gff_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/gff_writer.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/hts_path.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/hts_verbose.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/merge_variants.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/reader_base.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/reference.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/sam_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/sam_utils.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/sam_writer.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/tabix_indexer.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/text_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/text_writer.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/variant_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/vcf_concat.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/vcf_conversion.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/vcf_reader.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/vcf_writer.cc"
  # Patches replacing the TF-dependent originals:
  "${CMAKE_SOURCE_DIR}/patches/gfile_macos.cc"
  "${CMAKE_SOURCE_DIR}/patches/tfrecord_reader_macos.cc"
  "${CMAKE_SOURCE_DIR}/patches/tfrecord_writer_macos.cc"
)

# nucleus/core (only non-test .cc files)
set(NUCLEUS_CORE_SRCS
  "${CMAKE_SOURCE_DIR}/third_party/nucleus/core/status.cc"
)

# nucleus/util
set(NUCLEUS_UTIL_SRCS
  "${CMAKE_SOURCE_DIR}/third_party/nucleus/util/math.cc"
  "${CMAKE_SOURCE_DIR}/third_party/nucleus/util/utils.cc"
)

add_library(nucleus_io STATIC
  ${NUCLEUS_IO_SRCS}
  ${NUCLEUS_CORE_SRCS}
  ${NUCLEUS_UTIL_SRCS}
)

target_include_directories(nucleus_io PUBLIC
  "${CMAKE_SOURCE_DIR}"
  "${CMAKE_BINARY_DIR}/proto_gen"
  "${CMAKE_SOURCE_DIR}/cmake/tf_stubs"
  "${ABSL_PREFIX}/include"
  "${RE2_PREFIX}/include"
)
# Force-include the TF compat umbrella so CHECK / TF_PREDICT_FALSE are
# available without explicit #include in each nucleus .cc file.
target_compile_options(nucleus_io PRIVATE
  "-include${CMAKE_SOURCE_DIR}/cmake/tf_stubs/tensorflow/core/platform/tf_compat.h"
)

target_link_libraries(nucleus_io PUBLIC
  htslib::htslib
  re2::re2
  proto_nucleus
  proto_tf_example
  absl::strings
  absl::status
  absl::statusor
  absl::log
  absl::check
  absl::crc32c
  ZLIB::ZLIB
)

# Tests deferred to tests/native/ smoke tests (upstream tests depend on test_utils.h).
