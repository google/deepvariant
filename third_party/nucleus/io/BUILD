# Code for reading and writing genomics data.
# load("//devtools/python/blaze:python3.bzl", "py2and3_test")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

#
# Python API
#

nucleus_py_library(
    name = "io_py",
    deps = [
        ":bed",
        ":bedgraph",
        ":clif_postproc",
        ":fasta",
        ":fastq",
        ":genomics_reader",
        ":genomics_writer",
        ":gff",
        ":sam",
        ":vcf",
    ],
)

nucleus_py_library(
    name = "genomics_reader",
    srcs = ["genomics_reader.py"],
    deps = ["@absl_py//absl/logging"],
)

py_test(
    name = "genomics_reader_test",
    srcs = ["genomics_reader_test.py"],
    deps = [
        ":genomics_reader",
        "@absl_py//absl/testing:absltest",
    ],
)

nucleus_py_library(
    name = "genomics_writer",
    srcs = ["genomics_writer.py"],
    deps = ["@absl_py//absl/logging"],
)

nucleus_py_library(
    name = "vcf",
    srcs = ["vcf.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:vcf_reader",
        "//third_party/nucleus/io/python:vcf_writer",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:vcf_constants",
    ],
)

py_test(
    name = "vcf_test",
    size = "small",
    srcs = ["vcf_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":vcf",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "sam",
    srcs = ["sam.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:sam_reader",
        "//third_party/nucleus/io/python:sam_writer",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
    ],
)

py_test(
    name = "sam_test",
    size = "small",
    srcs = ["sam_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":sam",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "bed",
    srcs = ["bed.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:bed_reader",
        "//third_party/nucleus/io/python:bed_writer",
        "//third_party/nucleus/protos:bed_py_pb2",
    ],
)

py_test(
    name = "bed_test",
    size = "small",
    srcs = ["bed_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bed",
        "//third_party/nucleus/protos:bed_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "bedgraph",
    srcs = ["bedgraph.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:bedgraph_reader",
        "//third_party/nucleus/io/python:bedgraph_writer",
        "//third_party/nucleus/protos:bedgraph_py_pb2",
    ],
)

py_test(
    name = "bedgraph_test",
    size = "small",
    srcs = ["bedgraph_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bedgraph",
        "//third_party/nucleus/protos:bedgraph_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "fasta",
    srcs = ["fasta.py"],
    deps = [
        ":genomics_reader",
        "//third_party/nucleus/io/python:in_memory_fasta_reader",
        "//third_party/nucleus/io/python:indexed_fasta_reader",
        "//third_party/nucleus/io/python:unindexed_fasta_reader",
        "//third_party/nucleus/protos:fasta_py_pb2",
        "//third_party/nucleus/util:ranges",
    ],
)

py_test(
    name = "fasta_test",
    size = "small",
    srcs = ["fasta_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":fasta",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "fastq",
    srcs = ["fastq.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:fastq_reader",
        "//third_party/nucleus/io/python:fastq_writer",
        "//third_party/nucleus/protos:fastq_py_pb2",
    ],
)

py_test(
    name = "fastq_test",
    size = "small",
    srcs = ["fastq_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":fastq",
        "//third_party/nucleus/protos:fastq_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "gff",
    srcs = ["gff.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:gff_reader",
        "//third_party/nucleus/io/python:gff_writer",
        "//third_party/nucleus/protos:gff_py_pb2",
    ],
)

nucleus_py_library(
    name = "clif_postproc",
    srcs = ["clif_postproc.py"],
)

#
# C++ API
#

NUCLEUS_COPTS = [
    "-Wsign-compare",
]

cc_library(
    name = "io_cpp",
    deps = [
        ":bed_reader",
        ":bed_writer",
        ":bedgraph_reader",
        ":bedgraph_writer",
        ":fastq_reader",
        ":fastq_writer",
        ":gff_reader",
        ":gff_writer",
        ":hts_path",
        ":hts_verbose",
        ":in_memory_fasta_reader",
        ":indexed_fasta_reader",
        ":reader_base",
        ":reference",
        ":sam_reader",
        ":sam_writer",
        ":text_reader",
        ":text_writer",
        ":vcf_conversion",
        ":vcf_reader",
        ":vcf_writer",
    ],
)

test_suite(
    name = "smoke_tests",
    tests = ["hts_test"],
)

cc_library(
    name = "bed_reader",
    srcs = ["bed_reader.cc"],
    hdrs = ["bed_reader.h"],
    deps = [
        ":reader_base",
        ":text_reader",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:bed_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "bed_reader_test",
    size = "small",
    srcs = ["bed_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bed_reader",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "bed_writer",
    srcs = ["bed_writer.cc"],
    hdrs = ["bed_writer.h"],
    deps = [
        ":text_writer",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:bed_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "bed_writer_test",
    size = "small",
    srcs = ["bed_writer_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bed_writer",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "bedgraph_reader",
    srcs = ["bedgraph_reader.cc"],
    hdrs = ["bedgraph_reader.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":reader_base",
        ":text_reader",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:bedgraph_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "bedgraph_reader_test",
    size = "small",
    srcs = ["bedgraph_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bedgraph_reader",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "bedgraph_writer",
    srcs = ["bedgraph_writer.cc"],
    hdrs = ["bedgraph_writer.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":text_writer",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:bedgraph_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "bedgraph_writer_test",
    size = "small",
    srcs = ["bedgraph_writer_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bedgraph_reader",
        ":bedgraph_writer",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "fastq_reader",
    srcs = ["fastq_reader.cc"],
    hdrs = ["fastq_reader.h"],
    deps = [
        ":reader_base",
        ":text_reader",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:fastq_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "fastq_reader_test",
    size = "small",
    srcs = ["fastq_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":fastq_reader",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "fastq_writer",
    srcs = ["fastq_writer.cc"],
    hdrs = ["fastq_writer.h"],
    deps = [
        ":text_writer",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:fastq_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "fastq_writer_test",
    size = "small",
    srcs = ["fastq_writer_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":fastq_writer",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "sam_utils",
    srcs = ["sam_utils.cc"],
    hdrs = ["sam_utils.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":hts_path",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "@htslib",
    ],
)

cc_test(
    name = "sam_utils_test",
    size = "small",
    srcs = ["sam_utils_test.cc"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":sam_utils",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "@com_google_googletest//:gtest_main",
        "@htslib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "sam_reader",
    srcs = ["sam_reader.cc"],
    hdrs = ["sam_reader.h"],
    deps = [
        ":hts_path",
        ":reader_base",
        ":sam_utils",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:samplers",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "sam_reader_test",
    size = "small",
    srcs = ["sam_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":sam_reader",
        ":sam_writer",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "sam_writer",
    srcs = ["sam_writer.cc"],
    hdrs = ["sam_writer.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":hts_path",
        ":sam_utils",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "sam_writer_test",
    size = "small",
    srcs = ["sam_writer_test.cc"],
    copts = NUCLEUS_COPTS,
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":sam_reader",
        ":sam_writer",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/base:core_headers",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "vcf_reader",
    srcs = ["vcf_reader.cc"],
    hdrs = ["vcf_reader.h"],
    deps = [
        ":hts_path",
        ":reader_base",
        ":vcf_conversion",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "vcf_reader_test",
    size = "medium",
    srcs = ["vcf_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":vcf_reader",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
        "@protobuf_archive//:protobuf",
    ],
)

cc_library(
    name = "vcf_writer",
    srcs = ["vcf_writer.cc"],
    hdrs = ["vcf_writer.h"],
    deps = [
        ":hts_path",
        ":vcf_conversion",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "vcf_writer_test",
    size = "small",
    srcs = ["vcf_writer_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":vcf_writer",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
        "@protobuf_archive//:protobuf",
    ],
)

cc_library(
    name = "reference",
    srcs = ["reference.cc"],
    hdrs = ["reference.h"],
    deps = [
        ":reader_base",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/vendor:statusor",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

# Abstract value-parameterized tests:
# https://github.com/google/googletest/blob/master/googletest/docs/AdvancedGuide.md#creating-value-parameterized-abstract-tests
cc_library(
    name = "reference_test",
    testonly = True,
    srcs = ["reference_test.cc"],
    hdrs = ["reference_test.h"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":reference",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "indexed_fasta_reader",
    srcs = ["indexed_fasta_reader.cc"],
    hdrs = ["indexed_fasta_reader.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":hts_path",
        ":reader_base",
        ":reference",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:optional",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "indexed_fasta_reader_test",
    size = "small",
    srcs = ["indexed_fasta_reader_test.cc"],
    deps = [
        ":indexed_fasta_reader",
        ":reader_base",
        ":reference_test",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "unindexed_fasta_reader",
    srcs = ["unindexed_fasta_reader.cc"],
    hdrs = ["unindexed_fasta_reader.h"],
    deps = [
        ":reader_base",
        ":reference",
        ":text_reader",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:optional",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "unindexed_fasta_reader_test",
    size = "small",
    srcs = ["unindexed_fasta_reader_test.cc"],
    deps = [
        ":reference_test",
        ":unindexed_fasta_reader",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "in_memory_fasta_reader",
    srcs = ["in_memory_fasta_reader.cc"],
    hdrs = ["in_memory_fasta_reader.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":reader_base",
        ":reference",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "in_memory_fasta_reader_test",
    size = "small",
    srcs = ["in_memory_fasta_reader_test.cc"],
    deps = [
        ":in_memory_fasta_reader",
        ":reader_base",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "reader_base",
    srcs = ["reader_base.cc"],
    hdrs = ["reader_base.h"],
    deps = [
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/synchronization",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "reader_base_test",
    size = "small",
    srcs = ["reader_base_test.cc"],
    deps = [
        ":reader_base",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/vendor:status_matchers",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "vcf_conversion",
    srcs = ["vcf_conversion.cc"],
    hdrs = ["vcf_conversion.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_library(
    name = "hts_path",
    srcs = ["hts_path.cc"],
    hdrs = ["hts_path.h"],
    deps = [
        "//third_party/nucleus/platform:types",
        "@com_google_absl//absl/strings",
        "@htslib",
    ],
)

cc_library(
    name = "hts_verbose",
    srcs = ["hts_verbose.cc"],
    hdrs = ["hts_verbose.h"],
    deps = ["@htslib"],
)

cc_test(
    name = "hts_verbose_test",
    size = "small",
    srcs = ["hts_verbose_test.cc"],
    deps = [
        ":hts_verbose",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_test(
    name = "hts_test",
    size = "small",
    srcs = ["hts_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":hts_path",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:port",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@htslib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "text_writer",
    srcs = ["text_writer.cc"],
    hdrs = ["text_writer.h"],
    copts = NUCLEUS_COPTS,
    deps = [
        ":hts_path",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_library(
    name = "text_reader",
    srcs = ["text_reader.cc"],
    hdrs = ["text_reader.h"],
    deps = [
        ":hts_path",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "text_io_test",
    srcs = ["text_io_test.cc"],
    deps = [
        ":text_reader",
        ":text_writer",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

#
# Misc
#

py_test(
    name = "genomics_io_noplugin_test",
    size = "small",
    srcs = ["genomics_io_noplugin_test.py"],
    deps = [
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
    ],
)

cc_library(
    name = "gff_reader",
    srcs = ["gff_reader.cc"],
    hdrs = ["gff_reader.h"],
    deps = [
        ":reader_base",
        ":text_reader",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:gff_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:optional",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "gff_reader_test",
    srcs = ["gff_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":gff_reader",
        "//third_party/nucleus/protos:gff_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_library(
    name = "gff_writer",
    srcs = ["gff_writer.cc"],
    hdrs = ["gff_writer.h"],
    deps = [
        ":text_writer",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:gff_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "gff_writer_test",
    srcs = ["gff_writer_test.cc"],
    deps = [
        ":gff_writer",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:gff_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

py_test(
    name = "gff_test",
    srcs = ["gff_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":gff",
        "//third_party/nucleus/protos:gff_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

nucleus_py_binary(
    name = "converter",
    srcs = ["converter.py"],
    deps = [
        "//third_party/nucleus/io:io_py",
        "@absl_py//absl:app",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "converter_test",
    srcs = ["converter_test.py"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":converter",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)
