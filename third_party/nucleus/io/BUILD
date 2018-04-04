# Code for reading and writing genomics data.

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

py_library(
    name = "genomics_reader",
    srcs = ["genomics_reader.py"],
    deps = ["@com_google_absl_py//absl/logging"],
)

py_test(
    name = "genomics_reader_test",
    srcs = ["genomics_reader_test.py"],
    deps = [
        ":genomics_reader",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_library(
    name = "genomics_writer",
    srcs = ["genomics_writer.py"],
    deps = ["@com_google_absl_py//absl/logging"],
)

py_library(
    name = "vcf",
    srcs = ["vcf.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:vcf_reader",
        "//third_party/nucleus/io/python:vcf_writer",
        "//third_party/nucleus/protos:index_py_pb2",
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
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "sam",
    srcs = ["sam.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//third_party/nucleus/io/python:sam_reader",
        "//third_party/nucleus/protos:index_py_pb2",
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
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
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
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "fasta",
    srcs = ["fasta.py"],
    deps = [
        ":genomics_reader",
        "//third_party/nucleus/io/python:reference_fai",
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
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
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
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

#
# C++ API
#

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
        "//third_party/nucleus/protos:bed_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "//third_party/nucleus/vendor:zlib_compression_options",
        "//third_party/nucleus/vendor:zlib_inputstream",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "bed_reader_test",
    size = "small",
    srcs = ["bed_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bed_reader",
        "//third_party/nucleus/protos:bed_cc_pb2",
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
    name = "bed_writer",
    srcs = ["bed_writer.cc"],
    hdrs = ["bed_writer.h"],
    deps = [
        "//third_party/nucleus/protos:bed_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "//third_party/nucleus/vendor:zlib_compression_options",
        "//third_party/nucleus/vendor:zlib_outputbuffer",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "bed_writer_test",
    size = "small",
    srcs = ["bed_writer_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":bed_writer",
        "//third_party/nucleus/protos:bed_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
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
        "//third_party/nucleus/protos:fastq_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "//third_party/nucleus/vendor:zlib_compression_options",
        "//third_party/nucleus/vendor:zlib_inputstream",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "fastq_reader_test",
    size = "small",
    srcs = ["fastq_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":fastq_reader",
        "//third_party/nucleus/protos:fastq_cc_pb2",
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
    name = "fastq_writer",
    srcs = ["fastq_writer.cc"],
    hdrs = ["fastq_writer.h"],
    deps = [
        "//third_party/nucleus/protos:fastq_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "//third_party/nucleus/vendor:zlib_compression_options",
        "//third_party/nucleus/vendor:zlib_outputbuffer",
        "@com_google_absl//absl/memory",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "fastq_writer_test",
    size = "small",
    srcs = ["fastq_writer_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":fastq_writer",
        "//third_party/nucleus/protos:fastq_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
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
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:index_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:samplers",
        "//third_party/nucleus/vendor:statusor",
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
        "//third_party/nucleus/protos:index_cc_pb2",
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
    name = "vcf_reader",
    srcs = ["vcf_reader.cc"],
    hdrs = ["vcf_reader.h"],
    deps = [
        ":hts_path",
        ":reader_base",
        ":vcf_conversion",
        "//third_party/nucleus/protos:index_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "vcf_reader_test",
    size = "small",
    srcs = ["vcf_reader_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":vcf_reader",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "vcf_writer",
    srcs = ["vcf_writer.cc"],
    hdrs = ["vcf_writer.h"],
    deps = [
        ":hts_path",
        ":vcf_conversion",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "vcf_writer_test",
    size = "small",
    srcs = ["vcf_writer_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":vcf_writer",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "reference",
    srcs = ["reference.cc"],
    hdrs = ["reference.h"],
    deps = [
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
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "reference_fai",
    srcs = ["reference_fai.cc"],
    hdrs = ["reference_fai.h"],
    deps = [
        ":hts_path",
        ":reference",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "reference_fai_test",
    size = "small",
    srcs = ["reference_fai_test.cc"],
    deps = [
        ":reference_fai",
        ":reference_test",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "reader_base",
    srcs = ["reader_base.cc"],
    hdrs = ["reader_base.h"],
    deps = [
        "//third_party/nucleus/vendor:statusor",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "reader_base_test",
    size = "small",
    srcs = ["reader_base_test.cc"],
    deps = [
        ":reader_base",
        "//third_party/nucleus/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "vcf_conversion",
    srcs = ["vcf_conversion.cc"],
    hdrs = ["vcf_conversion.h"],
    deps = [
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

py_test(
    name = "genomics_io_noplugin_test",
    size = "small",
    srcs = ["genomics_io_noplugin_test.py"],
    deps = [
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/testing:py_test_utils",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

cc_library(
    name = "hts_path",
    srcs = ["hts_path.cc"],
    hdrs = ["hts_path.h"],
    deps = [
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
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
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:port",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@htslib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "clif_postproc",
    srcs = ["clif_postproc.py"],
)

# -----------------------------------------------------------
#
# These deps are purged by copybara and are not really part
# of the open-source version of genomics. They should
# probably be moved elsewhere.
#
# -----------------------------------------------------------

# -----------------------------------------------------------
# Done with non-OSS targets
# -----------------------------------------------------------
