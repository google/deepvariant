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
        "//deepvariant/util:ranges",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:vcf_constants",
        "//deepvariant/util/genomics:index_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/io/python:vcf_reader",
        "//deepvariant/util/io/python:vcf_writer",
    ],
)

py_test(
    name = "vcf_test",
    size = "small",
    srcs = ["vcf_test.py"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":vcf",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
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
        "//deepvariant/util/io/python:sam_reader",
        "//deepvariant/util/genomics:index_py_pb2",
        "//deepvariant/util/genomics:reads_py_pb2",
        # "//third_party/py/tfbam_lib",
    ],
)

py_test(
    name = "sam_test",
    size = "small",
    srcs = ["sam_test.py"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":sam",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
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
        "//deepvariant/util/genomics:bed_py_pb2",
        "//deepvariant/util/io/python:bed_reader",
        "//deepvariant/util/io/python:bed_writer",
    ],
)

py_test(
    name = "bed_test",
    size = "small",
    srcs = ["bed_test.py"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":bed",
        "//deepvariant/util/genomics:bed_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "fasta",
    srcs = ["fasta.py"],
    deps = [
        ":genomics_reader",
        "//deepvariant/util/io/python:reference_fai",
    ],
)

py_test(
    name = "fasta_test",
    size = "small",
    srcs = ["fasta_test.py"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":fasta",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
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
        "//deepvariant/util/genomics:fastq_py_pb2",
        "//deepvariant/util/io/python:fastq_reader",
        "//deepvariant/util/io/python:fastq_writer",
    ],
)

py_test(
    name = "fastq_test",
    size = "small",
    srcs = ["fastq_test.py"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":fastq",
        "//deepvariant/util/genomics:fastq_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:bed_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "//deepvariant/util/vendor:zlib_compression_options",
        "//deepvariant/util/vendor:zlib_inputstream",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "bed_reader_test",
    size = "small",
    srcs = ["bed_reader_test.cc"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":bed_reader",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:bed_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/testing:gunit_extras",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:bed_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "//deepvariant/util/vendor:zlib_compression_options",
        "//deepvariant/util/vendor:zlib_outputbuffer",
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
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":bed_writer",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:bed_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:fastq_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "//deepvariant/util/vendor:zlib_compression_options",
        "//deepvariant/util/vendor:zlib_inputstream",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "fastq_reader_test",
    size = "small",
    srcs = ["fastq_reader_test.cc"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":fastq_reader",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:fastq_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/testing:gunit_extras",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:fastq_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "//deepvariant/util/vendor:zlib_compression_options",
        "//deepvariant/util/vendor:zlib_outputbuffer",
        "@com_google_absl//absl/memory",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "fastq_writer_test",
    size = "small",
    srcs = ["fastq_writer_test.cc"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":fastq_writer",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:fastq_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util:samplers",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:index_cc_pb2",
        "//deepvariant/util/genomics:position_cc_pb2",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "sam_reader_test",
    size = "small",
    srcs = ["sam_reader_test.cc"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":sam_reader",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:index_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/testing:gunit_extras",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util:cpp_math",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:index_cc_pb2",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "vcf_reader_test",
    size = "small",
    srcs = ["vcf_reader_test.cc"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":vcf_reader",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/testing:gunit_extras",
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
        "//deepvariant/util:cpp_math",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/genomics:struct_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "vcf_writer_test",
    size = "small",
    srcs = ["vcf_writer_test.cc"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":vcf_writer",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/vendor:statusor",
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
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":reference",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/vendor:statusor",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/testing:cpp_test_utils",
        "//deepvariant/util/vendor:status_matchers",
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
        "//deepvariant/util/vendor:statusor",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "reader_base_test",
    size = "small",
    srcs = ["reader_base_test.cc"],
    deps = [
        ":reader_base",
        "//deepvariant/util/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "vcf_conversion",
    srcs = ["vcf_conversion.cc"],
    hdrs = ["vcf_conversion.h"],
    deps = [
        "//deepvariant/util:cpp_math",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/vendor:statusor",
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
        "//deepvariant/util/io:sam",
        "//deepvariant/util/testing:py_test_utils",
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
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":hts_path",
        "//deepvariant/util:port",
        "//deepvariant/util/testing:cpp_test_utils",
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
