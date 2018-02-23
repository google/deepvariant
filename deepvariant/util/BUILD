# Utility routines for working with genomics data and genomics protobufs.

load("//:third_party/clif.bzl", "pyclif_proto_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])  # Apache 2.0 License

exports_files(["LICENSE"])

filegroup(
    name = "testdata",
    srcs = glob(["testdata/*"]),
)

test_suite(
    name = "smoke_tests",
    tests = ["hts_test"],
)

# C++ header providing portabililty macros and functions.
cc_library(
    name = "port",
    hdrs = ["port.h"],
)

py_library(
    name = "core_py",
    srcs = [],
    deps = [
        ":genomics_io",
        ":io_utils",
        ":ranges",
        ":variantutils",
        "//deepvariant/util/python:reference_fai",
        "//deepvariant/util/python:sam_reader",
        "//deepvariant/util/python:vcf_reader",
        "//deepvariant/util/python:vcf_writer",
    ],
)

cc_library(
    name = "sam_reader",
    srcs = ["sam_reader.cc"],
    hdrs = ["sam_reader.h"],
    deps = [
        ":cpp_utils",
        ":hts_path",
        ":reader_base",
        ":samplers",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:position_cc_pb2",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
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
    data = [":testdata"],
    deps = [
        ":cpp_test_utils",
        ":cpp_utils",
        ":sam_reader",
        "//deepvariant/util/testing:gunit_extras",
        "//deepvariant/util/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "vcf_conversion",
    srcs = ["vcf_conversion.cc"],
    hdrs = ["vcf_conversion.h"],
    deps = [
        ":cpp_math",
        ":cpp_utils",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "@com_google_absl//absl/memory",
        "@htslib",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_library(
    name = "vcf_reader",
    srcs = ["vcf_reader.cc"],
    hdrs = ["vcf_reader.h"],
    deps = [
        ":cpp_math",
        ":cpp_utils",
        ":hts_path",
        ":reader_base",
        ":vcf_conversion",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
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
    data = [":testdata"],
    deps = [
        ":cpp_test_utils",
        ":cpp_utils",
        ":vcf_reader",
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
        ":cpp_math",
        ":cpp_utils",
        ":hts_path",
        ":vcf_conversion",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/genomics:struct_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
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
    data = [":testdata"],
    deps = [
        ":cpp_test_utils",
        ":cpp_utils",
        ":vcf_writer",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
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
        "//deepvariant/util/protos:core_cc_pb2",
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
    data = [
        "testdata/test.fasta",
        "testdata/test.fasta.fai",
    ],
    deps = [
        ":cpp_test_utils",
        ":cpp_utils",
        ":reference",
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
        ":cpp_utils",
        ":hts_path",
        ":reference",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
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
        ":cpp_test_utils",
        ":cpp_utils",
        ":reference_fai",
        ":reference_test",
        "//deepvariant/util/vendor:status_matchers",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "cpp_utils",
    srcs = ["utils.cc"],
    hdrs = ["utils.h"],
    deps = [
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:position_cc_pb2",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/genomics:struct_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "cpp_utils_test",
    size = "small",
    srcs = ["utils_test.cc"],
    deps = [
        ":cpp_test_utils",
        ":cpp_utils",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:struct_cc_pb2",
        "//deepvariant/util/testing:gunit_extras",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "py_utils",
    srcs = ["utils.py"],
    deps = [
        "//deepvariant/util:cigar",
        "//deepvariant/util:ranges",
    ],
)

py_library(
    name = "cigar",
    srcs = ["cigar.py"],
    deps = [
        "//deepvariant/util/genomics:cigar_py_pb2",
    ],
)

py_test(
    name = "cigar_test",
    size = "small",
    srcs = ["cigar_test.py"],
    deps = [
        ":cigar",
        "//deepvariant/util/genomics:cigar_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "py_utils_test",
    size = "small",
    srcs = ["utils_test.py"],
    main = "utils_test.py",
    deps = [
        ":py_test_utils",
        ":py_utils",
        ":ranges",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "cpp_test_utils",
    testonly = True,
    srcs = ["test_utils.cc"],
    hdrs = ["test_utils.h"],
    deps = [
        ":cpp_utils",
        ":reader_base",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
        "//deepvariant/util/vendor:statusor",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "py_test_utils",
    testonly = True,
    srcs = ["test_utils.py"],
    deps = [
        "//deepvariant/util:cigar",
        "//deepvariant/util/genomics:position_py_pb2",
        "//deepvariant/util/genomics:reads_py_pb2",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/flags",
    ],
)

py_library(
    name = "proto_utils",
    srcs = ["proto_utils.py"],
    deps = [
        "@protobuf_archive//:protobuf_python",
    ],
)

cc_library(
    name = "cpp_math",
    srcs = ["math.cc"],
    hdrs = ["math.h"],
    deps = [
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "cpp_math_test",
    size = "small",
    srcs = ["math_test.cc"],
    deps = [
        ":cpp_math",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "genomics_math",
    srcs = ["genomics_math.py"],
    deps = ["//deepvariant/util/python:math"],
)

py_test(
    name = "genomics_math_test",
    size = "small",
    srcs = ["genomics_math_test.py"],
    deps = [
        ":genomics_math",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "ranges",
    srcs = ["ranges.py"],
    deps = [
        "//deepvariant/util/genomics:position_py_pb2",
        "//deepvariant/util/genomics:range_py_pb2",
    ],
)

py_test(
    name = "ranges_test",
    size = "small",
    srcs = ["ranges_test.py"],
    data = [":testdata"],
    deps = [
        ":py_test_utils",
        ":ranges",
        "//deepvariant/util/genomics:position_py_pb2",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/protos:core_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "variantcall_utils",
    srcs = ["variantcall_utils.py"],
    deps = [
        ":vcf_constants",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
    ],
)

py_test(
    name = "variantcall_utils_test",
    size = "small",
    srcs = ["variantcall_utils_test.py"],
    deps = [
        ":variantcall_utils",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "variantutils",
    srcs = ["variantutils.py"],
    deps = [
        ":ranges",
        ":vcf_constants",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
    ],
)

py_test(
    name = "variantutils_test",
    size = "small",
    srcs = ["variantutils_test.py"],
    deps = [
        ":py_test_utils",
        ":ranges",
        ":variantutils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "struct_utils",
    srcs = ["struct_utils.py"],
    deps = [
        "//deepvariant/util/genomics:struct_py_pb2",
    ],
)

py_test(
    name = "struct_utils_test",
    srcs = ["struct_utils_test.py"],
    deps = [
        ":struct_utils",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "vcf_constants",
    srcs = ["vcf_constants.py"],
    deps = [
        ":struct_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
    ],
)

py_test(
    name = "vcf_constants_test",
    srcs = ["vcf_constants_test.py"],
    deps = [
        ":struct_utils",
        ":vcf_constants",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "genomics_io",
    srcs = ["genomics_io.py"],
    deps = [
        "//deepvariant/util/protos:core_py_pb2",
        "//deepvariant/util/python:reference_fai",
        "//deepvariant/util/python:sam_reader",
        "//deepvariant/util/python:vcf_reader",
        "//deepvariant/util/python:vcf_writer",
    ],
)

py_test(
    name = "genomics_io_test",
    size = "small",
    srcs = ["genomics_io_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    deps = [
        ":genomics_io",
        ":py_test_utils",
        ":ranges",
        "//deepvariant/util/genomics:struct_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "genomics_io_noplugin_test",
    size = "small",
    srcs = ["genomics_io_noplugin_test.py"],
    deps = [
        ":genomics_io",
        ":py_test_utils",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_test(
    name = "genomics_io_gcs_test",
    size = "small",
    srcs = ["genomics_io_gcs_test.py"],
    # GCS access doesn't work in the sanitized environment in our test system.
    tags = [
        "local",
    ],
    deps = [
        ":genomics_io",
        ":htslib_gcp_oauth_lib",
        ":ranges",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_library(
    name = "io_utils",
    srcs = ["io_utils.py"],
)

py_test(
    name = "io_utils_test",
    size = "small",
    srcs = ["io_utils_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    deps = [
        ":io_utils",
        ":py_test_utils",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/protos:core_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "errors",
    srcs = ["errors.py"],
    deps = [
        "@com_google_absl_py//absl/logging",
    ],
)

py_test(
    name = "errors_test",
    srcs = ["errors_test.py"],
    deps = [
        ":errors",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "samplers",
    hdrs = ["samplers.h"],
    deps = ["@org_tensorflow//tensorflow/core:lib"],
)

cc_test(
    name = "samplers_test",
    size = "small",
    srcs = ["samplers_test.cc"],
    deps = [
        ":samplers",
        "//deepvariant/util:cpp_test_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
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

cc_library(
    name = "tf_str",
    srcs = ["tf_str.cc"],
    deps = [
        "@org_tensorflow//tensorflow/core:lib",
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

py_binary(
    name = "cloud_utils",
    srcs = ["cloud_utils.py"],
    deps = ["@protobuf_archive//:protobuf_python"],
)

py_test(
    name = "cloud_utils_test",
    srcs = ["cloud_utils_test.py"],
    deps = [
        ":cloud_utils",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_library(
    name = "htslib_gcp_oauth_lib",
    srcs = ["htslib_gcp_oauth.py"],
    deps = [
        ":cloud_utils",
        "@com_google_absl_py//absl:app",
    ],
)

py_binary(
    name = "htslib_gcp_oauth",
    srcs = ["htslib_gcp_oauth.py"],
    main = "htslib_gcp_oauth.py",
    deps = [
        ":cloud_utils",
        "@com_google_absl_py//absl:app",
    ],
)

cc_test(
    name = "hts_test",
    size = "small",
    srcs = ["hts_test.cc"],
    data = [":testdata"],
    deps = [
        ":hts_path",
        ":port",
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
