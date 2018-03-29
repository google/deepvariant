# Utility routines for working with genomics data and genomics protobufs.

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

# C++ header providing portabililty macros and functions.
cc_library(
    name = "port",
    hdrs = ["port.h"],
)

# redacted
# /util.
py_library(
    name = "core_py",
    srcs = [],
    deps = [
        ":io_utils",
        ":ranges",
        ":variant_utils",
        "//deepvariant/util/io:reference_fai",
        "//deepvariant/util/io:sam_reader",
        "//deepvariant/util/io:vcf_reader",
        "//deepvariant/util/io:vcf_writer",
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
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "cpp_utils_test",
    size = "small",
    srcs = ["utils_test.cc"],
    deps = [
        ":cpp_utils",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
        "//deepvariant/util/genomics:struct_cc_pb2",
        "//deepvariant/util/testing:cpp_test_utils",
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
        ":py_utils",
        ":ranges",
        "//deepvariant/util/testing:py_test_utils",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
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
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":ranges",
        "//deepvariant/util/genomics:position_py_pb2",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
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
    name = "variant_utils",
    srcs = ["variant_utils.py"],
    deps = [
        ":ranges",
        ":vcf_constants",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
    ],
)

py_test(
    name = "variant_utils_test",
    size = "small",
    srcs = ["variant_utils_test.py"],
    deps = [
        ":ranges",
        ":variant_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
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
    name = "io_utils",
    srcs = ["io_utils.py"],
)

py_test(
    name = "io_utils_test",
    size = "small",
    srcs = ["io_utils_test.py"],
    data = ["//deepvariant/util/testdata"],
    deps = [
        ":io_utils",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/testing:py_test_utils",
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
        "//deepvariant/util/testing:cpp_test_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
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
