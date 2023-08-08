# Utility routines for working with genomics data and genomics protobufs.

# Placeholder: load py_library
# Placeholder: load py_test

package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

py_library(
    name = "util_py",
    deps = [
        ":cigar",
        ":errors",
        ":genomics_math",
        ":proto_utils",
        ":py_utils",
        ":ranges",
        ":sequence_utils",
        ":struct_utils",
        ":variant_utils",
        ":variantcall_utils",
        ":vcf_constants",
        ":vis",
    ],
)

cc_library(
    name = "util_cpp",
    deps = [
        ":cpp_math",
        ":cpp_utils",
        ":port",
        ":samplers",
    ],
)

# C++ header providing portabililty macros and functions.
cc_library(
    name = "port",
    hdrs = ["port.h"],
)

cc_library(
    name = "cpp_utils",
    srcs = ["utils.cc"],
    hdrs = ["utils.h"],
    deps = [
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "cpp_utils_test",
    size = "small",
    srcs = ["utils_test.cc"],
    deps = [
        ":cpp_utils",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "py_utils",
    srcs = ["utils.py"],
    deps = ["//third_party/nucleus/util/python:utils"],
)

py_library(
    name = "cigar",
    srcs = ["cigar.py"],
    deps = ["//third_party/nucleus/protos:cigar_py_pb2"],
)

py_test(
    name = "cigar_test",
    size = "small",
    srcs = ["cigar_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":cigar",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "py_utils_test",
    size = "small",
    srcs = ["utils_test.py"],
    main = "utils_test.py",
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_utils",
        ":ranges",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "proto_utils",
    srcs = ["proto_utils.py"],
    deps = [
        "@com_google_protobuf//:protobuf_python",
    ],
)

cc_library(
    name = "cpp_math",
    srcs = ["math.cc"],
    hdrs = ["math.h"],
    deps = [
        "@com_google_absl//absl/strings:cord",
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
    deps = ["//third_party/nucleus/util/python:math"],
)

py_test(
    name = "genomics_math_test",
    size = "small",
    srcs = ["genomics_math_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":genomics_math",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "ranges",
    srcs = ["ranges.py"],
    deps = [
        "//third_party/nucleus/io:bed",
        "//third_party/nucleus/protos:position_py_pb2",
        "//third_party/nucleus/protos:range_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "ranges_test",
    size = "small",
    srcs = ["ranges_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":ranges",
        "//third_party/nucleus/protos:position_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "variantcall_utils",
    srcs = ["variantcall_utils.py"],
    deps = [
        ":vcf_constants",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
    ],
)

py_test(
    name = "variantcall_utils_test",
    size = "small",
    srcs = ["variantcall_utils_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":variantcall_utils",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "variant_utils",
    srcs = ["variant_utils.py"],
    deps = [
        ":ranges",
        ":variantcall_utils",
        ":vcf_constants",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
    ],
)

py_test(
    name = "variant_utils_test",
    size = "small",
    srcs = ["variant_utils_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":ranges",
        ":variant_utils",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "struct_utils",
    srcs = ["struct_utils.py"],
    deps = [
        "//third_party/nucleus/protos:struct_py_pb2",
    ],
)

py_test(
    name = "struct_utils_test",
    srcs = ["struct_utils_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":struct_utils",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "vcf_constants",
    srcs = ["vcf_constants.py"],
    deps = [
        ":struct_utils",
        "//third_party/nucleus/protos:variants_py_pb2",
    ],
)

py_test(
    name = "vcf_constants_test",
    srcs = ["vcf_constants_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":struct_utils",
        ":vcf_constants",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "errors",
    srcs = ["errors.py"],
    deps = [
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "errors_test",
    srcs = ["errors_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":errors",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "samplers",
    hdrs = ["samplers.h"],
    deps = [
        "//third_party/nucleus/platform:types",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "samplers_test",
    size = "small",
    srcs = ["samplers_test.cc"],
    deps = [
        ":samplers",
        "//third_party/nucleus/testing:cpp_test_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "sequence_utils",
    srcs = ["sequence_utils.py"],
)

py_test(
    name = "sequence_utils_test",
    srcs = ["sequence_utils_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":sequence_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "proto_ptr",
    hdrs = [
        "proto_ptr.h",
    ],
)

cc_library(
    name = "proto_clif_converter",
    srcs = ["proto_clif_converter.cc"],
    hdrs = ["proto_clif_converter.h"],
    deps = [
        ":proto_ptr",
        "@clif//:cpp_runtime",
        "@com_google_protobuf//:proto_api",
        "@com_google_protobuf//:protobuf",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

py_library(
    name = "vis",
    srcs = ["vis.py"],
    deps = ["//third_party/nucleus/protos:variants_py_pb2"],
)

py_test(
    name = "vis_test",
    size = "small",
    srcs = ["vis_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":vis",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)
