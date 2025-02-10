# Placeholder: load py_test
load("@org_tensorflow//tensorflow:tensorflow.default.bzl", "pybind_extension", "pybind_library")

package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

pybind_extension(
    name = "math",
    srcs = ["math_pybind.cc"],
    deps = [
        "//third_party/nucleus/util:cpp_math",
        "@pybind11",
    ],
)

pybind_library(
    name = "math_cclib",
    srcs = ["math_pybind.cc"],
    deps = [
        "//third_party/nucleus/util:cpp_math",
    ],
)

py_test(
    name = "math_wrap_test",
    size = "small",
    srcs = ["math_wrap_test.py"],
    data = [],
    deps = [
        ":math",
        "@absl_py//absl/testing:absltest",
    ],
)

pybind_library(
    name = "type_caster_nucleus_proto_ptr",
    hdrs = ["type_caster_nucleus_proto_ptr.h"],
    deps = [
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/log:absl_check",
        "@com_google_absl//absl/log:absl_log",
        "@com_google_protobuf//:proto_api",
        "@com_google_protobuf//:protobuf",
    ],
)

pybind_extension(
    name = "utils",
    srcs = ["utils_pybind.cc"],
    deps = [
        ":type_caster_nucleus_proto_ptr",
        "//third_party/nucleus/util:cpp_utils",
    ],
)

pybind_library(
    name = "utils_cclib",
    srcs = ["utils_pybind.cc"],
    deps = [
        ":type_caster_nucleus_proto_ptr",
        "//third_party/nucleus/util:cpp_utils",
    ],
)

py_test(
    name = "utils_wrap_test",
    size = "small",
    srcs = ["utils_wrap_test.py"],
    deps = [
        ":utils",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/protos:range_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@com_google_protobuf//:protobuf_python",
    ],
)
