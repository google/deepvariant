# Testing utilities that are temporarily located here until they can
# be separately open source.
# load("//devtools/python/blaze:python3.bzl", "py2and3_test")

load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(default_visibility = ["//visibility:public"])

licenses(["notice"])

exports_files(["LICENSE"])

test_suite(
    name = "smoke_tests",
    tests = [
        "gunit_extras_test",
    ],
)

cc_library(
    name = "testing_cpp",
    testonly = True,
    deps = [
        ":cpp_test_utils",
        ":gunit_extras",
    ],
)

nucleus_py_library(
    name = "testing_py",
    testonly = True,
    deps = [
        ":py_test_utils",
    ],
)

cc_library(
    name = "gunit_extras",
    testonly = True,
    srcs = [
        "protocol-buffer-matchers.cc",
    ],
    hdrs = [
        "protocol-buffer-matchers.h",
    ],
    deps = [
        "//third_party/nucleus/platform:types",
        "@com_google_absl//absl/strings",
        "@com_googlesource_code_re2//:re2",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "gunit_extras_test",
    size = "small",
    srcs = ["gunit_extras_test.cc"],
    deps = [
        ":gunit_extras",
        "//net/proto2/proto:empty_cc_pb2",
        "@com_google_googletest//:gtest_main",
    ],
)

cc_library(
    name = "cpp_test_utils",
    testonly = True,
    srcs = ["test_utils.cc"],
    hdrs = ["test_utils.h"],
    deps = [
        "//third_party/nucleus/io:reader_base",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

nucleus_py_library(
    name = "py_test_utils",
    testonly = True,
    srcs = ["test_utils.py"],
    deps = [
        "//third_party/nucleus/protos:position_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:cigar",
        "@io_abseil_py//absl/flags",
    ],
)

py_test(
    name = "py_test_utils_test",
    size = "small",
    srcs = ["test_utils_test.py"],
    main = "test_utils_test.py",
    deps = [
        ":py_test_utils",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "@io_abseil_py//absl/testing:absltest",
    ],
)

py_test(
    name = "protobuf_implementation_test",
    size = "small",
    srcs = ["protobuf_implementation_test.py"],
    deps = [
        "//third_party/nucleus/io:sam",
        "@io_abseil_py//absl/testing:absltest",
        "@protobuf_archive//:protobuf_python",
    ],
)
