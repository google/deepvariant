# Testing utilities that are temporarily located here until they can
# be separately open source.
load("//devtools/python/blaze:python3.bzl", "py2and3_test")

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
    name = "gunit_extras",
    testonly = True,
    srcs = [
        "protocol-buffer-matchers.cc",
    ],
    hdrs = [
        "protocol-buffer-matchers.h",
    ],
    deps = [
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
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "py_test_utils",
    testonly = True,
    srcs = ["test_utils.py"],
    deps = [
        "//third_party/nucleus/protos:position_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:cigar",
        "@com_google_absl_py//absl/flags",
    ],
)

py2and3_test(
    name = "py_test_utils_test",
    size = "small",
    srcs = ["test_utils_test.py"],
    main = "test_utils_test.py",
    deps = [
        ":py_test_utils",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py2and3_test(
    name = "protobuf_implementation_test",
    size = "small",
    srcs = ["protobuf_implementation_test.py"],
    deps = [
        "//net/proto2/python/public",
        "//third_party/nucleus/io:sam",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)
