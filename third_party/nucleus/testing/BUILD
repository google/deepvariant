# Testing utilities that are temporarily located here until they can
# be separately open source.

# Placeholder: load py_library
# Placeholder: load py_test

package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

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

py_library(
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
        "@com_google_protobuf//:protobuf",
        "@com_googlesource_code_re2//:re2",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
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
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
    ],
)

py_test(
    name = "py_test_utils_test",
    size = "small",
    srcs = ["test_utils_test.py"],
    main = "test_utils_test.py",
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_test_utils",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "@absl_py//absl/testing:absltest",
    ],
)

py_test(
    name = "protobuf_implementation_test",
    size = "small",
    srcs = ["protobuf_implementation_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        "//third_party/nucleus/io:sam",
        "@absl_py//absl/testing:absltest",
        "@com_google_protobuf//:protobuf_python",
    ],
)
