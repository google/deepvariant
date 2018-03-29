# Testing utilities that are temporarily located here until they can
# be separately open source.

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
        "@com_google_googletest//:gtest_main",
        "@protobuf_archive//:protobuf",
    ],
)

cc_library(
    name = "cpp_test_utils",
    testonly = True,
    srcs = ["test_utils.cc"],
    hdrs = ["test_utils.h"],
    deps = [
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/io:reader_base",
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
