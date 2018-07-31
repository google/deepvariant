load("//:third_party/clif.bzl", "pyclif_proto_library")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

nucleus_py_library(
    name = "realigner",
    srcs = ["realigner.py"],
    deps = [
        ":aligner",
        ":window_selector",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/realigner/python:debruijn_graph",
        "//deepvariant/realigner/python:fast_pass_aligner",
        "//deepvariant/vendor:timer",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/flags",
    ],
)

py_test(
    name = "realigner_test",
    size = "small",
    srcs = ["realigner_test.py"],
    data = ["//deepvariant:testdata"],
    deps = [
        ":realigner",
        ":utils",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "cpp_window_selector",
    srcs = ["window_selector.cc"],
    hdrs = ["window_selector.h"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant/protos:realigner_cc_pb2",
    ],
)

nucleus_py_library(
    name = "window_selector",
    srcs = ["window_selector.py"],
    deps = [
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/realigner/python:window_selector",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:ranges",
    ],
)

py_test(
    name = "window_selector_test",
    size = "small",
    srcs = ["window_selector_test.py"],
    deps = [
        ":window_selector",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:realigner_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "debruijn_graph",
    srcs = ["debruijn_graph.cc"],
    hdrs = ["debruijn_graph.h"],
    deps = [
        "//deepvariant/protos:realigner_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

nucleus_py_library(
    name = "aligner",
    srcs = ["aligner.py"],
    deps = [
        ":utils",
        "//deepvariant/realigner/python:ssw",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "@io_abseil_py//absl/logging",
    ],
)

py_test(
    name = "aligner_test",
    size = "small",
    srcs = ["aligner_test.py"],
    deps = [
        ":aligner",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:realigner_py_pb2",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "//third_party/nucleus/util:cigar",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "utils",
    srcs = ["utils.py"],
    deps = [
        "//third_party/nucleus/protos:cigar_py_pb2",
    ],
)

cc_library(
    name = "ssw",
    srcs = ["ssw.cc"],
    hdrs = ["ssw.h"],
    deps = [
        "@libssw//:ssw_cpp",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "ssw_test",
    size = "small",
    srcs = ["ssw_test.cc"],
    deps = [
        ":ssw",
        "@libssw//:ssw_cpp",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_library(
    name = "fast_pass_aligner",
    srcs = [
        "fast_pass_aligner.cc",
    ],
    hdrs = [
        "fast_pass_aligner.h",
    ],
    deps = [
        ":ssw",
        "//deepvariant/protos:realigner_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@com_googlesource_code_re2//:re2",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "fast_pass_aligner_test",
    size = "small",
    srcs = ["fast_pass_aligner_test.cc"],
    data = ["//deepvariant:testdata"],
    deps = [
        "fast_pass_aligner",
        "//deepvariant/protos:realigner_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
        "@protobuf_archive//:protobuf",
    ],
)
