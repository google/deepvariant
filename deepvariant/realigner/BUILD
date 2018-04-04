load("//:third_party/clif.bzl", "pyclif_proto_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_library(
    name = "realigner",
    srcs = ["realigner.py"],
    deps = [
        ":aligner",
        ":window_selector",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/realigner/python:debruijn_graph",
        "//deepvariant/vendor:timer",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
        "@com_google_absl_py//absl/flags",
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
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:ranges",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "window_selector",
    srcs = ["window_selector.py"],
    deps = [
        ":utils",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "//third_party/nucleus/protos:range_py_pb2",
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
        "//third_party/nucleus/protos:range_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
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
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

py_library(
    name = "aligner",
    srcs = ["aligner.py"],
    deps = [
        ":utils",
        "//deepvariant/realigner/python:ssw",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "@com_google_absl_py//absl/logging",
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
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
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
