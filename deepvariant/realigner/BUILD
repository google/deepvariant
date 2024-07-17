# Placeholder: load py_library
# Placeholder: load py_test

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_library(
    name = "realigner",
    srcs = ["realigner.py"],
    deps = [
        ":window_selector",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/realigner/python:debruijn_graph",
        "//deepvariant/realigner/python:fast_pass_aligner",
        "//deepvariant/vendor:timer",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:cigar",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/flags",
    ],
)

py_test(
    name = "realigner_test",
    size = "small",
    srcs = ["realigner_test.py"],
    data = ["//deepvariant:testdata"],
    python_version = "PY3",
    deps = [
        ":realigner",
        ":utils",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:realigner_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:cigar",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "cpp_window_selector",
    srcs = ["window_selector.cc"],
    hdrs = ["window_selector.h"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//deepvariant/protos:realigner_cc_pb2",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
    ],
)

py_library(
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
    python_version = "PY3",
    deps = [
        ":window_selector",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:realigner_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "debruijn_graph",
    srcs = ["debruijn_graph.cc"],
    hdrs = ["debruijn_graph.h"],
    deps = [
        "//deepvariant/protos:realigner_cc_pb2",
        "//third_party/nucleus/platform:types",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/container:node_hash_set",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
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
        "@com_google_absl//absl/strings:cord",
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
        "@com_google_absl//absl/strings:str_format",
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
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@com_googlesource_code_re2//:re2",
    ],
)

cc_test(
    name = "fast_pass_aligner_test",
    size = "small",
    srcs = ["fast_pass_aligner_test.cc"],
    data = ["//deepvariant:testdata"],
    deps = [
        "fast_pass_aligner",
        ":ssw",
        "//deepvariant/protos:realigner_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_protobuf//:protobuf",
        "@org_tensorflow//tensorflow/core:test",
    ],
)
