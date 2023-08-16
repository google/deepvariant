# Placeholder: load py_test
load("//third_party:clif.bzl", "py_clif_cc")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_clif_cc(
    name = "fast_pass_aligner",
    srcs = ["fast_pass_aligner.clif"],
    py_deps = [],
    pyclif_deps = [
        "//third_party/nucleus/protos:reads_pyclif",
        "//deepvariant/protos:realigner_pyclif",
    ],
    deps = ["//deepvariant/realigner:fast_pass_aligner"],
)

py_clif_cc(
    name = "debruijn_graph",
    srcs = ["debruijn_graph.clif"],
    py_deps = [],
    pyclif_deps = [
        "//deepvariant/protos:realigner_pyclif",
        "//third_party/nucleus/protos:reads_pyclif",
    ],
    deps = [
        "//deepvariant/realigner:debruijn_graph",
        "//third_party/nucleus/util:proto_clif_converter",
    ],
)

py_test(
    name = "debruijn_graph_wrap_test",
    size = "small",
    srcs = ["debruijn_graph_wrap_test.py"],
    data = ["//deepvariant:testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":debruijn_graph",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:realigner_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

# CLIF wrap for the SSW C++ interface.
py_clif_cc(
    name = "ssw",
    srcs = ["ssw.clif"],
    py_deps = ["//third_party/nucleus/io:clif_postproc"],
    deps = ["//deepvariant/realigner:ssw"],
)

py_test(
    name = "ssw_wrap_test",
    size = "small",
    srcs = ["ssw_wrap_test.py"],
    python_version = "PY3",
    deps = [
        ":ssw",
        "@absl_py//absl/testing:absltest",
    ],
)

py_test(
    name = "ssw_misc_test",
    size = "small",
    srcs = ["ssw_misc_test.py"],
    python_version = "PY3",
    deps = [
        ":ssw",
        "@absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "window_selector",
    srcs = ["window_selector.clif"],
    clif_deps = [
        "//deepvariant/python:allelecounter",
    ],
    py_deps = [],
    pyclif_deps = [
        "//deepvariant/protos:deepvariant_pyclif",
        "//deepvariant/protos:realigner_pyclif",
    ],
    deps = ["//deepvariant/realigner:cpp_window_selector"],
)
