load("//:third_party/clif.bzl", "py_clif_cc")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_clif_cc(
    name = "fast_pass_aligner",
    srcs = ["fast_pass_aligner.clif"],
    clif_deps = [
        "//third_party/nucleus/io/python:reference_fai",  # other py_clif_cc rules
    ],
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
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [
        "//deepvariant/protos:realigner_pyclif",
        "//third_party/nucleus/protos:reads_pyclif",
    ],
    deps = ["//deepvariant/realigner:debruijn_graph"],
)

py_test(
    name = "debruijn_graph_wrap_test",
    size = "small",
    srcs = ["debruijn_graph_wrap_test.py"],
    data = ["//deepvariant:testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":debruijn_graph",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:realigner_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
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
    deps = [
        ":ssw",
        "@io_abseil_py//absl/testing:absltest",
    ],
)

py_test(
    name = "ssw_misc_test",
    size = "small",
    srcs = ["ssw_misc_test.py"],
    deps = [
        ":ssw",
        "@io_abseil_py//absl/testing:absltest",
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
