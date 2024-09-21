# Placeholder: load py_test
load("@org_tensorflow//tensorflow:tensorflow.default.bzl", "pybind_extension", "pybind_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

pybind_extension(
    name = "fast_pass_aligner",
    srcs = ["fast_pass_aligner_pybind.cc"],
    deps = [
        "//deepvariant/realigner:fast_pass_aligner",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "fast_pass_aligner_cclib",
    srcs = ["fast_pass_aligner_pybind.cc"],
    deps = [
        "//deepvariant/realigner:fast_pass_aligner",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "debruijn_graph",
    srcs = ["debruijn_graph_pybind.cc"],
    deps = [
        "//deepvariant/realigner:debruijn_graph",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "debruijn_graph_cclib",
    srcs = ["debruijn_graph_pybind.cc"],
    deps = [
        "//deepvariant/realigner:debruijn_graph",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
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
pybind_extension(
    name = "ssw",
    srcs = ["ssw_pybind.cc"],
    deps = [
        "//deepvariant/realigner:ssw",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "ssw_cclib",
    srcs = ["ssw_pybind.cc"],
    deps = [
        "//deepvariant/realigner:ssw",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "ssw_wrap_test",
    size = "small",
    srcs = ["ssw_wrap_test.py"],
    python_version = "PY3",
    deps = [
        ":ssw",
        "//third_party/nucleus/io:clif_postproc",
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
        "//third_party/nucleus/io:clif_postproc",
        "@absl_py//absl/testing:absltest",
    ],
)

pybind_extension(
    name = "window_selector",
    srcs = ["window_selector_pybind.cc"],
    deps = [
        "//deepvariant/realigner:cpp_window_selector",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "window_selector_cclib",
    srcs = ["window_selector_pybind.cc"],
    deps = [
        "//deepvariant/realigner:cpp_window_selector",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)
