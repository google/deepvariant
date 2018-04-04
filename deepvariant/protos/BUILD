load("//:third_party/clif.bzl", "py_clif_cc", "pyclif_proto_library")
load("@protobuf_archive//:protobuf.bzl", "py_proto_library")
load("@protobuf_archive//:protobuf.bzl", "cc_proto_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

proto_library(
    name = "deepvariant_proto",
    srcs = ["deepvariant.proto"],
    deps = [
        ":realigner_proto",  # NO COPYBARA
        "//third_party/nucleus/protos:position_proto",  # NO COPYBARA
        "//third_party/nucleus/protos:reads_proto",  # NO COPYBARA
        "//third_party/nucleus/protos:variants_proto",  # NO COPYBARA
    ],
)

cc_proto_library(
    name = "deepvariant_cc_pb2",
    srcs = ["deepvariant.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [
        ":realigner_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

py_proto_library(
    name = "deepvariant_py_pb2",
    srcs = ["deepvariant.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
    deps = [
        ":realigner_py_pb2",
        "//third_party/nucleus/protos:position_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
    ],
)

pyclif_proto_library(
    name = "deepvariant_pyclif",
    proto_lib = ":deepvariant_proto",  # NO COPYBARA
)

py_proto_library(
    name = "realigner_py_pb2",
    srcs = ["realigner.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
    deps = ["//third_party/nucleus/protos:range_py_pb2"],
)

cc_proto_library(
    name = "realigner_cc_pb2",
    srcs = ["realigner.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = ["//third_party/nucleus/protos:range_cc_pb2"],
)

proto_library(
    name = "realigner_proto",
    srcs = ["realigner.proto"],
    deps = ["//third_party/nucleus/protos:range_proto"],  # NO COPYBARA
)

pyclif_proto_library(
    name = "realigner_pyclif",
    proto_lib = ":realigner_proto",  # NO COPYBARA
)
