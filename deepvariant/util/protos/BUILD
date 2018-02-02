load("//:third_party/clif.bzl", "py_clif_cc", "pyclif_proto_library")
load("@protobuf_archive//:protobuf.bzl", "py_proto_library")
load("@protobuf_archive//:protobuf.bzl", "cc_proto_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

proto_library(
    name = "core_proto",
    srcs = ["core.proto"],
)

cc_proto_library(
    name = "core_cc_pb2",
    srcs = ["core.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [],
)

py_proto_library(
    name = "core_py_pb2",
    srcs = ["core.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
)

pyclif_proto_library(
    name = "core_pyclif",
    proto_lib = ":core_proto",
)
