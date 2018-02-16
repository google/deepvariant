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
    deps = [
        "//deepvariant/util/genomics:reference_proto",  # NO COPYBARA
        "//deepvariant/util/genomics:variants_proto",  # NO COPYBARA
    ],
)

cc_proto_library(
    name = "core_cc_pb2",
    srcs = ["core.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
    ],
)

py_proto_library(
    name = "core_py_pb2",
    srcs = ["core.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
    ],
)

pyclif_proto_library(
    name = "core_pyclif",
    proto_lib = ":core_proto",  # NO COPYBARA
)
