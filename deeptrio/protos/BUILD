load("@com_google_protobuf//:protobuf.bzl", "cc_proto_library", "py_proto_library")
load("//third_party:clif.bzl", "pyclif_proto_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

# ------------------------------------------------------------------------------
# deeptrio proto
# ------------------------------------------------------------------------------

proto_library(
    name = "deeptrio_proto",
    srcs = ["deeptrio.proto"],
    deps = [
        "//deepvariant/protos:deepvariant_proto",  # NO COPYBARA
        "//deepvariant/protos:realigner_proto",  # NO COPYBARA
        "//deepvariant/protos:resources_proto",  # NO COPYBARA
        "//third_party/nucleus/protos:reads_proto",  # NO COPYBARA
    ],
)

cc_proto_library(
    name = "deeptrio_cc_pb2",
    srcs = ["deeptrio.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [
        "//deepvariant/protos:realigner_cc_pb2",
        "//deepvariant/protos:resources_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

# buildozer: disable=LANG_proto_library-colocation
# buildozer: disable=LANG_proto_library-single-deps
py_proto_library(
    name = "deeptrio_py_pb2",
    srcs = ["deeptrio.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/protos:resources_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
    ],
)

pyclif_proto_library(
    name = "deeptrio_pyclif",
    proto_lib = ":deeptrio_proto",  # NO COPYBARA
)
