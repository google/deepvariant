load("@com_google_protobuf//:protobuf.bzl", "cc_proto_library", "py_proto_library")
load("//third_party:clif.bzl", "pyclif_proto_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

# ------------------------------------------------------------------------------
# deepvariant proto
# ------------------------------------------------------------------------------

proto_library(
    name = "deepvariant_proto",
    srcs = ["deepvariant.proto"],
    deps = [
        ":realigner_proto",  # NO COPYBARA
        ":resources_proto",  # NO COPYBARA
        "//third_party/nucleus/protos:position_proto",  # NO COPYBARA
        "//third_party/nucleus/protos:reads_proto",  # NO COPYBARA
        "//third_party/nucleus/protos:variants_proto",  # NO COPYBARA
    ],
)

cc_proto_library(
    name = "deepvariant_cc_pb2",
    srcs = ["deepvariant.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [
        ":realigner_cc_pb2",
        ":resources_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

py_proto_library(
    name = "deepvariant_py_pb2",
    srcs = ["deepvariant.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [
        ":realigner_py_pb2",
        ":resources_py_pb2",
        "//third_party/nucleus/protos:position_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
    ],
)

pyclif_proto_library(
    name = "deepvariant_pyclif",
    proto_lib = ":deepvariant_proto",  # NO COPYBARA
)

# ------------------------------------------------------------------------------
# realigner proto
# ------------------------------------------------------------------------------

py_proto_library(
    name = "realigner_py_pb2",
    srcs = ["realigner.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = ["//third_party/nucleus/protos:range_py_pb2"],
)

cc_proto_library(
    name = "realigner_cc_pb2",
    srcs = ["realigner.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
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

# ------------------------------------------------------------------------------
# resources proto
# ------------------------------------------------------------------------------

py_proto_library(
    name = "resources_py_pb2",
    srcs = ["resources.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

cc_proto_library(
    name = "resources_cc_pb2",
    srcs = ["resources.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

proto_library(
    name = "resources_proto",
    srcs = ["resources.proto"],
)
