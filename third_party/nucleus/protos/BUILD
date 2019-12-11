load("//third_party:clif.bzl", "pyclif_proto_library")
load("@com_google_protobuf//:protobuf.bzl", "cc_proto_library", "py_proto_library")

package(default_visibility = [
    "//visibility:public",
])

licenses(["notice"])

exports_files(["LICENSE"])

# ----------------------------------------------------------------------
# proto_library's
# ----------------------------------------------------------------------

proto_library(
    name = "bed_proto",
    srcs = [
        "bed.proto",
    ],
)

proto_library(
    name = "bedgraph_proto",
    srcs = [
        "bedgraph.proto",
    ],
)

proto_library(
    name = "cigar_proto",
    srcs = ["cigar.proto"],
)

proto_library(
    name = "fasta_proto",
    srcs = [
        "fasta.proto",
    ],
    deps = [
        ":range_proto",  # NO COPYBARA
        ":reference_proto",  # NO COPYBARA
    ],
)

proto_library(
    name = "fastq_proto",
    srcs = [
        "fastq.proto",
    ],
)

proto_library(
    name = "gff_proto",
    srcs = [
        "gff.proto",
    ],
    deps = [":range_proto"],  # NO COPYBARA
)

proto_library(
    name = "position_proto",
    srcs = ["position.proto"],
)

proto_library(
    name = "range_proto",
    srcs = ["range.proto"],
)

proto_library(
    name = "reads_proto",
    srcs = [
        "reads.proto",
    ],
    deps = [
        ":cigar_proto",  # NO COPYBARA
        ":position_proto",  # NO COPYBARA
        ":reference_proto",  # NO COPYBARA
        ":struct_proto",  # NO COPYBARA
    ],
)

proto_library(
    name = "reference_proto",
    srcs = [
        "reference.proto",
    ],
    deps = [
        ":range_proto",  # NO COPYBARA
    ],
)

proto_library(
    name = "struct_proto",
    srcs = ["struct.proto"],
)

proto_library(
    name = "variants_proto",
    srcs = [
        "variants.proto",
    ],
    deps = [
        ":reference_proto",  # NO COPYBARA
        ":struct_proto",  # NO COPYBARA
    ],
)

proto_library(
    name = "feature_proto",
    srcs = ["feature.proto"],
)

proto_library(
    name = "example_proto",
    srcs = ["feature.proto"],
    deps = [
        ":feature_proto",  # NO COPYBARA
    ],
)

# ----------------------------------------------------------------------
# cc_proto_library's
# ----------------------------------------------------------------------

# redacted
# boilerplate in these cc_proto_library rules and in this file in general.

cc_proto_library(
    name = "bed_cc_pb2",
    srcs = ["bed.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "bedgraph_cc_pb2",
    srcs = ["bedgraph.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "cigar_cc_pb2",
    srcs = ["cigar.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "fasta_cc_pb2",
    srcs = ["fasta.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [
        ":range_cc_pb2",
        ":reference_cc_pb2",
    ],
)

cc_proto_library(
    name = "fastq_cc_pb2",
    srcs = ["fastq.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "gff_cc_pb2",
    srcs = ["gff.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [":range_cc_pb2"],
)

cc_proto_library(
    name = "position_cc_pb2",
    srcs = ["position.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "range_cc_pb2",
    srcs = ["range.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "reads_cc_pb2",
    srcs = ["reads.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [
        ":cigar_cc_pb2",
        ":position_cc_pb2",
        ":reference_cc_pb2",
        ":struct_cc_pb2",
    ],
)

cc_proto_library(
    name = "reference_cc_pb2",
    srcs = ["reference.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [
        ":range_cc_pb2",
    ],
)

cc_proto_library(
    name = "struct_cc_pb2",
    srcs = ["struct.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "variants_cc_pb2",
    srcs = ["variants.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [
        ":reference_cc_pb2",
        ":struct_cc_pb2",
    ],
)

cc_proto_library(
    name = "feature_cc_pb2",
    srcs = ["feature.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
)

cc_proto_library(
    name = "example_cc_pb2",
    srcs = ["example.proto"],
    default_runtime = "@com_google_protobuf//:protobuf",
    protoc = "@com_google_protobuf//:protoc",
    deps = [
        ":feature_cc_pb2",
    ],
)

cc_library(
    name = "all_nucleus_protos_cc",
    deps = [
        ":bed_cc_pb2",
        ":bedgraph_cc_pb2",
        ":cigar_cc_pb2",
        ":fasta_cc_pb2",
        ":fastq_cc_pb2",
        ":gff_cc_pb2",
        ":position_cc_pb2",
        ":range_cc_pb2",
        ":reads_cc_pb2",
        ":reference_cc_pb2",
        ":struct_cc_pb2",
        ":variants_cc_pb2",
    ],
)

# ----------------------------------------------------------------------
# py_proto_library's
# ----------------------------------------------------------------------

py_proto_library(
    name = "bed_py_pb2",
    srcs = ["bed.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "bedgraph_py_pb2",
    srcs = ["bedgraph.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "cigar_py_pb2",
    srcs = ["cigar.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "gff_py_pb2",
    srcs = ["gff.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [":range_py_pb2"],
)

py_proto_library(
    name = "fasta_py_pb2",
    srcs = ["fasta.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [
        ":range_py_pb2",
        ":reference_py_pb2",
    ],
)

py_proto_library(
    name = "fastq_py_pb2",
    srcs = ["fastq.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "position_py_pb2",
    srcs = ["position.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "range_py_pb2",
    srcs = ["range.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "reads_py_pb2",
    srcs = ["reads.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [
        ":cigar_py_pb2",
        ":position_py_pb2",
        ":reference_py_pb2",
        ":struct_py_pb2",
    ],
)

py_proto_library(
    name = "reference_py_pb2",
    srcs = ["reference.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [
        ":range_py_pb2",
    ],
)

py_proto_library(
    name = "struct_py_pb2",
    srcs = ["struct.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "variants_py_pb2",
    srcs = ["variants.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [
        ":reference_py_pb2",
        ":struct_py_pb2",
    ],
)

py_proto_library(
    name = "feature_py_pb2",
    srcs = ["feature.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
)

py_proto_library(
    name = "example_py_pb2",
    srcs = ["example.proto"],
    default_runtime = "@com_google_protobuf//:protobuf_python",
    protoc = "@com_google_protobuf//:protoc",
    py_libs = ["//third_party/nucleus:__init__py"],
    deps = [
        ":feature_py_pb2",
    ],
)

# ----------------------------------------------------------------------
# pyclif targets
# ----------------------------------------------------------------------

pyclif_proto_library(
    name = "bed_pyclif",
    proto_lib = ":bed_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "bedgraph_pyclif",
    proto_lib = ":bedgraph_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "gff_pyclif",
    proto_lib = ":gff_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "fasta_pyclif",
    proto_lib = ":fasta_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "fastq_pyclif",
    proto_lib = ":fastq_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "position_pyclif",
    proto_lib = ":position_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "range_pyclif",
    proto_lib = ":range_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "reads_pyclif",
    proto_lib = ":reads_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "reference_pyclif",
    proto_lib = ":reference_proto",  # NO COPYBARA
)

pyclif_proto_library(
    name = "variants_pyclif",
    proto_lib = ":variants_proto",  # NO COPYBARA
)
