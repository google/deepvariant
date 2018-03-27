load("//:third_party/clif.bzl", "py_clif_cc", "pyclif_proto_library")
load("@protobuf_archive//:protobuf.bzl", "py_proto_library")
load("@protobuf_archive//:protobuf.bzl", "cc_proto_library")

package(default_visibility = [
    "//visibility:public",
])

# ----------------------------------------------------------------------
# Vendor'ed ListValue from google/protobuf/struct.proto
# ----------------------------------------------------------------------

proto_library(
    name = "struct_proto",
    srcs = ["struct.proto"],
)

cc_proto_library(
    name = "struct_cc_pb2",
    srcs = ["struct.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

py_proto_library(
    name = "struct_py_pb2",
    srcs = ["struct.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
)

# DO we need a pyclif rule?

# ----------------------------------------------------------------------
# Cigar
# ----------------------------------------------------------------------

proto_library(
    name = "cigar_proto",
    srcs = ["cigar.proto"],
)

cc_proto_library(
    name = "cigar_cc_pb2",
    srcs = ["cigar.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

py_proto_library(
    name = "cigar_py_pb2",
    srcs = ["cigar.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
)

pyclif_proto_library(
    name = "cigar_pyclif",
    proto_lib = ":cigar_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# Position
# ----------------------------------------------------------------------

proto_library(
    name = "position_proto",
    srcs = ["position.proto"],
)

cc_proto_library(
    name = "position_cc_pb2",
    srcs = ["position.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

py_proto_library(
    name = "position_py_pb2",
    srcs = ["position.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
)

pyclif_proto_library(
    name = "position_pyclif",
    proto_lib = ":position_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# Range
# ----------------------------------------------------------------------

proto_library(
    name = "range_proto",
    srcs = ["range.proto"],
)

# redacted
# redacted
cc_proto_library(
    name = "range_cc_pb2",
    srcs = ["range.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

py_proto_library(
    name = "range_py_pb2",
    srcs = ["range.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
)

pyclif_proto_library(
    name = "range_pyclif",
    proto_lib = ":range_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# Variants (VCF)
# ----------------------------------------------------------------------

proto_library(
    name = "variants_proto",
    srcs = ["variants.proto"],
    deps = [
        ":index_proto",  # NO COPYBARA
        ":reference_proto",  # NO COPYBARA
        ":struct_proto",  # NO COPYBARA
    ],
)

cc_proto_library(
    name = "variants_cc_pb2",
    srcs = ["variants.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [
        ":index_cc_pb2",
        ":reference_cc_pb2",
        ":struct_cc_pb2",
    ],
)

py_proto_library(
    name = "variants_py_pb2",
    srcs = ["variants.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
    deps = [
        ":index_py_pb2",
        ":reference_py_pb2",
        ":struct_py_pb2",
    ],
)

pyclif_proto_library(
    name = "variants_pyclif",
    proto_lib = ":variants_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# Reads (SAM/BAM)
# ----------------------------------------------------------------------

proto_library(
    name = "reads_proto",
    srcs = ["reads.proto"],
    deps = [
        ":cigar_proto",  # NO COPYBARA
        ":index_proto",  # NO COPYBARA
        ":position_proto",  # NO COPYBARA
        ":reference_proto",  # NO COPYBARA
        ":struct_proto",  # NO COPYBARA
    ],
)

cc_proto_library(
    name = "reads_cc_pb2",
    srcs = ["reads.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [
        ":cigar_cc_pb2",
        ":index_cc_pb2",
        ":position_cc_pb2",
        ":reference_cc_pb2",
        ":struct_cc_pb2",
    ],
)

py_proto_library(
    name = "reads_py_pb2",
    srcs = ["reads.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
    deps = [
        ":cigar_py_pb2",
        ":index_py_pb2",
        ":position_py_pb2",
        ":reference_py_pb2",
        ":struct_py_pb2",
    ],
)

pyclif_proto_library(
    name = "reads_pyclif",
    proto_lib = ":reads_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# Reference
# ----------------------------------------------------------------------

proto_library(
    name = "reference_proto",
    srcs = ["reference.proto"],
)

cc_proto_library(
    name = "reference_cc_pb2",
    srcs = ["reference.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

py_proto_library(
    name = "reference_py_pb2",
    srcs = ["reference.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
)

pyclif_proto_library(
    name = "reference_pyclif",
    proto_lib = ":reference_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# Index
# ----------------------------------------------------------------------

proto_library(
    name = "index_proto",
    srcs = ["index.proto"],
)

cc_proto_library(
    name = "index_cc_pb2",
    srcs = ["index.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

py_proto_library(
    name = "index_py_pb2",
    srcs = ["index.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
)

pyclif_proto_library(
    name = "index_pyclif",
    proto_lib = ":index_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# FASTQ records
# ----------------------------------------------------------------------

proto_library(
    name = "fastq_proto",
    srcs = ["fastq.proto"],
    deps = [],
)

cc_proto_library(
    name = "fastq_cc_pb2",
    srcs = ["fastq.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [],
)

py_proto_library(
    name = "fastq_py_pb2",
    srcs = ["fastq.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
    deps = [],
)

pyclif_proto_library(
    name = "fastq_pyclif",
    proto_lib = ":fastq_proto",  # NO COPYBARA
)

# ----------------------------------------------------------------------
# BED records
# ----------------------------------------------------------------------

proto_library(
    name = "bed_proto",
    srcs = ["bed.proto"],
    deps = [],
)

cc_proto_library(
    name = "bed_cc_pb2",
    srcs = ["bed.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [],
)

py_proto_library(
    name = "bed_py_pb2",
    srcs = ["bed.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    srcs_version = "PY2AND3",
    deps = [],
)

pyclif_proto_library(
    name = "bed_pyclif",
    proto_lib = ":bed_proto",  # NO COPYBARA
)
