load("//:third_party/clif.bzl", "py_clif_cc", "pyclif_proto_library")
load("@protobuf_archive//:protobuf.bzl", "py_proto_library")
load("@protobuf_archive//:protobuf.bzl", "cc_proto_library")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

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

# ----------------------------------------------------------------------
# cc_proto_library's
# ----------------------------------------------------------------------

# redacted
# boilerplate in these cc_proto_library rules and in this file in general.

cc_proto_library(
    name = "bed_cc_pb2",
    srcs = ["bed.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

cc_proto_library(
    name = "cigar_cc_pb2",
    srcs = ["cigar.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

cc_proto_library(
    name = "fasta_cc_pb2",
    srcs = ["fasta.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [
        ":range_cc_pb2",
        ":reference_cc_pb2",
    ],
)

cc_proto_library(
    name = "fastq_cc_pb2",
    srcs = ["fastq.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

cc_proto_library(
    name = "gff_cc_pb2",
    srcs = ["gff.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [":range_cc_pb2"],
)

cc_proto_library(
    name = "position_cc_pb2",
    srcs = ["position.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

cc_proto_library(
    name = "range_cc_pb2",
    srcs = ["range.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

cc_proto_library(
    name = "reads_cc_pb2",
    srcs = ["reads.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
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
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [
        ":range_cc_pb2",
    ],
)

cc_proto_library(
    name = "struct_cc_pb2",
    srcs = ["struct.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
)

cc_proto_library(
    name = "variants_cc_pb2",
    srcs = ["variants.proto"],
    default_runtime = "@protobuf_archive//:protobuf",
    protoc = "@protobuf_archive//:protoc",
    deps = [
        ":reference_cc_pb2",
        ":struct_cc_pb2",
    ],
)

# ----------------------------------------------------------------------
# fast_proto
# ----------------------------------------------------------------------

cc_binary(
    name = "load_descriptors.so",
    srcs = ["load_descriptors.cc"],
    linkshared = 1,
    linkstatic = 0,
    deps = [
        ":bed_cc_pb2",
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
        "@local_config_python//:python_headers",
    ],
)

# fast_protos is a single target py_libraries can depend on to make sure
# they get the fastest Python protocol buffer implementation.
nucleus_py_library(
    name = "fast_protos",
    srcs = ["load_descriptors.py"],
    data = [":load_descriptors.so"],
)

# ----------------------------------------------------------------------
# py_proto_library's
# ----------------------------------------------------------------------

py_proto_library(
    name = "bed_py_pb2",
    srcs = ["bed.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
)

py_proto_library(
    name = "cigar_py_pb2",
    srcs = ["cigar.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
)

py_proto_library(
    name = "gff_py_pb2",
    srcs = ["gff.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
    deps = [":range_py_pb2"],
)

py_proto_library(
    name = "fasta_py_pb2",
    srcs = ["fasta.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
    deps = [
        ":range_py_pb2",
        ":reference_py_pb2",
    ],
)

py_proto_library(
    name = "fastq_py_pb2",
    srcs = ["fastq.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
)

py_proto_library(
    name = "position_py_pb2",
    srcs = ["position.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
)

py_proto_library(
    name = "range_py_pb2",
    srcs = ["range.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
)

py_proto_library(
    name = "reads_py_pb2",
    srcs = ["reads.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
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
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
    deps = [
        ":range_py_pb2",
    ],
)

py_proto_library(
    name = "struct_py_pb2",
    srcs = ["struct.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
)

py_proto_library(
    name = "variants_py_pb2",
    srcs = ["variants.proto"],
    default_runtime = "@protobuf_archive//:protobuf_python",
    protoc = "@protobuf_archive//:protoc",
    py_libs = [":fast_protos"],
    srcs_version = "PY2AND3",
    deps = [
        ":reference_py_pb2",
        ":struct_py_pb2",
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
