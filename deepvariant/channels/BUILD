package(
    default_visibility = [
        "//visibility:public",
    ],
)

cc_library(
    name = "channel",
    srcs = ["channel.cc"],
    hdrs = ["channel.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "mapping_quality_channel",
    srcs = ["mapping_quality_channel.cc"],
    hdrs = ["mapping_quality_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "avg_base_quality_channel",
    srcs = ["avg_base_quality_channel.cc"],
    hdrs = ["avg_base_quality_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/log",
    ],
)

cc_library(
    name = "gap_compressed_identity_channel",
    srcs = ["gap_compressed_identity_channel.cc"],
    hdrs = ["gap_compressed_identity_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "gc_content_channel",
    srcs = ["gc_content_channel.cc"],
    hdrs = ["gc_content_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "homopolymer_weighted_channel",
    srcs = ["homopolymer_weighted_channel.cc"],
    hdrs = ["homopolymer_weighted_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "read_supports_variant_channel",
    srcs = ["read_supports_variant_channel.cc"],
    hdrs = ["read_supports_variant_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "identity_channel",
    srcs = ["identity_channel.cc"],
    hdrs = ["identity_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "read_mapping_percent_channel",
    srcs = ["read_mapping_percent_channel.cc"],
    hdrs = ["read_mapping_percent_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "allele_frequency_channel",
    srcs = ["allele_frequency_channel.cc"],
    hdrs = ["allele_frequency_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "is_homopolymer_channel",
    srcs = ["is_homopolymer_channel.cc"],
    hdrs = ["is_homopolymer_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "strand_channel",
    srcs = ["strand_channel.cc"],
    hdrs = ["strand_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "insert_size_channel",
    srcs = ["insert_size_channel.cc"],
    hdrs = ["insert_size_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "haplotype_tag_channel",
    srcs = ["haplotype_tag_channel.cc"],
    hdrs = ["haplotype_tag_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/log",
    ],
)

cc_library(
    name = "blank_channel",
    srcs = ["blank_channel.cc"],
    hdrs = ["blank_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)
