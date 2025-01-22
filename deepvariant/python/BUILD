# Placeholder: load py_test
load("@org_tensorflow//tensorflow:tensorflow.default.bzl", "pybind_extension", "pybind_library")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

pybind_extension(
    name = "allelecounter",
    srcs = ["allelecounter_pybind.cc"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@com_google_absl//absl/strings",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "allelecounter_cclib",
    srcs = ["allelecounter_pybind.cc"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@com_google_absl//absl/strings",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "allelecounter_wrap_test",
    size = "small",
    srcs = ["allelecounter_wrap_test.py"],
    data = [
        "//deepvariant:testdata",
        "//third_party/nucleus/testdata",
    ],
    deps = [
        ":allelecounter",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

pybind_extension(
    name = "postprocess_variants",
    srcs = ["postprocess_variants_pybind.cc"],
    deps = [
        "//deepvariant:postprocess_variants_lib",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "postprocess_variants_cclib",
    srcs = ["postprocess_variants_pybind.cc"],
    deps = [
        "//deepvariant:postprocess_variants_lib",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "variant_calling",
    srcs = ["variant_calling_pybind.cc"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant:variant_calling",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "variant_calling_cclib",
    srcs = ["variant_calling_pybind.cc"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant:variant_calling",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "variant_calling_wrap_test",
    size = "small",
    srcs = ["variant_calling_wrap_test.py"],
    data = ["//deepvariant:testdata"],
    deps = [
        ":allelecounter",
        ":variant_calling",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

pybind_extension(
    name = "direct_phasing",
    srcs = ["direct_phasing_pybind.cc"],
    deps = [
        "//deepvariant:direct_phasing",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "direct_phasing_cclib",
    srcs = ["direct_phasing_pybind.cc"],
    deps = [
        "//deepvariant:direct_phasing",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "pileup_image_native",
    srcs = ["pileup_image_native_pybind.cc"],
    deps = [
        "//deepvariant:pileup_image_native",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@com_google_absl//absl/strings",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "pileup_image_native_cclib",
    srcs = ["pileup_image_native_pybind.cc"],
    deps = [
        "//deepvariant:pileup_image_native",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@com_google_absl//absl/strings",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "variant_calling_multisample",
    srcs = ["variant_calling_multisample_pybind.cc"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant:variant_calling_multisample",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "variant_calling_multisample_cclib",
    srcs = ["variant_calling_multisample_pybind.cc"],
    deps = [
        "//deepvariant:allelecounter",
        "//deepvariant:variant_calling_multisample",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "variant_calling_multisample_wrap_test",
    size = "small",
    srcs = ["variant_calling_multisample_wrap_test.py"],
    data = ["//deepvariant:testdata"],
    deps = [
        ":allelecounter",
        ":variant_calling_multisample",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

pybind_extension(
    name = "make_examples_native",
    srcs = ["make_examples_native_pybind.cc"],
    deps = [
        "//deepvariant:make_examples_native",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "make_examples_native_cclib",
    srcs = ["make_examples_native_pybind.cc"],
    deps = [
        "//deepvariant:make_examples_native",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)
