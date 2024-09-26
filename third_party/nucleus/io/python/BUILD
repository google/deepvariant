# Placeholder: load py_test
load("@org_tensorflow//tensorflow:tensorflow.default.bzl", "pybind_extension", "pybind_library")

package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = ["//visibility:public"],
)

licenses(["notice"])

cc_library(
    name = "tstring_clif_conversions",
    hdrs = ["tstring_clif_conversions.h"],
    deps = [
        "@clif//:cpp_runtime",
        "@org_tensorflow//tensorflow/tsl/platform:tstring",
    ],
)

pybind_extension(
    name = "vcf_reader",
    srcs = ["vcf_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "vcf_reader_cclib",
    srcs = ["vcf_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "vcf_reader_wrap_test",
    size = "small",
    srcs = ["vcf_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":vcf_reader",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
    ],
)

pybind_extension(
    name = "vcf_writer",
    srcs = ["vcf_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "vcf_writer_cclib",
    srcs = ["vcf_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "vcf_concat",
    srcs = ["vcf_concat_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_concat",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "vcf_concat_cclib",
    srcs = ["vcf_concat_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:vcf_concat",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "vcf_writer_wrap_test",
    size = "small",
    srcs = ["vcf_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":vcf_writer",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "gbz_reader",
    srcs = ["gbz_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gbz_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "gbz_reader_cclib",
    srcs = ["gbz_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gbz_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "sam_reader",
    srcs = ["sam_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:sam_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "sam_reader_cclib",
    srcs = ["sam_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:sam_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "sam_reader_wrap_test",
    size = "small",
    srcs = ["sam_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":sam_reader",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "sam_writer",
    srcs = ["sam_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:sam_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "sam_writer_cclib",
    srcs = ["sam_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:sam_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "reference",
    srcs = ["reference_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:reader_base",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "reference_cclib",
    srcs = ["reference_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:reader_base",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "reference_wrap_test",
    size = "small",
    srcs = ["reference_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":reference",
        "//third_party/nucleus/protos:fasta_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "bed_reader",
    srcs = ["bed_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bed_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "bed_reader_cclib",
    srcs = ["bed_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bed_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "bed_reader_wrap_test",
    size = "small",
    srcs = ["bed_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":bed_reader",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/protos:bed_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "bed_writer",
    srcs = ["bed_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bed_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "bed_writer_cclib",
    srcs = ["bed_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bed_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "bed_writer_wrap_test",
    size = "small",
    srcs = ["bed_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":bed_writer",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:bed_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "bedgraph_reader",
    srcs = ["bedgraph_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bedgraph_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "bedgraph_reader_cclib",
    srcs = ["bedgraph_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bedgraph_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "bedgraph_writer",
    srcs = ["bedgraph_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bedgraph_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "bedgraph_writer_cclib",
    srcs = ["bedgraph_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:bedgraph_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "fastq_reader",
    srcs = ["fastq_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:fastq_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "fastq_reader_cclib",
    srcs = ["fastq_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:fastq_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "fastq_reader_wrap_test",
    size = "small",
    srcs = ["fastq_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":fastq_reader",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/protos:fastq_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "fastq_writer",
    srcs = ["fastq_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:fastq_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "fastq_writer_cclib",
    srcs = ["fastq_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:fastq_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "fastq_writer_wrap_test",
    size = "small",
    srcs = ["fastq_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":fastq_writer",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/io:fastq",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:fastq_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "hts_verbose",
    srcs = ["hts_verbose_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:hts_verbose",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "hts_verbose_cclib",
    srcs = ["hts_verbose_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:hts_verbose",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "hts_verbose_test",
    size = "small",
    srcs = ["hts_verbose_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":hts_verbose",
        "@absl_py//absl/testing:absltest",
    ],
)

pybind_extension(
    name = "gff_reader",
    srcs = ["gff_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gff_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "gff_reader_cclib",
    srcs = ["gff_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gff_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "gff_reader_wrap_test",
    size = "small",
    srcs = ["gff_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":gff_reader",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/protos:gff_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "gff_writer",
    srcs = ["gff_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gff_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "gff_writer_cclib",
    srcs = ["gff_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gff_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

py_test(
    name = "gff_writer_wrap_test",
    size = "small",
    srcs = ["gff_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":gff_writer",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:gff_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

pybind_extension(
    name = "tabix_indexer",
    srcs = ["tabix_indexer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:tabix_indexer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "tabix_indexer_cclib",
    srcs = ["tabix_indexer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:tabix_indexer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "tfrecord_writer",
    srcs = ["tfrecord_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:tfrecord_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "tfrecord_writer_cclib",
    srcs = ["tfrecord_writer_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:tfrecord_writer",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "tfrecord_reader",
    srcs = ["tfrecord_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:tfrecord_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "tfrecord_reader_cclib",
    srcs = ["tfrecord_reader_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:tfrecord_reader",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "gfile",
    srcs = ["gfile_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gfile_cc",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "gfile_cclib",
    srcs = ["gfile_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:gfile_cc",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_extension(
    name = "merge_variants",
    srcs = ["merge_variants_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:merge_variants",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)

pybind_library(
    name = "merge_variants_cclib",
    srcs = ["merge_variants_pybind.cc"],
    deps = [
        "//third_party/nucleus/core/python:type_caster_nucleus_status",
        "//third_party/nucleus/core/python:type_caster_nucleus_statusor",
        "//third_party/nucleus/io:merge_variants",
        "//third_party/nucleus/util/python:type_caster_nucleus_proto_ptr",
        "@pybind11_protobuf//pybind11_protobuf:native_proto_caster",
    ],
)
