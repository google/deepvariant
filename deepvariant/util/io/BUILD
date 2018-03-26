# Code for reading and writing genomics data.

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

py_library(
    name = "genomics_reader",
    srcs = ["genomics_reader.py"],
    deps = ["@com_google_absl_py//absl/logging"],
)

py_test(
    name = "genomics_reader_test",
    srcs = ["genomics_reader_test.py"],
    deps = [
        ":genomics_reader",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_library(
    name = "genomics_writer",
    srcs = ["genomics_writer.py"],
    deps = ["@com_google_absl_py//absl/logging"],
)

py_library(
    name = "vcf",
    srcs = ["vcf.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//deepvariant/util:vcf_constants",
        "//deepvariant/util/genomics:index_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/genomics:vcf_py_pb2",
        "//deepvariant/util/python:vcf_reader",
        "//deepvariant/util/python:vcf_writer",
    ],
)

py_test(
    name = "vcf_test",
    size = "small",
    srcs = ["vcf_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    deps = [
        ":vcf",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "sam",
    srcs = ["sam.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//deepvariant/util/genomics:index_py_pb2",
        "//deepvariant/util/genomics:reads_py_pb2",
        "//deepvariant/util/python:sam_reader",
        # "//third_party/py/tfbam_lib",
    ],
)

py_test(
    name = "sam_test",
    size = "small",
    srcs = ["sam_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    deps = [
        ":sam",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:struct_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "fasta",
    srcs = ["fasta.py"],
    deps = [
        ":genomics_reader",
        "//deepvariant/util/python:reference_fai",
    ],
)

py_test(
    name = "fasta_test",
    size = "small",
    srcs = ["fasta_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    deps = [
        ":fasta",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:struct_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "fastq",
    srcs = ["fastq.py"],
    deps = [
        ":genomics_reader",
        ":genomics_writer",
        "//deepvariant/util/genomics:fastq_py_pb2",
        "//deepvariant/util/python:fastq_reader",
        "//deepvariant/util/python:fastq_writer",
    ],
)

py_test(
    name = "fastq_test",
    size = "small",
    srcs = ["fastq_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    deps = [
        ":fastq",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util/genomics:fastq_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)
