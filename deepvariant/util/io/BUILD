# Code for reading and writing genomics data.

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])  # Apache 2.0 License

exports_files(["LICENSE"])

py_library(
    name = "genomics_writer",
    srcs = ["genomics_writer.py"],
)

py_library(
    name = "vcf",
    srcs = ["vcf.py"],
    deps = [
        ":genomics_writer",
        "//deepvariant/util/protos:core_py_pb2",
        "//deepvariant/util/python:vcf_writer",
    ],
)

py_test(
    name = "vcf_test",
    size = "small",
    srcs = ["vcf_test.py"],
    deps = [
        ":vcf",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/protos:core_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "sam",
    srcs = ["sam.py"],
    deps = [
        ":genomics_writer",
    ],
)

py_test(
    name = "sam_test",
    size = "small",
    srcs = ["sam_test.py"],
    deps = [
        ":sam",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/protos:core_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)
