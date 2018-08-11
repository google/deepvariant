load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

nucleus_py_binary(
    name = "generate_trained_model",
    srcs = ["generate_trained_model.py"],
    main = "generate_trained_model.py",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/python:allelecounter",
        "//third_party/nucleus/io:io_py",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:util_py",
    ],
)

py_test(
    name = "generate_trained_model_test",
    srcs = ["generate_trained_model_test.py"],
    data = ["//deepvariant:testdata"],
    deps = [
        ":generate_trained_model",
        "//deepvariant:py_testdata",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_binary(
    name = "model_evaluation",
    srcs = ["model_evaluation.py"],
    main = "model_evaluation.py",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/python:allelecounter",
        "//third_party/nucleus/io:io_py",
        "//third_party/nucleus/protos:cigar_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:util_py",
    ],
)

py_test(
    name = "model_evaluation_test",
    srcs = ["model_evaluation_test.py"],
    data = ["//deepvariant:testdata"],
    deps = [
        ":model_evaluation",
        "//deepvariant:py_testdata",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)
