# Placeholder: load py_test

package(
    default_visibility = [
        "//deepvariant:__subpackages__",
    ],
)

py_library(
    name = "make_small_model_examples",
    srcs = ["make_small_model_examples.py"],
    srcs_version = "PY3",
    deps = [
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:variant_utils",
    ],
)

py_test(
    name = "make_small_model_examples_test",
    srcs = ["make_small_model_examples_test.py"],
    python_version = "PY3",
    deps = [
        ":make_small_model_examples",
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "inference",
    srcs = ["inference.py"],
    srcs_version = "PY3",
    deps = [
        ":make_small_model_examples",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:variantcall_utils",
    ],
)

py_test(
    name = "inference_test",
    srcs = ["inference_test.py"],
    python_version = "PY3",
    deps = [
        ":inference",
        ":make_small_model_examples",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)
