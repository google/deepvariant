# Placeholder: load py_test

package(
    default_visibility = [
        "//deepvariant:__subpackages__",
    ],
)

py_library(
    name = "make_small_model_examples",
    srcs = ["make_small_model_examples.py"],
    deps = [
        "//deepvariant/labeler:customized_classes_labeler",
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:variant_utils",
    ],
)

py_library(
    name = "small_model_json",
    srcs = ["small_model_json.py"],
    deps = ["//deepvariant:dv_vcf_constants"],
)

py_test(
    name = "small_model_json_test",
    srcs = ["small_model_json_test.py"],
    deps = [
        ":small_model_json",
        "//deepvariant:dv_vcf_constants",
        "@absl_py//absl/testing:absltest",
    ],
)

py_test(
    name = "make_small_model_examples_test",
    srcs = ["make_small_model_examples_test.py"],
    deps = [
        ":make_small_model_examples",
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "keras_config",
    srcs = ["keras_config.py"],
    deps = [
        ":make_small_model_examples",
    ],
)

py_test(
    name = "keras_config_reproduce_test",
    srcs = ["keras_config_reproduce_test.py"],
    deps = [
        ":keras_config",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "inference",
    srcs = ["inference.py"],
    deps = [
        ":keras_config",
        ":make_small_model_examples",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
    ],
)

py_test(
    name = "inference_test",
    srcs = ["inference_test.py"],
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
