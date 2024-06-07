# Placeholder: load py_binary
# Placeholder: load py_library
# Placeholder: load py_test

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_binary(
    name = "labeled_examples_to_vcf",
    srcs = ["labeled_examples_to_vcf.py"],
    python_version = "PY3",
    srcs_version = "PY2AND3",
    deps = [
        ":labeled_examples_to_vcf_main_lib",
        "@absl_py//absl:app",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "labeled_examples_to_vcf_main_lib",
    srcs = ["labeled_examples_to_vcf.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant:dv_utils",
        "//deepvariant:dv_vcf_constants",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "labeled_examples_to_vcf_test",
    srcs = ["labeled_examples_to_vcf_test.py"],
    data = ["//deepvariant:testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":labeled_examples_to_vcf_main_lib",
        "//deepvariant:py_testdata",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "variant_labeler",
    srcs = ["variant_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant:dv_utils",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "variant_labeler_test",
    size = "small",
    srcs = ["variant_labeler_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":variant_labeler",
        "//deepvariant:py_testdata",
        "//third_party/nucleus/io:vcf",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "positional_labeler",
    srcs = ["positional_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":variant_labeler",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "positional_labeler_test",
    size = "small",
    srcs = ["positional_labeler_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":positional_labeler",
        "//deepvariant:py_testdata",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "haplotype_labeler",
    srcs = ["haplotype_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
    ],
)

py_library(
    name = "customized_classes_labeler",
    srcs = ["customized_classes_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":positional_labeler",
        ":variant_labeler",
        "//third_party/nucleus/util:struct_utils",
    ],
)

py_test(
    name = "customized_classes_labeler_test",
    size = "small",
    srcs = ["customized_classes_labeler_test.py"],
    data = ["//deepvariant:testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":customized_classes_labeler",
        "//deepvariant:py_testdata",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "haplotype_labeler_test",
    srcs = ["haplotype_labeler_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":haplotype_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)
