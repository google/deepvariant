# KEEP THIS LINE to import nucleus.bzl externally

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_binary(
    name = "labeled_examples_to_vcf",
    srcs = ["labeled_examples_to_vcf.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant:dv_vcf_constants",
        "//deepvariant:tf_utils",
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
    srcs_version = "PY2AND3",
    deps = [
        ":labeled_examples_to_vcf",
        "//deepvariant:py_testdata",
        "//deepvariant/testing:flagsaver",
        "//third_party/nucleus/io:vcf",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "variant_labeler",
    srcs = ["variant_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant:tf_utils",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
    ],
)

py_test(
    name = "variant_labeler_test",
    size = "small",
    srcs = ["variant_labeler_test.py"],
    srcs_version = "PY2AND3",
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
    srcs_version = "PY2AND3",
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
    srcs_version = "PY2AND3",
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
    srcs_version = "PY2AND3",
    deps = [
        ":haplotype_labeler",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)
