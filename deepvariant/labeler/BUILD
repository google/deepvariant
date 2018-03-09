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
        "//deepvariant:tf_utils",
        "//deepvariant/util:genomics_io",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:variantcall_utils",
        "//deepvariant/util/io:vcf",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
    ],
)

py_test(
    name = "labeled_examples_to_vcf_test",
    srcs = ["labeled_examples_to_vcf_test.py"],
    data = ["//deepvariant:testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":labeled_examples_to_vcf",
        "//deepvariant:py_test_utils",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/util/io:vcf",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "variant_labeler",
    srcs = ["variant_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:variantcall_utils",
    ],
)

py_test(
    name = "variant_labeler_test",
    size = "small",
    srcs = ["variant_labeler_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":variant_labeler",
        "//deepvariant:py_test_utils",
        "//deepvariant/util:in_memory_vcf_reader",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "positional_labeler",
    srcs = ["positional_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":variant_labeler",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:variantcall_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/logging",
    ],
)

py_test(
    name = "positional_labeler_test",
    size = "small",
    srcs = ["positional_labeler_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":positional_labeler",
        "//deepvariant:py_test_utils",
        "//deepvariant/util:in_memory_vcf_reader",
        "//deepvariant/util:ranges",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "haplotype_labeler",
    srcs = ["haplotype_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":variant_labeler",
        "//deepvariant/util:ranges",
        "//deepvariant/util:variant_utils",
    ],
)

py_test(
    name = "haplotype_labeler_test",
    srcs = ["haplotype_labeler_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":haplotype_labeler",
        "//deepvariant/util:in_memory_vcf_reader",
        "//deepvariant/util:ranges",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)
