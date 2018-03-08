package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_library(
    name = "variant_labeler",
    srcs = ["variant_labeler.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/util:variant_utils",
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
        "//deepvariant/util:ranges",
        "//deepvariant/util:variant_utils",
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
        ":positional_labeler_test",
        "//deepvariant/util:ranges",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)
