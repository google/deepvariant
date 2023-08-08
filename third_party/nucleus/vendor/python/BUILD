# Placeholder: load py_test
load("//third_party:clif.bzl", "py_clif_cc")

package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

py_clif_cc(
    name = "statusor_examples",
    srcs = ["statusor_examples.clif"],
    deps = [
        "//third_party/nucleus/vendor:statusor_clif_converters",
        "//third_party/nucleus/vendor:statusor_examples",
    ],
)

py_test(
    name = "statusor_examples_test",
    size = "small",
    srcs = ["statusor_examples_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":statusor_examples",
        "@absl_py//absl/testing:absltest",
    ],
)
