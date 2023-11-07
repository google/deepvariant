load("//third_party:clif.bzl", "py_clif_cc")
# Placeholder: load py_test

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
        "//third_party/nucleus/core:statusor_clif_converters",
        "//third_party/nucleus/core:statusor_examples",
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
        "//third_party/py/google/protobuf:use_fast_cpp_protos",
        "@absl_py//absl/testing:absltest",
    ],
)
