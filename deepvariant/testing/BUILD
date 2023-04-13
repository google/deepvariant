# Testing utilities that are temporarily located here until they can
# be separately open source.

load("//third_party/bazel_rules/rules_python/python:py_library.bzl", "py_library")

package(default_visibility = ["//:__subpackages__"])

py_library(
    name = "tf_test_utils",
    srcs = ["tf_test_utils.py"],
    srcs_version = "PY3",
    deps = [
        "//deepvariant:dv_constants",
        "//deepvariant:modeling",
        "@absl_py//absl/flags",
    ],
)
