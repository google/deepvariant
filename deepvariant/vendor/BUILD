load("//third_party/bazel_rules/rules_python/python:py_library.bzl", "py_library")
load("//third_party/bazel_rules/rules_python/python:py_test.bzl", "py_test")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_library(
    name = "timer",
    srcs = ["timer.py"],
)

py_test(
    name = "timer_test",
    size = "small",
    srcs = ["timer_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":timer",
        "@absl_py//absl/testing:absltest",
    ],
)
