load("//devtools/python/blaze:pytype.bzl", "pytype_strict_binary")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_binary(
    name = "preprocess_truth",
    srcs = ["preprocess_truth.py"],
    deps = [
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)
