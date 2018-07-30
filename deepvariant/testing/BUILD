# Testing utilities that are temporarily located here until they can
# be separately open source.
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(default_visibility = ["//:__subpackages__"])

nucleus_py_library(
    name = "flagsaver",
    srcs = ["flagsaver.py"],
    srcs_version = "PY2AND3",
    deps = ["@io_abseil_py//absl/flags"],
)

nucleus_py_library(
    name = "tf_test_utils",
    srcs = ["tf_test_utils.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant:dv_constants",
        "//deepvariant:modeling",
        "@io_abseil_py//absl/flags",
    ],
)
