# Testing utilities that are temporarily located here until they can
# be separately open source.

package(default_visibility = ["//:__subpackages__"])

py_library(
    name = "flagsaver",
    srcs = ["flagsaver.py"],
    srcs_version = "PY2AND3",
    deps = ["@io_abseil_py//absl/flags"],
)

py_library(
    name = "tf_test_utils",
    srcs = ["tf_test_utils.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant:dv_constants",
        "//deepvariant:modeling",
        "@io_abseil_py//absl/flags",
    ],
)
