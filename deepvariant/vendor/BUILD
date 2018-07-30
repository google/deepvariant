load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

nucleus_py_library(
    name = "timer",
    srcs = ["timer.py"],
)

py_test(
    name = "timer_test",
    size = "small",
    srcs = ["timer_test.py"],
    deps = [
        ":timer",
        "@io_abseil_py//absl/testing:absltest",
    ],
)
