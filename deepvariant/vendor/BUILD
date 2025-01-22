# Placeholder: load py_library
# Placeholder: load py_test

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
    deps = [
        ":timer",
        "@absl_py//absl/testing:absltest",
    ],
)
