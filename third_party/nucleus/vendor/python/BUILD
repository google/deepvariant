load("//:third_party/clif.bzl", "py_clif_cc")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

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
    deps = [
        ":statusor_examples",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)
