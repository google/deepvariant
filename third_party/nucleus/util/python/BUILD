load("//:third_party/clif.bzl", "py_clif_cc")
# load("//devtools/python/blaze:python3.bzl", "py2and3_test")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

py_clif_cc(
    name = "math",
    srcs = ["math.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [],
    deps = ["//third_party/nucleus/util:cpp_math"],
)

py_test(
    name = "math_wrap_test",
    size = "small",
    srcs = ["math_wrap_test.py"],
    data = [],
    srcs_version = "PY2AND3",
    deps = [
        ":math",
        "@io_abseil_py//absl/testing:absltest",
    ],
)
