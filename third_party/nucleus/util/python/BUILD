load("//:third_party/clif.bzl", "py_clif_cc")

# load("//devtools/python/blaze:python3.bzl", "py2and3_test")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

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
        "@absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "utils",
    srcs = ["utils.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [
        "//third_party/nucleus/protos:range_pyclif",
        "//third_party/nucleus/protos:reads_pyclif",
    ],
    deps = [
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:proto_clif_converter",
    ],
)
