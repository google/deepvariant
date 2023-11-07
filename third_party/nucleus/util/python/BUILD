load("//third_party:clif.bzl", "py_clif_cc")
# Placeholder: load py_test

package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

py_clif_cc(
    name = "math",
    srcs = ["math.clif"],
    py_deps = [],
    pyclif_deps = [],
    deps = ["//third_party/nucleus/util:cpp_math"],
)

py_test(
    name = "math_wrap_test",
    size = "small",
    srcs = ["math_wrap_test.py"],
    data = [],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":math",
        "//third_party/py/google/protobuf:use_fast_cpp_protos",
        "@absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "utils",
    srcs = ["utils.clif"],
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
