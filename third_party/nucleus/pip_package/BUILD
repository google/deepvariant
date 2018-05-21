# Description:
#   Tools for building the Nucleus pip package.

load("@org_tensorflow//tensorflow:tensorflow.bzl", "transitive_hdrs")

licenses(["notice"])

transitive_hdrs(
    name = "included_headers",
    deps = [
        "//third_party/nucleus/io:io_cpp",
        # redacted
        # dependencies:
        # "//third_party/nucleus/testing:testing_cpp",
        # "//third_party/nucleus/testing:testing_py",
        "//third_party/nucleus/util:util_cpp",
        "//third_party/nucleus/vendor:vendor_cpp",
    ],
)

filegroup(
    name = "licenses",
    data = [
        # redacted
        # by Nucleus dependencies, including:
        # "//third_party/nucleus/opensource_only/third_party/abseil_cpp.LICENSE",
        # "//third_party/nucleus/LICENSE",
        # "//third_party/htslib/LICENSE",
        # "@clif//:LICENSE",
        # "@com_google_absl_py//:LICENSE",
        # "@com_google_googletest//:LICENSE",
        "@org_tensorflow//:LICENSE",
        "@protobuf_archive//:LICENSE",
    ],
)

PIP_DEPS = [
    ":licenses",
    ":included_headers",
    "setup.py",
    "//third_party/nucleus/io:io_py",
    "//third_party/nucleus/util:util_py",
]

sh_binary(
    name = "build_pip_package",
    srcs = ["build_pip_package.sh"],
    data = PIP_DEPS,
)
