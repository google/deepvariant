# Description:
#   Tools for building the Nucleus pip package.

load("@org_tensorflow//tensorflow:tensorflow.bzl", "transitive_hdrs")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

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
        "//third_party:abseil_cpp.LICENSE",
        "//third_party:bazel.LICENSE",
        "//third_party:clif.LICENSE",
        "//:LICENSE",
        "@htslib//:LICENSE",
        # The following dependencies contain LICENSE files, but
        # don't define proper BUILD files and/or don't export the LICENSE files.
        # redacted
        # write pull requests for exporting them.
        # "@com_google_absl_py//:LICENSE",
        # "@com_google_googletest//:LICENSE",
        # "@com_googlesource_code_cctz//:LICENSE.txt",
        # "@io_bazel_rules_closure//:LICENSE",
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
