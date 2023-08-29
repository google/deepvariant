# Description:
#   Tools for building the Nucleus pip package.

licenses(["notice"])

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
        # TODO: Copy these LICENSES to third_party, or
        # write pull requests for exporting them.
        # "@com_google_absl_py//:LICENSE",
        # "@com_google_googletest//:LICENSE",
        # "@com_googlesource_code_cctz//:LICENSE.txt",
        # "@io_bazel_rules_closure//:LICENSE",
        "@org_tensorflow//:LICENSE",
        "@com_google_protobuf//:LICENSE",
    ],
)

# This is a fake target, designed to populate build_pip_package.runfiles
# with the right dependency files (including the __init__.py's).  It is testonly
# so that it may depend on testonly targets in nucleus/testing.
py_binary(
    name = "setup",
    testonly = 1,
    srcs = ["setup.py"],
    data = [
        "setup.cfg",
    ],
    legacy_create_init = 1,
    deps = [
        "//third_party/nucleus/core/python:statusor_examples",
        "//third_party/nucleus/examples:examples_py",
        "//third_party/nucleus/io:io_py",
        "//third_party/nucleus/testing:testing_py",
        "//third_party/nucleus/util:util_py",
    ],
)

filegroup(
    name = "egg_files",
    data = [
        "//third_party/nucleus/pip_package:egg_files/PKG-INFO",
        "//third_party/nucleus/pip_package:egg_files/dependency_links.txt",
        "//third_party/nucleus/pip_package:egg_files/not-zip-safe",
        "//third_party/nucleus/pip_package:egg_files/requires.txt",
        "//third_party/nucleus/pip_package:egg_files/top_level.txt",
    ],
)

PIP_DEPS = [
    ":egg_files",
    ":licenses",
    ":setup",
    "//:all_extensions",
]

sh_binary(
    name = "build_pip_package",
    testonly = 1,
    srcs = ["build_pip_package.sh"],
    data = PIP_DEPS,
)
