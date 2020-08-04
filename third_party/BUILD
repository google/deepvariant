load("//tools/build_defs/testing:bzl_library.bzl", "bzl_library")

exports_files([
    "abseil_cpp.LICENSE",
    "boost.LICENSE",
    "tensorflow_models.LICENSE",
])

bzl_library(
    name = "clif_bzl",
    srcs = ["clif.bzl"],
    parse_tests = False,
    visibility = ["//visibility:private"],
)
