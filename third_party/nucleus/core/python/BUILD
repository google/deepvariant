# Placeholder: load py_test
load("@org_tensorflow//tensorflow:tensorflow.default.bzl", "pybind_extension", "pybind_library")

package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

pybind_library(
    name = "type_caster_nucleus_status",
    hdrs = ["type_caster_nucleus_status.h"],
    deps = [
        "//third_party/nucleus/core:status",
    ],
)

pybind_library(
    name = "type_caster_nucleus_statusor",
    hdrs = ["type_caster_nucleus_statusor.h"],
    deps = [
        ":type_caster_nucleus_status",
        "//third_party/nucleus/core:status",
        "//third_party/nucleus/core:statusor",
    ],
)

pybind_extension(
    name = "statusor_examples",
    srcs = ["statusor_examples_pybind.cc"],
    deps = [
        ":type_caster_nucleus_status",
        ":type_caster_nucleus_statusor",
        "//third_party/nucleus/core:statusor_examples",
    ],
)

pybind_library(
    name = "statusor_examples_cclib",
    srcs = ["statusor_examples_pybind.cc"],
    deps = [
        ":type_caster_nucleus_status",
        ":type_caster_nucleus_statusor",
        "//third_party/nucleus/core:statusor_examples",
    ],
)

py_test(
    name = "statusor_examples_test",
    size = "small",
    srcs = ["statusor_examples_test.py"],
    deps = [
        ":statusor_examples",
        "@absl_py//absl/testing:absltest",
    ],
)
