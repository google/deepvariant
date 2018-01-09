load("//tools:zip_dir.bzl", "zip_dir")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

test_suite(
    name = "smoke_tests",
    tests = [
        "//deepvariant/core:smoke_tests",
        "//deepvariant/environment_tests:smoke_tests",
        "//deepvariant/testing:smoke_tests",
    ],
)

filegroup(
    name = "binaries",
    srcs = [
        "//deepvariant:binaries",
    ],
)

exports_files(["LICENSE"])

filegroup(
    name = "licenses",
    srcs = [
        ":LICENSE",
        "//third_party:abseil_cpp.LICENSE",  # redacted
        "//third_party:boost.LICENSE",
        "//third_party:tensorflow_models.LICENSE",  # redacted
        "@com_google_protobuf//:LICENSE",
        "@com_googlesource_code_re2//:LICENSE",
        "@htslib//:LICENSE",
        "@libssw//:README.md",  # SSW license embedded in the README.
        "@org_tensorflow//:LICENSE",
    ],
)

zip_dir(
    name = "licenses_zip",
    srcs = [":licenses"],
    zipname = "licenses.zip",
)
