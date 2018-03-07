load("//tools:zip_dir.bzl", "zip_dir")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

# redacted

filegroup(
    name = "binaries",
    srcs = [
        "//nucleus:binaries",
    ],
)

exports_files(["LICENSE"])

filegroup(
    name = "licenses",
    srcs = [
        ":LICENSE",
        "//third_party:abseil_cpp.LICENSE",  # redacted
        "@com_google_protobuf//:LICENSE",
        "@com_googlesource_code_re2//:LICENSE",
        "@htslib//:LICENSE",
        "@org_tensorflow//:LICENSE",
    ],
)

zip_dir(
    name = "licenses_zip",
    srcs = [":licenses"],
    zipname = "licenses.zip",
)
