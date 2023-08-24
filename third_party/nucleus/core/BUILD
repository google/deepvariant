package(
    #default_applicable_licenses = ["//third_party/nucleus:license"],
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

cc_library(
    name = "status",
    srcs = ["status.cc"],
    hdrs = [
        "status.h",
    ],
    deps = [
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:optional",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)
