# Description:
# (https://github.com/vgteam/libhandlegraph)
# libhandlegraph: A handle-based abstraction for graph access

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # MIT

exports_files(["LICENSE"])

includes = [
    "src/include",
    ".",
]

LIBRARY_COPTS = [
    "-w",  # turn off all warnings
]

SRCS = glob(
    ["src/*cpp"],
    exclude = [],
)

PRIVATE_HDRS = glob(["src/include/**/*"])

cc_library(
    name = "libhandlegraph",
    srcs = PRIVATE_HDRS + SRCS,
    hdrs = PRIVATE_HDRS,
    copts = LIBRARY_COPTS,
    includes = includes,
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:absl_log",
    ],
)
