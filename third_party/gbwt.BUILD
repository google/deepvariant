# Description:
# (https://github.com/jltsiren/gbwt)
# Graph BWT is an independent implementation of the graph extension (gPBWT) of the positional
# Burrows-Wheeler transform (PBWT).

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # MIT

exports_files(["LICENSE"])

includes = [
    "include",
    ".",
]

LIBRARY_COPTS = [
    "-w",  # turn off all warnings
]

SRCS = glob(
    ["src/*cpp"],
    exclude = [
        "src/remove_seq.cpp",
        "src/metadata_tool.cpp",
        "src/merge_gbwt.cpp",
        "src/build_ri.cpp",
        "src/build_gbwt.cpp",
        "src/benchmark.cpp",
    ],
)

PRIVATE_HDRS = glob(["include/gbwt/*.h"])

cc_library(
    name = "gbwt",
    srcs = PRIVATE_HDRS + SRCS,
    hdrs = PRIVATE_HDRS,
    copts = LIBRARY_COPTS,
    includes = includes,
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:absl_log",
        "@sdsl_lite",
    ],
)

cc_binary(
    name = "build_gbwt",
    srcs = ["src/build_gbwt.cpp"],
    copts = LIBRARY_COPTS,
    visibility = ["//visibility:public"],
    deps = [
        ":gbwt",
        "@com_google_absl//absl/log",
        "@sdsl_lite",
    ],
)

cc_binary(
    name = "build_ri",
    srcs = ["src/build_ri.cpp"],
    copts = LIBRARY_COPTS,
    visibility = ["//visibility:public"],
    deps = [
        ":gbwt",
        "@com_google_absl//absl/log",
        "@sdsl_lite",
    ],
)

cc_binary(
    name = "merge_gbwt",
    srcs = ["src/merge_gbwt.cpp"],
    copts = LIBRARY_COPTS,
    visibility = ["//visibility:public"],
    deps = [
        ":gbwt",
        "@com_google_absl//absl/log",
        "@sdsl_lite",
    ],
)

cc_binary(
    name = "metadata_tool",
    srcs = ["src/metadata_tool.cpp"],
    copts = LIBRARY_COPTS,
    visibility = ["//visibility:public"],
    deps = [
        ":gbwt",
        "@com_google_absl//absl/log",
        "@sdsl_lite",
    ],
)

cc_binary(
    name = "remove_seq",
    srcs = ["src/remove_seq.cpp"],
    copts = LIBRARY_COPTS,
    visibility = ["//visibility:public"],
    deps = [
        ":gbwt",
        "@com_google_absl//absl/log",
        "@sdsl_lite",
    ],
)
