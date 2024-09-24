# Description:
# (https://github.com/jltsiren/gbwtgraph)
# GBWTGraph is a handle graph based on the GBWT. Its data model is based on the graph as an
# alignment of haplotypes.

package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # MIT

exports_files(["LICENSE"])

include_gbwtgraph = "include"

includes = [
    include_gbwtgraph,
    ".",
]

LIBRARY_COPTS = [
    "-w",  # turn off all warnings
]

SRCS = glob(
    ["src/*cpp"],
    exclude = [
        "src/gbz_stats.cpp",
        "src/gfa2gbwt.cpp",
        "src/kmer_freq.cpp",
        "src/subgraph_query.cpp",
    ],
)

PRIVATE_HDRS = glob(["include/gbwtgraph/*.h"])

cc_library(
    name = "gbwtgraph",
    srcs = PRIVATE_HDRS + SRCS,
    hdrs = PRIVATE_HDRS,
    copts = LIBRARY_COPTS,
    includes = includes,
    visibility = ["//visibility:public"],
    deps = [
        "@com_google_absl//absl/log",
        "@gbwt",
        "@libhandlegraph",
        "@sdsl_lite",
    ],
)

cc_binary(
    name = "subgraph_query",
    srcs = ["src/subgraph_query.cpp"],
    copts = LIBRARY_COPTS,
    deps = [
        ":gbwtgraph",
        "@com_google_absl//absl/log",
    ],
)

cc_binary(
    name = "gbz_stats",
    srcs = ["src/gbz_stats.cpp"],
    copts = LIBRARY_COPTS,
    deps = [":gbwtgraph"],
)
