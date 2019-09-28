# Description:
# (https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
#
package(default_visibility = ["//visibility:public"])

licenses(["notice"])  # MIT

# The license for libssw is embedded in README.md.
exports_files(["README.md"])

cc_library(
    name = "ssw",
    srcs = ["src/ssw.c"],
    hdrs = ["src/ssw.h"],
    copts = ["-fno-inline"],  # gcc-5.4 bug
)

cc_binary(
    name = "_c_example",
    srcs = ["src/example.c"],
    deps = [":ssw"],
)

cc_library(
    name = "ssw_cpp",
    srcs = ["src/ssw_cpp.cpp"],
    hdrs = ["src/ssw_cpp.h"],
    deps = [":ssw"],
)

cc_binary(
    name = "_cpp_example",
    srcs = ["src/example.cpp"],
    deps = [":ssw_cpp"],
)
