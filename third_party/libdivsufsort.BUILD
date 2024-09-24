# Description:
# (https://github.com/simongog/libdivsufsort)
# libdivsufsort is a software library that implements a lightweight suffix array construction
# algorithm.

package(
    default_visibility = ["//visibility:public"],
    features = [
        "-layering_check",
        "-parse_headers",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

# libdivsufsort defines two different libraries compiled from the same
# codebase. The libraries export different symbols and differ on the integer
# type used in the suffix array index. "libdivsufsort" uses 32-bit integers
# allowing input sizes up to 2GiB, while "libdivsufsort64" uses 64-bit integers
# allowing larger input sizes (8 EiB) but also requiring twice as much memory
# for the smaller inputs.

LIBDIVSUFSORT_SRCS = [
    "lib/divsufsort.c",
    "lib/sssort.c",
    "lib/trsort.c",
    "lib/utils.c",
]

# This list includes both generated headers and existing headers in the include/
# directory.
LIBDIVSUFSORT_PRIVATE_HDRS = ["include/config.h"] + glob(["include/*.h"])

LIBRARY_COPTS = [
    "-Wall",
    "-Wextra",
    "-DHAVE_CONFIG_H=1",
]

include_libdivsufsort = "include"

includes = [
    include_libdivsufsort,
    ".",
]

# libdivsufsort using 32-bit integers for the suffix array.
cc_library(
    name = "libdivsufsort",
    srcs = LIBDIVSUFSORT_SRCS + LIBDIVSUFSORT_PRIVATE_HDRS,
    hdrs = ["include/divsufsort.h"],
    copts = ["-Werror"] + LIBRARY_COPTS,
    includes = includes,
)

# libdivsufsort using 64-bit integers for the suffix array.
cc_library(
    name = "libdivsufsort64",
    srcs = LIBDIVSUFSORT_SRCS + LIBDIVSUFSORT_PRIVATE_HDRS,
    hdrs = ["include/divsufsort64.h"],
    copts = ["-Werror"] + (LIBRARY_COPTS + ["-DBUILD_DIVSUFSORT64"]),
    includes = includes,
)

cc_library(
    name = "libdivsufsort_system",
    hdrs = [
        "include/divsufsort.h",
        "include/divsufsort64.h",
    ],
    copts = ["-Werror"],
    includes = ["include"],
    visibility = ["//third_party:__subpackages__"],
    deps = [
        ":libdivsufsort",
        ":libdivsufsort64",
    ],
)

cc_library(
    name = "pydivsufsort",
    hdrs = ["pydivsufsort.h"],
    copts = ["-Werror"],
    deps = [
        ":libdivsufsort",
    ],
)

cc_test(
    name = "libdivsufsort_test",
    size = "small",
    srcs = ["tests/libdivsufsort_test.cc"],
    copts = ["-Werror"],
    deps = [
        ":libdivsufsort",
        ":libdivsufsort64",
        "//testing/base/public:gunit_main",
    ],
)

# Genrules.
# libdivsufsort uses cmake to generate some headers, including the public
# header. We use genrules here to generate those headers instead.
commom_awk_replaces = (
    "gsub(/#cmakedefine/, \"#define\"); " +
    "gsub(/@DIVSUFSORT_EXPORT@/, \"\"); " +
    "gsub(/@DIVSUFSORT_IMPORT@/, \"\"); " +
    "gsub(/@INLINE@/, \"inline\"); " +
    "gsub(/@INCFILE@/, \"#include <inttypes.h>\"); " +
    "gsub(/@SAUCHAR_TYPE@/, \"uint8_t\"); " +
    "gsub(/@SAINT32_TYPE@/, \"int32_t\"); " +
    "gsub(/@SAINT_PRId@/, \"PRId32\"); "
)

genrule(
    name = "config_h",
    srcs = ["include/config.h.cmake"],
    outs = ["include/config.h"],
    cmd = ("awk '{ " +
           "gsub(/@HAVE_IO_H 1@/, \"HAVE_IO_H 0\"); " +
           commom_awk_replaces +
           "print; }' $(<) > $(@)"),
)

genrule(
    name = "divsufsort_h",
    srcs = ["include/divsufsort.h.cmake"],
    outs = ["include/divsufsort.h"],
    cmd = ("awk '{ " +
           "gsub(/@W64BIT@/, \"\"); " +
           "gsub(/@SAINDEX_TYPE@/, \"int32_t\"); " +
           "gsub(/@SAINDEX_PRId@/, \"PRId32\"); " +
           commom_awk_replaces +
           "print; }' $(<) > $(@)"),
)

genrule(
    name = "divsufsort64_h",
    srcs = ["include/divsufsort.h.cmake"],
    outs = ["include/divsufsort64.h"],
    cmd = ("awk '{ " +
           "gsub(/@W64BIT@/, \"64\"); " +
           "gsub(/@SAINDEX_TYPE@/, \"int64_t\"); " +
           "gsub(/@SAINDEX_PRId@/, \"PRId64\"); " +
           commom_awk_replaces +
           "print; }' $(<) > $(@)"),
)
