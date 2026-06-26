# Description:
#   libdeflate - fast DEFLATE/zlib/gzip compression and decompression.
#   Used by htslib (when HAVE_LIBDEFLATE is set) to accelerate bgzf (BAM)
#   (de)compression.

licenses(["notice"])  # MIT

exports_files(["COPYING"])

# All architecture-specific sources are guarded by ARCH_* #ifdefs and compile to
# nothing on non-matching targets, so it is safe to compile every lib/**/*.c.
cc_library(
    name = "libdeflate",
    srcs = glob([
        "lib/**/*.c",
        "lib/**/*.h",
        "common_defs.h",
    ]),
    hdrs = ["libdeflate.h"],
    copts = ["-O3"],
    includes = ["."],
    visibility = ["//visibility:public"],
)
