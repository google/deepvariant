load("//:third_party/clif.bzl", "py_clif_cc")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

cc_library(
    name = "statusor",
    hdrs = [
        "statusor.h",
    ],
    deps = [
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_library(
    name = "status_matchers",
    testonly = True,
    hdrs = ["status_matchers.h"],
    deps = [
        ":statusor",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_test(
    name = "statusor_test",
    size = "small",
    srcs = ["statusor_test.cc"],
    deps = [
        ":status_matchers",
        ":statusor",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "statusor_clif_converters",
    srcs = ["statusor_clif_converters.cc"],
    hdrs = ["statusor_clif_converters.h"],
    deps = [
        ":statusor",
        "@clif//:cpp_runtime",
        "@local_config_python//:python_headers",
    ],
)

cc_library(
    name = "statusor_examples",
    hdrs = ["statusor_examples.h"],
    deps = [
        ":statusor",
    ],
)

# The zlib_* headers come from tensorflow/core/lib/io/.
# The corresponding implementations are not included as that causes multiple
# definition conflicts in the CLIF/Python readers.
cc_library(
    name = "zlib_compression_options",
    hdrs = ["zlib_compression_options.h"],
    deps = [
        "@org_tensorflow//tensorflow/core:lib",
        "@zlib_archive//:zlib",
    ],
)

cc_library(
    name = "zlib_inputstream",
    hdrs = ["zlib_inputstream.h"],
    deps = [
        ":zlib_compression_options",
        "@org_tensorflow//tensorflow/core:lib",
        "@zlib_archive//:zlib",
    ],
)

cc_library(
    name = "zlib_outputbuffer",
    hdrs = ["zlib_outputbuffer.h"],
    deps = [
        ":zlib_compression_options",
        "@org_tensorflow//tensorflow/core:lib",
        "@zlib_archive//:zlib",
    ],
)
