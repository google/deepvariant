load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

package(default_visibility = ["//:__subpackages__"])

test_suite(
    name = "smoke_tests",
    tests = [
        "protobuf_implementation_test",
        ":env_smoke_test",
    ],
)

py_test(
    name = "env_smoke_test",
    size = "small",
    srcs = ["env_smoke_test.py"],
    srcs_version = "PY2AND3",
)

py_test(
    name = "protobuf_implementation_test",
    size = "small",
    srcs = ["protobuf_implementation_test.py"],
    srcs_version = "PY2AND3",
    deps = ["@protobuf_archive//:protobuf_python"],
)
