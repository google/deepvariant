# KEEP THIS LINE to import nucleus.bzl externally
# load("//devtools/python/blaze:python3.bzl", "py2and3_test")

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
    python_version = "PY2",
    srcs_version = "PY2AND3",
)

py_test(
    name = "protobuf_implementation_test",
    size = "small",
    srcs = ["protobuf_implementation_test.py"],
    python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = ["@com_google_protobuf//:protobuf_python"],
)
