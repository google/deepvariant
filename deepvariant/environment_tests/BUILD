# KEEP THIS LINE to import nucleus.bzl externally

# Placeholder: load py_test

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
)

py_test(
    name = "protobuf_implementation_test",
    size = "small",
    srcs = ["protobuf_implementation_test.py"],
    deps = ["@com_google_protobuf//:protobuf_python"],
)
