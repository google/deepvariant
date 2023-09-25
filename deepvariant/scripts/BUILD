# Build run_deepvariant.py and run_deeptrio.py.

# Placeholder: load py_binary
# Placeholder: load py_test

filegroup(
    name = "binaries",
    srcs = [
        "run_deepsomatic",
        "run_deeptrio",
        "run_deepvariant",
    ],
)

py_library(
    name = "run_deepvariant_lib",
    srcs = ["run_deepvariant.py"],
    srcs_version = "PY2AND3",
    deps = [
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@com_google_protobuf//:protobuf_python",
    ],
)

py_library(
    name = "run_deeptrio_lib",
    srcs = ["run_deeptrio.py"],
    srcs_version = "PY2AND3",
    deps = [
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@com_google_protobuf//:protobuf_python",
    ],
)

py_library(
    name = "run_deepsomatic_lib",
    srcs = ["run_deepsomatic.py"],
    srcs_version = "PY2AND3",
    deps = [
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@com_google_protobuf//:protobuf_python",
    ],
)

py_test(
    name = "run_deepvariant_test",
    size = "small",
    srcs = ["run_deepvariant_test.py"],
    python_version = "PY3",
    srcs_version = "PY2AND3",
    deps = [
        ":run_deepvariant_lib",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "run_deeptrio_test",
    size = "small",
    srcs = ["run_deeptrio_test.py"],
    python_version = "PY3",
    srcs_version = "PY2AND3",
    deps = [
        ":run_deeptrio_lib",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "run_deepsomatic_test",
    size = "small",
    srcs = ["run_deepsomatic_test.py"],
    python_version = "PY3",
    srcs_version = "PY2AND3",
    deps = [
        ":run_deepsomatic_lib",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

# These binaries are for local dry-run development. They are not used in production.
py_binary(
    name = "run_deepvariant",
    srcs = ["run_deepvariant.py"],
    srcs_version = "PY3",
    deps = [
        ":run_deepvariant_lib",
    ],
)

py_binary(
    name = "run_deeptrio",
    srcs = ["run_deeptrio.py"],
    srcs_version = "PY3",
    deps = [
        ":run_deeptrio_lib",
    ],
)

py_binary(
    name = "run_deepsomatic",
    srcs = ["run_deepsomatic.py"],
    srcs_version = "PY3",
    deps = [
        ":run_deepsomatic_lib",
    ],
)
