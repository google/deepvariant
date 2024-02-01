# Placeholder: load py_binary
# Placeholder: load py_library
# Placeholder: load py_test

package(
    default_visibility = [
        "//visibility:public",
    ],
)

filegroup(
    name = "binaries",
    srcs = [
        "make_examples",
        "//deepvariant:call_variants",
        "//deepvariant:postprocess_variants",
        "//deepvariant:vcf_stats_report",
    ],
)

py_library(
    name = "dt_constants",
    srcs = ["dt_constants.py"],
    srcs_version = "PY2AND3",
)

py_library(
    name = "py_testdata",
    testonly = True,
    srcs = ["testdata.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//third_party/nucleus/testing:py_test_utils",
    ],
)

py_library(
    name = "make_examples_lib",
    srcs = ["make_examples.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":dt_constants",
        "//deepvariant:logging_level",
        "//deepvariant:make_examples_core",
        "//deepvariant:make_examples_options",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
    ],
)

py_binary(
    name = "make_examples",
    srcs = [
        "make_examples.py",
    ],
    main = "make_examples.py",
    python_version = "PY3",
    deps = [
        ":dt_constants",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant:logging_level",
        "//deepvariant:make_examples_core",
        "//deepvariant:make_examples_options",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl/flags",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@absl_py//absl:app",
    ],
)

py_test(
    name = "make_examples_test",
    size = "large",
    srcs = ["make_examples_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    shard_count = 2,
    srcs_version = "PY2AND3",
    deps = [
        ":make_examples_lib",
        ":py_testdata",
        "//deepvariant:dv_constants",
        "//deepvariant:dv_utils",
        "//deepvariant:make_examples_core",
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)
