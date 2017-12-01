load("//:third_party/clif.bzl", "py_clif_cc")

package(
    default_visibility = [
        "//visibility:public",
        # Uncomment and comment out the above line to replace the public visibilily
        # of these low-level IO libraries with more restricted access. May become the
        # default in the future, but is useful now to flag uses of these libraries.
        # "//deepvariant/core:__subpackages__",
        # "//deepvariant/python:__subpackages__",
    ],
)

py_clif_cc(
    name = "vcf_reader",
    srcs = ["vcf_reader.clif"],
    clif_deps = [],
    py_deps = [
        "//deepvariant/core:clif_postproc",
    ],
    pyclif_deps = [
        "//deepvariant/core/protos:core_pyclif",
        "//deepvariant/core/genomics:range_pyclif",
        "//deepvariant/core/genomics:variants_pyclif",
    ],
    visibility = ["//deepvariant/core:__subpackages__"],
    deps = [
        "//deepvariant/core:vcf_reader",
        "//deepvariant/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "vcf_reader_wrap_test",
    size = "small",
    srcs = ["vcf_reader_wrap_test.py"],
    data = [
        "//deepvariant/core:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_reader",
        "//deepvariant/core:py_test_utils",
        "//deepvariant/core:ranges",
        "//deepvariant/core/protos:core_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "vcf_writer",
    srcs = ["vcf_writer.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [
        "//deepvariant/core/genomics:variants_pyclif",
        "//deepvariant/core/protos:core_pyclif",
    ],
    deps = [
        "//deepvariant/core:vcf_writer",
        "//deepvariant/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "vcf_writer_wrap_test",
    size = "small",
    srcs = ["vcf_writer_wrap_test.py"],
    data = [
        "//deepvariant/core:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_writer",
        "//deepvariant/core:genomics_io",
        "//deepvariant/core:io_utils",
        "//deepvariant/core:py_test_utils",
        "//deepvariant/core/genomics:variants_py_pb2",
        "//deepvariant/core/protos:core_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "sam_reader",
    srcs = ["sam_reader.clif"],
    clif_deps = [],  # other py_clif_cc rules
    py_deps = [
        "//deepvariant/core:clif_postproc",
    ],
    pyclif_deps = [
        "//deepvariant/core/protos:core_pyclif",
        "//deepvariant/core/genomics:range_pyclif",
        "//deepvariant/core/genomics:reads_pyclif",
    ],
    visibility = ["//deepvariant/core:__subpackages__"],
    deps = [
        "//deepvariant/core:sam_reader",
        "//deepvariant/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "sam_reader_wrap_test",
    size = "small",
    srcs = ["sam_reader_wrap_test.py"],
    data = [
        "//deepvariant/core:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":sam_reader",
        "//deepvariant/core:clif_postproc",
        "//deepvariant/core:py_test_utils",
        "//deepvariant/core:ranges",
        "//deepvariant/core/protos:core_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "reference_fai",
    srcs = ["reference_fai.clif"],
    pyclif_deps = [
        "//deepvariant/core/genomics:range_pyclif",
        "//deepvariant/core/protos:core_pyclif",
    ],
    deps = [
        "//deepvariant/core:reference_fai",
        "//deepvariant/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "reference_wrap_test",
    size = "small",
    srcs = ["reference_wrap_test.py"],
    data = [
        "//deepvariant/core:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":reference_fai",
        "//deepvariant/core:py_test_utils",
        "//deepvariant/core:ranges",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "math",
    srcs = ["math.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [],
    visibility = ["//deepvariant/core:__subpackages__"],
    deps = ["//deepvariant/core:cpp_math"],
)

py_test(
    name = "math_wrap_test",
    size = "small",
    srcs = ["math_wrap_test.py"],
    data = [],
    srcs_version = "PY2AND3",
    deps = [
        ":math",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "hts_verbose",
    srcs = ["hts_verbose.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [],
    visibility = [
        "//deepvariant:__subpackages__",
        "//deepvariant/core:__subpackages__",
        "//internal:__subpackages__",
    ],
    deps = [
        "//deepvariant/core:hts_verbose",
    ],
)

py_test(
    name = "hts_verbose_test",
    size = "small",
    srcs = ["hts_verbose_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":hts_verbose",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)
