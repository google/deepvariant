load("//:third_party/clif.bzl", "py_clif_cc")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

py_clif_cc(
    name = "vcf_reader",
    srcs = ["vcf_reader.clif"],
    clif_deps = [],
    py_deps = [
        "//deepvariant/util:clif_postproc",
    ],
    pyclif_deps = [
        "//deepvariant/util/genomics:range_pyclif",
        "//deepvariant/util/genomics:reference_pyclif",
        "//deepvariant/util/genomics:variants_pyclif",
        "//deepvariant/util/genomics:vcf_pyclif",
    ],
    deps = [
        "//deepvariant/util:vcf_reader",
        "//deepvariant/util/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "vcf_reader_wrap_test",
    size = "small",
    srcs = ["vcf_reader_wrap_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_reader",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/genomics:vcf_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "vcf_writer",
    srcs = ["vcf_writer.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [
        "//deepvariant/util/genomics:variants_pyclif",
        "//deepvariant/util/genomics:vcf_pyclif",
    ],
    deps = [
        "//deepvariant/util:vcf_writer",
        "//deepvariant/util/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "vcf_writer_wrap_test",
    size = "small",
    srcs = ["vcf_writer_wrap_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_writer",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/genomics:vcf_py_pb2",
        "//deepvariant/util/io:vcf",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "sam_reader",
    srcs = ["sam_reader.clif"],
    clif_deps = [],  # other py_clif_cc rules
    py_deps = [
        "//deepvariant/util:clif_postproc",
    ],
    pyclif_deps = [
        "//deepvariant/util/genomics:range_pyclif",
        "//deepvariant/util/genomics:reads_pyclif",
        "//deepvariant/util/genomics:reference_pyclif",
        "//deepvariant/util/genomics:sam_pyclif",
    ],
    deps = [
        "//deepvariant/util:sam_reader",
        "//deepvariant/util/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "sam_reader_wrap_test",
    size = "small",
    srcs = ["sam_reader_wrap_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":sam_reader",
        "//deepvariant/util:clif_postproc",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:index_py_pb2",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:sam_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "reference_fai",
    srcs = ["reference_fai.clif"],
    pyclif_deps = [
        "//deepvariant/util/genomics:range_pyclif",
        "//deepvariant/util/genomics:reference_pyclif",
    ],
    deps = [
        "//deepvariant/util:reference_fai",
        "//deepvariant/util/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "reference_wrap_test",
    size = "small",
    srcs = ["reference_wrap_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":reference_fai",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util:ranges",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "fastq_reader",
    srcs = ["fastq_reader.clif"],
    clif_deps = [],
    py_deps = [
        "//deepvariant/util:clif_postproc",
    ],
    pyclif_deps = [
        "//deepvariant/util/genomics:fastq_pyclif",
    ],
    deps = [
        "//deepvariant/util:fastq_reader",
        "//deepvariant/util/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "fastq_reader_wrap_test",
    size = "small",
    srcs = ["fastq_reader_wrap_test.py"],
    data = [
        "//deepvariant/util:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":fastq_reader",
        "//deepvariant/util:py_test_utils",
        "//deepvariant/util/genomics:fastq_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "math",
    srcs = ["math.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [],
    deps = ["//deepvariant/util:cpp_math"],
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
    deps = [
        "//deepvariant/util:hts_verbose",
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
