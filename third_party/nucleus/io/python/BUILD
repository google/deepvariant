load("//:third_party/clif.bzl", "py_clif_cc")

# load("//devtools/python/blaze:python3.bzl", "py2and3_test")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

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
    py_deps = [
        "//third_party/nucleus/io:clif_postproc",
    ],
    pyclif_deps = [
        "//third_party/nucleus/protos:range_pyclif",
        "//third_party/nucleus/protos:reference_pyclif",
        "//third_party/nucleus/protos:variants_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "vcf_reader_wrap_test",
    size = "small",
    srcs = ["vcf_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_reader",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "vcf_writer",
    srcs = ["vcf_writer.clif"],
    pyclif_deps = [
        "//third_party/nucleus/protos:variants_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:vcf_writer",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "vcf_writer_wrap_test",
    size = "small",
    srcs = ["vcf_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_writer",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:io_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "sam_reader",
    srcs = ["sam_reader.clif"],
    py_deps = [
        "//third_party/nucleus/io:clif_postproc",
    ],
    pyclif_deps = [
        "//third_party/nucleus/protos:range_pyclif",
        "//third_party/nucleus/protos:reads_pyclif",
        "//third_party/nucleus/protos:reference_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:sam_reader",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "sam_reader_wrap_test",
    size = "small",
    srcs = ["sam_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":sam_reader",
        "//third_party/nucleus/io:clif_postproc",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "reference",
    srcs = ["reference.clif"],
    pyclif_deps = [
        "//third_party/nucleus/protos:range_pyclif",
        "//third_party/nucleus/protos:reference_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_clif_cc(
    name = "reference_fai",
    srcs = ["reference_fai.clif"],
    clif_deps = [
        ":reference",
    ],
    deps = [
        "//third_party/nucleus/io:reference_fai",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "reference_wrap_test",
    size = "small",
    srcs = ["reference_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":reference_fai",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "fasta_reader",
    srcs = ["fasta_reader.clif"],
    pyclif_deps = [
        "//third_party/nucleus/protos:reference_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:fasta_reader",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_clif_cc(
    name = "bed_reader",
    srcs = ["bed_reader.clif"],
    py_deps = [
        "//third_party/nucleus/io:clif_postproc",
    ],
    pyclif_deps = [
        "//third_party/nucleus/protos:bed_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:bed_reader",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "bed_reader_wrap_test",
    size = "small",
    srcs = ["bed_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":bed_reader",
        "//third_party/nucleus/protos:bed_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "bed_writer",
    srcs = ["bed_writer.clif"],
    pyclif_deps = [
        "//third_party/nucleus/protos:bed_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:bed_writer",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "bed_writer_wrap_test",
    size = "small",
    srcs = ["bed_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":bed_writer",
        "//third_party/nucleus/protos:bed_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:io_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "fastq_reader",
    srcs = ["fastq_reader.clif"],
    py_deps = [
        "//third_party/nucleus/io:clif_postproc",
    ],
    pyclif_deps = [
        "//third_party/nucleus/protos:fastq_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:fastq_reader",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "fastq_reader_wrap_test",
    size = "small",
    srcs = ["fastq_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":fastq_reader",
        "//third_party/nucleus/protos:fastq_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "fastq_writer",
    srcs = ["fastq_writer.clif"],
    pyclif_deps = [
        "//third_party/nucleus/protos:fastq_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:fastq_writer",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "fastq_writer_wrap_test",
    size = "small",
    srcs = ["fastq_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":fastq_writer",
        "//third_party/nucleus/io:fastq",
        "//third_party/nucleus/protos:fastq_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:io_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "hts_verbose",
    srcs = ["hts_verbose.clif"],
    deps = [
        "//third_party/nucleus/io:hts_verbose",
    ],
)

py_test(
    name = "hts_verbose_test",
    size = "small",
    srcs = ["hts_verbose_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":hts_verbose",
        "@io_abseil_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "gff_reader",
    srcs = ["gff_reader.clif"],
    py_deps = [
        "//third_party/nucleus/io:clif_postproc",
    ],
    pyclif_deps = [
        "//third_party/nucleus/protos:gff_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:gff_reader",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "gff_reader_wrap_test",
    size = "small",
    srcs = ["gff_reader_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":gff_reader",
        "//third_party/nucleus/protos:gff_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_clif_cc(
    name = "gff_writer",
    srcs = ["gff_writer.clif"],
    py_deps = [
        "//third_party/nucleus/io:clif_postproc",
    ],
    pyclif_deps = [
        "//third_party/nucleus/protos:gff_pyclif",
    ],
    deps = [
        "//third_party/nucleus/io:gff_writer",
        "//third_party/nucleus/vendor:statusor_clif_converters",
    ],
)

py_test(
    name = "gff_writer_wrap_test",
    size = "small",
    srcs = ["gff_writer_wrap_test.py"],
    data = ["//third_party/nucleus/testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":gff_writer",
        "//third_party/nucleus/protos:gff_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)
