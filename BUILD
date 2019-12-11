load("//third_party/nucleus/tools:zip_dir.bzl", "zip_dir")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

test_suite(
    name = "smoke_tests",
    tests = [
        "//deepvariant/core:smoke_tests",
        "//deepvariant/environment_tests:smoke_tests",
        "//deepvariant/testing:smoke_tests",
    ],
)

filegroup(
    name = "binaries",
    srcs = [
        "//deepvariant:binaries",
    ],
)

exports_files(["LICENSE"])

filegroup(
    name = "licenses",
    srcs = [
        ":LICENSE",
        "//third_party:abseil_cpp.LICENSE",  # redacted
        "//third_party:boost.LICENSE",
        "//third_party:tensorflow_models.LICENSE",  # redacted
        "@com_google_protobuf//:LICENSE",
        "@com_googlesource_code_re2//:LICENSE",
        "@htslib//:LICENSE",
        "@libssw//:README.md",  # SSW license embedded in the README.
        "@org_tensorflow//:LICENSE",
    ],
)

zip_dir(
    name = "licenses_zip",
    srcs = [":licenses"],
    zipname = "licenses.zip",
)

cc_library(
    name = "all_extensions",
    srcs = [],
    deps = [
        "//deepvariant/python:allelecounter_cclib",
        "//deepvariant/python:pileup_image_native_cclib",
        "//deepvariant/python:postprocess_variants_cclib",
        "//deepvariant/python:variant_calling_cclib",
        "//deepvariant/realigner/python:debruijn_graph_cclib",
        "//deepvariant/realigner/python:fast_pass_aligner_cclib",
        "//deepvariant/realigner/python:ssw_cclib",
        "//deepvariant/realigner/python:window_selector_cclib",
        "//third_party/nucleus/io/python:bed_reader_cclib",
        "//third_party/nucleus/io/python:bed_writer_cclib",
        "//third_party/nucleus/io/python:bedgraph_reader_cclib",
        "//third_party/nucleus/io/python:bedgraph_writer_cclib",
        "//third_party/nucleus/io/python:fastq_reader_cclib",
        "//third_party/nucleus/io/python:fastq_writer_cclib",
        "//third_party/nucleus/io/python:gff_reader_cclib",
        "//third_party/nucleus/io/python:gff_writer_cclib",
        "//third_party/nucleus/io/python:gfile_cclib",
        "//third_party/nucleus/io/python:hts_verbose_cclib",
        "//third_party/nucleus/io/python:reference_cclib",
        "//third_party/nucleus/io/python:sam_reader_cclib",
        "//third_party/nucleus/io/python:sam_writer_cclib",
        "//third_party/nucleus/io/python:tabix_indexer_cclib",
        "//third_party/nucleus/io/python:tfrecord_reader_cclib",
        "//third_party/nucleus/io/python:tfrecord_writer_cclib",
        "//third_party/nucleus/io/python:vcf_reader_cclib",
        "//third_party/nucleus/io/python:vcf_writer_cclib",
        "//third_party/nucleus/util/python:math_cclib",
        "//third_party/nucleus/util/python:utils_cclib",
        "//third_party/nucleus/vendor/python:statusor_examples_cclib",
    ],
)

# Until https://github.com/bazelbuild/bazel/issues/4815 is fixed,
# we need to specify and force the use of a py_runtime.
# Python3.6 (at least 3.5.3) is needed because of altair. See this bug:
# https://github.com/altair-viz/altair/issues/972
# Ubuntu16.04 default is Python 3.5.2, which is why we need to specifically
# update the Python version.
py_runtime(
    name = "deepvariant_python_runtime",
    files = [],
    interpreter_path = select({
        "@bazel_tools//tools/python:PY3": "/usr/bin/python3.6",
    }),
)
