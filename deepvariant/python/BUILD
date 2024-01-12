load("//third_party:clif.bzl", "py_clif_cc")
# Placeholder: load py_test

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_clif_cc(
    name = "allelecounter",
    srcs = ["allelecounter.clif"],
    clif_deps = [
        "//third_party/nucleus/io/python:reference",
    ],
    pyclif_deps = [
        "//third_party/nucleus/protos:reads_pyclif",
        "//third_party/nucleus/protos:range_pyclif",
        "//third_party/nucleus/protos:cigar_pyclif",
        "//deepvariant/protos:deepvariant_pyclif",
    ],
    deps = [
        "//deepvariant:allelecounter",
        "//third_party/nucleus/util:proto_clif_converter",
    ],
)

py_test(
    name = "allelecounter_wrap_test",
    size = "small",
    srcs = ["allelecounter_wrap_test.py"],
    data = [
        "//deepvariant:testdata",
        "//third_party/nucleus/testdata",
    ],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":allelecounter",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_clif_cc(
    name = "postprocess_variants",
    srcs = ["postprocess_variants.clif"],
    pyclif_deps = [
        "//third_party/nucleus/protos:reference_pyclif",
        "//deepvariant/protos:deepvariant_pyclif",
    ],
    deps = ["//deepvariant:postprocess_variants_lib"],
)

py_clif_cc(
    name = "variant_calling",
    srcs = ["variant_calling.clif"],
    clif_deps = [
        ":allelecounter",
        "//third_party/nucleus/io/python:reference",  # other py_clif_cc rules
        "//third_party/nucleus/io/python:vcf_reader",
    ],
    pyclif_deps = [
        "//deepvariant/protos:deepvariant_pyclif",
    ],
    deps = ["//deepvariant:variant_calling"],
)

py_test(
    name = "variant_calling_wrap_test",
    size = "small",
    srcs = ["variant_calling_wrap_test.py"],
    data = ["//deepvariant:testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":allelecounter",
        ":variant_calling",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_clif_cc(
    name = "direct_phasing",
    srcs = ["direct_phasing.clif"],
    deps = [
        "//deepvariant:direct_phasing",
        "//deepvariant/protos:deepvariant_pyclif",
        "//third_party/nucleus/protos:reads_pyclif",
        "//third_party/nucleus/util:proto_clif_converter",
    ],
)

py_clif_cc(
    name = "pileup_image_native",
    srcs = ["pileup_image_native.clif"],
    deps = [
        ":clif_converters",
        "//deepvariant:pileup_image_native",
        "//deepvariant/protos:deepvariant_pyclif",
        "//third_party/nucleus/protos:reads_pyclif",
        "//third_party/nucleus/util:proto_clif_converter",
    ],
)

cc_library(
    name = "clif_converters",
    srcs = ["clif_converters.cc"],
    hdrs = ["clif_converters.h"],
    deps = [
        "//deepvariant:pileup_image_native",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "@clif//:cpp_runtime",
        "@com_google_absl//absl/log:check",
        "@local_config_python//:numpy_headers",
    ],
)

py_clif_cc(
    name = "variant_calling_multisample",
    srcs = ["variant_calling_multisample.clif"],
    clif_deps = [
        ":allelecounter",
        "//third_party/nucleus/io/python:reference",  # other py_clif_cc rules
        "//third_party/nucleus/io/python:vcf_reader",
    ],
    pyclif_deps = [
        "//deepvariant/protos:deepvariant_pyclif",
    ],
    deps = [
        "//deepvariant:variant_calling_multisample",
        "//third_party/nucleus/util:proto_clif_converter",
    ],
)

py_test(
    name = "variant_calling_multisample_wrap_test",
    size = "small",
    srcs = ["variant_calling_multisample_wrap_test.py"],
    data = ["//deepvariant:testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":allelecounter",
        ":variant_calling_multisample",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)
