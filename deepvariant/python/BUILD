load("//:third_party/clif.bzl", "py_clif_cc")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

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
        "//deepvariant/protos:deepvariant_pyclif",
    ],
    deps = ["//deepvariant:allelecounter"],
)

py_test(
    name = "allelecounter_wrap_test",
    size = "small",
    srcs = ["allelecounter_wrap_test.py"],
    data = [
        "//deepvariant:testdata",
        "//third_party/nucleus/testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":allelecounter",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
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
        "//third_party/nucleus/io/python:reference_fai",  # other py_clif_cc rules
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
    srcs_version = "PY2AND3",
    deps = [
        ":allelecounter",
        ":variant_calling",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "pileup_image_native",
    srcs = ["pileup_image_native.clif"],
    deps = [
        "//deepvariant:pileup_image_native",
        "//deepvariant/protos:deepvariant_pyclif",
        "//deepvariant/python:clif_converters",
        "//third_party/nucleus/protos:reads_pyclif",
    ],
)

cc_library(
    name = "clif_converters",
    srcs = ["clif_converters.cc"],
    hdrs = ["clif_converters.h"],
    deps = [
        "//deepvariant:pileup_image_native",
        "@clif//:cpp_runtime",
        "@local_config_python//:numpy_headers",
        "@local_config_python//:python_headers",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)
