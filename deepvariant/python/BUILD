load("//:third_party/clif.bzl", "py_clif_cc")

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_clif_cc(
    name = "allelecounter",
    srcs = ["allelecounter.clif"],
    clif_deps = [
        "//deepvariant/util/python:reference_fai",  # other py_clif_cc rules
    ],
    py_deps = [],
    pyclif_deps = [
        "//deepvariant/util/genomics:reads_pyclif",
        "//deepvariant/util/genomics:range_pyclif",
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
        "//deepvariant/util:testdata",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":allelecounter",
        "//deepvariant:py_test_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/util:genomics_io",
        "//deepvariant/util:ranges",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "postprocess_variants",
    srcs = ["postprocess_variants.clif"],
    clif_deps = [
    ],
    py_deps = [],
    pyclif_deps = [
        "//deepvariant/util/genomics:variants_pyclif",
        "//deepvariant/util/protos:core_pyclif",
        "//deepvariant/protos:deepvariant_pyclif",
    ],
    deps = ["//deepvariant:postprocess_variants_lib"],
)

py_clif_cc(
    name = "variant_calling",
    srcs = ["variant_calling.clif"],
    clif_deps = [
        ":allelecounter",
        "//deepvariant/util/python:reference_fai",  # other py_clif_cc rules
    ],
    py_deps = [],
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
        "//deepvariant:py_test_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/util:genomics_io",
        "//deepvariant/util:ranges",
        "@com_google_absl_py//absl/testing:absltest",
    ],
)

py_clif_cc(
    name = "pileup_image_native",
    srcs = ["pileup_image_native.clif"],
    clif_deps = [],
    py_deps = [],
    pyclif_deps = [],
    deps = [
        "//deepvariant:pileup_image_native",
        "//deepvariant/protos:deepvariant_pyclif",
        "//deepvariant/python:clif_converters",
        "//deepvariant/util/genomics:reads_pyclif",
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
