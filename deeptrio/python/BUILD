load("//third_party:clif.bzl", "py_clif_cc")
# KEEP THIS LINE to import nucleus.bzl externally

package(
    default_visibility = [
        "//visibility:public",
    ],
)

py_clif_cc(
    name = "variant_calling_deeptrio",
    srcs = ["variant_calling_deeptrio.clif"],
    clif_deps = [
        "//deepvariant/python:allelecounter",
        "//third_party/nucleus/io/python:reference",  # other py_clif_cc rules
        "//third_party/nucleus/io/python:vcf_reader",
    ],
    pyclif_deps = [
        "//deepvariant/protos:deepvariant_pyclif",
    ],
    deps = [
        "//deeptrio:variant_calling_deeptrio",
        "//third_party/nucleus/util:proto_clif_converter",
    ],
)

py_test(
    name = "variant_calling_deeptrio_wrap_test",
    size = "small",
    srcs = ["variant_calling_deeptrio_wrap_test.py"],
    data = ["//deepvariant:testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":variant_calling_deeptrio",
        "//deepvariant:py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
    ],
)
