# Example programs for demonstrating Nucleus's APIs.

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

py_binary(
    name = "add_ad_to_vcf",
    srcs = ["add_ad_to_vcf.py"],
    deps = [
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:variantcall_utils",
        "//deepvariant/util:vcf_constants",
        "//deepvariant/util/io:vcf",
    ],
)

# redacted

py_binary(
    name = "ascii_pileup",
    srcs = ["ascii_pileup.py"],
    deps = [
        "//deepvariant/util:ranges",
        "//deepvariant/util/io:sam",
    ],
)

py_binary(
    name = "count_variants",
    srcs = ["count_variants.py"],
    deps = [
        "//deepvariant/util:variant_utils",
        "//deepvariant/util/io:vcf",
    ],
)

py_binary(
    name = "validate_vcf",
    srcs = ["validate_vcf.py"],
    deps = [
        "//deepvariant/util:variant_utils",
        "//deepvariant/util/io:fasta",
        "//deepvariant/util/io:vcf",
    ],
)

# redacted

py_binary(
    name = "ngs_errors",
    srcs = ["ngs_errors.py"],
    deps = [
        "//deepvariant/util:cigar",
        "//deepvariant/util:py_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util/genomics:sam_py_pb2",
        "//deepvariant/util/io:fasta",
        "//deepvariant/util/io:genomics_writer",
        "//deepvariant/util/io:sam",
        "//deepvariant/util/io:vcf",
    ],
)
