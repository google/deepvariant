# Example programs for demonstrating Nucleus's APIs.

package(
    default_visibility = [
        "//visibility:public",
    ],
)

licenses(["notice"])

exports_files(["LICENSE"])

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
