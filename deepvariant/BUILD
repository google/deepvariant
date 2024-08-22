# Placeholder: load py_binary
# Placeholder: load py_library
# Placeholder: load py_test

package(
    default_visibility = [
        "//visibility:public",
    ],
)

filegroup(
    name = "binaries",
    srcs = [
        "call_variants",
        "make_examples",
        "make_examples_pangenome_aware_dv",
        "make_examples_somatic",
        "multisample_make_examples",
        "runtime_by_region_vis",
        "show_examples",
        "train",
        "vcf_stats_report",
    ],
)

filegroup(
    name = "testdata",
    testonly = True,
    srcs = glob([
        "testdata/input/*",
        "testdata/obsolete/*",
        "testdata/*",
        "multiallelic_model/**",
    ]),
)

test_suite(
    name = "gpu_tests",
    tags = [
        "gpu_test",
        "manual",
    ],
    tests = [
        ":call_variants_accelerator_test",
    ],
)

# Umbrella library for all the Python under DeepVariant.
py_library(
    name = "deepvariant_py",
    srcs = [],
    deps = [
        ":allele_frequency",
        ":call_variants_lib",
        ":data_providers",
        ":dv_utils",
        ":dv_vcf_constants",
        ":haplotypes",
        ":keras_modeling",
        ":make_examples_core",
        ":make_examples_lib",
        ":pileup_image",
        "//deepvariant/python:allelecounter",
        "//deepvariant/python:direct_phasing",
        "//deepvariant/python:variant_calling",
        "//deepvariant/realigner",
        "//third_party/nucleus/io/python:merge_variants",
    ],
)

cc_library(
    name = "allelecounter",
    srcs = ["allelecounter.cc"],
    hdrs = ["allelecounter.h"],
    deps = [
        ":utils",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "allelecounter_test",
    size = "small",
    srcs = ["allelecounter_test.cc"],
    data = ["//third_party/nucleus/testdata"],
    deps = [
        ":allelecounter",
        ":testing_utils",
        ":utils",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core:statusor",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "postprocess_variants_lib",
    srcs = ["postprocess_variants.cc"],
    hdrs = ["postprocess_variants.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/strings:cord",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core/platform/cloud:gcs_file_system",
    ],
)

cc_test(
    name = "postprocess_variants_lib_test",
    size = "small",
    srcs = ["postprocess_variants_test.cc"],
    data = [
        ":testdata",
    ],
    deps = [
        ":postprocess_variants_lib",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "variant_calling",
    srcs = ["variant_calling.cc"],
    hdrs = ["variant_calling.h"],
    deps = [
        ":allelecounter",
        ":utils",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core:statusor",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:samplers",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "variant_calling_test",
    size = "small",
    srcs = ["variant_calling_test.cc"],
    data = [":testdata"],
    deps = [
        ":utils",
        ":variant_calling",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest_main",
        "@com_google_protobuf//:protobuf",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "variant_calling_multisample",
    srcs = ["variant_calling_multisample.cc"],
    hdrs = ["variant_calling_multisample.h"],
    deps = [
        ":allelecounter",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:samplers",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_test(
    name = "variant_calling_multisample_test",
    size = "small",
    srcs = ["variant_calling_multisample_test.cc"],
    deps = [
        ":utils",
        ":variant_calling_multisample",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_test(
    name = "variant_calling_multisample_somatic_test",
    size = "small",
    srcs = ["variant_calling_multisample_somatic_test.cc"],
    deps = [
        ":utils",
        ":variant_calling_multisample",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_test(
    name = "variant_calling_multisample_trio_test",
    size = "small",
    srcs = ["variant_calling_multisample_trio_test.cc"],
    deps = [
        ":allelecounter",
        ":utils",
        ":variant_calling_multisample",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "variant_caller",
    srcs = ["variant_caller.py"],
    srcs_version = "PY3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/python:variant_calling",
        "//deepvariant/python:variant_calling_multisample",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
    ],
)

py_test(
    name = "variant_caller_test",
    size = "small",
    srcs = ["variant_caller_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_testdata",
        ":variant_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "very_sensitive_caller",
    srcs = ["very_sensitive_caller.py"],
    srcs_version = "PY3",
    deps = [
        ":variant_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
    ],
)

py_test(
    name = "very_sensitive_caller_test",
    size = "small",
    srcs = ["very_sensitive_caller_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_testdata",
        ":very_sensitive_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "vcf_candidate_importer",
    srcs = ["vcf_candidate_importer.py"],
    srcs_version = "PY3",
    deps = [
        ":variant_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//third_party/nucleus/io:vcf",
    ],
)

py_test(
    name = "vcf_candidate_importer_test",
    size = "small",
    srcs = ["vcf_candidate_importer_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_testdata",
        ":vcf_candidate_importer",
        "//deepvariant/labeler:labeled_examples_to_vcf_main_lib",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "variant_caller_trio_test",
    size = "small",
    srcs = ["variant_caller_trio_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":variant_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "very_sensitive_caller_trio_test",
    size = "small",
    srcs = ["very_sensitive_caller_trio_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":very_sensitive_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "utils",
    srcs = ["utils.cc"],
    hdrs = ["utils.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "utils_test",
    size = "small",
    srcs = ["utils_test.cc"],
    deps = [
        ":utils",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "make_examples_lib",
    srcs = ["make_examples.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":logging_level",
        ":make_examples_core",
        ":make_examples_options",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/python:direct_phasing",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
    ],
)

cc_library(
    name = "direct_phasing",
    srcs = ["direct_phasing.cc"],
    hdrs = ["direct_phasing.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/core:statusor",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_test(
    name = "direct_phasing_test",
    size = "small",
    srcs = ["direct_phasing_test.cc"],
    deps = [
        ":direct_phasing",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "make_examples_core",
    srcs = ["make_examples_core.py"],
    srcs_version = "PY3",
    deps = [
        ":allele_frequency",
        ":calling_regions_utils",
        ":dv_constants",
        ":dv_utils",
        ":dv_utils_using_clif",
        ":dv_vcf_constants",
        ":pileup_image",
        ":resources_main_lib",
        ":sample",
        ":variant_caller",
        ":vcf_candidate_importer",
        ":very_sensitive_caller",
        "//deepvariant/labeler:customized_classes_labeler",
        "//deepvariant/labeler:haplotype_labeler",
        "//deepvariant/labeler:positional_labeler",
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/python:direct_phasing",
        "//deepvariant/python:make_examples_native",
        "//deepvariant/realigner",
        "//deepvariant/small_model:inference",
        "//deepvariant/small_model:make_small_model_examples",
        "//deepvariant/vendor:timer",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/protos:range_py_pb2",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:struct_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/logging",
        "@com_google_protobuf//:protobuf_python",
    ],
)

py_test(
    name = "make_examples_core_test",
    size = "large",
    srcs = ["make_examples_core_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    shard_count = 2,
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":make_examples_core",
        ":make_examples_lib",
        ":py_testdata",
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/protos:realigner_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "make_examples_options",
    srcs = ["make_examples_options.py"],
    srcs_version = "PY3",
    deps = [
        ":calling_regions_utils",
        ":dv_constants",
        ":dv_utils",
        ":exclude_contigs",
        ":make_examples_core",
        ":pileup_image",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/python:direct_phasing",
        "//deepvariant/realigner",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:errors",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "make_examples",
    srcs = [
        "make_examples.py",
    ],
    main = "make_examples.py",
    python_version = "PY3",
    deps = [
        ":dv_constants",
        ":logging_level",
        ":make_examples_core",
        ":make_examples_options",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl/flags",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@absl_py//absl:app",
    ],
)

cc_library(
    name = "merge_phased_reads_lib",
    srcs = [
        "merge_phased_reads.cc",
    ],
    hdrs = [
        "merge_phased_reads.h",
    ],
    deps = [
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_google_absl//absl/types:span",
        "@com_googlesource_code_re2//:re2",
    ],
)

cc_binary(
    name = "merge_phased_reads",
    srcs = [
        "merge_phased_reads_main.cc",
    ],
    deps = [
        ":merge_phased_reads_lib",
        "@com_google_absl//absl/flags:flag",
        "@com_google_absl//absl/log:check",
    ],
)

cc_test(
    name = "merge_phased_reads_test",
    size = "small",
    srcs = [
        "merge_phased_reads_test.cc",
    ],
    deps = [
        ":merge_phased_reads_lib",
        "@com_google_absl//absl/hash:hash_testing",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "make_examples_somatic_lib",
    srcs = [
        "make_examples_somatic.py",
    ],
    srcs_version = "PY2AND3",
    deps = [
        ":dv_constants",
        ":logging_level",
        ":make_examples_core",
        ":make_examples_options",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl/flags",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@absl_py//absl:app",
    ],
)

py_binary(
    name = "make_examples_somatic",
    srcs = [
        "make_examples_somatic.py",
    ],
    main = "make_examples_somatic.py",
    python_version = "PY3",
    deps = [
        ":dv_constants",
        ":logging_level",
        ":make_examples_core",
        ":make_examples_options",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl/flags",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@absl_py//absl:app",
    ],
)

py_binary(
    name = "make_examples_pangenome_aware_dv",
    srcs = [
        "make_examples_pangenome_aware_dv.py",
    ],
    main = "make_examples_pangenome_aware_dv.py",
    python_version = "PY3",
    deps = [
        ":dv_constants",
        ":logging_level",
        ":make_examples_core",
        ":make_examples_options",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl/flags",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@absl_py//absl:app",
    ],
)

py_test(
    name = "make_examples_somatic_test",
    srcs = ["make_examples_somatic_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":make_examples_somatic_lib",
        ":py_testdata",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "allele_frequency",
    srcs = ["allele_frequency.py"],
    srcs_version = "PY3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "allele_frequency_test",
    srcs = ["allele_frequency_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    deps = [
        ":allele_frequency",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "calling_regions_utils",
    srcs = ["calling_regions_utils.py"],
    srcs_version = "PY3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/protos:range_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "calling_regions_utils_test",
    size = "small",
    srcs = ["calling_regions_utils_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":calling_regions_utils",
        "//third_party/nucleus/protos:range_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "make_examples_test",
    size = "large",
    srcs = ["make_examples_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    shard_count = 2,
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":dv_utils",
        ":make_examples_core",
        ":make_examples_lib",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "//third_party/nucleus/util:vis",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "keras_modeling",
    srcs = ["keras_modeling.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":dv_utils",
        ":dv_utils_using_clif",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "keras_modeling_test",
    size = "large",
    srcs = ["keras_modeling_test.py"],
    deps = [
        ":dv_config",
        ":dv_constants",
        ":keras_modeling",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_binary(
    name = "call_variants",
    srcs = ["call_variants.py"],
    data = [
       "examples_from_stream.so",
    ],
    python_version = "PY3",
    deps = [
        "call_variants_lib",
        ":dv_config",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_binary(
    name = "train",
    srcs = [
        "train.py",
    ],
    main = "train.py",
    python_version = "PY3",
    deps = [
        ":data_providers",
        ":dv_constants",
        ":dv_utils",
        ":keras_modeling",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

# py_test(
#     name = "train_test",
#     size = "large",
#     srcs = ["train_test.py"],
#     data = [":testdata"],
#     python_version = "PY3",
#     srcs_version = "PY3",
#     deps = [
#         ":dv_config",
#         ":py_testdata",
#         ":train",
#         "//third_party/nucleus/testing:py_test_utils",
#         "@absl_py//absl/testing:flagsaver",
#         "@absl_py//absl/testing:parameterized",
#         "//third_party/py/tensorflow:tensorflow_no_contrib",
#     ],
# )

py_binary(
    name = "train_inceptionv3",
    srcs = [
        "train_inceptionv3.py",
    ],
    main = "train_inceptionv3.py",
    deps = [
        ":average_model_checkpoint_patched",
        ":data_providers",
        ":dv_constants",
        ":dv_utils",
        ":keras_modeling",
        "//third_party/nucleus/io:sharded_file_utils",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "call_variants_lib",
    srcs = ["call_variants.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_utils",
        ":keras_modeling",
        ":logging_level",
        "@bazel_tools//tools/python/runfiles",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:vis",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "call_variants_test",
    size = "large",
    srcs = ["call_variants_test.py"],
    data = [":testdata"],
    deps = [
        ":call_variants_lib",
        ":dv_utils",
        ":keras_modeling",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tfrecord",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",  # buildcleaner: keep
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "pileup_image",
    srcs = ["pileup_image.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":sample",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:pileup_image_native",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "pileup_image_test",
    size = "small",
    srcs = ["pileup_image_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":make_examples_core",
        ":pileup_image",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:pileup_image_native",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "pileup_image_native",
    srcs = ["pileup_image_native.cc"],
    hdrs = ["pileup_image_native.h"],
    deps = [
        ":pileup_channel_lib",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "pileup_channel_lib",
    srcs = ["pileup_channel_lib.cc"],
    hdrs = ["pileup_channel_lib.h"],
    deps = [
        ":allele_frequency_channel",
        ":avg_base_quality_channel",
        ":blank_channel",
        ":channel",
        ":gap_compressed_identity_channel",
        ":gc_content_channel",
        ":haplotype_tag_channel",
        ":homopolymer_weighted_channel",
        ":identity_channel",
        ":insert_size_channel",
        ":is_homopolymer_channel",
        ":mapping_quality_channel",
        ":read_mapping_percent_channel",
        ":read_supports_variant_channel",
        ":strand_channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
    ],
)

cc_test(
    name = "pileup_channel_lib_test",
    srcs = ["pileup_channel_lib_test.cc"],
    deps = [
        ":allele_frequency_channel",
        ":avg_base_quality_channel",
        ":blank_channel",
        ":channel",
        ":gap_compressed_identity_channel",
        ":gc_content_channel",
        ":haplotype_tag_channel",
        ":homopolymer_weighted_channel",
        ":identity_channel",
        ":insert_size_channel",
        ":is_homopolymer_channel",
        ":mapping_quality_channel",
        ":pileup_channel_lib",
        ":read_mapping_percent_channel",
        ":read_supports_variant_channel",
        ":strand_channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "postprocess_variants_py_lib",
    srcs = ["postprocess_variants.py"],
    data = glob([
        "multiallelic_model/**",
    ]),
    srcs_version = "PY3",
    deps = [
        # END_INTERNAL
        ":dv_vcf_constants",
        ":logging_level",
        ":dv_utils",
        ":dv_constants",
        ":haplotypes",
        ":calling_regions_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:postprocess_variants",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tabix",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/io/python:merge_variants",
        "//third_party/nucleus/io/python:vcf_concat",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:struct_utils",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "postprocess_variants",
    srcs = [
        "postprocess_variants.py",
    ],
    main = "postprocess_variants.py",
    python_version = "PY3",
    deps = [
        ":postprocess_variants_py_lib",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_test(
    name = "postprocess_variants_test",
    size = "medium",
    srcs = ["postprocess_variants_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    tags = ["requires-net:external"],
    deps = [
        ":dv_constants",
        ":dv_vcf_constants",
        ":postprocess_variants_py_lib",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:range_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:struct_utils",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "py_testdata",
    testonly = True,
    srcs = ["testdata.py"],
    srcs_version = "PY3",
    deps = [
        "//third_party/nucleus/testing:py_test_utils",
    ],
)

py_library(
    name = "data_providers",
    srcs = ["data_providers.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:sharded_file_utils",
        "@com_google_protobuf//:protobuf_python",
    ],
)

py_library(
    name = "data_providers_test_lib",
    testonly = True,
    srcs = ["data_providers_test.py"],
    srcs_version = "PY3",
    deps = [
        ":data_providers",
        ":dv_config",
        ":dv_constants",
        ":dv_utils",
        ":py_testdata",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "data_providers_test",
    size = "medium",
    srcs = ["data_providers_test.py"],
    data = [":testdata"],
    main = "data_providers_test.py",
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":data_providers_test_lib",
    ],
)

py_library(
    name = "haplotypes",
    srcs = ["haplotypes.py"],
    srcs_version = "PY3",
    deps = [
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "haplotypes_test",
    size = "small",
    srcs = ["haplotypes_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":haplotypes",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "dv_utils",
    srcs = ["dv_utils.py"],
    srcs_version = "PY3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "dv_utils_test",
    size = "small",
    srcs = ["dv_utils_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":dv_utils",
        ":py_testdata",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "dv_utils_using_clif",
    srcs = ["dv_utils_using_clif.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
    ],
)

py_test(
    name = "dv_utils_using_clif_test",
    size = "small",
    srcs = ["dv_utils_using_clif_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":dv_utils_using_clif",
        ":py_testdata",
        "//third_party/nucleus/protos:variants_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "logging_level",
    srcs = ["logging_level.py"],
    deps = ["@absl_py//absl/flags"],
)

py_library(
    name = "exclude_contigs",
    srcs = ["exclude_contigs.py"],
)

py_test(
    name = "exclude_contigs_test",
    size = "small",
    srcs = ["exclude_contigs_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":exclude_contigs",
        "@absl_py//absl/testing:absltest",
    ],
)

py_library(
    name = "dv_constants",
    srcs = ["dv_constants.py"],
    srcs_version = "PY3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
    ],
)

py_library(
    name = "dv_vcf_constants",
    srcs = ["dv_vcf_constants.py"],
    srcs_version = "PY3",
    deps = [
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:vcf_constants",
    ],
)

py_test(
    name = "dv_vcf_constants_test",
    size = "small",
    srcs = ["dv_vcf_constants_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":dv_vcf_constants",
        "//third_party/nucleus/protos:reference_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_binary(
    name = "resources",
    srcs = ["resources.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [":resources_main_lib"],
)

py_library(
    name = "resources_main_lib",
    srcs = ["resources.py"],
    srcs_version = "PY3",
    deps = ["//deepvariant/protos:resources_py_pb2"],
)

py_test(
    name = "resources_test",
    srcs = ["resources_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":resources_main_lib",
        "@absl_py//absl/testing:absltest",
    ],
)

py_binary(
    name = "vcf_stats_report",
    srcs = ["vcf_stats_report.py"],
    main = "vcf_stats_report.py",
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":vcf_stats_report_lib",
        "@absl_py//absl/flags",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "vcf_stats_report_lib",
    srcs = ["vcf_stats_report.py"],
    srcs_version = "PY3",
    deps = [
        ":vcf_stats",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/util:errors",
    ],
)

py_library(
    name = "vcf_stats",
    srcs = ["vcf_stats.py"],
    srcs_version = "PY3",
    deps = [
        ":vcf_stats_vis",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "vcf_stats_test",
    srcs = ["vcf_stats_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_testdata",
        ":vcf_stats",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "vcf_stats_vis",
    srcs = ["vcf_stats_vis.py"],
    srcs_version = "PY3",
)

py_test(
    name = "vcf_stats_vis_test",
    srcs = ["vcf_stats_vis_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":vcf_stats_vis",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_binary(
    name = "show_examples",
    srcs = ["show_examples.py"],
    python_version = "PY3",
    deps = [
        ":show_examples_lib",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "show_examples_lib",
    srcs = ["show_examples.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:vis",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "show_examples_test",
    srcs = ["show_examples_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_testdata",
        ":show_examples_lib",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "average_model_checkpoint_patched",
    srcs = ["average_model_checkpoint_patched.py"],
)

py_library(
    name = "dv_config",
    srcs = ["dv_config.py"],
)

py_library(
    name = "dashboard_utils",
    srcs = ["dashboard_utils.py"],
    srcs_version = "PY3",
)

py_test(
    name = "dashboard_utils_test",
    srcs = ["dashboard_utils_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":dashboard_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_binary(
    name = "runtime_by_region_vis",
    srcs = ["runtime_by_region_vis.py"],
    python_version = "PY3",
    deps = [
        ":runtime_by_region_vis_lib",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "runtime_by_region_vis_lib",
    srcs = ["runtime_by_region_vis.py"],
    srcs_version = "PY3",
    deps = [
        ":dashboard_utils",
        "//third_party/nucleus/io:sharded_file_utils",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
    ],
)

py_test(
    name = "runtime_by_region_vis_test",
    srcs = ["runtime_by_region_vis_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":py_testdata",
        ":runtime_by_region_vis_lib",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_binary(
    name = "multisample_make_examples",
    srcs = [
        "multisample_make_examples.py",
    ],
    main = "multisample_make_examples.py",
    python_version = "PY3",
    deps = [
        ":logging_level",
        ":make_examples_core",
        ":make_examples_options",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl/flags",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
        "@absl_py//absl:app",
    ],
)

py_library(
    name = "sample",
    srcs = ["sample.py"],
    srcs_version = "PY3",
    deps = [
        ":variant_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/protos:reads_py_pb2",
    ],
)

py_library(
    name = "convert_to_saved_model_lib",
    srcs = ["convert_to_saved_model.py"],
    deps = [
        ":dv_utils",
        ":keras_modeling",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "convert_to_saved_model",
    srcs = ["convert_to_saved_model.py"],
    python_version = "PY3",
    deps = [":convert_to_saved_model_lib"],
)

cc_library(
    name = "make_examples_native",
    srcs = ["make_examples_native.cc"],
    hdrs = ["make_examples_native.h"],
    deps = [
        ":alt_aligned_pileup_lib",
        ":fast_pipeline_utils",
        ":pileup_image_native",
        ":stream_examples",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/io:tfrecord_writer",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_protobuf//:protobuf_lite",
        "@com_googlesource_code_re2//:re2",
        "@org_tensorflow//tensorflow/core/example:example_protos_cc",
    ],
)

cc_test(
    name = "make_examples_native_test",
    srcs = ["make_examples_native_test.cc"],
    deps = [
        ":make_examples_native",
        ":pileup_image_native",
        ":testing_utils",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/strings",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_test(
    name = "pileup_image_native_test",
    srcs = ["pileup_image_native_test.cc"],
    deps = [
        ":pileup_image_native",
        ":testing_utils",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "@com_google_absl//absl/log",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "testing_utils",
    testonly = True,
    srcs = ["testing_utils.cc"],
    hdrs = ["testing_utils.h"],
    deps = [
        ":pileup_image_native",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
        "@com_google_protobuf//:protobuf",
    ],
)

cc_library(
    name = "alt_aligned_pileup_lib",
    srcs = ["alt_aligned_pileup_lib.cc"],
    hdrs = ["alt_aligned_pileup_lib.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//deepvariant/realigner:fast_pass_aligner",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "mapping_quality_channel",
    srcs = ["mapping_quality_channel.cc"],
    hdrs = ["mapping_quality_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "avg_base_quality_channel",
    srcs = ["avg_base_quality_channel.cc"],
    hdrs = ["avg_base_quality_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/log",
    ],
)

cc_library(
    name = "gap_compressed_identity_channel",
    srcs = ["gap_compressed_identity_channel.cc"],
    hdrs = ["gap_compressed_identity_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "gc_content_channel",
    srcs = ["gc_content_channel.cc"],
    hdrs = ["gc_content_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "homopolymer_weighted_channel",
    srcs = ["homopolymer_weighted_channel.cc"],
    hdrs = ["homopolymer_weighted_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "read_supports_variant_channel",
    srcs = ["read_supports_variant_channel.cc"],
    hdrs = ["read_supports_variant_channel.h"],
    deps = [
        ":channel",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/log:check",
    ],
)

cc_library(
    name = "identity_channel",
    srcs = ["identity_channel.cc"],
    hdrs = ["identity_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "channel",
    srcs = ["channel.cc"],
    hdrs = ["channel.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "read_mapping_percent_channel",
    srcs = ["read_mapping_percent_channel.cc"],
    hdrs = ["read_mapping_percent_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "allele_frequency_channel",
    srcs = ["allele_frequency_channel.cc"],
    hdrs = ["allele_frequency_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/types:span",
    ],
)

cc_library(
    name = "is_homopolymer_channel",
    srcs = ["is_homopolymer_channel.cc"],
    hdrs = ["is_homopolymer_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "strand_channel",
    srcs = ["strand_channel.cc"],
    hdrs = ["strand_channel.h"],
    deps = [
        ":channel",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "insert_size_channel",
    srcs = ["insert_size_channel.cc"],
    hdrs = ["insert_size_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_library(
    name = "haplotype_tag_channel",
    srcs = ["haplotype_tag_channel.cc"],
    hdrs = ["haplotype_tag_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/log",
    ],
)

cc_library(
    name = "blank_channel",
    srcs = ["blank_channel.cc"],
    hdrs = ["blank_channel.h"],
    deps = [
        ":channel",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
    ],
)

cc_test(
    name = "alt_aligned_pileup_lib_test",
    srcs = ["alt_aligned_pileup_lib_test.cc"],
    deps = [
        ":alt_aligned_pileup_lib",
        ":testing_utils",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_absl//absl/types:span",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

cc_library(
    name = "stream_examples",
    srcs = ["stream_examples.cc"],
    hdrs = ["stream_examples.h"],
    deps = [
        ":fast_pipeline_utils",
        ":pileup_image_native",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/types:span",
    ],
)

cc_binary(
    name = "fast_pipeline",
    srcs = [
        "fast_pipeline.cc",
        "fast_pipeline.h",
    ],
    deps = [
        ":fast_pipeline_utils",
        "@com_google_absl//absl/flags:flag",
        "@com_google_absl//absl/flags:parse",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/strings",
    ],
)

cc_library(
    name = "fast_pipeline_utils",
    srcs = ["fast_pipeline_utils.h"],
    hdrs = ["fast_pipeline_utils.h"],
    deps = [
        "@com_google_absl//absl/strings",
    ],
)
