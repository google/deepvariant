# Placeholder: load py_library
# Placeholder: load py_test
# Placeholder: load py_binary

package(
    default_visibility = [
        "//visibility:public",
    ],
)

filegroup(
    name = "binaries",
    srcs = [
        "call_variants",
        "call_variants_keras",
        "freeze_graph",
        "make_examples",
        "make_examples_somatic",
        "model_eval",
        "model_train",
        "multisample_make_examples",
        "postprocess_variants_nucleus",
        "runtime_by_region_vis",
        "show_examples",
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
        ":call_variants_keras_lib",
        ":call_variants_lib",
        ":data_providers",
        ":dv_utils",
        ":dv_vcf_constants",
        ":haplotypes",
        ":keras_modeling",
        ":make_examples_core",
        ":make_examples_lib",
        ":model_eval_lib",
        ":model_train_lib",
        ":modeling",
        ":pileup_image",
        ":postprocess_variants_nucleus_py_lib",
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
        ":utils",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/memory",
        "@com_google_absl//absl/strings",
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
        "@com_google_absl//absl/log",
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
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
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
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:samplers",
        "//third_party/nucleus/vendor:statusor",
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
        "@com_google_absl//absl/container:node_hash_map",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/strings",
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
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "//third_party/nucleus/vendor:statusor",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/hash",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
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
        "//third_party/nucleus/vendor:statusor",
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
        ":dv_constants",
        ":dv_utils",
        ":dv_vcf_constants",
        ":pileup_image",
        ":resources_main_lib",
        ":variant_caller",
        ":vcf_candidate_importer",
        ":very_sensitive_caller",
        "//deepvariant/labeler:customized_classes_labeler",
        "//deepvariant/labeler:haplotype_labeler",
        "//deepvariant/labeler:positional_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/python:direct_phasing",
        "//deepvariant/realigner",
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
        ":dv_utils",
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
        ":dv_constants",
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
    name = "merge_phased_reads_cpp",
    srcs = [
        "merge_phased_reads.cc",
    ],
    hdrs = [
        "merge_phased_reads.h",
    ],
    deps = [
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/flags:flag",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_googlesource_code_re2//:re2",
    ],
)

cc_binary(
    name = "merge_phased_reads",
    srcs = [
        "merge_phased_reads.cc",
        "merge_phased_reads.h",
    ],
    deps = [
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
        "@com_google_absl//absl/flags:flag",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@com_google_absl//absl/status",
        "@com_google_absl//absl/status:statusor",
        "@com_google_absl//absl/strings",
        "@com_google_absl//absl/strings:str_format",
        "@com_googlesource_code_re2//:re2",
    ],
)

cc_test(
    name = "merge_phased_reads_test",
    size = "small",
    srcs = [
        "merge_phased_reads_test.cc",
    ],
    deps = [
        ":merge_phased_reads_cpp",
        "@com_google_absl//absl/container:flat_hash_map",
        "@com_google_absl//absl/container:flat_hash_set",
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

py_test(
    name = "make_examples_test",
    size = "large",
    srcs = ["make_examples_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    shard_count = 2,
    srcs_version = "PY3",
    deps = [
        ":dv_utils",
        ":make_examples_core",
        ":make_examples_lib",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
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
    name = "call_variants_lib",
    srcs = ["call_variants.py"],
    srcs_version = "PY3",
    deps = [
        ":data_providers",
        ":dv_utils",
        ":logging_level",
        ":modeling",
        ":openvino_estimator",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "call_variants",
    srcs = ["call_variants.py"],
    python_version = "PY3",
    deps = [
        ":call_variants_main_lib",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "keras_modeling",
    srcs = ["keras_modeling.py"],
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":dv_utils",
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
    name = "call_variants_keras",
    srcs = ["call_variants_keras.py"],
    python_version = "PY3",
    deps = [
        "call_variants_keras_lib",
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
    deps = [
        ":data_providers",
        ":dv_constants",
        ":dv_utils",
        ":keras_modeling",
        "//third_party/nucleus/io:sharded_file_utils",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

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
    name = "call_variants_keras_lib",
    srcs = ["call_variants_keras.py"],
    srcs_version = "PY3",
    deps = [
        ":call_variants_main_lib",
        ":keras_modeling",
    ],
)

py_library(
    name = "call_variants_main_lib",
    srcs = [
        "call_variants.py",
    ],
    deps = [
        ":call_variants_lib",
        "@absl_py//absl/flags",
    ],
)

py_test(
    name = "call_variants_test",
    size = "large",
    srcs = ["call_variants_test.py"],
    args = ["--use_tpu=false"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":call_variants_main_lib",
        ":dv_utils",
        ":modeling",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "call_variants_accelerator_test",
    size = "large",
    srcs = ["call_variants_accelerator_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    tags = [
        "gpu_test",
        "manual",
        "notap",
    ],
    deps = [
        ":call_variants_main_lib",
        ":modeling",
        ":py_testdata",
        "@absl_py//absl/testing:absltest",
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
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_library(
    name = "pileup_channel_lib",
    hdrs = ["pileup_channel_lib.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "@com_google_absl//absl/container:btree",
        "@com_google_absl//absl/log",
        "@com_google_absl//absl/log:check",
    ],
)

cc_test(
    name = "pileup_channel_lib_test",
    srcs = ["pileup_channel_lib_test.cc"],
    deps = [
        ":pileup_channel_lib",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:struct_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:test",
    ],
)

py_library(
    name = "postprocess_variants_nucleus_py_lib",
    srcs = ["postprocess_variants_nucleus.py"],
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
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:postprocess_variants",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tabix",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/io/python:merge_variants",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:struct_utils",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "postprocess_variants_nucleus",
    srcs = [
        "postprocess_variants_nucleus.py",
    ],
    main = "postprocess_variants_nucleus.py",
    python_version = "PY3",
    deps = [
        ":postprocess_variants_nucleus_py_lib",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_test(
    name = "postprocess_variants_nucleus_test",
    size = "medium",
    srcs = ["postprocess_variants_nucleus_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    srcs_version = "PY3",
    tags = ["requires-net:external"],
    deps = [
        ":dv_constants",
        ":dv_vcf_constants",
        ":postprocess_variants_nucleus_py_lib",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
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
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:postprocess_variants",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tabix",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/io/python:merge_variants",
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
        ":dv_constants",
        ":dv_utils",
        ":keras_modeling",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:sharded_file_utils",
        "@absl_py//absl/logging",
        "@com_google_protobuf//:protobuf_python",
    ],
)

py_library(
    name = "openvino_estimator",
    srcs = ["openvino_estimator.py"],
    srcs_version = "PY3",
    deps = [
        "@absl_py//absl/logging",
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
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:variant_utils",
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
    name = "modeling",
    srcs = ["modeling.py"],
    srcs_version = "PY3",
    deps = [
        ":attention_inception_v3",
        ":attention_module",
        ":data_providers",
        ":dv_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "modeling_test",
    size = "large",
    srcs = ["modeling_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":dv_constants",
        ":dv_utils",
        ":modeling",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "model_train_lib",
    srcs = ["model_train.py"],
    srcs_version = "PY3",
    deps = [
        ":data_providers",
        ":logging_level",
        ":modeling",
        "//third_party/nucleus/util:proto_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "model_train",
    srcs = ["model_train.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":model_train_main_lib",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "model_train_main_lib",
    srcs = ["model_train.py"],
    srcs_version = "PY3",
    deps = [
        ":model_train_lib",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "model_train_test",
    size = "medium",
    srcs = ["model_train_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    shard_count = 10,
    srcs_version = "PY3",
    tags = ["nozapfhahn"],  # timeout
    deps = [
        ":data_providers_test_lib",
        ":dv_constants",
        ":model_train_main_lib",
        ":modeling",
        ":py_testdata",
        "//deepvariant/testing:tf_test_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:flagsaver",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "model_eval_lib",
    srcs = ["model_eval.py"],
    srcs_version = "PY3",
    deps = [
        ":data_providers",
        ":logging_level",
        ":modeling",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "model_eval",
    srcs = ["model_eval.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":model_eval_main_lib",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
    ],
)

py_library(
    name = "model_eval_main_lib",
    srcs = ["model_eval.py"],
    srcs_version = "PY3",
    deps = [
        ":model_eval_lib",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_test(
    name = "model_eval_test",
    size = "large",
    srcs = ["model_eval_test.py"],
    data = [":testdata"],
    python_version = "PY3",
    shard_count = 10,
    srcs_version = "PY3",
    deps = [
        ":data_providers_test_lib",
        ":model_eval_main_lib",
        ":py_testdata",
        "//deepvariant/testing:tf_test_utils",
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
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
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
        "//deepvariant/protos:deepvariant_py_pb2",
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
    name = "attention_module",
    srcs = ["attention_module.py"],
    srcs_version = "PY3",
)

py_test(
    name = "attention_module_test",
    srcs = ["attention_module_test.py"],
    python_version = "PY3",
    srcs_version = "PY3",
    deps = [
        ":attention_module",
        ":dv_constants",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "attention_inception_v3",
    srcs = ["attention_inception_v3.py"],
    srcs_version = "PY3",
    deps = [":attention_module"],
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
    name = "freeze_graph",
    srcs = [
        "freeze_graph.py",
    ],
    main = "freeze_graph.py",
    python_version = "PY3",
    deps = [
        ":modeling",
        "@absl_py//absl:app",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        # TODO: Figure out why we need this dep.
        "@absl_py//absl/testing:parameterized",  # buildcleaner: keep
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
