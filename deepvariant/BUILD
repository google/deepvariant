# load("//devtools/python/blaze:python3.bzl", "py2and3_test")

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
        "model_eval",
        "model_train",
        "postprocess_variants",
        "vcf_stats_report",
    ],
)

filegroup(
    name = "testdata",
    testonly = True,
    srcs = glob([
        "testdata/*.bam*",
        "testdata/*.bed*",
        "testdata/*.cram*",
        "testdata/*.tfrecord*",
        "testdata/*.vcf*",
        "testdata/*.fasta*",
        "testdata/*.pbtxt*",
        "testdata/*.pckl*",
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
        ":call_variants_lib",
        ":data_providers",
        ":dv_vcf_constants",
        ":haplotypes",
        ":make_examples_lib",
        ":model_eval_lib",
        ":model_train_lib",
        ":modeling",
        ":pileup_image",
        ":postprocess_variants_py_lib",
        ":tf_utils",
        "//deepvariant/python:allelecounter",
        "//deepvariant/python:variant_calling",
        "//deepvariant/realigner",
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
        "//third_party/nucleus/io:reference",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
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
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/protos:reference_cc_pb2",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "@org_tensorflow//tensorflow/core:lib",
        # redacted
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
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/io:vcf_reader",
        "//third_party/nucleus/protos:range_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:cpp_math",
        "//third_party/nucleus/util:cpp_utils",
        "//third_party/nucleus/util:samplers",
        "@org_tensorflow//tensorflow/core:lib",
        "@protobuf_archive//:protobuf",
    ],
)

cc_test(
    name = "variant_calling_test",
    size = "small",
    srcs = ["variant_calling_test.cc"],
    deps = [
        ":utils",
        ":variant_calling",
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/testing:cpp_test_utils",
        "//third_party/nucleus/testing:gunit_extras",
        "//third_party/nucleus/util:cpp_utils",
        "@com_google_googletest//:gtest_main",
        "@org_tensorflow//tensorflow/core:lib",
        "@org_tensorflow//tensorflow/core:test",
        "@protobuf_archive//:protobuf",
    ],
)

py_library(
    name = "variant_caller",
    srcs = ["variant_caller.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/python:variant_calling",
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
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
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
    srcs_version = "PY2AND3",
    deps = [
        ":variant_caller",
    ],
)

py_test(
    name = "very_sensitive_caller_test",
    size = "small",
    srcs = ["very_sensitive_caller_test.py"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":py_testdata",
        ":very_sensitive_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "vcf_caller",
    srcs = ["vcf_caller.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":variant_caller",
        "//third_party/nucleus/io:vcf",
    ],
)

py_test(
    name = "vcf_caller_test",
    size = "small",
    srcs = ["vcf_caller_test.py"],
    data = [":testdata"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":py_testdata",
        ":vcf_caller",
        "//deepvariant/labeler:labeled_examples_to_vcf_main_lib",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:variant_utils",
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
    srcs_version = "PY2AND3",
    deps = [
        ":exclude_contigs",
        ":logging_level",
        ":pileup_image",
        ":resources_main_lib",
        ":tf_utils",
        ":vcf_caller",
        ":very_sensitive_caller",
        "//deepvariant:dv_constants",
        "//deepvariant/labeler:customized_classes_labeler",
        "//deepvariant/labeler:haplotype_labeler",
        "//deepvariant/labeler:positional_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/realigner",
        "//deepvariant/vendor:timer",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "make_examples",
    srcs = [
        "make_examples.py",
    ],
    default_python_version = "PY2",
    main = "make_examples.py",
    deps = [
        ":make_examples_lib",
    ],
)

py_test(
    name = "make_examples_test",
    size = "large",
    srcs = ["make_examples_test.py"],
    data = [":testdata"],
    default_python_version = "PY2",
    shard_count = 2,
    srcs_version = "PY2AND3",
    deps = [
        ":make_examples_lib",
        ":py_testdata",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "call_variants_lib",
    srcs = ["call_variants.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers",
        ":logging_level",
        ":modeling",
        ":tf_utils",
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
    default_python_version = "PY2",
    deps = [":call_variants_main_lib"],
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
    size = "medium",
    srcs = ["call_variants_test.py"],
    args = ["--use_tpu=false"],
    data = [":testdata"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":call_variants_main_lib",
        ":modeling",
        ":py_testdata",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/testing:tf_test_utils",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "call_variants_accelerator_test",
    size = "large",
    srcs = ["call_variants_accelerator_test.py"],
    data = [":testdata"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
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
    srcs_version = "PY2AND3",
    deps = [
        ":dv_constants",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:pileup_image_native",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
    ],
)

py_test(
    name = "pileup_image_test",
    size = "small",
    srcs = ["pileup_image_test.py"],
    data = [":testdata"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
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
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//third_party/nucleus/protos:cigar_cc_pb2",
        "//third_party/nucleus/protos:position_cc_pb2",
        "//third_party/nucleus/protos:reads_cc_pb2",
        "//third_party/nucleus/protos:variants_cc_pb2",
        "//third_party/nucleus/util:proto_ptr",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

py_library(
    name = "postprocess_variants_py_lib",
    srcs = ["postprocess_variants.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":dv_vcf_constants",
        ":logging_level",
        ":tf_utils",
        ":vcf_stats",
        "//deepvariant:dv_constants",
        "//deepvariant:haplotypes",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:postprocess_variants",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sharded_file_utils",
        "//third_party/nucleus/io:tabix",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
    ],
)

py_binary(
    name = "postprocess_variants",
    srcs = [
        "postprocess_variants.py",
    ],
    default_python_version = "PY2",
    main = "postprocess_variants.py",
    deps = [
        ":postprocess_variants_py_lib",
        "@absl_py//absl/flags",
    ],
)

py_test(
    name = "postprocess_variants_test",
    size = "small",
    srcs = ["postprocess_variants_test.py"],
    data = [":testdata"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":postprocess_variants_py_lib",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:tfrecord",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:vcf_constants",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "py_testdata",
    testonly = True,
    srcs = ["testdata.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//third_party/nucleus/testing:py_test_utils",
    ],
)

py_library(
    name = "data_providers",
    srcs = ["data_providers.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":dv_constants",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:sharded_file_utils",
        "@protobuf_archive//:protobuf_python",
    ],
)

py_library(
    name = "data_providers_test_lib",
    testonly = True,
    srcs = ["data_providers_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers",
        ":dv_constants",
        ":modeling",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/io:tfrecord",
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
    default_python_version = "PY2",
    main = "data_providers_test.py",
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers_test_lib",
    ],
)

py_library(
    name = "haplotypes",
    srcs = ["haplotypes.py"],
    srcs_version = "PY2AND3",
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
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":haplotypes",
        ":py_testdata",
        "//deepvariant/testing:flagsaver",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "modeling",
    srcs = ["modeling.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers",
        ":tf_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/logging",
        "@org_tensorflow_slim//:inception",
        "@org_tensorflow_slim//:mobilenet_v1",
        "@org_tensorflow_slim//:resnet_v2",
    ],
)

py_test(
    name = "modeling_test",
    size = "medium",
    srcs = ["modeling_test.py"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":dv_constants",
        ":modeling",
        ":tf_utils",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "model_train_lib",
    srcs = ["model_train.py"],
    srcs_version = "PY2AND3",
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
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [":model_train_main_lib"],
)

py_library(
    name = "model_train_main_lib",
    srcs = ["model_train.py"],
    srcs_version = "PY2AND3",
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
    default_python_version = "PY2",
    shard_count = 10,
    srcs_version = "PY2AND3",
    tags = ["nozapfhahn"],  # timeout
    deps = [
        ":data_providers_test_lib",
        ":model_train_main_lib",
        ":modeling",
        ":py_testdata",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/testing:tf_test_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "model_eval_lib",
    srcs = ["model_eval.py"],
    srcs_version = "PY2AND3",
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
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [":model_eval_main_lib"],
)

py_library(
    name = "model_eval_main_lib",
    srcs = ["model_eval.py"],
    srcs_version = "PY2AND3",
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
    default_python_version = "PY2",
    shard_count = 10,
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers_test_lib",
        ":model_eval_main_lib",
        ":py_testdata",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/testing:tf_test_utils",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:variant_utils",
        "@absl_py//absl/flags",
        "@absl_py//absl/testing:absltest",
        "@absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "tf_utils",
    srcs = ["tf_utils.py"],
    srcs_version = "PY2AND3",
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
    name = "tf_utils_test",
    size = "small",
    srcs = ["tf_utils_test.py"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":py_testdata",
        ":tf_utils",
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
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":exclude_contigs",
        "@absl_py//absl/testing:absltest",
    ],
)

py_library(
    name = "dv_constants",
    srcs = ["dv_constants.py"],
    srcs_version = "PY2AND3",
    deps = [],
)

py_library(
    name = "dv_vcf_constants",
    srcs = ["dv_vcf_constants.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:vcf_constants",
    ],
)

py_test(
    name = "dv_vcf_constants_test",
    size = "small",
    srcs = ["dv_vcf_constants_test.py"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
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
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [":resources_main_lib"],
)

py_library(
    name = "resources_main_lib",
    srcs = ["resources.py"],
    srcs_version = "PY2AND3",
    deps = ["//deepvariant/protos:resources_py_pb2"],
)

py_test(
    name = "resources_test",
    srcs = ["resources_test.py"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":resources_main_lib",
        "@absl_py//absl/testing:absltest",
    ],
)

py_binary(
    name = "vcf_stats_report",
    srcs = ["vcf_stats_report.py"],
    default_python_version = "PY2",
    main = "vcf_stats_report.py",
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_stats_report_lib",
        "@absl_py//absl/flags",
    ],
)

py_library(
    name = "vcf_stats_report_lib",
    srcs = ["vcf_stats_report.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_stats",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/util:errors",
    ],
)

py_library(
    name = "vcf_stats",
    srcs = ["vcf_stats.py"],
    srcs_version = "PY2AND3",
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
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
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
    srcs_version = "PY2AND3",
)

py_test(
    name = "vcf_stats_vis_test",
    srcs = ["vcf_stats_vis_test.py"],
    default_python_version = "PY2",
    srcs_version = "PY2AND3",
    deps = [
        ":vcf_stats_vis",
        "//third_party/nucleus/testing:py_test_utils",
        "@absl_py//absl/testing:absltest",
    ],
)
