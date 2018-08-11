load("//:third_party/clif.bzl", "pyclif_proto_library")
load("//:tools/nucleus.bzl", "nucleus_py_library", "nucleus_py_extension", "nucleus_py_binary")

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
        "//deepvariant/realigner/allele_count_linear:generate_trained_model",
    ],
)

filegroup(
    name = "testdata",
    testonly = True,
    srcs = glob([
        "testdata/*.bam*",
        "testdata/*.bed*",
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
nucleus_py_library(
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
        "@com_google_absl//absl/strings",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "allelecounter_test",
    size = "small",
    srcs = ["allelecounter_test.cc"],
    deps = [
        ":allelecounter",
        ":utils",
        "//third_party/nucleus/io:reference_fai",
        "//third_party/nucleus/io:reference_test",
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

nucleus_py_library(
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
    srcs_version = "PY2AND3",
    deps = [
        ":py_testdata",
        ":variant_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
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

nucleus_py_library(
    name = "make_examples_lib",
    srcs = ["make_examples.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":exclude_contigs",
        ":logging_level",
        ":pileup_image",
        ":resources",
        ":tf_utils",
        ":variant_caller",
        "//deepvariant/labeler:customized_classes_labeler",
        "//deepvariant/labeler:haplotype_labeler",
        "//deepvariant/labeler:positional_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/realigner",
        "//deepvariant/vendor:timer",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:sam",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/io/python:hts_verbose",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:py_utils",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

nucleus_py_binary(
    name = "make_examples",
    srcs = [
        "make_examples.py",
    ],
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
    shard_count = 2,
    srcs_version = "PY2AND3",
    deps = [
        ":make_examples",
        ":py_testdata",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:reads_py_pb2",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "call_variants_lib",
    srcs = ["call_variants.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers",
        ":logging_level",
        ":modeling",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

nucleus_py_binary(
    name = "call_variants",
    srcs = [
        "call_variants.py",
    ],
    main = "call_variants.py",
    deps = [
        ":call_variants_lib",
        "@io_abseil_py//absl/flags",
    ],
)

py_test(
    name = "call_variants_test",
    size = "medium",
    srcs = ["call_variants_test.py"],
    args = ["--use_tpu=false"],
    data = [":testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":call_variants",
        ":modeling",
        ":py_testdata",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/testing:tf_test_utils",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "call_variants_accelerator_test",
    size = "large",
    srcs = ["call_variants_accelerator_test.py"],
    data = [":testdata"],
    srcs_version = "PY2AND3",
    tags = [
        "gpu_test",
        "manual",
        "notap",
    ],
    deps = [
        ":call_variants",
        ":modeling",
        ":py_testdata",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
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
    srcs_version = "PY2AND3",
    deps = [
        ":pileup_image",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:pileup_image_native",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:ranges",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
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
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

nucleus_py_library(
    name = "postprocess_variants_py_lib",
    srcs = ["postprocess_variants.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":dv_vcf_constants",
        ":logging_level",
        "//deepvariant:haplotypes",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:postprocess_variants",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/io:vcf",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:ranges",
        "//third_party/nucleus/util:variant_utils",
        "//third_party/nucleus/util:variantcall_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

nucleus_py_binary(
    name = "postprocess_variants",
    srcs = [
        "postprocess_variants.py",
    ],
    main = "postprocess_variants.py",
    deps = [
        ":postprocess_variants_py_lib",
        "@io_abseil_py//absl/flags",
    ],
)

py_test(
    name = "postprocess_variants_test",
    size = "small",
    srcs = ["postprocess_variants_test.py"],
    data = [":testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":postprocess_variants_py_lib",
        ":py_testdata",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//third_party/nucleus/io:fasta",
        "//third_party/nucleus/protos:reference_py_pb2",
        "//third_party/nucleus/protos:struct_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:errors",
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:vcf_constants",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "py_testdata",
    testonly = True,
    srcs = ["testdata.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//third_party/nucleus/testing:py_test_utils",
    ],
)

nucleus_py_library(
    name = "data_providers",
    srcs = ["data_providers.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":dv_constants",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/util:io_utils",
        "@protobuf_archive//:protobuf_python",
    ],
)

nucleus_py_library(
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
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "data_providers_test",
    size = "medium",
    srcs = ["data_providers_test.py"],
    data = [":testdata"],
    main = "data_providers_test.py",
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers_test_lib",
    ],
)

nucleus_py_library(
    name = "haplotypes",
    srcs = ["haplotypes.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//third_party/nucleus/util:genomics_math",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

py_test(
    name = "haplotypes_test",
    size = "small",
    srcs = ["haplotypes_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":haplotypes",
        ":py_testdata",
        "//deepvariant/testing:flagsaver",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "modeling",
    srcs = ["modeling.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers",
        ":tf_utils",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
        "@org_tensorflow_slim//:inception",
        "@org_tensorflow_slim//:mobilenet_v1",
        "@org_tensorflow_slim//:resnet_v2",
    ],
)

py_test(
    name = "modeling_test",
    size = "medium",
    srcs = ["modeling_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":dv_constants",
        ":modeling",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "model_train_lib",
    srcs = ["model_train.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers",
        ":logging_level",
        ":modeling",
        "//third_party/nucleus/util:proto_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

nucleus_py_binary(
    name = "model_train",
    srcs = ["model_train.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":model_train_lib",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

py_test(
    name = "model_train_test",
    size = "medium",
    srcs = ["model_train_test.py"],
    data = [":testdata"],
    shard_count = 10,
    srcs_version = "PY2AND3",
    tags = ["nozapfhahn"],  # timeout
    deps = [
        ":data_providers_test_lib",
        ":model_train",
        ":modeling",
        ":py_testdata",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/testing:tf_test_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "model_eval_lib",
    srcs = ["model_eval.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers",
        ":logging_level",
        ":modeling",
        "//third_party/nucleus/util:proto_utils",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

nucleus_py_binary(
    name = "model_eval",
    srcs = ["model_eval.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":model_eval_lib",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/logging",
    ],
)

py_test(
    name = "model_eval_test",
    size = "large",
    srcs = ["model_eval_test.py"],
    data = [":testdata"],
    shard_count = 10,
    srcs_version = "PY2AND3",
    deps = [
        ":data_providers_test_lib",
        ":model_eval",
        ":py_testdata",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/testing:tf_test_utils",
        "//third_party/nucleus/testing:py_test_utils",
        "//third_party/nucleus/util:variant_utils",
        "@io_abseil_py//absl/flags",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "tf_utils",
    srcs = ["tf_utils.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:io_utils",
        "//third_party/nucleus/util:ranges",
    ],
)

py_test(
    name = "tf_utils_test",
    size = "small",
    srcs = ["tf_utils_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":py_testdata",
        ":tf_utils",
        "//third_party/nucleus/protos:variants_py_pb2",
        "//third_party/nucleus/util:io_utils",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_library(
    name = "logging_level",
    srcs = ["logging_level.py"],
    deps = ["@io_abseil_py//absl/flags"],
)

nucleus_py_library(
    name = "exclude_contigs",
    srcs = ["exclude_contigs.py"],
)

py_test(
    name = "exclude_contigs_test",
    size = "small",
    srcs = ["exclude_contigs_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":exclude_contigs",
        "@io_abseil_py//absl/testing:absltest",
    ],
)

nucleus_py_library(
    name = "dv_constants",
    srcs = ["dv_constants.py"],
    srcs_version = "PY2AND3",
    deps = [],
)

nucleus_py_library(
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
    srcs_version = "PY2AND3",
    deps = [
        ":dv_vcf_constants",
        "//third_party/nucleus/protos:reference_py_pb2",
        "@io_abseil_py//absl/testing:absltest",
        "@io_abseil_py//absl/testing:parameterized",
    ],
)

nucleus_py_binary(
    name = "resources",
    srcs = ["resources.py"],
    srcs_version = "PY2AND3",
    deps = ["//deepvariant/protos:resources_py_pb2"],
)

py_test(
    name = "resources_test",
    srcs = ["resources_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":resources",
        "@io_abseil_py//absl/testing:absltest",
    ],
)
