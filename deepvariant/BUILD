load("//:third_party/clif.bzl", "pyclif_proto_library")

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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util:reference",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:position_cc_pb2",
        "//deepvariant/util/genomics:range_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
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
        "//deepvariant/util:cpp_test_utils",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util:reference_fai",
        "//deepvariant/util:reference_test",
        "//deepvariant/util/genomics:position_cc_pb2",
        "//deepvariant/util/testing:gunit_extras",
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
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/genomics:reference_cc_pb2",
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
        "//deepvariant/util:cpp_test_utils",
        "//deepvariant/util/genomics:reference_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
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
        "//deepvariant/util:cpp_math",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util:samplers",
        "//deepvariant/util/genomics:variants_cc_pb2",
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
        "//deepvariant/util:cpp_test_utils",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/testing:gunit_extras",
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
        "//deepvariant/util:genomics_math",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:vcf_constants",
        "//deepvariant/util/genomics:variants_py_pb2",
    ],
)

py_test(
    name = "variant_caller_test",
    size = "small",
    srcs = ["variant_caller_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":py_test_utils",
        ":variant_caller",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/util:variant_utils",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "utils",
    srcs = ["utils.cc"],
    hdrs = ["utils.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//deepvariant/util:cpp_utils",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

cc_test(
    name = "utils_test",
    size = "small",
    srcs = ["utils_test.cc"],
    deps = [
        ":utils",
        "//deepvariant/util:cpp_test_utils",
        "//deepvariant/util/testing:gunit_extras",
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
        ":tf_utils",
        ":variant_caller",
        "//deepvariant/labeler:variant_labeler",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:allelecounter",
        "//deepvariant/realigner",
        "//deepvariant/util:errors",
        "//deepvariant/util:genomics_io",
        "//deepvariant/util:htslib_gcp_oauth",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:proto_utils",
        "//deepvariant/util:py_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/io:sam",
        "//deepvariant/util/io:vcf",
        "//deepvariant/util/protos:core_py_pb2",
        "//deepvariant/util/python:hts_verbose",
        "//deepvariant/vendor:timer",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
    ],
)

py_binary(
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
    size = "medium",
    srcs = ["make_examples_test.py"],
    data = [":testdata"],
    shard_count = 2,
    srcs_version = "PY2AND3",
    deps = [
        ":make_examples",
        ":py_test_utils",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/protos:realigner_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/util:errors",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:vcf_constants",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/io:vcf",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "call_variants_lib",
    srcs = ["call_variants.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":logging_level",
        ":modeling",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/util:errors",
        "//deepvariant/util:htslib_gcp_oauth",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:proto_utils",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
    ],
)

py_binary(
    name = "call_variants",
    srcs = [
        "call_variants.py",
    ],
    main = "call_variants.py",
    deps = [
        ":call_variants_lib",
        "@com_google_absl_py//absl/flags",
    ],
)

py_test(
    name = "call_variants_test",
    size = "medium",
    srcs = ["call_variants_test.py"],
    data = [":testdata"],
    srcs_version = "PY2AND3",
    deps = [
        ":call_variants",
        ":modeling",
        ":py_test_utils",
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/util:errors",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:variant_utils",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
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
        ":py_test_utils",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "pileup_image",
    srcs = ["pileup_image.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:pileup_image_native",
        "//deepvariant/util:py_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/protos:core_py_pb2",
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
        ":py_test_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:pileup_image_native",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/python:reference_fai",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

cc_library(
    name = "pileup_image_native",
    srcs = ["pileup_image_native.cc"],
    hdrs = ["pileup_image_native.h"],
    deps = [
        "//deepvariant/protos:deepvariant_cc_pb2",
        "//deepvariant/util/genomics:cigar_cc_pb2",
        "//deepvariant/util/genomics:position_cc_pb2",
        "//deepvariant/util/genomics:reads_cc_pb2",
        "//deepvariant/util/genomics:variants_cc_pb2",
        "//deepvariant/util/protos:core_cc_pb2",
        "@org_tensorflow//tensorflow/core:lib",
    ],
)

py_library(
    name = "postprocess_variants_py_lib",
    srcs = ["postprocess_variants.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":logging_level",
        "//deepvariant:haplotypes",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/python:postprocess_variants",
        "//deepvariant/util:errors",
        "//deepvariant/util:genomics_io",
        "//deepvariant/util:genomics_math",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:proto_utils",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util:vcf_constants",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "//deepvariant/util/io:vcf",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
    ],
)

py_binary(
    name = "postprocess_variants",
    srcs = [
        "postprocess_variants.py",
    ],
    main = "postprocess_variants.py",
    deps = [
        ":postprocess_variants_py_lib",
        "@com_google_absl_py//absl/flags",
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
        ":py_test_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/util:errors",
        "//deepvariant/util:genomics_math",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:vcf_constants",
        "//deepvariant/util/genomics:reference_py_pb2",
        "//deepvariant/util/genomics:struct_py_pb2",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "py_test_utils",
    testonly = True,
    srcs = ["test_utils.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/util:py_test_utils",
    ],
)

py_library(
    name = "data_providers",
    srcs = ["data_providers.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":tf_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/util:io_utils",
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
        ":modeling",
        ":pileup_image",
        ":py_test_utils",
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:variant_utils",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_test(
    name = "data_providers_test",
    size = "small",
    srcs = ["data_providers_test.py"],
    data = [":testdata"],
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
        "//deepvariant/util:genomics_math",
        "//deepvariant/util:variant_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
    ],
)

py_test(
    name = "haplotypes_test",
    size = "small",
    srcs = ["haplotypes_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":haplotypes",
        ":py_test_utils",
        "//deepvariant/testing:flagsaver",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "modeling",
    srcs = ["modeling.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":pileup_image",
        ":tf_utils",
        "@com_google_absl_py//absl/logging",
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
        ":modeling",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
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
        "//deepvariant/util:proto_utils",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
    ],
)

py_binary(
    name = "model_train",
    srcs = ["model_train.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":model_train_lib",
        "@com_google_absl_py//absl/flags",
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
        ":py_test_utils",
        "//deepvariant/testing:flagsaver",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
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
        "//deepvariant/util:proto_utils",
        "//deepvariant/util:variant_utils",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/logging",
    ],
)

py_binary(
    name = "model_eval",
    srcs = ["model_eval.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":model_eval_lib",
        "@com_google_absl_py//absl/flags",
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
        ":modeling",
        ":pileup_image",
        ":py_test_utils",
        "//deepvariant/testing:flagsaver",
        "//deepvariant/util:variant_utils",
        "@com_google_absl_py//absl/flags",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "tf_utils",
    srcs = ["tf_utils.py"],
    srcs_version = "PY2AND3",
    deps = [
        "//deepvariant/protos:deepvariant_py_pb2",
        "//deepvariant/util:io_utils",
        "//deepvariant/util:ranges",
        "//deepvariant/util/genomics:variants_py_pb2",
    ],
)

py_test(
    name = "tf_utils_test",
    size = "small",
    srcs = ["tf_utils_test.py"],
    srcs_version = "PY2AND3",
    deps = [
        ":py_test_utils",
        ":tf_utils",
        "//deepvariant/util:io_utils",
        "//deepvariant/util/genomics:variants_py_pb2",
        "@com_google_absl_py//absl/testing:absltest",
        "@com_google_absl_py//absl/testing:parameterized",
    ],
)

py_library(
    name = "logging_level",
    srcs = ["logging_level.py"],
    deps = ["@com_google_absl_py//absl/flags"],
)

py_library(
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
        "@com_google_absl_py//absl/testing:absltest",
    ],
)
