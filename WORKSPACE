# The workspace name appears at the top of the runfiles tree,
# and in paths to tests, so to keep python happy it is best
# if it is unique.
workspace(name = "com_google_deepvariant")

load("@bazel_tools//tools/build_defs/repo:http.bzl", "http_archive")

# Note: absl_py and com_google_absl (the Python and C++ abseil libraries) are
# provided by TensorFlow.

# CCTZ (Time-zone framework).
# redacted
# work in bazel, so we need to include this to enable nucleus's use of
# //absl/{time,synchronization}
http_archive(
    name = "com_googlesource_code_cctz",
    strip_prefix = "cctz-master",
    urls = ["https://github.com/google/cctz/archive/master.zip"],
)

# This is the 1.10.2 release of htslib.
http_archive(
    name = "htslib",
    build_file = "//:third_party/htslib.BUILD",
    sha256 = "f7994e9636f8a4032dea477a8613f5f73b330c23b5538e45666ce7306240ac14",
    strip_prefix = "htslib-1.10.2",
    urls = [
        "https://github.com/samtools/htslib/archive/1.10.2.zip",
    ],
)

http_archive(
    name = "libssw",
    build_file = "//:third_party/libssw.BUILD",
    sha256 = "10b9305e5a580ee5319f736d3581916f6c873ef4475bd0c0e564c2934334732c",
    strip_prefix = "Complete-Striped-Smith-Waterman-Library-1.0",
    urls = [
        "https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library/archive/v1.0.tar.gz",
    ],
)

# Import tensorflow.  Note path.
local_repository(
    name = "org_tensorflow",
    path = "../tensorflow",
)

# Required boilerplate for tf_workspace().
# This is copied from https://github.com/tensorflow/tensorflow/blob/v2.0.0/WORKSPACE.
http_archive(
    name = "io_bazel_rules_closure",
    sha256 = "5b00383d08dd71f28503736db0500b6fb4dda47489ff5fc6bed42557c07c6ba9",
    strip_prefix = "rules_closure-308b05b2419edb5c8ee0471b67a40403df940149",
    urls = [
        "https://storage.googleapis.com/mirror.tensorflow.org/github.com/bazelbuild/rules_closure/archive/308b05b2419edb5c8ee0471b67a40403df940149.tar.gz",
        "https://github.com/bazelbuild/rules_closure/archive/308b05b2419edb5c8ee0471b67a40403df940149.tar.gz",  # 2019-06-13
    ],
)

# This needs to be in sync with the version of protobuf used by TensorFlow,
# which is currently defined in @tensorflow/tensorflow/workspace.bzl.
# We supply our # own BUILD file, though, so we can prevent ODR violations by
# putting all of Nucleus's C++ binary dependencies into a single library.
# That BUILD file must be kept in sync with the version of protobuf used.
http_archive(
    name = "com_google_protobuf",
    build_file = "//:third_party/protobuf.BUILD",
    sha256 = "b9e92f9af8819bbbc514e2902aec860415b70209f31dfc8c4fa72515a5df9d59",
    # This protobuf release is based on protobuf 3.8.0.
    strip_prefix = "protobuf-310ba5ee72661c081129eb878c1bbcec936b20f0",
    urls = [
        "https://storage.googleapis.com/mirror.tensorflow.org/github.com/protocolbuffers/protobuf/archive/310ba5ee72661c081129eb878c1bbcec936b20f0.tar.gz",
        "https://github.com/protocolbuffers/protobuf/archive/310ba5ee72661c081129eb878c1bbcec936b20f0.tar.gz",
    ],
)

# bazel_skylib is now a required dependency of protobuf_archive.
http_archive(
    name = "bazel_skylib",
    sha256 = "bbccf674aa441c266df9894182d80de104cabd19be98be002f6d478aaa31574d",
    strip_prefix = "bazel-skylib-2169ae1c374aab4a09aa90e65efe1a3aad4e279b",
    urls = ["https://github.com/bazelbuild/bazel-skylib/archive/2169ae1c374aab4a09aa90e65efe1a3aad4e279b.tar.gz"],
)

# Import all of the tensorflow dependencies.
load("@org_tensorflow//tensorflow:workspace.bzl", "tf_workspace")

tf_workspace(tf_repo_name = "org_tensorflow")

# Pull in slim.
# slim is  located inside the tensorflow/models repository.
# The slim subdirectory in the tensorflow/models repository has its own
# WORKSPACE file so we need to strip a prefix to make it the root of the
# repository.
# The prefix is "models-<commit>/slim"
# where commit is the full commit.
# Pin to the lastest version that builds for now. See internal#comment4.
http_archive(
    name = "org_tensorflow_slim",
    strip_prefix = "models-6d140f139cf02ceb87afa76024c4b502a556a3e5/slim",
    urls = ["https://github.com/tensorflow/models/archive/6d140f1.tar.gz"],
)

new_local_repository(
    name = "clif",
    build_file = "third_party/clif.BUILD",
    path = "/usr/local",
)
