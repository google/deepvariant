# The workspace name appears at the top of the runfiles tree,
# and in paths to tests, so to keep python happy it is best
# if it is unique.
workspace(name = "com_google_deepvariant")

# Abseil libraries
git_repository(
    name = "com_google_absl_py",
    # redacted
    commit = "5e343642d987268df199b4c851b7dd3d687ac316",
    remote = "https://github.com/abseil/abseil-py.git",
)
# Note: com_google_absl (the C++ abseil library) is provided by TensorFlow.

# Note: we are using a post-1.6 build release that fixes a double-free.
new_http_archive(
    name = "htslib",
    build_file = "third_party/htslib.BUILD",
    sha256 = "7743e379fa27fdbaa81d4efc97adc5e0b2c5ade3cd09a93e311ea0c6b3a4ddf6",
    strip_prefix = "htslib-57fa9be5255475b2cf9331db32848590a8ea8eb9",
    urls = [
        "https://github.com/samtools/htslib/archive/57fa9be5255475b2cf9331db32848590a8ea8eb9.zip"
    ],
)

new_http_archive(
    name = "libssw",
    build_file = "third_party/libssw.BUILD",
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

# Required boilerplate for tf_workspace(), apparently.
# This is copied from https://github.com/tensorflow/tensorflow/blob/v1.7.0/WORKSPACE.
# Because in build-prereq.sh we are pinning by "git checkout v1.7.0".
# Note: Make sure to update both in future versions.
http_archive(
    name = "io_bazel_rules_closure",
    sha256 = "6691c58a2cd30a86776dd9bb34898b041e37136f2dc7e24cadaeaf599c95c657",
    strip_prefix = "rules_closure-08039ba8ca59f64248bb3b6ae016460fe9c9914f",
    urls = [
        "https://mirror.bazel.build/github.com/bazelbuild/rules_closure/archive/08039ba8ca59f64248bb3b6ae016460fe9c9914f.tar.gz",
        "https://github.com/bazelbuild/rules_closure/archive/08039ba8ca59f64248bb3b6ae016460fe9c9914f.tar.gz",  # 2018-01-16
    ],
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
# Pin to the lastest version that builds for now. See b/68431494#comment4.
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
