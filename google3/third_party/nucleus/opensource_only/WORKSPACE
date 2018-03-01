workspace(name = "nucleus")

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

# Import tensorflow.  Note path.
local_repository(
    name = "org_tensorflow",
    path = "../tensorflow",
)

# Required boilerplate for tf_workspace(), apparently.
http_archive(
    name = "io_bazel_rules_closure",
    sha256 = "25f5399f18d8bf9ce435f85c6bbf671ec4820bc4396b3022cc5dc4bc66303609",
    strip_prefix = "rules_closure-0.4.2",
    urls = [
        "http://mirror.bazel.build/github.com/bazelbuild/rules_closure/archive/0.4.2.tar.gz",  # 2017-08-29
        "https://github.com/bazelbuild/rules_closure/archive/0.4.2.tar.gz",
    ],
)

# Import all of the tensorflow dependencies.
load("@org_tensorflow//tensorflow:workspace.bzl", "tf_workspace")

tf_workspace(tf_repo_name = "org_tensorflow")

new_local_repository(
    name = "clif",
    build_file = "third_party/clif.BUILD",
    path = "/usr/local",
)
