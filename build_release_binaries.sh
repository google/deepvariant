#!/bin/bash
# Copyright 2017 Google LLC.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from this
#    software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

#
# Build release binaries using our standard compiler flags.
# Note we use the lowest common denominator compiler optimization options (For
# Google Cloud Engine chipsets lowest possible is Sandy Bridge).

# NOLINT
source settings.sh

set -e

# Bazel's --build_python_zip replaces our carefully engineered symbolic links
# with copies.  This function puts the symbolic links back.
function fix_zip_file {
  orig_zip_file=$1

  # Step 1:  Copy the zip file to a temporary place.
  TMPDIR=$(mktemp -d -t tmp.XXXXXXXXXXX)
  # The .zip version of the binary doesn't have the header that makes it
  # self-executable.  We use that version because otherwise unzip would
  # complain and raise an error code.
  cp "${orig_zip_file}.zip" "${TMPDIR}"

  # Step 2: Unzip it.
  pushd "${TMPDIR}" > /dev/null
  BN=$(basename "${orig_zip_file}")
  unzip -qq "${BN}.zip"

  # Step 3: Restore the symbolic links.
  find "runfiles/com_google_deepvariant" -name '*.so' ! -name 'examples_from_stream.so' -exec ln --force -s --relative "runfiles/com_google_protobuf/python/google/protobuf/pyext/_message.so" {} \;

  # Step 4: Fix the __main__.py's use of zipfile, which can't handle
  # symbolic links.  Replace it with an invocation of unzip, which can.
  # The lines we replace are
  # with zipfile.ZipFile(zip_path) as zf:
  #   for info in zf.infolist():
  #     zf.extract(info, dest_dir)
  #     # UNC-prefixed paths must be absolute/normalized. See

  sed -i 's/  with zipfile.ZipFile(zip_path) as zf:/  if True:/' __main__.py
  sed -i 's/  for info in zf.infolist():/  if True:/' __main__.py
  sed -i 's/  zf.extract(info, dest_dir)/  os.system("unzip -qq " + zip_path + " -d " + dest_dir)/' __main__.py
  sed -i 's/  # UNC-prefixed paths must be absolute\/normalized. See/  return/' __main__.py

  # Step 5: Zip it back up, with zip --symbolic
  rm -f "${BN}.zip"
  ZIP_OUT="/tmp/${BN}.zip"
  rm -f "${ZIP_OUT}"
  zip -q --symlinks -r "${ZIP_OUT}" *

  # Step 6: Make the zip file self-executable
  SELF_ZIP="/tmp/${BN}"
  # If the Python interpreter discovers it is being run from part of a zip
  # file, it will uncompress and run the __main__.py.  This is the trick that
  # bazel uses to make a self-executable zip, see for example
  # https://github.com/bazelbuild/bazel/blob/558b717e906156477b1c6bd29d049a0fb8e18b27/src/main/java/com/google/devtools/build/lib/bazel/rules/python/BazelPythonSemantics.java#L193
  echo '#!/usr/bin/env python3' | cat - "${ZIP_OUT}" > "${SELF_ZIP}"

  # Step 7: Copy it back and make it executable.
  popd > /dev/null
  rm -f "${orig_zip_file}"
  mv "${SELF_ZIP}" "${orig_zip_file}"
  chmod +x "${orig_zip_file}"

  # Step 8: DeepVariant also uses "${orig_zip_file}.zip" in many of its
  # instructions, so make sure that we also copy that.
  rm -f "${orig_zip_file}.zip"
  mv "${ZIP_OUT}" "${orig_zip_file}.zip"
  # No executable bit because the .zip version is not self-executing and
  # must be invoked as
  #   python3 ${orig_zip_file}.zip
}

# Building examples_from_stream.so C++ library. It cannot be built correctly
# with the default bazel setup, so we build it manually.
# examples_from_stream.so is used by call_variants target therefore it has to
# be built before :binaries.
TF_CFLAGS=( $(python3 -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_compile_flags()))') )
TF_LFLAGS=( $(python3 -c 'import tensorflow as tf; print(" ".join(tf.sysconfig.get_link_flags()))') )

# shellcheck disable=SC2068
g++ -std=c++14 -shared \
        deepvariant/stream_examples_kernel.cc  \
        deepvariant/stream_examples_ops.cc \
        -o deepvariant/examples_from_stream.so \
        -fPIC \
        -l:libtensorflow_framework.so.2  \
        -I. \
        ${TF_CFLAGS[@]} \
        ${TF_LFLAGS[@]} \
        -D_GLIBCXX_USE_CXX11_ABI=1 \
        --std=c++17 \
        -DEIGEN_MAX_ALIGN_BYTES=64 \
        -O2

# shellcheck disable=SC2086
bazel build -c opt \
  //deepvariant:fast_pipeline

# shellcheck disable=SC2086
bazel build -c opt \
  --output_filter=DONT_MATCH_ANYTHING \
  --noshow_loading_progress \
  --show_result=0 \
  ${DV_COPT_FLAGS} \
  --build_python_zip \
  :binaries

# shellcheck disable=SC2086
bazel build -c opt \
  --output_filter=DONT_MATCH_ANYTHING \
  --noshow_loading_progress \
  --show_result=0 \
  ${DV_COPT_FLAGS} \
  --build_python_zip \
  //deepvariant/labeler:labeled_examples_to_vcf

# shellcheck disable=SC2086
bazel build -c opt \
  --output_filter=DONT_MATCH_ANYTHING \
  --noshow_loading_progress \
  --show_result=0 \
  ${DV_COPT_FLAGS} \
  --build_python_zip \
  //deepvariant:convert_to_saved_model

# shellcheck disable=SC2086
bazel build -c opt \
  --output_filter=DONT_MATCH_ANYTHING \
  --noshow_loading_progress \
  --show_result=0 \
  ${DV_COPT_FLAGS} \
  --build_python_zip \
  :binaries-deeptrio

# shellcheck disable=SC2086
bazel build  -c opt \
  --output_filter=DONT_MATCH_ANYTHING \
  --noshow_loading_progress \
  --show_result=0 \
  --noshow_progress \
  ${DV_COPT_FLAGS} \
  :licenses_zip

# Bazel understandably doesn't like it when its output files are edited, so
# make sure all the builds are done before we fix things.

# TODO: Replace this hand-made list with a find command.
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/train"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/call_variants"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/load_gbz_into_shared_memory"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/make_examples"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/make_examples_pangenome_aware_dv"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/make_examples_somatic"
fix_zip_file "bazel-out/k8-opt/bin/deeptrio/make_examples"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/postprocess_variants"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/vcf_stats_report"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/show_examples"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/runtime_by_region_vis"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/convert_to_saved_model"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/multisample_make_examples"
fix_zip_file "bazel-out/k8-opt/bin/deepvariant/labeler/labeled_examples_to_vcf"
