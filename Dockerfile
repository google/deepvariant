# Copyright 2019 Google LLC.
# This is used to build the DeepVariant release docker image.
# It can also be used to build local images, especially if you've made changes
# to the code.
# Example command:
# $ git clone https://github.com/google/deepvariant.git
# $ cd deepvariant
# $ sudo docker build -t deepvariant .
#
# To build for GPU, use a command like:
# $ sudo docker build --build-arg=FROM_IMAGE=nvidia/cuda:12.3.2-cudnn9-devel-ubuntu22.04 --build-arg=DV_GPU_BUILD=1 -t deepvariant_gpu .


ARG FROM_IMAGE=ubuntu:22.04
# PYTHON_VERSION is also set in settings.sh.
ARG PYTHON_VERSION=3.10
ARG DV_GPU_BUILD=0
ARG VERSION=1.10.0-rc1
ENV VERSION=${VERSION}
ARG TF_ENABLE_ONEDNN_OPTS=1

#======================================#
# Stage 1: Install samtools + bcftools #
#======================================#
FROM condaforge/miniforge3:24.9.2-0 AS hts_utils
RUN conda config --add channels bioconda
RUN conda create -n bio \
                    bioconda::bcftools=1.15 \
                    bioconda::samtools=1.15 \
    && conda clean -a

#==========================#
# Stage 2: Download Models #
#==========================#
FROM alpine:latest AS download_models

RUN apk add --no-cache wget parallel

# Copy models
# The hybrid and ont models mix naming conventions:
# hybrid_pacbio_illumina --> hybrid
# ont_r104 --> ont
RUN parallel --halt now,fail=1 --verbose --jobs 10 \
  "mkdir -p /opt/models/{1}/variables && wget -O /opt/models/{1}/{2} https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.{=1 s/_.*// =}.savedmodel/{2}" ::: \
  wgs wes pacbio masseq ont_r104 hybrid_pacbio_illumina ::: \
  fingerprint.pb saved_model.pb model.example_info.json variables/variables.data-00000-of-00001 variables/variables.index && \
  chmod -R +r /opt/models/

# Download small models
RUN parallel --halt now,fail=1 --verbose --jobs 10 \
  "mkdir -p /opt/smallmodels/{1}/variables && wget -O /opt/smallmodels/{1}/{2} https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.{=1 s/_.*// =}.smallmodel/{2}" ::: \
  wgs pacbio ont_r104 ::: \
  fingerprint.pb saved_model.pb keras_metadata.pb variables/variables.data-00000-of-00001 variables/variables.index && \
  chmod -R +r /opt/smallmodels/

#===================#
# Stage 3: Build DV #
#===================#
FROM ${FROM_IMAGE} AS prereq
LABEL maintainer="https://github.com/google/deepvariant/issues"

# DV_GPU_BUILD, PYTHON_VERSION, and TF_ENABLE_ONEDNN_OPTS are used by
# ./build-prereq.sh and by tensorflow during the build.
ARG DV_GPU_BUILD
ENV DV_GPU_BUILD=${DV_GPU_BUILD}

ARG PYTHON_VERSION
ENV PYTHON_VERSION ${PYTHON_VERSION}

ARG TF_ENABLE_ONEDNN_OPTS
ENV TF_ENABLE_ONEDNN_OPTS ${TF_ENABLE_ONEDNN_OPTS}

# Copy over just ./build-prereq.sh, ./run-prereq.sh, and ./tools/build_absl.sh
# so we can cache these build steps.

WORKDIR /opt/deepvariant
COPY ./build-prereq.sh \
     ./run-prereq.sh \
     ./settings.sh \
     /opt/deepvariant/
COPY ./tools/build_absl.sh /opt/deepvariant/tools/
COPY ./third_party/tensorflow.bzl.patch /opt/deepvariant/third_party/
RUN ./build-prereq.sh

#=====================================#
# Stage 4: Build DeepVariant Binaries #
#=====================================#
FROM prereq AS builder

COPY . /opt/deepvariant
RUN PATH="${HOME}/bin:${PATH}" ./build_release_binaries.sh # PATH for bazel

#===============================#
# Stage 5: Integrate everything #
#===============================#

FROM ${FROM_IMAGE}

ENV DV_BIN_PATH=/opt/deepvariant/bin

# Install libraries
RUN apt-get -y update && \
  apt-get install -y parallel python3-pip unzip && \
  PATH="${HOME}/.local/bin:$PATH" python3 -m pip install absl-py==0.13.0 && \
  apt-get clean autoclean && \
  apt-get autoremove -y --purge && \
  rm -rf /var/lib/apt/lists/*

# Since samtools/bcftools is relatively static, we copy them first.
# Copy over samtools and bcftools
COPY --from=hts_utils /opt/conda/envs/bio/bin /opt/conda/envs/bio/bin
COPY --from=hts_utils /opt/conda/envs/bio/lib /opt/conda/envs/bio/lib

# Integrate everything.
RUN echo "Acquire::http::proxy \"$http_proxy\";\n" \
         "Acquire::https::proxy \"$https_proxy\";" > "/etc/apt/apt.conf"

WORKDIR /opt/
COPY --from=builder /usr/local/lib/python3.10/dist-packages /usr/local/lib/python3.10/dist-packages
COPY --from=builder /opt/deepvariant/bazel-bin/licenses.zip .

# Copy over zip binaries.
COPY --from=builder \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/make_examples.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/call_variants.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/postprocess_variants.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/vcf_stats_report.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/show_examples.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/runtime_by_region_vis.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/multisample_make_examples.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/labeler/labeled_examples_to_vcf.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/convert_to_saved_model.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/make_examples_somatic.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/train.zip \
      /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/fast_pipeline \
      /opt/deepvariant/scripts/run_deepvariant.py \
      /opt/deepvariant/bin/

# Create shell wrappers for python zip files for easier use.
RUN \
  BASH_HEADER='#!/bin/bash' && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/make_examples.zip "$@"' > \
    /opt/deepvariant/bin/make_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/call_variants.zip "$@"' > \
    /opt/deepvariant/bin/call_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/postprocess_variants.zip "$@"' > \
    /opt/deepvariant/bin/postprocess_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/vcf_stats_report.zip "$@"' > \
    /opt/deepvariant/bin/vcf_stats_report && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/show_examples.zip "$@"' > \
    /opt/deepvariant/bin/show_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/runtime_by_region_vis.zip "$@"' > \
    /opt/deepvariant/bin/runtime_by_region_vis && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/multisample_make_examples.zip "$@"' > \
    /opt/deepvariant/bin/multisample_make_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 -u /opt/deepvariant/bin/labeled_examples_to_vcf.zip "$@"' > \
    /opt/deepvariant/bin/labeled_examples_to_vcf && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/convert_to_saved_model.zip "$@"' > \
    /opt/deepvariant/bin/convert_to_saved_model && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 -u /opt/deepvariant/bin/make_examples_somatic.zip "$@"' > \
    /opt/deepvariant/bin/make_examples_somatic && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 -u /opt/deepvariant/bin/run_deepvariant.py "$@"' > \
    /opt/deepvariant/bin/run_deepvariant && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    '/usr/bin/python3 /opt/deepvariant/bin/train.zip "$@"' > \
    /opt/deepvariant/bin/train && \
  chmod -R +x /opt/deepvariant/bin

# Copy over models
COPY --from=download_models /opt/models /opt/models
COPY --from=download_models /opt/smallmodels /opt/smallmodels

ENV PATH="${PATH}":${DV_BIN_PATH}:/opt/conda/envs/bio/bin
# This to use Keras 2.x with TF 2.16.1 or higher.
ENV TF_USE_LEGACY_KERAS=1
WORKDIR /opt/deepvariant

CMD ["/opt/deepvariant/bin/run_deepvariant", "--help"]
