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
# $ sudo docker build --build-arg=FROM_IMAGE=nvidia/cuda:11.3.1-cudnn8-devel-ubuntu20.04 --build-arg=DV_GPU_BUILD=1 -t deepvariant_gpu .


ARG FROM_IMAGE=ubuntu:20.04
# PYTHON_VERSION is also set in settings.sh.
ARG PYTHON_VERSION=3.8
ARG DV_GPU_BUILD=0
ARG VERSION=1.5.0

FROM continuumio/miniconda3 as conda_setup
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge
RUN conda create -n bio \
                    bioconda::bcftools=1.10 \
                    bioconda::samtools=1.10 \
    && conda clean -a

FROM ${FROM_IMAGE} as builder
COPY --from=conda_setup /opt/conda /opt/conda
LABEL maintainer="https://github.com/google/deepvariant/issues"

ARG DV_GPU_BUILD
ENV DV_GPU_BUILD=${DV_GPU_BUILD}

# Copying DeepVariant source code
COPY . /opt/deepvariant

ARG VERSION
ENV VERSION=${VERSION}

WORKDIR /opt/deepvariant

RUN echo "Acquire::http::proxy \"$http_proxy\";\n" \
         "Acquire::https::proxy \"$https_proxy\";" > "/etc/apt/apt.conf"

RUN ./build-prereq.sh \
  && PATH="${HOME}/bin:${PATH}" ./build_release_binaries.sh  # PATH for bazel

FROM ${FROM_IMAGE}
ARG DV_GPU_BUILD
ARG VERSION
ARG PYTHON_VERSION
ENV DV_GPU_BUILD=${DV_GPU_BUILD}
ENV VERSION ${VERSION}
ENV PYTHON_VERSION ${PYTHON_VERSION}

RUN echo "Acquire::http::proxy \"$http_proxy\";\n" \
         "Acquire::https::proxy \"$https_proxy\";" > "/etc/apt/apt.conf"

WORKDIR /opt/
COPY --from=builder /opt/deepvariant/bazel-bin/licenses.zip .

WORKDIR /opt/deepvariant/bin/
COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /opt/deepvariant/run-prereq.sh .
COPY --from=builder /opt/deepvariant/settings.sh .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/make_examples.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/call_variants.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/call_variants_slim.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/postprocess_variants.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/vcf_stats_report.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/show_examples.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/runtime_by_region_vis.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/multisample_make_examples.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/labeler/labeled_examples_to_vcf.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/make_examples_somatic.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/scripts/run_deepvariant.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/scripts/run_deepsomatic.zip  .

RUN ./run-prereq.sh

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python${PYTHON_VERSION} 0 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python${PYTHON_VERSION} 0

# Create shell wrappers for python zip files for easier use.
RUN \
  BASH_HEADER='#!/bin/bash' && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/make_examples.zip "$@"' > \
    /opt/deepvariant/bin/make_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/call_variants.zip "$@"' > \
    /opt/deepvariant/bin/call_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/call_variants_slim.zip "$@"' > \
    /opt/deepvariant/bin/call_variants_slim && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/postprocess_variants.zip "$@"' > \
    /opt/deepvariant/bin/postprocess_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/vcf_stats_report.zip "$@"' > \
    /opt/deepvariant/bin/vcf_stats_report && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/show_examples.zip "$@"' > \
    /opt/deepvariant/bin/show_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/runtime_by_region_vis.zip "$@"' > \
    /opt/deepvariant/bin/runtime_by_region_vis && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/multisample_make_examples.zip "$@"' > \
    /opt/deepvariant/bin/multisample_make_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 -u /opt/deepvariant/bin/labeled_examples_to_vcf.zip "$@"' > \
    /opt/deepvariant/bin/labeled_examples_to_vcf && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 -u /opt/deepvariant/bin/make_examples_somatic.zip "$@"' > \
    /opt/deepvariant/bin/make_examples_somatic && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 -u /opt/deepvariant/bin/run_deepvariant.zip "$@"' > \
    /opt/deepvariant/bin/run_deepvariant && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 -u /opt/deepvariant/bin/run_deepsomatic.zip "$@"' > \
    /opt/deepvariant/bin/run_deepsomatic && \
  chmod +x /opt/deepvariant/bin/make_examples \
    /opt/deepvariant/bin/call_variants \
    /opt/deepvariant/bin/call_variants_slim \
    /opt/deepvariant/bin/postprocess_variants \
    /opt/deepvariant/bin/vcf_stats_report \
    /opt/deepvariant/bin/show_examples \
    /opt/deepvariant/bin/runtime_by_region_vis \
    /opt/deepvariant/bin/multisample_make_examples \
    /opt/deepvariant/bin/run_deepvariant \
    /opt/deepvariant/bin/run_deepsomatic \
    /opt/deepvariant/bin/labeled_examples_to_vcf \
    /opt/deepvariant/bin/make_examples_somatic

# Copy models
WORKDIR /opt/models/wgs
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.example_info.json .
RUN chmod +r /opt/models/wgs/model.ckpt*

WORKDIR /opt/models/wes
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.example_info.json .
RUN chmod +r /opt/models/wes/model.ckpt*

WORKDIR /opt/models/pacbio
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.example_info.json .
RUN chmod +r /opt/models/pacbio/model.ckpt*

WORKDIR /opt/models/hybrid_pacbio_illumina
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-hybrid_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-hybrid_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-hybrid_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-hybrid_standard/model.ckpt.example_info.json .
RUN chmod +r /opt/models/hybrid_pacbio_illumina/model.ckpt*

WORKDIR /opt/models/ont_r104
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-ont_r104/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-ont_r104/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-ont_r104/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-ont_r104/model.ckpt.example_info.json .

RUN chmod +r /opt/models/ont_r104/model.ckpt*

ENV PATH="${PATH}":/opt/conda/bin:/opt/conda/envs/bio/bin:/opt/deepvariant/bin

RUN apt-get -y update && \
  apt-get install -y parallel python3-pip && \
  PATH="${HOME}/.local/bin:$PATH" python3 -m pip install absl-py==0.13.0 && \
  apt-get clean autoclean && \
  apt-get autoremove -y --purge && \
  rm -rf /var/lib/apt/lists/*

WORKDIR /opt/deepvariant

CMD ["/opt/deepvariant/bin/run_deepvariant", "--help"]
