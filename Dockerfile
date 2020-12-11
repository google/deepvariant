# Copyright 2019 Google LLC.
# This is used to build the DeepVariant release docker image.
# It can also be used to build local images, especially if you've made changes
# to the code.
# Example command:
# $ git clone https://github.com/google/deepvariant.git
# $ cd deepvariant
# $ sudo docker build -t deepvariant . --build-arg DV_OPENVINO_BUILD=1
#
# To build for GPU, use a command like:
# $ sudo docker build --build-arg=FROM_IMAGE=nvidia/cuda:10.1-cudnn7-devel-ubuntu18.04 --build-arg=DV_GPU_BUILD=1 -t deepvariant_gpu .


ARG FROM_IMAGE=ubuntu:18.04
ARG DV_GPU_BUILD=0
ARG DV_OPENVINO_BUILD=0
ARG VERSION=1.1.0

FROM ${FROM_IMAGE} as builder
LABEL maintainer="https://github.com/google/deepvariant/issues"

ARG DV_GPU_BUILD
ENV DV_GPU_BUILD=${DV_GPU_BUILD}

# Copying DeepVariant source code
COPY . /opt/deepvariant

ARG DV_OPENVINO_BUILD
ENV DV_OPENVINO_BUILD=${DV_OPENVINO_BUILD}

ARG VERSION
ENV VERSION=${VERSION}

WORKDIR /opt/deepvariant

RUN echo "Acquire::http::proxy \"$http_proxy\";\n" \
         "Acquire::https::proxy \"$https_proxy\";" > "/etc/apt/apt.conf"

RUN ./build-prereq.sh \
  && PATH="${HOME}/bin:${PATH}" ./build_release_binaries.sh  # PATH for bazel

FROM ${FROM_IMAGE}
ARG DV_GPU_BUILD
ARG DV_OPENVINO_BUILD
ARG VERSION
ENV DV_GPU_BUILD=${DV_GPU_BUILD}
ENV DV_OPENVINO_BUILD=${DV_OPENVINO_BUILD}
ENV VERSION ${VERSION}

RUN echo "Acquire::http::proxy \"$http_proxy\";\n" \
         "Acquire::https::proxy \"$https_proxy\";" > "/etc/apt/apt.conf"

WORKDIR /opt/
COPY --from=builder /opt/deepvariant/bazel-bin/licenses.zip .

WORKDIR /opt/deepvariant/bin/
COPY --from=builder /opt/deepvariant/run-prereq.sh .
COPY --from=builder /opt/deepvariant/settings.sh .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/make_examples.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/call_variants.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/postprocess_variants.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/vcf_stats_report.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/show_examples.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/runtime_by_region_vis.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/model_train.zip .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/model_eval.zip  .
COPY --from=builder /opt/deepvariant/scripts/run_deepvariant.py .

RUN ./run-prereq.sh

# Create shell wrappers for python zip files for easier use.
RUN \
  BASH_HEADER='#!/bin/bash' && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/make_examples.zip "$@"' > \
    /opt/deepvariant/bin/make_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/call_variants.zip "$@"' > \
    /opt/deepvariant/bin/call_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/postprocess_variants.zip "$@"' > \
    /opt/deepvariant/bin/postprocess_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/model_train.zip "$@"' > \
    /opt/deepvariant/bin/model_train && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/model_eval.zip "$@"' > \
    /opt/deepvariant/bin/model_eval && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/vcf_stats_report.zip "$@"' > \
    /opt/deepvariant/bin/vcf_stats_report && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/show_examples.zip "$@"' > \
    /opt/deepvariant/bin/show_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python /opt/deepvariant/bin/runtime_by_region_vis.zip "$@"' > \
    /opt/deepvariant/bin/runtime_by_region_vis && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python -u /opt/deepvariant/bin/run_deepvariant.py "$@"' > \
    /opt/deepvariant/bin/run_deepvariant && \
  chmod +x /opt/deepvariant/bin/make_examples \
    /opt/deepvariant/bin/call_variants \
    /opt/deepvariant/bin/postprocess_variants \
    /opt/deepvariant/bin/vcf_stats_report \
    /opt/deepvariant/bin/show_examples \
    /opt/deepvariant/bin/runtime_by_region_vis \
    /opt/deepvariant/bin/model_train \
    /opt/deepvariant/bin/model_eval \
    /opt/deepvariant/bin/run_deepvariant

# Copy models
WORKDIR /opt/models/wgs
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wgs_standard/model.ckpt.meta .
RUN chmod +r /opt/models/wgs/model.ckpt*

WORKDIR /opt/models/wes
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-wes_standard/model.ckpt.meta .
RUN chmod +r /opt/models/wes/model.ckpt*

WORKDIR /opt/models/pacbio
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-pacbio_standard/model.ckpt.meta .
RUN chmod +r /opt/models/pacbio/model.ckpt*

WORKDIR /opt/models/hybrid_pacbio_illumina
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-hybrid_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-hybrid_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/DeepVariant-inception_v3-${VERSION}+data-hybrid_standard/model.ckpt.meta .
RUN chmod +r /opt/models/hybrid_pacbio_illumina/model.ckpt*

RUN apt-get -y update && \
  apt-get install -y parallel python3-pip && \
  PATH="${HOME}/.local/bin:$PATH" python3 -m pip install absl-py==0.8.1 && \
  apt-get clean autoclean && \
  apt-get autoremove -y --purge && \
  rm -rf /var/lib/apt/lists/*

ENV PATH="/opt/deepvariant/bin:${PATH}"

# Setup OpenVINO environment.
ENV LD_LIBRARY_PATH=/opt/intel/openvino/inference_engine/lib/intel64/:/opt/intel/openvino/inference_engine/external/tbb/lib/:/opt/intel/openvino/deployment_tools/ngraph/lib/:$LD_LIBRARY_PATH \
    PYTHONPATH=/opt/intel/openvino/python/python3.6/:/opt/intel/openvino/deployment_tools/model_optimizer/:$PYTHONPATH

WORKDIR /opt/deepvariant

CMD ["/opt/deepvariant/bin/run_deepvariant", "--help"]
