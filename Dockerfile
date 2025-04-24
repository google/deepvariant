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
# $ sudo docker build --build-arg=FROM_IMAGE=nvidia/cuda:11.8.0-cudnn8-devel-ubuntu22.04 --build-arg=DV_GPU_BUILD=1 -t deepvariant_gpu .


ARG FROM_IMAGE=ubuntu:22.04
# PYTHON_VERSION is also set in settings.sh.
ARG PYTHON_VERSION=3.10
ARG DV_GPU_BUILD=0
ARG VERSION=1.9.0
ARG TF_ENABLE_ONEDNN_OPTS=1

FROM condaforge/miniforge3:24.9.2-0 as conda_setup
RUN conda config --add channels bioconda
RUN conda create -n bio \
                    bioconda::bcftools=1.15 \
                    bioconda::samtools=1.15 \
    && conda clean -a

FROM ${FROM_IMAGE} as builder
COPY --from=conda_setup /opt/conda /opt/conda
LABEL maintainer="https://github.com/google/deepvariant/issues"

ARG DV_GPU_BUILD
ENV DV_GPU_BUILD=${DV_GPU_BUILD}
ENV DV_BIN_PATH=/opt/deepvariant/bin

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
ARG TF_ENABLE_ONEDNN_OPTS
ENV DV_GPU_BUILD=${DV_GPU_BUILD}
ENV VERSION ${VERSION}
ENV PYTHON_VERSION ${PYTHON_VERSION}
ENV TF_ENABLE_ONEDNN_OPTS ${TF_ENABLE_ONEDNN_OPTS}

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
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/postprocess_variants.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/vcf_stats_report.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/show_examples.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/runtime_by_region_vis.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/multisample_make_examples.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/labeler/labeled_examples_to_vcf.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/convert_to_saved_model.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/make_examples_somatic.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/train.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/fast_pipeline .
COPY --from=builder /opt/deepvariant/scripts/run_deepvariant.py .

RUN ./run-prereq.sh

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python${PYTHON_VERSION} 0 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python${PYTHON_VERSION} 0

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
  chmod +x /opt/deepvariant/bin/make_examples \
    /opt/deepvariant/bin/call_variants \
    /opt/deepvariant/bin/postprocess_variants \
    /opt/deepvariant/bin/vcf_stats_report \
    /opt/deepvariant/bin/show_examples \
    /opt/deepvariant/bin/runtime_by_region_vis \
    /opt/deepvariant/bin/multisample_make_examples \
    /opt/deepvariant/bin/run_deepvariant \
    /opt/deepvariant/bin/labeled_examples_to_vcf \
    /opt/deepvariant/bin/convert_to_saved_model \
    /opt/deepvariant/bin/make_examples_somatic \
    /opt/deepvariant/bin/train

# Copy models
WORKDIR /opt/models/wgs
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wgs.savedmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wgs.savedmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wgs.savedmodel/example_info.json .
WORKDIR /opt/models/wgs/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wgs.savedmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wgs.savedmodel/variables/variables.index .
RUN chmod -R +r /opt/models/wgs/*

WORKDIR /opt/models/wes
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wes.savedmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wes.savedmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wes.savedmodel/example_info.json .
WORKDIR /opt/models/wes/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wes.savedmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.wes.savedmodel/variables/variables.index .
RUN chmod -R +r /opt/models/wes/*

WORKDIR /opt/models/pacbio
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.pacbio.savedmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.pacbio.savedmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.pacbio.savedmodel/example_info.json .
WORKDIR /opt/models/pacbio/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.pacbio.savedmodel/variables/variables.index .
RUN chmod -R +r /opt/models/pacbio/*

WORKDIR /opt/models/hybrid_pacbio_illumina
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.hybrid.savedmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.hybrid.savedmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.hybrid.savedmodel/example_info.json .
WORKDIR /opt/models/hybrid_pacbio_illumina/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.hybrid.savedmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.hybrid.savedmodel/variables/variables.index .
RUN chmod -R +r /opt/models/hybrid_pacbio_illumina/*

WORKDIR /opt/models/ont_r104
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.ont.savedmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.ont.savedmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.ont.savedmodel/example_info.json .
WORKDIR /opt/models/ont_r104/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.ont.savedmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.ont.savedmodel/variables/variables.index .
RUN chmod -R +r /opt/models/ont_r104/*

WORKDIR /opt/models/masseq
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.masseq.savedmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.masseq.savedmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.masseq.savedmodel/example_info.json .
WORKDIR /opt/models/masseq/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.masseq.savedmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/savedmodels/deepvariant.masseq.savedmodel/variables/variables.index .
RUN chmod -R +r /opt/models/masseq/*

# Copy small models
WORKDIR /opt/smallmodels/wgs
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.wgs.smallmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.wgs.smallmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.wgs.smallmodel/keras_metadata.pb .
WORKDIR /opt/smallmodels/wgs/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.wgs.smallmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.wgs.smallmodel/variables/variables.index .
RUN chmod -R +r /opt/smallmodels/wgs/*

WORKDIR /opt/smallmodels/pacbio
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.pacbio.smallmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.pacbio.smallmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.pacbio.smallmodel/keras_metadata.pb .
WORKDIR /opt/smallmodels/pacbio/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.pacbio.smallmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.pacbio.smallmodel/variables/variables.index .
RUN chmod -R +r /opt/smallmodels/pacbio/*

WORKDIR /opt/smallmodels/ont_r104
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.ont.smallmodel/fingerprint.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.ont.smallmodel/saved_model.pb .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.ont.smallmodel/keras_metadata.pb .
WORKDIR /opt/smallmodels/ont_r104/variables
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.ont.smallmodel/variables/variables.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepVariant/${VERSION}/smallmodels/deepvariant.ont.smallmodel/variables/variables.index .
RUN chmod -R +r /opt/smallmodels/ont_r104/*

ENV PATH="${PATH}":/opt/conda/bin:/opt/conda/envs/bio/bin:/opt/deepvariant/bin

RUN apt-get -y update && \
  apt-get install -y parallel python3-pip && \
  PATH="${HOME}/.local/bin:$PATH" python3 -m pip install absl-py==0.13.0 && \
  apt-get clean autoclean && \
  apt-get autoremove -y --purge && \
  rm -rf /var/lib/apt/lists/*

WORKDIR /opt/deepvariant

CMD ["/opt/deepvariant/bin/run_deepvariant", "--help"]
