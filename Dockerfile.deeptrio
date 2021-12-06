# Copyright 2019 Google LLC.
# This is used to build the DeepTrio release docker image.
# It can also be used to build local images, especially if you've made changes
# to the code.
# Example command:
# $ git clone https://github.com/google/deepvariant.git
# $ cd deepvariant
# $ sudo docker build -f Dockerfile.deeptrio -t deeptrio .
#
# To build for GPU, use a command like:
# $ sudo docker build -f Dockerfile.deeptrio --build-arg=FROM_IMAGE=nvidia/cuda:11.3.0-cudnn8-devel-ubuntu20.04 --build-arg=DV_GPU_BUILD=1 -t deeptrio_gpu .


ARG FROM_IMAGE=ubuntu:20.04
# PYTHON_VERSION is also set in settings.sh.
ARG PYTHON_VERSION=3.8
ARG DV_GPU_BUILD=0
ARG VERSION_DEEPTRIO=1.3.0

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

WORKDIR /opt/deepvariant

RUN ./build-prereq.sh \
  && PATH="${HOME}/bin:${PATH}" ./build_release_binaries.sh  # PATH for bazel

FROM ${FROM_IMAGE}
ARG DV_GPU_BUILD
ARG VERSION_DEEPTRIO
ARG PYTHON_VERSION
ENV DV_GPU_BUILD=${DV_GPU_BUILD}
ENV VERSION_DEEPTRIO ${VERSION_DEEPTRIO}
ENV PYTHON_VERSION ${PYTHON_VERSION}

WORKDIR /opt/
COPY --from=builder /opt/deepvariant/bazel-bin/licenses.zip .

WORKDIR /opt/deepvariant/bin/
COPY --from=builder /opt/conda /opt/conda
COPY --from=builder /opt/deepvariant/run-prereq.sh .
COPY --from=builder /opt/deepvariant/settings.sh .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/call_variants.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/postprocess_variants.zip  .
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deepvariant/runtime_by_region_vis.zip  .
COPY --from=builder /opt/deepvariant/scripts/run_deeptrio.py ./deeptrio/
COPY --from=builder /opt/deepvariant/bazel-out/k8-opt/bin/deeptrio/make_examples.zip  ./deeptrio/
RUN ./run-prereq.sh

RUN update-alternatives --install /usr/bin/python3 python3 /usr/bin/python${PYTHON_VERSION} 0 && \
    update-alternatives --install /usr/bin/python python /usr/bin/python${PYTHON_VERSION} 0

# Create shell wrappers for python zip files for easier use.
RUN \
  BASH_HEADER='#!/bin/bash' && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/deeptrio/make_examples.zip "$@"' > \
    /opt/deepvariant/bin/deeptrio/make_examples && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/call_variants.zip "$@"' > \
    /opt/deepvariant/bin/call_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/postprocess_variants.zip "$@"' > \
    /opt/deepvariant/bin/postprocess_variants && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 /opt/deepvariant/bin/runtime_by_region_vis.zip "$@"' > \
    /opt/deepvariant/bin/runtime_by_region_vis && \
  printf "%s\n%s\n" \
    "${BASH_HEADER}" \
    'python3 -u /opt/deepvariant/bin/deeptrio/run_deeptrio.py "$@"' > \
    /opt/deepvariant/bin/deeptrio/run_deeptrio && \
  chmod +x /opt/deepvariant/bin/deeptrio/make_examples \
    /opt/deepvariant/bin/call_variants \
    /opt/deepvariant/bin/postprocess_variants \
    /opt/deepvariant/bin/runtime_by_region_vis \
    /opt/deepvariant/bin/deeptrio/run_deeptrio

# Copy models
WORKDIR /opt/models/deeptrio/wgs/child
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_child_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_child_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_child_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_child_standard/model.ckpt.input_shape .
RUN chmod +r /opt/models/deeptrio/wgs/child/model.ckpt*

WORKDIR /opt/models/deeptrio/wgs/parent
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_parent_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_parent_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_parent_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wgs_parent_standard/model.ckpt.input_shape .
RUN chmod +r /opt/models/deeptrio/wgs/parent/model.ckpt*

WORKDIR /opt/models/deeptrio/pacbio/child
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_child_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_child_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_child_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_child_standard/model.ckpt.input_shape .
RUN chmod +r /opt/models/deeptrio/pacbio/child/model.ckpt*

WORKDIR /opt/models/deeptrio/pacbio/parent
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_parent_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_parent_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_parent_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-pacbio_parent_standard/model.ckpt.input_shape .
RUN chmod +r /opt/models/deeptrio/pacbio/parent/model.ckpt*

WORKDIR /opt/models/deeptrio/wes/child
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_child_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_child_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_child_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_child_standard/model.ckpt.input_shape .
RUN chmod +r /opt/models/deeptrio/wes/child/model.ckpt*

WORKDIR /opt/models/deeptrio/wes/parent
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_parent_standard/model.ckpt.data-00000-of-00001 .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_parent_standard/model.ckpt.index .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_parent_standard/model.ckpt.meta .
ADD https://storage.googleapis.com/deepvariant/models/DeepTrio/${VERSION_DEEPTRIO}/DeepTrio-inception_v3-${VERSION_DEEPTRIO}+data-wes_parent_standard/model.ckpt.input_shape .
RUN chmod +r /opt/models/deeptrio/wes/parent/model.ckpt*

ENV PATH="${PATH}":/opt/conda/bin:/opt/conda/envs/bio/bin:/opt/deepvariant/bin/deeptrio:/opt/deepvariant/bin

RUN apt-get -y update && \
  apt-get install -y parallel python3-pip && \
  PATH="${HOME}/.local/bin:$PATH" python3 -m pip install absl-py==0.13.0 && \
  apt-get clean autoclean && \
  apt-get autoremove -y --purge && \
  rm -rf /var/lib/apt/lists/*


WORKDIR /opt/deepvariant

CMD ["/opt/deepvariant/bin/deeptrio/run_deeptrio", "--help"]
