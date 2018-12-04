#!/bin/bash
# Copyright 2018 Google LLC.
set -euo pipefail

SCRIPT_DIR="$(dirname "$0")"
DV_GPU_BUILD=1 DV_INSTALL_GPU_DRIVERS=1 "${SCRIPT_DIR}"/run_wes_case_study_binaries.sh
