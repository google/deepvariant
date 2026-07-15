#!/usr/bin/env bash
# macOS arm64 build prerequisites for the DeepVariant Apple Silicon native port.
# Replaces upstream build-prereq.sh (which is Linux/Ubuntu-only).
#
# Idempotent: safe to re-run.
set -euo pipefail

if [[ "$(uname)" != "Darwin" || "$(uname -m)" != "arm64" ]]; then
  echo "error: this script targets macOS arm64 only" >&2
  exit 1
fi

if (( $(sw_vers -productVersion | cut -d. -f1) < 14 )); then
  echo "error: macOS 14 (Sonoma) or newer required" >&2
  exit 1
fi

echo "==> Xcode Command Line Tools"
if ! xcode-select -p >/dev/null 2>&1; then
  echo "    not installed; running xcode-select --install"
  xcode-select --install
  echo "    re-run this script after the CLT installer finishes"
  exit 1
fi

echo "==> Homebrew"
if ! command -v brew >/dev/null 2>&1; then
  echo "    Homebrew not found; install from https://brew.sh and re-run" >&2
  exit 1
fi

echo "==> brew dependencies"
# Build-time only; none of these end up in the shipped binary.
BREW_DEPS=(
  cmake
  ninja
  pkg-config
  pyenv
  git-lfs
  bash       # /usr/bin/bash on macOS is too old (3.2) for some scripts
)
for dep in "${BREW_DEPS[@]}"; do
  if brew list "${dep}" >/dev/null 2>&1; then
    echo "    ${dep}: ok"
  else
    echo "    installing ${dep}"
    brew install "${dep}"
  fi
done

echo "==> git-lfs hook"
git lfs install --skip-repo

echo "==> environment"
echo "    cmake:   $(cmake --version | head -1)"
echo "    ninja:   $(ninja --version)"
echo "    clang:   $(clang --version | head -1)"
echo "    pyenv:   $(pyenv --version)"
echo "    brew:    $(brew --version | head -1)"

echo "==> ready. next:"
# NB: use 'build-macos', not 'build' — the repo has a Bazel 'BUILD' file at the
# root and macOS's default case-insensitive filesystem treats a 'build/'
# directory as the same name, clobbering it. README/docs/release already
# standardize on 'build-macos'.
echo "    cmake -S . -B build-macos -G Ninja"
echo "    cmake --build build-macos --parallel"
echo "    ctest --test-dir build-macos --output-on-failure"
