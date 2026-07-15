#!/usr/bin/env bash
# One-shot release build: clean + cmake + sign + notarize.
# Produces a release-ready ./build-macos/bin/deepvariant.

set -euo pipefail
cd "$(dirname "$0")/.."

echo "==> Clean build dir"
rm -rf build-macos

echo "==> Configure (Release)"
cmake -B build-macos -DCMAKE_BUILD_TYPE=Release

echo "==> Build (parallel)"
cmake --build build-macos --target deepvariant -j

echo "==> ctest"
ctest --test-dir build-macos --output-on-failure

if [[ -n "${DEVELOPER_ID:-}" ]]; then
  ./release/sign.sh build-macos/bin/deepvariant
  if [[ "${NOTARIZE:-no}" == "yes" ]]; then
    ./release/notarize.sh build-macos/bin/deepvariant
  else
    echo "==> NOTARIZE=yes not set; skipping Apple notary submission"
  fi
else
  echo "==> DEVELOPER_ID not set; skipping codesign"
fi

echo "==> Release artefact: $(pwd)/build-macos/bin/deepvariant"
ls -la build-macos/bin/deepvariant
otool -L build-macos/bin/deepvariant | head -8
