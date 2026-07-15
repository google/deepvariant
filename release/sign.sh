#!/usr/bin/env bash
# Code-sign the deepvariant binary with the user's Developer ID.
#
# Requires:
#   - Apple Developer ID Application certificate installed in keychain
#   - $DEVELOPER_ID set to the cert's Common Name, e.g.
#       "Developer ID Application: Benjamin Demaille (TEAMID)"
#
# Usage: ./release/sign.sh path/to/deepvariant

set -euo pipefail
BIN="${1:?usage: $0 <path/to/deepvariant>}"
ID="${DEVELOPER_ID:?error: set DEVELOPER_ID env var to the certificate Common Name}"

echo "==> codesign --force --options=runtime --timestamp ${BIN}"
codesign \
  --force \
  --options=runtime \
  --timestamp \
  --sign "${ID}" \
  "${BIN}"

echo "==> verify"
codesign --verify --deep --strict --verbose=2 "${BIN}"

echo "==> spctl assess (Gatekeeper)"
spctl --assess --verbose "${BIN}" || {
  echo "  (spctl will fail until notarisation completes — expected at this stage)"
}
echo "==> done — signed in place: ${BIN}"
