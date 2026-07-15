#!/usr/bin/env bash
# Submit the signed deepvariant binary to Apple's notary service and
# staple the ticket. Run AFTER ./release/sign.sh.
#
# Requires:
#   - Apple Developer account credentials stored as a notarytool keychain
#     profile named "deepvariant-notary":
#       xcrun notarytool store-credentials deepvariant-notary \
#         --apple-id <APPLE_ID> --team-id <TEAMID> --password <APP_SPECIFIC_PWD>
#   - Tools: xcrun notarytool, ditto, stapler
#
# Usage: ./release/notarize.sh path/to/deepvariant

set -euo pipefail
BIN="${1:?usage: $0 <path/to/deepvariant>}"
PROFILE="${NOTARY_PROFILE:-deepvariant-notary}"

WORK="$(mktemp -d)"
trap 'rm -rf "${WORK}"' EXIT

echo "==> Packaging ${BIN} into ${WORK}/deepvariant.zip"
ditto -c -k --keepParent "${BIN}" "${WORK}/deepvariant.zip"

echo "==> Submitting to Apple notary (xcrun notarytool, profile=${PROFILE})"
xcrun notarytool submit "${WORK}/deepvariant.zip" \
  --keychain-profile "${PROFILE}" \
  --wait

echo "==> Stapling the ticket onto ${BIN}"
xcrun stapler staple "${BIN}"

echo "==> Final Gatekeeper check"
spctl --assess --verbose "${BIN}"
echo "==> done — notarised + stapled"
