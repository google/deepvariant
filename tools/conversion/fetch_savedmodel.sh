#!/usr/bin/env bash
# Pull a DeepVariant / DeepTrio / DeepSomatic SavedModel from
# gs://deepvariant/models/.../<NAME>.savedmodel/
#
# Usage:
#   ./fetch_savedmodel.sh <NAME>
#
# Supported NAMEs:
#   DeepVariant  : wgs, wes, pacbio, ont_r104, hybrid_pacbio_illumina,
#                  masseq, rnaseq, pangenome_aware_deepvariant
#   DeepTrio     : deeptrio.<wgs|wes|pacbio|ont>.<child|parent>
#                  (some variants ship parent1/parent2 separately —
#                   the script tries the listed pattern first)
#   DeepSomatic  : deepsomatic.<wgs|wes|pacbio|ont|ffpe_wgs|ffpe_wes>

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p models

NAME="${1:?usage: $0 <model_name>}"
DST="models/${NAME}"
DV_VERSION="${DV_VERSION:-1.10.0}"

if [[ -f "${DST}/saved_model.pb" ]]; then
  echo "==> ${DST} already populated, skipping"
  exit 0
fi
mkdir -p "${DST}/variables"

# Pick the upstream bucket path from the NAME prefix. Trio uses underscore
# between subtype and member (e.g. wgs_child), somatic uses underscore
# between subtype and tumor_only (e.g. wgs_tumor_only).
case "${NAME}" in
  deeptrio.*)
    SUITE="DeepTrio"
    SHORT="${NAME#deeptrio.}"
    SAVEDMODEL="deeptrio.${SHORT}.savedmodel"
    ;;
  deepsomatic.*)
    SUITE="DeepSomatic"
    SHORT="${NAME#deepsomatic.}"
    SAVEDMODEL="deepsomatic.${SHORT}.savedmodel"
    ;;
  *)
    SUITE="DeepVariant"
    SAVEDMODEL="deepvariant.${NAME}.savedmodel"
    ;;
esac

BASE="https://storage.googleapis.com/deepvariant/models/${SUITE}/${DV_VERSION}/savedmodels/${SAVEDMODEL}"

FILES=(
  "saved_model.pb"
  "fingerprint.pb"
  "model.example_info.json"
  "variables/variables.data-00000-of-00001"
  "variables/variables.index"
)

for f in "${FILES[@]}"; do
  url="${BASE}/${f}"
  out="${DST}/${f}"
  echo "==> fetching ${url}"
  if ! curl -fL --retry 3 --connect-timeout 15 -o "${out}" "${url}"; then
    case "${f}" in
      fingerprint.pb|model.example_info.json)
        echo "    (${f} optional, skipping)"
        rm -f "${out}"
        ;;
      *)
        echo "error: failed to fetch ${url}" >&2
        exit 1
        ;;
    esac
  fi
done

echo "==> ${NAME} SavedModel ready at ${DST}"
ls "${DST}/" "${DST}/variables/" 2>/dev/null
