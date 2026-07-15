#!/usr/bin/env bash
# GBZ → BAM extractor for pangenome-aware DV testing.
#
# v2's native binary supports BAM-only --pangenome (GBZ requires the
# upstream load_gbz_into_shared_memory preprocessor + gbwt/gbwtgraph
# deps not in Homebrew). To test against Docker output without GBZ
# runtime support, pre-extract the synthetic reads from the GBZ once.
#
# Usage:
#   GBZ=/path/to/hprc.gbz REF=/path/to/GRCh38.fa REFS_BAM=/path/to/HG003.bam \
#     OUT=/path/to/pangenome.bam REGION=chr20:10000000-10100000 \
#     ./dump_gbz_to_bam.sh
set -euo pipefail
GBZ="${GBZ:?required: GBZ}"
REGION="${REGION:?required: REGION (chrN:start-end)}"
OUT="${OUT:?required: OUT (output BAM path)}"
REFS_BAM="${REFS_BAM:?required: REFS_BAM (template BAM for header)}"
DV_VERSION="${DV_VERSION:-pangenome_aware_deepvariant-1.10.0}"

GBZ_DIR="$(cd "$(dirname "${GBZ}")" && pwd)"
GBZ_BASE="$(basename "${GBZ}")"
OUT_DIR="$(cd "$(dirname "${OUT}")" && pwd)"
OUT_BASE="$(basename "${OUT}")"
REFS_DIR="$(cd "$(dirname "${REFS_BAM}")" && pwd)"
REFS_BASE="$(basename "${REFS_BAM}")"

# Convert region "chr:start-end" to start/end ints.
CHR="${REGION%:*}"; SE="${REGION#*:}"; START="${SE%-*}"; END="${SE#*-}"

PY="$(dirname "$0")/dump_gbz_to_bam.py"

docker run --rm \
  -v "${GBZ_DIR}:/gbz:ro" \
  -v "${REFS_DIR}:/refs:ro" \
  -v "${OUT_DIR}:/out" \
  -v "$(dirname "${PY}"):/scripts:ro" \
  -v /tmp:/tmp \
  -e GBZ="/gbz/${GBZ_BASE}" \
  -e REFS_BAM="/refs/${REFS_BASE}" \
  -e OUT="/out/${OUT_BASE}" \
  -e CHR="${CHR}" -e START="${START}" -e END="${END}" \
  "google/deepvariant:${DV_VERSION}" bash -c '
  cd /tmp && rm -rf dv_src && mkdir dv_src && cd dv_src && \
    unzip -o -q /opt/deepvariant/bin/make_examples_pangenome_aware_dv.zip
  pip install --quiet pysam 2>&1 | tail -1
  PYTHONPATH=/tmp/dv_src/runfiles/com_google_protobuf/python:\
/tmp/dv_src/runfiles/com_google_deepvariant:\
/tmp/dv_src/runfiles \
    GBZ="$GBZ" REFS_BAM="$REFS_BAM" OUT="$OUT" \
    CHR="$CHR" START="$START" END="$END" \
    python3 /scripts/dump_gbz_to_bam.py
'
