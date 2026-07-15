#!/usr/bin/env bash
# Wrap a command with `powermetrics` GPU/ANE residency capture.
#
# Usage:
#     sudo ./tools/verify_gpu.sh path/to/log.json -- ./build/deepvariant call_variants ...
#
# powermetrics requires root, so this is intended to be invoked under sudo
# during validation runs. It writes a JSON summary to the first arg, then runs
# the rest of the command line. Pure shell + awk — no Python.
set -euo pipefail

if [[ $# -lt 3 ]]; then
  echo "usage: $0 <out.json> -- <command...>" >&2
  exit 2
fi

OUT="$1"; shift
[[ "$1" == "--" ]] || { echo "expected '--' separator" >&2; exit 2; }
shift

TMP_PM="$(mktemp -t verify_gpu.pm.XXXXXX)"
trap 'rm -f "${TMP_PM}"' EXIT

if [[ "$(id -u)" -ne 0 ]]; then
  echo "warning: not running as root; powermetrics will fail. re-run with sudo." >&2
fi

# Start powermetrics in background.
powermetrics --samplers gpu_power,ane_power -i 500 -f text > "${TMP_PM}" 2>/dev/null &
PM_PID=$!

# Run the workload.
RC=0
"$@" || RC=$?

# Stop powermetrics.
kill -INT "${PM_PID}" 2>/dev/null || true
wait "${PM_PID}" 2>/dev/null || true

# Summarise with awk — no interpreter required beyond /usr/bin/awk.
awk '
  /^GPU Power:/ {
    n_gpu++; sum_gpu += $3;
    if ($3 > max_gpu) max_gpu = $3;
    if ($3 > 50) act_gpu++;
    next
  }
  /^ANE Power:/ {
    n_ane++; sum_ane += $3;
    if ($3 > max_ane) max_ane = $3;
    if ($3 > 50) act_ane++;
    next
  }
  END {
    function blk(label, n, sum, mx, act,    mean, pct) {
      mean = (n > 0) ? sum / n : 0
      pct  = (n > 0) ? 100.0 * act / n : 0
      printf "  \"%s\": {\"n\": %d, \"mean_mw\": %.1f, \"max_mw\": %d, \"active_pct\": %.1f}", \
        label, n, mean, mx, pct
    }
    print "{"
    blk("gpu", n_gpu, sum_gpu, max_gpu, act_gpu); print ","
    blk("ane", n_ane, sum_ane, max_ane, act_ane); print ""
    print "}"
  }
' "${TMP_PM}" | tee "${OUT}"

exit "${RC}"
