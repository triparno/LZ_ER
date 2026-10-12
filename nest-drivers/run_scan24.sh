#!/usr/bin/env bash
# 24-point energy scan of the ER leakage, LZ EFT2024 configuration.
#   ./run_scan24.sh [n_parallel] [path-to-lz_scan]     default 8, ./lz_scan
#   or:  bash run_scan24.sh 8 /path/to/nest/build/lz_scan
#   or:  LZ_SCAN_BIN=/path/to/lz_scan ./run_scan24.sh 8
#
# Writes to <repo>/scan/scan24_results.txt regardless of the working directory.
#
# Uses lz_scan (one energy per process). Each job gets a distinct seed --
# NEST's RandomGen is a singleton seeded deterministically, so without this
# every process would return identical numbers.
#
# ~1.1e10 trials total, about 2x the 8-point run.
# Trial counts target ~50 hits on P(Ne<409) below 85 keV, ~25 above,
# capped at 6e8 per point. The 79 keV point is a window edge: its
# P(window) is deliberately under-sampled, the edge is set by kinematics
# (Nq - 5435 < Ne < 409), not measured.

set -u
# Path to the compiled driver. Override with LZ_SCAN_BIN or the second argument.
BIN="${LZ_SCAN_BIN:-${2:-./lz_scan}}"
JOBS=${1:-8}
# Anchor the output at the repository root, not the working directory.
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/../scan/scan24_results.txt"
TMP="$(mktemp "${OUT}.XXXXXX")"
trap 'rm -f "$TMP"' EXIT

read -r -d '' JOBLIST << 'EOJ'
45      307000000  201
50      405000000  202
55      400000000  203
60      367000000  204
63      352000000  205
65      346000000  206
67      345000000  207
69      349000000  208
71      357000000  209
73      371000000  210
75      390000000  211
77      416000000  212
79      450000000  213
81      491000000  214
83      543000000  215
85      600000000  216
88      363000000  217
91      442000000  218
94      547000000  219
97      600000000  220
100     600000000  221
104     600000000  222
107     600000000  223
110     600000000  224
EOJ

if [ ! -x "$BIN" ]; then echo "error: $BIN not found or not executable" >&2; exit 1; fi

printf '#%6s %8s %7s %9s %8s %7s %12s %7s %11s %10s %6s %11s %10s\n' \
  E Nq NexNi meanNe sdNe skewNe N n_ceil P_ceil dP_ceil n_win P_win dP_win > "$TMP"

echo "$JOBLIST" | awk 'NF==3' \
  | xargs -P "$JOBS" -L1 bash -c '"'"$BIN"'" $0 $1 $2' >> "$TMP"

rows=$(grep -vc '^#' "$TMP")
if [ "$rows" -ne 24 ]; then
  echo "ERROR: $rows data rows, expected 24. $OUT left unchanged; partial output in $TMP" >&2
  trap - EXIT
  exit 1
fi
# sort the data rows by energy, keeping the header first (on TMP, not OUT)
{ head -1 "$TMP"; tail -n +2 "$TMP" | sort -g -k1,1; } > "$TMP.s" && mv "$TMP.s" "$TMP"
mv "$TMP" "$OUT"; trap - EXIT
echo "wrote $OUT"
if command -v column >/dev/null 2>&1; then column -t < "$OUT"; else cat "$OUT"; fi
