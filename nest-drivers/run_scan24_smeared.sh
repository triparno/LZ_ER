#!/usr/bin/env bash
# The 24-point leakage scan carried through to S1c and S2c, so that the
# detection fluctuations omitted by run_scan24.sh are included.
#
#   ./run_scan24_smeared.sh [n_parallel] [path-to-lz_scan_s1s2] [s1c_idx] [s2c_idx]
#   or:  LZ_S1S2_BIN=/path/to/lz_scan_s1s2 ./run_scan24_smeared.sh 8 '' 3 5
#
# Writes <repo>/scan/scan24_smeared.txt. The existing scan24_results.txt is
# never touched, so the two can be compared.
#
# BEFORE THE FIRST RUN, get the two vector indices from the driver's probe
# mode -- the layout of GetS1/GetS2 differs between NEST versions and the
# driver will not guess:
#
#     ./lz_scan_s1s2 probe 70
#
# and pass the two indices it marks. Set them here once you know them.
#
# Same energies, trial counts and seeds as run_scan24.sh. The quanta-level
# columns will agree with scan24_results.txt statistically but not
# bit-for-bit: GetS1 and GetS2 draw from the same RandomGen singleton, so the
# GetQuanta sequence differs. The comparison that matters is paired within
# this file -- the same events classified both ways -- and that one is exact.
#
# Runtime is close to run_scan24.sh despite the extra work: the driver only
# calls GetS1/GetS2 for events near a cut boundary, and at a leakage of ~1e-7
# almost none are. The n_slow column reports how many took the slow path.

set -u
BIN="${LZ_S1S2_BIN:-${2:-./lz_scan_s1s2}}"
JOBS=${1:-8}
S1C_IDX=${3:-}
S2C_IDX=${4:-}

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="$HERE/../scan/scan24_smeared.txt"
TMP="$(mktemp "${OUT}.XXXXXX")"
trap 'rm -f "$TMP"' EXIT

if [ ! -x "$BIN" ]; then
  echo "error: $BIN not found or not executable" >&2; exit 1
fi
if [ -z "$S1C_IDX" ] || [ -z "$S2C_IDX" ]; then
  cat >&2 <<'MSG'
error: the S1c and S2c vector indices are required.

  Run the probe first:      ./lz_scan_s1s2 probe 70
  then pass what it marks:  ./run_scan24_smeared.sh 8 ./lz_scan_s1s2 <s1c> <s2c>

They are not defaulted on purpose. A wrong index would not crash; it would
quietly produce a scan of the wrong quantity.
MSG
  exit 1
fi

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

printf '#%6s %8s %7s %9s %8s %7s %12s %7s %11s %10s %6s %11s %10s %9s %11s %10s %8s %11s %10s %10s\n' \
  E Nq NexNi meanNe sdNe skewNe N \
  n_ceil P_ceil dP_ceil n_win P_win dP_win \
  n_ceil_s P_ceil_s dP_ceil_s n_win_s P_win_s dP_win_s n_slow > "$TMP"

echo "$JOBLIST" | awk 'NF==3' \
  | xargs -P "$JOBS" -L1 bash -c \
      '"'"$BIN"'" $0 $1 $2 '"$S1C_IDX $S2C_IDX" >> "$TMP"

rows=$(grep -vc '^#' "$TMP")
if [ "$rows" -ne 24 ]; then
  echo "ERROR: $rows data rows, expected 24. $OUT left unchanged;" \
       "partial output in $TMP" >&2
  trap - EXIT
  exit 1
fi

# sort the data rows by energy, keeping the header first (on TMP, not OUT)
{ head -1 "$TMP"; tail -n +2 "$TMP" | sort -g -k1,1; } > "$TMP.s" && mv "$TMP.s" "$TMP"
mv "$TMP" "$OUT"; trap - EXIT

echo "wrote $OUT"
if command -v column >/dev/null 2>&1; then column -t < "$OUT"; else cat "$OUT"; fi
