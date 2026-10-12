#!/usr/bin/env bash
# Runs every internal-consistency check in the archive and reports a summary.
#
#   ./verify.sh            the checks that need only Python
#   ./verify.sh --figures  those, plus regenerate the three figures
#
# None of this needs NEST: the scan output is committed. Rebuilding the scan
# itself is nest-drivers/HOWTO.md.

set -u
cd "$(dirname "${BASH_SOURCE[0]}")/analysis" || exit 1
FIGS=${1:-}
PASS=0; FAILED=""

# ---------------------------------------------------------------- preflight
# A partial or inconsistent install is the failure mode worth catching first.
# response.py holds the one definition of the response; a tree missing it, or
# carrying a copy of any script that predates it, will still run and will
# quietly answer on the wrong convention. Check before running anything.
NEED="response.py check_manuscript.py family_table.py figfam.py gap.py
      johnson.py smearing.py leakage.py fit_K.py printed_omega.py wvalue.py
      read_smeared_scan.py discriminate.py double_vacancy.py families.py
      window_ratio.py fig_ratio.py fig_leakage_scan.py"
IMPORTERS="check_manuscript.py family_table.py figfam.py gap.py johnson.py
           smearing.py discriminate.py leakage.py printed_omega.py fit_K.py
           fig_ratio.py fig_leakage_scan.py double_vacancy.py families.py
           window_ratio.py"
MISS=""; STALE=""; PRIVATE=""
for f in $NEED; do [ -f "$f" ] || MISS="$MISS $f"; done
for f in $IMPORTERS; do
  [ -f "$f" ] || continue
  grep -q "import response" "$f" || STALE="$STALE $f"
done
# Importing is not enough. gap.py imported response.py and still carried its
# own KNORM = 0.9566, and fit_K.py named its own ceiling 409.0 while everything
# else used another value -- a 1.75% error in every absolute number that the
# import check above passed straight over. So also refuse any private copy of
# the two constants that scale results: the normalisation and the ceiling.
for f in *.py; do
  case "$f" in response.py) continue ;; esac
  if grep -Eq '(^|[^._[:alnum:]])(K|KNORM|NE_CEIL|NLEAK|NLEAK0|DNDE|TOTBKG|WININT)[[:space:]]*=[[:space:]]*[0-9]' "$f" \
     || grep -q '10\*\*4\.15/g2' "$f"; then
    PRIVATE="$PRIVATE $f"
  fi
done
if [ -n "$MISS" ] || [ -n "$STALE" ] || [ -n "$PRIVATE" ]; then
  echo
  [ -n "$MISS" ]    && echo "MISSING from analysis/:$MISS"
  [ -n "$STALE" ]   && echo "STALE (does not import response.py):$STALE"
  [ -n "$PRIVATE" ] && echo "PRIVATE COPY of K or the ceiling:$PRIVATE"
  cat <<'MSG'

These scripts share one definition of the response, in response.py: omega from
RecombOmegaER unrounded, gamma_ER evaluated at the event-conditioned energy
rather than rounded to 1.9, the ceiling at the continuity-corrected 409.5, and
K = 0.9374 fit AT that ceiling. A script that does not import it carries its own
constants and will answer on a convention the paper no longer uses -- it will
run, and it will be wrong by a few per cent in the ratios and 11% in the
absolute tail probability. A private K or ceiling is refused even in a script
that does import response.py, because an import does not prevent a stale
constant sitting beside it.

Install the full analysis/ directory and re-run. Do not run the checks against
a partial tree: the ones that pass are not evidence that the rest agree.
MSG
  exit 1
fi

run () {                     # run <label> <script> [args...]
  local lbl="$1"; shift
  printf '  %-26s ' "$lbl"
  if out=$(timeout 900 python3 -I "$@" 2>&1); then
    printf 'ok\n'; PASS=$((PASS+1))
  else
    printf 'FAILED\n'; FAILED="$FAILED $lbl"
    printf '%s\n' "$out" | tail -15 | sed 's/^/      | /'
  fi
}

echo
echo "the two checks that can fail:"
run "manuscript vs code"  check_manuscript.py
run "quoted values"       leakage.py

echo
echo "the pieces they rest on:"
run "response reference"  response.py
run "normalisation K"     fit_K.py
run "printed vs impl Eq.4" printed_omega.py
run "OldW13eV channel"    wvalue.py
run "family table"        family_table.py
run "three families"      families.py
run "window weighting"    window_ratio.py
run "empty interval"      gap.py
run "Johnson S_U moments" johnson.py
run "smearing model"      smearing.py
run "discrimination floor" discriminate.py --selftest
run "calibration size"    discriminate.py
run "double-vacancy modes" double_vacancy.py
[ -f ../scan/scan24_smeared.txt ] \
  && run "measured smearing" read_smeared_scan.py \
  || echo "  measured smearing          skipped (no scan/scan24_smeared.txt)"

if [ "$FIGS" = "--figures" ]; then
  echo
  echo "figures:"
  run "lz_tail_family"       figfam.py
  run "lz_ratio_vs_deviation" fig_ratio.py
  run "lz_leakage_scan"      fig_leakage_scan.py
fi

echo
if [ -z "$FAILED" ]; then
  echo "$PASS checks passed; the manuscript and the code agree."
else
  echo "$PASS passed, FAILED:$FAILED"
  exit 1
fi
