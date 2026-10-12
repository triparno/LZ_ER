#!/usr/bin/env bash
# Does run_scan24.sh refuse a bad scan?
#
#   ./test_runner_guards.sh
#
# The scan is the one artefact in this archive that costs hours and is
# committed rather than regenerated, so the runner must never replace a good
# file with a bad one. This drives it with stub drivers that fail in each way
# a real one could, and asserts two things every time: the runner exits
# non-zero, and the existing scan file is untouched. The last case is a stub
# that behaves correctly, to show the guards do not reject a good run.
#
# Everything happens in a temporary tree; the committed scan is never at risk.

set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PASS=0; FAIL=0

# One energy per line, as the runner's own job list has them.
ENERGIES="45 50 55 60 63 65 67 69 71 73 75 77 79 81 83 85 88 91 94 97 100 104 107 110"

make_tree() {                      # make_tree -> echoes the temp root
    local root; root="$(mktemp -d)"
    mkdir -p "$root/nest-drivers" "$root/scan"
    cp "$HERE/run_scan24.sh" "$HERE/validate_scan.sh" "$root/nest-drivers/"
    printf 'SENTINEL -- a good scan file\n' > "$root/scan/scan24_results.txt"
    printf '%s\n' "$root"
}

# A stub driver takes "E N seed" and prints one 13-field row. MODE decides how
# it misbehaves. The trial count it echoes back is $2, which is what the
# runner's job list asked for, except where MODE says otherwise.
write_stub() {                     # write_stub ROOT MODE
    cat > "$1/nest-drivers/stub" <<STUB
#!/usr/bin/env bash
E=\$1; N=\$2
case "$2" in
  exit1)      printf '%7s %8.1f %7.4f %9.2f %8.2f %7.4f %12s %7d %11.5e %10.2e %6d %11.5e %10.2e\n' \\
                  "\$E" 4000 0.1812 1000 180 0.5 "\$N" 50 1.0e-7 1e-8 10 5e-9 1e-9; exit 1 ;;
  short)      [ "\$E" = "45" ] && exit 0
              printf '%7s %8.1f %7.4f %9.2f %8.2f %7.4f %12s %7d %11.5e %10.2e %6d %11.5e %10.2e\n' \\
                  "\$E" 4000 0.1812 1000 180 0.5 "\$N" 50 1.0e-7 1e-8 10 5e-9 1e-9 ;;
  wrongE)     [ "\$E" = "45" ] && E=46
              printf '%7s %8.1f %7.4f %9.2f %8.2f %7.4f %12s %7d %11.5e %10.2e %6d %11.5e %10.2e\n' \\
                  "\$E" 4000 0.1812 1000 180 0.5 "\$N" 50 1.0e-7 1e-8 10 5e-9 1e-9 ;;
  nan)        V=1.0e-7; [ "\$E" = "63" ] && V=nan
              printf '%7s %8.1f %7.4f %9.2f %8.2f %7.4f %12s %7d %11s %10.2e %6d %11.5e %10.2e\n' \\
                  "\$E" 4000 0.1812 1000 180 0.5 "\$N" 50 "\$V" 1e-8 10 5e-9 1e-9 ;;
  fields)     [ "\$E" = "71" ] && { printf '%7s %8.1f %7.4f\n' "\$E" 4000 0.1812; exit 0; }
              printf '%7s %8.1f %7.4f %9.2f %8.2f %7.4f %12s %7d %11.5e %10.2e %6d %11.5e %10.2e\n' \\
                  "\$E" 4000 0.1812 1000 180 0.5 "\$N" 50 1.0e-7 1e-8 10 5e-9 1e-9 ;;
  trials)     [ "\$E" = "85" ] && N=1000
              printf '%7s %8.1f %7.4f %9.2f %8.2f %7.4f %12s %7d %11.5e %10.2e %6d %11.5e %10.2e\n' \\
                  "\$E" 4000 0.1812 1000 180 0.5 "\$N" 50 1.0e-7 1e-8 10 5e-9 1e-9 ;;
  good)       printf '%7s %8.1f %7.4f %9.2f %8.2f %7.4f %12s %7d %11.5e %10.2e %6d %11.5e %10.2e\n' \\
                  "\$E" 4000 0.1812 1000 180 0.5 "\$N" 50 1.0e-7 1e-8 10 5e-9 1e-9 ;;
esac
STUB
    chmod +x "$1/nest-drivers/stub"
}

expect_refusal() {                 # expect_refusal MODE DESCRIPTION
    local mode="$1" what="$2" root rc
    root="$(make_tree)"; write_stub "$root" "$mode"
    ( cd "$root/nest-drivers" && bash run_scan24.sh 4 ./stub ) >/dev/null 2>&1
    rc=$?
    if [ "$rc" -eq 0 ]; then
        printf '  FAIL  %-34s runner exited 0\n' "$what"; FAIL=$((FAIL+1))
    elif ! grep -q SENTINEL "$root/scan/scan24_results.txt" 2>/dev/null; then
        printf '  FAIL  %-34s scan file was replaced\n' "$what"; FAIL=$((FAIL+1))
    else
        printf '  ok    %-34s refused, scan intact\n' "$what"; PASS=$((PASS+1))
    fi
    rm -rf "$root"
}

echo
echo "the runner must refuse:"
expect_refusal exit1  "a driver that exits non-zero"
expect_refusal short  "a missing energy"
expect_refusal wrongE "an energy not in the job list"
expect_refusal nan    "a non-finite field"
expect_refusal fields "a row with too few fields"
expect_refusal trials "fewer trials than were asked for"

echo
echo "and must accept:"
root="$(make_tree)"; write_stub "$root" good
if ( cd "$root/nest-drivers" && bash run_scan24.sh 4 ./stub ) >/dev/null 2>&1 \
   && [ "$(grep -vc '^#' "$root/scan/scan24_results.txt")" -eq 24 ] \
   && ! grep -q SENTINEL "$root/scan/scan24_results.txt"; then
    # and the rows must come out sorted by energy
    if [ "$(tail -n +2 "$root/scan/scan24_results.txt" | awk '{print $1}' | tr '\n' ' ')" \
         = "$ENERGIES " ]; then
        printf '  ok    %-34s accepted, 24 rows, sorted\n' "a correct run"; PASS=$((PASS+1))
    else
        printf '  FAIL  %-34s accepted but not sorted by energy\n' "a correct run"; FAIL=$((FAIL+1))
    fi
else
    printf '  FAIL  %-34s a correct run was refused\n' "a correct run"; FAIL=$((FAIL+1))
fi
rm -rf "$root"

echo
if [ "$FAIL" -eq 0 ]; then
    echo "$PASS guards hold."
    exit 0
fi
echo "$PASS passed, $FAIL FAILED."
exit 1
