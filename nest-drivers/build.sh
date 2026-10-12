#!/usr/bin/env bash
# Builds the three drivers against an existing NEST checkout, discovering the
# include directories and the library rather than assuming a layout.
#
#   ./build.sh /path/to/NEST [/path/to/NEST/build]
#   or:  NEST_SRC=... NEST_BUILD=... ./build.sh
#
# Use this in preference to adding CMake targets. The drivers are single
# translation units with one dependency, and NEST's own CMake does not
# necessarily export its include directories to targets that merely link the
# library -- which is the usual cause of "NEST.hh: No such file or directory".

set -u
SRC="${NEST_SRC:-${1:-}}"
BLD="${NEST_BUILD:-${2:-${SRC:-}/build}}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

if [ -z "$SRC" ] || [ ! -d "$SRC" ]; then
  echo "usage: $0 /path/to/NEST [/path/to/NEST/build]" >&2
  echo "       the first argument is the git checkout, the one containing" >&2
  echo "       include/ and src/, not the build directory." >&2
  exit 1
fi

# --- include directories ----------------------------------------------------
# Preferred route: reuse the flags CMake worked out when it compiled NEST. They
# already cover the vendored dependencies (gcem, and anything else NEST's own
# headers include), which a hand-written list of directories does not.
INCS=""
FLAGS=$(find "$BLD" -path '*.dir/flags.make' -print 2>/dev/null \
        | grep -i 'nest' | head -1)
[ -z "$FLAGS" ] && FLAGS=$(find "$BLD" -path '*.dir/flags.make' -print -quit 2>/dev/null || true)
if [ -n "$FLAGS" ]; then
  INCS=$(sed -n 's/^CXX_INCLUDES *= *//p' "$FLAGS" | head -1)
  [ -n "$INCS" ] && echo "include flags taken from $FLAGS"
fi

# Fallback: discover them. Covers the NEST headers, the detector headers, and
# every FetchContent dependency CMake unpacked under the build tree.
if [ -z "$INCS" ]; then
  HEADERS="NEST.hh RandomGen.hh VDetector.hh LZ_WS2024.hh gcem.hpp"
  MISSING=""
  for h in $HEADERS; do
    f=$(find "$SRC" "$BLD" -name "$h" -print -quit 2>/dev/null || true)
    if [ -z "$f" ]; then MISSING="$MISSING $h"; else
      d=$(dirname "$f")
      # a header included as <gcem.hpp> wants its own directory; one included
      # as <gcem/gcem.hpp> wants the parent. Add both and let the compiler pick.
      INCS="$INCS -I$d -I$(dirname "$d")"
    fi
  done
  for d in "$BLD"/_deps/*/include "$BLD"/_deps/*-src/include; do
    [ -d "$d" ] && INCS="$INCS -I$d"
  done
  if [ -n "$MISSING" ]; then
    echo "error: header(s) not found under $SRC or $BLD:$MISSING" >&2
    case "$MISSING" in
      *LZ_WS2024.hh*) cat >&2 <<'MSG'

LZ_WS2024.hh missing means the checkout predates the LZ configuration. The
numbers in this archive need commit 541b643; see BUILD.md.
MSG
        ;;
    esac
    case "$MISSING" in
      *gcem.hpp*) cat >&2 <<'MSG'

gcem.hpp is a header-only dependency that NEST's CMake fetches. If it is
nowhere under the build tree, NEST has not been configured yet: run
  mkdir -p BUILDDIR && cd BUILDDIR && cmake .. && make
and try again. cmake alone is enough to unpack it.
MSG
        ;;
    esac
    exit 1
  fi
  INCS=$(printf '%s\n' $INCS | sort -u | tr '\n' ' ')
fi

# --- library ----------------------------------------------------------------
LIB=$(find "$BLD" \( -name 'libNEST*.a' -o -name 'libNEST*.so' \
                     -o -name 'libNEST*.dylib' \) -print -quit 2>/dev/null || true)
if [ -z "$LIB" ]; then
  echo "error: no libNEST* under $BLD" >&2
  echo "       build NEST first:  mkdir -p $BLD && cd $BLD && cmake .. && make" >&2
  echo "       at commit 541b643 the library is libNESTCore.a" >&2
  exit 1
fi

echo "NEST source : $SRC"
echo "includes    :$INCS"
echo "library     : $LIB"
echo

STD="${CXXSTD:-c++17}"
EXTRA="${EXTRA_LIBS:-}"
rc=0
for d in lz_tail_test lz_scan lz_scan_s1s2; do
  [ -f "$HERE/$d.cpp" ] || { echo "skip $d (no $d.cpp here)"; continue; }
  if g++ -std="$STD" -O2 "$HERE/$d.cpp" -o "$HERE/$d" $INCS "$LIB" $EXTRA; then
    echo "built $HERE/$d"
  else
    echo "FAILED $d" >&2
    echo "  if the errors are undefined gsl_* symbols, re-run with" >&2
    echo "    EXTRA_LIBS=\"\$(gsl-config --libs)\" $0 $SRC $BLD" >&2
    rc=1
  fi
done

if [ $rc -eq 0 ]; then
  cat <<'MSG'

Next:
  ./lz_tail_test 1000000                 # verification, see HOWTO.md section 5
  ./lz_scan_s1s2 probe 70                # the two vector indices, section 8.1
MSG
fi
exit $rc
