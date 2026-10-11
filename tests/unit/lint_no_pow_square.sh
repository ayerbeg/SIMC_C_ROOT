#!/usr/bin/env bash
# tests/unit/lint_no_pow_square.sh - no std::pow(x, 2) in the library/tool sources (task P10.4, D-044).
# gfortran computes x**2 as x*x. gcc <= 15 folded pow(x, 2.) into x*x as well, but gcc 16 keeps the libm pow() call
# under -fmath-errno (the default in C++), and pow() can differ from x*x by 1 ulp: write squares as x*x (or sq(x)).
# Usage: lint_no_pow_square.sh SOURCE_DIR
set -uo pipefail
cd "$1" || exit 2
RE='pow[[:space:]]*\([^;]*,[[:space:]]*2(\.0*)?[[:space:]]*\)'
hits=$(grep -rnE "$RE" --include='*.cpp' --include='*.h' --include='*.hpp' src include tools 2>/dev/null |
       while IFS= read -r l; do code="${l#*:*:}"; code="${code%%//*}"; grep -qE "$RE" <<<"$code" && echo "$l"; done)
if [[ -n "$hits" ]]; then echo "std::pow(x, 2) found (write x*x, see the header of this script):"; echo "$hits"; exit 1; fi
echo "no pow(x, 2) in src/ include/ tools/"
