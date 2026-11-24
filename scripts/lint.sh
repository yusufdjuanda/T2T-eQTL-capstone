#!/usr/bin/env bash
# Static checks, run by CI (.github/workflows/tests.yml) and by `pixi run -e lint lint`:
#   1. every Python file compiles (syntax);
#   2. every shell script parses (bash -n) and ShellCheck finds no errors (warnings are not enforced);
#   3. no committed file hard-codes a cluster path (those belong in config/site.config, which is not committed);
#   4. the tables in metadata/ are well formed (same number of columns on every line).
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.."
fail=0

echo "1. Python syntax"
python -m compileall -q scripts tests workflow || fail=1

echo "2. shell syntax and ShellCheck (errors only)"
mapfile -t sh < <(git ls-files '*.sh')
for f in "${sh[@]}"; do bash -n "$f" || { echo "   $f"; fail=1; }; done
shellcheck --shell=bash --severity=error "${sh[@]}" || fail=1

echo "3. hard-coded cluster paths"
if git grep -nE '/storage/|/scratch/|/home/[a-z]' -- scripts workflow config tests reports ':!scripts/lint.sh'; then
  echo "   put site paths in config/site.config instead"
  fail=1
fi

echo "4. metadata tables"
for f in metadata/*.tsv; do
  awk -F'\t' -v f="$f" 'NR == 1 {n = NF} NF != n {print "   " f " line " NR ": " NF " columns, expected " n; bad = 1}
    END {exit bad}' "$f" || fail=1
done

if [[ $fail -eq 0 ]]; then echo "all checks passed"; else echo "some checks failed"; exit 1; fi
