#!/bin/bash
# Stage 04 verification of module experiments, Amendment D96, check O27 through the command line: `bin/chaoslab.jl run --round-log` (the runner, not patched)
# by rounds, with one count_round (verify/experiments-d96.jl count) between the runs, the closure run and a resume run; then the CLI log against the library log.
# Usage, from the repo root, after `verify/experiments-d96.jl lib <fixture dir> <out dir>`:  bash verify/experiments-d96-cli.sh <fixture dir> <out dir>
set -u
FIX="$1"; OUT="$2"
PRE="$FIX/inputs/p3-abstention-preregistration.sealed.md"; PANEL="$FIX/p3-abstention-test-panel.toml"
S="$OUT/store-cli"; M="$OUT/manifest-cli.json"; L="$OUT/round-log-cli.txt"
case "$OUT" in *"/data/"*|*"\\data\\"*) echo "refused: out dir inside data/"; exit 1;; esac
[ -e "$S" ] || [ -e "$M" ] || [ -e "$L" ] && { echo "refused: CLI outputs exist (use a fresh out dir)"; exit 1; }
mkdir -p "$S"; : > "$L"
runcli() {
  julia --project -t 4 bin/chaoslab.jl run "$PRE" "$PANEL" --store "$S" --manifest "$M" --round-log "$L" --threads 4 --allow-dirty > "$OUT/cli-last.txt" 2>&1; local e=$?
  tail -3 "$OUT/cli-last.txt"
  echo "O27 | CLI run $1 | exit $e | exit 0 | $([ "$e" = 0 ] && echo PASS || echo FAIL)"
}
count() { julia --project -t 4 verify/experiments-d96.jl count "$FIX" "$OUT" cli 2>&1 | grep "^O27"; }
for r in 1 2 3; do runcli "round $r"; count; done
runcli "closure"
runcli "closure again (resume)"
echo "O27 | chunk files CLI, library | $(ls "$S" | grep -c '^chunk-'), $(ls "$OUT/store-lib" | grep -c '^chunk-') | equal | -"
if cmp -s "$L" "$OUT/round-log-lib.txt"; then r=PASS; else r=FAIL; fi
echo "O27 | CLI round log equals the library round log (bytes) | $(sha256sum "$L" | cut -c1-16), $(sha256sum "$OUT/round-log-lib.txt" | cut -c1-16) | identical | $r"
