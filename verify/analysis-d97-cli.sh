#!/bin/bash
# Stage 04 verification of module analysis, Amendment D97, checks O28 and O29 through the command line, with no scratch Analysis patch:
# `bin/chaoslab.jl run --round-log` by rounds with `bin/chaoslab.jl count-round` between the runs, the closure run, a resume run, then `bin/chaoslab.jl analyze`; exit codes, the CLI round log
# and the CLI outputs against the library ones (verify/experiments-d96.jl lib / analyze, same fixture).
# Usage, from the repo root, after `julia --project -t 4 verify/experiments-d96.jl lib <fixture dir> <out dir> [tag]` and `... analyze <fixture dir> <out dir> [tag]`:
#   bash verify/analysis-d97-cli.sh <fixture dir> <out dir> [lib tag]         (FAULTS=1 in the environment for the second campaign: runs through the same fault hook, see FAULT NOTE)
# FAULT NOTE: the fault hook is a library hook (RunHooks), the CLI has none; with FAULTS=1 this script only runs `count-round` and `analyze` on the store and log of the library tag.
set -u
FIX="$1"; OUT="$2"; LT="${3:-lib}"
IN="$FIX/inputs"; PRE="$IN/p3-abstention-preregistration.sealed.md"; PANEL="$FIX/p3-abstention-test-panel.toml"
case "$OUT" in *"/data/"*|*"\\data\\"*) echo "refused: out dir inside data/"; exit 1;; esac
CRA=("$PRE" "$IN/p3-abstention-calibration-parameters.toml" "$IN/p3-abstention-test-parameters.toml")
P2=(--p2-manifest "$FIX/p2/p2-indicators-manifest.json" --p2-preregistration "$FIX/p2/p2-indicators-preregistration.md" --p2-parameters "$FIX/p2/p2-indicators-analysis-parameters.toml" --calibration-record "$IN/p3-abstention-calibration-record.json")
chk() { echo "D97 | $1 | exit $2 | exit 0 | $([ "$2" = 0 ] && echo PASS || echo FAIL)"; }
if [ "${FAULTS:-0}" = 1 ]; then
  S="$OUT/store-$LT"; M="$OUT/manifest-$LT.json"; L="$OUT/round-log-$LT.txt"
else
  S="$OUT/store-cli2"; M="$OUT/manifest-cli2.json"; L="$OUT/round-log-cli2.txt"
  [ -e "$S" ] || [ -e "$M" ] || [ -e "$L" ] && { echo "refused: CLI outputs exist (use a fresh out dir)"; exit 1; }
  mkdir -p "$S"; : > "$L"
  runcli() { julia --project -t 4 bin/chaoslab.jl run "$PRE" "$PANEL" --store "$S" --manifest "$M" --round-log "$L" --threads 4 --allow-dirty > "$OUT/cli2-last.txt" 2>&1; chk "CLI run $1" $?; }
  count() { julia --project -t 4 bin/chaoslab.jl count-round "${CRA[@]}" --store "$S" --panel "$PANEL" --round-log "$L" "${P2[@]}" --threads 4 > "$OUT/cli2-count-$1.txt" 2>&1; local e=$?; tail -2 "$OUT/cli2-count-$1.txt"; chk "CLI count-round after round $1" $e; }
  for r in 1 2 3; do runcli "round $r"; count $r; done
  runcli "closure"; runcli "closure again (resume)"
  julia --project -t 4 bin/chaoslab.jl count-round "${CRA[@]}" --store "$S" --panel "$PANEL" --round-log "$L" "${P2[@]}" --threads 4 > "$OUT/cli2-count-closed.txt" 2>&1; chk "CLI count-round on the closed log" $?; tail -2 "$OUT/cli2-count-closed.txt"
  if cmp -s "$L" "$OUT/round-log-$LT.txt"; then r=PASS; else r=FAIL; fi
  echo "D97 | CLI round log (count-round written) equals the library round log (bytes) | $(sha256sum "$L" | cut -c1-16), $(sha256sum "$OUT/round-log-$LT.txt" | cut -c1-16) | identical | $r"
fi
an() {   # an <store> <manifest> <round log> <out folder> <label>: `bin/chaoslab.jl analyze`
  julia --project -t 4 bin/chaoslab.jl analyze "${CRA[@]}" --store "$1" --panel "$PANEL" --round-log "$3" --manifest "$2" "${P2[@]}" --p1-preregistration "$FIX/p1/p1-le-benchmark-preregistration.md"     --p1-parameters "$FIX/p1/p1-le-benchmark-parameters.toml" --out "$4" --threads 4 > "$OUT/cli2-analyze-$5.txt" 2>&1; local e=$?; tail -4 "$OUT/cli2-analyze-$5.txt"; chk "CLI analyze ($5)" $e; ls "$4"
}
cmpout() {   # cmpout <CLI out folder> <library out folder>
  for f in p3-abstention-tables.md p3-abstention-test-results.json; do
    cmp -s "$1/$f" "$2/$f" && echo "D97 | CLI $f equals the library output on the same store (bytes) | - | identical | PASS" || echo "D97 | CLI $f equals the library output on the same store (bytes) | - | identical | FAIL"
  done
}
F="${FAULTS:-0}"
if [ "$F" = 1 ]; then
  an "$S" "$M" "$L" "$OUT/analysis-cli-libstore-1" "libstore-faults"; cmpout "$OUT/analysis-cli-libstore-1" "$OUT/analysis-$LT"
else
  an "$S" "$M" "$L" "$OUT/analysis-cli2" "own-store"
  an "$OUT/store-$LT" "$OUT/manifest-$LT.json" "$OUT/round-log-$LT.txt" "$OUT/analysis-cli-libstore-0" "libstore"; cmpout "$OUT/analysis-cli-libstore-0" "$OUT/analysis-$LT"
fi
