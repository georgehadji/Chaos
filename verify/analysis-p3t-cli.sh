#!/bin/bash
# D95q, D95o, D95d (V): the module end to end through the command line at 1, 4 and 23 threads; time, peak memory, determinism, read-only.
# Run from the repo root after verify/analysis-p3t.jl:  bash verify/analysis-p3t-cli.sh
REPO="C:/Users/tesse/Documents/Vibe Coding/Chaos"
D95="C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/d95v"
JL='C:\Users\tesse\.julia\juliaup\julia-1.13.2+0.x64.w64.mingw32\bin\julia.exe'
RUN() { python -I "$REPO/verify/analysis-p3t-run.py" "$@"; }
cd "$REPO"
listing() { (cd "$1" && find inputs p1 p2 store -type f -print0 | sort -z | xargs -0 sha256sum); }
for T in 1 4 23; do
  W="$D95/run-t$T"; rm -rf "$W"; mkdir -p "$W"
  cp -r "$D95/campaign/inputs" "$D95/campaign/p1" "$D95/campaign/p2" "$D95/campaign/store" "$W/"
  listing "$W" > "$W/in-before.txt"
  for r in 1 2 3; do
    (cd "$D95/stage/r$r" && find . -type f -print0 | sort -z | xargs -0 sha256sum) > "$W/stage-r$r-before.txt"
    ( cd "$W" && RUN "count-round r$r t=$T" "$W/cr$r.out" -- "$JL" -t $T --project="$REPO" "$REPO/bin/chaoslab.jl" count-round inputs/p3-abstention-preregistration.md \
      inputs/p3-abstention-calibration-parameters.toml inputs/p3-abstention-test-parameters.toml --store "$D95/stage/r$r" --panel inputs/p3-abstention-test-panel.toml \
      --round-log round-log.txt --p2-manifest p2/p2-manifest.json --p2-preregistration p2/p2-preregistration.md --p2-parameters p2/p2-parameters.toml \
      --calibration-record inputs/p3-abstention-calibration-record.json --threads $T )
    (cd "$D95/stage/r$r" && find . -type f -print0 | sort -z | xargs -0 sha256sum) | cmp - "$W/stage-r$r-before.txt" && echo "stage r$r unchanged"
    listing "$W" | cmp - "$W/in-before.txt" && echo "inputs unchanged after count-round $r; round-log bytes: $(wc -c < "$W/round-log.txt")"
  done
  for O in module-out $( [ $T = 4 ] && echo module-out-b ); do
    ( cd "$W" && RUN "analyze t=$T out=$O" "$W/an-$O.out" -- "$JL" -t $T --project="$REPO" "$REPO/bin/chaoslab.jl" analyze inputs/p3-abstention-preregistration.md \
      inputs/p3-abstention-calibration-parameters.toml inputs/p3-abstention-test-parameters.toml --store store --panel inputs/p3-abstention-test-panel.toml \
      --round-log round-log.txt --manifest inputs/p3-abstention-test-manifest.json --p2-manifest p2/p2-manifest.json --p2-preregistration p2/p2-preregistration.md \
      --p2-parameters p2/p2-parameters.toml --calibration-record inputs/p3-abstention-calibration-record.json --p1-preregistration p1/p1-preregistration.md \
      --p1-parameters p1/p1-parameters.toml --out $O --threads $T )
    listing "$W" | cmp - "$W/in-before.txt" && echo "inputs unchanged after analyze ($O)"
    echo "files in $O: $(ls "$W/$O" | tr '\n' ' ')"
  done
done
echo "--- identity of the outputs across runs"
for T in 4 23; do
  cmp "$D95/run-t1/round-log.txt" "$D95/run-t$T/round-log.txt" && echo "round-log t1 == t$T"
  for f in p3-abstention-tables.md p3-abstention-test-results.json; do cmp "$D95/run-t1/module-out/$f" "$D95/run-t$T/module-out/$f" && echo "$f t1 == t$T"; done
done
for f in p3-abstention-tables.md p3-abstention-test-results.json; do cmp "$D95/run-t4/module-out/$f" "$D95/run-t4/module-out-b/$f" && echo "$f out == out-b (t4)"; done
cmp "$D95/campaign/round-log.txt" "$D95/run-t1/round-log.txt" && echo "round-log == part-1 library log"
for f in p3-abstention-tables.md p3-abstention-test-results.json; do cmp "$D95/campaign/module-out/$f" "$D95/run-t1/module-out/$f" && echo "$f == part-1 library output"; done
echo "--- absolute paths in the outputs:"; grep -l -i "Users\|AppData\|Temp" "$D95/run-t1/module-out/"* || echo none
