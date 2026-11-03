#!/bin/bash
# D95a (V): `calibrate` and the p1 and p2 `analyze` on the real inputs are unchanged by Amendment D95. Golden outputs: the module at commit 263717c (the last commit before D95, a git worktree
# under the scratch root); new outputs: the module at HEAD. Both read data/ (never write it) and write only under <root>/regress/. Figures are skipped by both (root project has no CairoMakie).
# Run from the repo root:  bash verify/analysis-p3t-regress.sh
REPO="C:/Users/tesse/Documents/Vibe Coding/Chaos"
D95="C:/Users/tesse/AppData/Local/Temp/claude/C--Users-tesse-Documents-Vibe-Coding-Chaos/5df15ee1-47df-4875-85bd-1bd4d7c665a6/scratchpad/d95v"
JL='C:\Users\tesse\.julia\juliaup\julia-1.13.2+0.x64.w64.mingw32\bin\julia.exe'
GOLD="$D95/wt-d94"
ST="$REPO/workspaces/study/stages"
cd "$REPO"
if [ ! -d "$GOLD" ]; then git worktree add --detach "$GOLD" 263717c && cp "$REPO/Manifest.toml" "$GOLD/Manifest.toml"; fi
one() {  # $1 project dir, $2 label
  P="$1"; L="$2"; O="$D95/regress/$L"; rm -rf "$O"; mkdir -p "$O"
  for cmd in calibrate p2 p1; do
    case $cmd in
      calibrate) A=(calibrate "$ST/01-design/output/p3-abstention-preregistration.md" "$ST/01-design/references/p3-abstention-calibration-parameters.toml" --out "$O/calibrate"
                    --p2-store "$REPO/data/p2-indicators" --p2-manifest "$ST/03-run/output/p2-indicators-manifest.json" --p2-preregistration "$ST/01-design/output/p2-indicators-preregistration.md"
                    --p2-parameters "$ST/04-analyze/references/p2-indicators-analysis-parameters.toml" --panel "$ST/02-panel/output/p1-le-benchmark-panel.toml"
                    --p1-store "$REPO/data/p1-le-benchmark" --p1-manifest "$ST/03-run/output/p1-le-benchmark-manifest.json" --p1-preregistration "$ST/01-design/output/p1-le-benchmark-preregistration.md"
                    --p1-parameters "$ST/04-analyze/references/p1-le-benchmark-analysis-parameters.toml") ;;
      p2) A=(analyze "$ST/01-design/output/p2-indicators-preregistration.md" "$ST/03-run/output/p2-indicators-manifest.json" "$ST/02-panel/output/p1-le-benchmark-panel.toml"
                 "$ST/04-analyze/references/p2-indicators-analysis-parameters.toml" --store "$REPO/data/p2-indicators" --out "$O/p2"
                 --p1-store "$REPO/data/p1-le-benchmark" --p1-preregistration "$ST/01-design/output/p1-le-benchmark-preregistration.md" --p1-manifest "$ST/03-run/output/p1-le-benchmark-manifest.json"
                 --p1-parameters "$ST/04-analyze/references/p1-le-benchmark-analysis-parameters.toml") ;;
      p1) A=(analyze "$ST/01-design/output/p1-le-benchmark-preregistration.md" "$ST/03-run/output/p1-le-benchmark-manifest.json" "$ST/02-panel/output/p1-le-benchmark-panel.toml"
                 "$ST/04-analyze/references/p1-le-benchmark-analysis-parameters.toml" --store "$REPO/data/p1-le-benchmark" --out "$O/p1") ;;
    esac
    python -I "$REPO/verify/analysis-p3t-run.py" "$L $cmd" "$O/$cmd.out" -- "$JL" -t 23 --project="$P" "$P/bin/chaoslab.jl" "${A[@]}"
    tail -c 300 "$O/$cmd.out" | tr '\n' ' '; echo
  done
}
one "$GOLD" golden
one "$REPO" head
python -I "$REPO/verify/analysis-p3t-regress.py" "$D95/regress"
