# Implementation notes: experiments, Amendment D96 (module-build stage 03)

Julia 1.13.2, 4 threads. Spec, `schema/`, `_config/`, `CHANGELOG.md`, `src/Experiments/CONTEXT.md`, the decision log, `data/` and `verify/` untouched. No dependency added. This file is read by stage 04.

## Result

| Run | Pass | Fail | Error | Broken | Time |
|-----|------|------|-------|--------|------|
| RED (stage 02, whole package suite) | 1 404 601 | 0 | 1 (D96a) | 1 | 28 min 13 s |
| GREEN, D96 script (helpers, `ex_subsets.jl`, the five D94 files, `ex_d96.jl`) | 2 477 | 0 | 0 | 0 | 65 s |
| of which `ex_d96.jl` | 20 | 0 | 0 | 0 | 7.7 s |
| GREEN, whole package suite (`Pkg.test()`, default environment, 4 threads) | 1 404 618 | 0 | 0 | 1 | 25 min 16 s |

The broken result is the existing one (weighted-birkhoff). Quality 12 of 12 (Aqua and JET clean), layering 2 of 2, experiments 9 264 pass, analysis 2 476 pass. D86b passes in the default environment (see the stage 02 note on `--project=verify`).

Logs (raw, scratchpad): `d96-green.log` (first GREEN run, 2 fails of the test below), `d96-green2.log` (after the test correction), `d96-pkg-test.log` (package suite).

## Change

`src/Experiments/run_indicator_report.jl`, `_indicator_manifest`: the first item of a spec is found with `findfirst`; a spec with no item in the plan is skipped (R28, P116: omitted, not listed empty). 3 lines added, 1 removed; the file is 76 lines, the function 18. No other file of `src/` changed. A Sonnet subagent was not used: the change is the guard of the discovery run (scratchpad `d96/runner_patch.jl`) in three lines.

## Test corrected (back to stage 02, owner decision 2026-11-03)

D96a compared the manifest budget of each spec with the raw fixture budget (Float64 and Double64). The block keeps only the precisions of the arms holding the spec (D86 item 1, `run_indicator.jl`, behaviour before D96), and the fixture's arms use Float64 only. The expectation was wrong, not the code: `ex_d96.jl:78` now compares with the Float64 line alone. Recorded in the test plan, section "Amendment D96".

## Not done here

Coverage was not measured for the three lines: D96a reaches the skip, and the closure run reaches the listing. No benchmark change: the table is built once per plan. O27 (fixture campaign through the library and the CLI without patch) is stage 04.
