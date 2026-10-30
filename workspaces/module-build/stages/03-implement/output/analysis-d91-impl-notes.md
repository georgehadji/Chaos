# Implementation notes: analysis, Amendment D91 (stage 03, module side)

Spec: `analysis-spec.md`, "Amendment D91" (E45, I22 to I26, I29, D91b to D91e, D91n, D91o). Tests: `test/analysis/an_p2_figures.jl` (not edited). The drawers in `ext/` are not part of this note.

## Result

Folder suite (`julia --project -t 4 --check-bounds=yes`, `test/analysis/runtests.jl`): 1287 pass, 0 fail, 0 error (1214 before; the new file adds 73). The first run had 2 failures (a drawer-package name inside a string of `cli.jl` tripped the layering tests D91n and AC20); fixed, second run green. Not run: the 17-minute full suite, Aqua, JET, coverage, `bench/`.

## What changed

- `src/Analysis/pipeline.jl`: new `_draw_figures(res, out, wanted)` returns `(status, error)`; `error` is the exception message (`sprint(showerror, e)`) when the drawer raised, else empty. `_figures` is now `_draw_figures(...).status`, so P1 (`analyze_p1`, its result fields and its warning) is unchanged.
- `src/Analysis/p2_pipeline.jl`: `analyze_p2` builds the result with `status = :ok`, writes the tables file (as before, after both input re-checks), then adds `figures` (a `Dict{String, Any}` "F1".."F6" holding the items `res.items.Fk`) and calls `_draw_figures`; the result gets `status` (`:ok`, `:figures_skipped`, `:figures_failed`) and `figures_error`. The tables file is written from the result before `figures` is added, so its bytes are unchanged (I23). The result has all it had plus `tables`, `figures`, `figures_error`. Docstring of `analyze_p2` updated.
- `src/Analysis/cli.jl`: for study p2-indicators a failed draw prints "figures failed: <message>" (exit 2); the P1 line is byte-unchanged. The skipped line is now the one P1 line for both studies ("figures skipped: the plotting package is not loaded (run under the verify environment)", exit 2); "not built yet" is gone.

## Deviations from the spec

- The drawer receives the result that already holds `tables` (the path) and `figures`; for P1 the drawer still receives the result without `tables` (unchanged).
- The skipped line does not name CairoMakie: the tests D91n and AC20 forbid plotting-package names in the code of `src/Analysis` (strings included). It says "the plotting package" and "verify", which D91d accepts.

## Open point for the owner

The reading P11 in `P2_READINGS` (`p2_pipeline.jl`) still says the drawers are a follow-up and that the command exits with code 2 (figures skipped). It is now outdated, but the text is part of the tables file, which stage 04 accepted (`workspaces/study/stages/04-analyze/output/p2-indicators-tables.md`, Readings section); changing it changes those bytes, so it is left unchanged. No test pins the text; the decision (change it and regenerate the tables of stage 04, or keep it) is the owner's.

## Limits

Real pictures (D91a, D91f to D91m) are checked by `verify/analysis-p2-figures.jl`, not here. A multi-line exception message would print over several lines in the "figures failed" line.

## Full suite and one test fix (main session)

The full `Pkg.test()` first stopped on `an_p2_figures.jl` line 8: `import Logging` fails there because the Logging stdlib is not a test dependency (the folder suite ran under `--project` and did not see it). The line was changed to `const Logging = Base.CoreLogging` (`with_logger` and `NullLogger` live there); no test body was edited, so no test result changed. This is the one touch of `test/` after stage 02, a defect of the test file found at stage 03. Second full run: passed (the analysis folder 1287 tests, the rest as before).
