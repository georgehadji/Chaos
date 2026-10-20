# Implementation notes: analysis

Stage 03, branch `module/analysis`, 2026-10-20. Spec: `analysis-spec.md` (decision D60). Tests: `test/analysis/` and `analysis-test-plan.md`. Read by stage 04.

## Result

- `test/analysis/`: 416 pass, 0 fail, 0 error (416). Full `Pkg.test` (default threads): 54898 pass, 0 fail, 0 error, 5 min 59 s. Quality testset 12/12 (Aqua 11, JET 1, no ignores), layering 2/2, every earlier module green (diagnostics 245, experiments 4861, lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, precision 150, standard-map 474, tangent-core 808).
- Suite times: analysis 39 to 58 s (the five synthetic studies of 55 chunks take about 1.5 s each, 9 s for the first one with compilation), budget 3 min per slug (D35). The experiments suite took 3 min 20 s in this run, over its own 3 min ceiling, which stays an open item from stage 05 of `experiments`; nothing here changed it.
- Coverage of `src/Analysis`: 661 of 668 executable lines, 99.0 percent (analysis suite alone, `--code-coverage=user`, `.cov` files counted then deleted). Uncovered: the invalid-TOML branch of `read_params`; the error constructor `_missing` for an absent required item; the threaded branch of `_pmap` (the suite runs on one thread); the counting of `breakdown_tau` items in F7; the success message and return of the command (needs CairoMakie).
- Size: longest function 25 lines (`_decide`), longest file 183 lines (`pipeline.jl`); no `@fastmath`; no accumulated logarithms (the module averages stored values only); no dashes.
- Figures: `ext/ChaoslabCairoMakieExt.jl` is not reachable from the unit suite (CairoMakie is not a test dependency). Smoke test under `julia --project=verify` on the synthetic study: `analyze(...; figures = true)` returned status `:ok` and wrote the seven files `an-fixture-fig-F1.pdf` to `F7.pdf`, each starting `%PDF-` (6.9 to 13.7 kB). The real figures are stage 04.

## Files

`src/Analysis/`: `Analysis.jl` (module, exports `analyze`, `read_params`, `read_store`), `params.jl` (E0 table reader), `stats.jl` (stream, bootstrap interval, ratio, band), `h1.jl` (H1 and controls), `sets.jl` (sets E7, R(n), n_conv, H2, H3), `h4.jl`, `labels.jl` (labels E14, H5), `store.jl` (reader), `checks.jl` (E17 and I4), `orbits.jl` (per-orbit quantities), `figdata.jl` (F1 to F7 tables), `pipeline.jl` (`analyze`, decisions), `tables.jl` (tables file), `cli.jl`. Outside the module: `src/Chaoslab.jl` (include and re-export), `Project.toml` (weak dependency CairoMakie, extension entry, compat `0.15`), `ext/ChaoslabCairoMakieExt.jl`, `bin/chaoslab.jl` (routes `analyze`, loads CairoMakie when present), `bench/analysis.jl`, `verify/Manifest.toml` (re-resolved: the entry of Chaoslab gained the extension), `test/layering.jl` (layer list, stage 02), `_config/architecture.md`, `_config/build-order.md`, `_config/decision-log.md` (D60, stage 02).

## Deviations from the spec and readings the code fixes

1. **Spec gap: the K values of S.** E7 defines S as the class C orbits at K = 2 pi and K = 10, but the parameter table E0 has no key for these two values and the approved tests (`an_sets.jl`) fix S by them. The code names them as a constant (`STRONG_K` in `sets.jl`, match to 1e-9), which is a study value in code against "the module holds no study value". It is the only one. Remedy if the owner wants it generic: a spec amendment adding a key (for example `strong_K`) and a change of `an_params.jl`, `an_helpers.jl` and `an_sets.jl` in stage 02.
2. W and E_ch take class C only at K > 0 (the preregistration says "at K <= K_c" and "K > 0"; class C does not exist at K = 0, the tests include a class C orbit at K = 0 that must stay out).
3. The half point of the rotation numbers (E14, P5) is the recorded step nearest N/2, exactly N/2 where it is recorded. The ladder items at N_B come from items of horizon N_1 (Float64, Double64) whose grid has 5.01e5 and not 5e5 as nearest recorded point, so the exact N/2 does not exist there; with the nearest step the second-half number is (x_N - x_h) / (2 pi (N - h)). For the main items (N/2 recorded) the two rules agree.
4. The drift flag check (I4) recomputes the flag of every item as the runner does: lifted |p| above the bound, a wrapped item with `2 pi` times the winding added. The test plan first said "wrapped items are not recomputed"; corrected in the plan.
5. Tables file: figure data are printed as `| col | col |` tables with numbers to 12 significant digits; the digest of a figure is the SHA-256 of that table text (LF line ends, lines joined by `\n`, no final newline).
6. `analyze` takes `--store` (default: the default store of the Experiments module for the preregistration) and `--out` (default: `workspaces/study/stages/04-analyze/output` beside the manifest's stage folder); the spec's Compute section shows only the four positional arguments.
7. The CLI run always requests the figures; without CairoMakie it returns exit code 2 after writing the tables (the test environment has no CairoMakie, so exit code 0 is a stage 04 check).

## Tests changed in this stage

The audit "Tests untouched" compares `test/analysis/` with its end of stage 02: three assertions were corrected (listed in `analysis-test-plan.md`, "Corrections found in stage 03": two convergence-time boundaries that sat exactly on a grid step, and the discrete end points in the stream-difference test). One line of an older suite changed: `test/experiments/ex_inputs.jl` AC28 ("nothing below Experiments imports it") now skips the new layer Analysis, which sits above Experiments by D60.

## Performance (bench/analysis.jl, this machine, 5 samples)

| Case | Time | Allocation |
|------|------|------------|
| bootstrap interval, 200 orbits, 1e4 resamples | 40 ms | 16 MiB |
| bootstrap band, 40 orbits, 201 steps, 1e3 resamples | 100 ms | 4 MiB |
| orbit convergence time, 8 methods, 201 steps | about 1 us | 0 |
| relative spread, 8 methods, 201 steps | about 3 us | 1 KiB |

Projection for the real study (derived, not measured): the F2 band at 1e4 resamples costs about 1 s per cell, 27 cells about 30 s; the H2, H3 and per-cell intervals under 1 s in all; reading 1363 chunk files (475 MB) is the unmeasured part and decides the 10 minute and 8 GB limits of acceptance criteria 21 (stage 04).

## Known limits and for stage 04

- **Not run on the real store, on purpose.** The real manifest and store were not given to `analyze` in this stage, so that no hypothesis outcome was seen before the independent verification (O3). Acceptance criteria 3, 19 and 21 (closed form at K = 0 on the real store, independent recomputation, wall time and memory) are stage 04. The real run also needs the study's parameter file (the transcription of E0, spec open point P1) in the study stage 04 folder; it does not exist yet. The final manifest is on branch `study/p1-le-benchmark-run` (commit e35ee5a), not on this branch.
- The reader keeps every item in memory (about 15569 series of 201 to 221 steps with their running maxima, estimated well below 1 GB) and reads the files one after the other; only the F2 bands run in parallel (the `threads` keyword).
- A required item or recorded step that the store lacks (for example N_B/10 in a ladder item) refuses the run with an error naming the orbit and the step; nothing is guessed.
- The extension was drawn only on the synthetic study; the styling of the real figures is for the stage 04 checkpoint.
- `verify/Manifest.toml` changed by two lines (the extension entry); `Manifest.toml` of Julia 1.10 is still unresolved (open item from the earlier stage).

## Audit

| Check | Result |
|-------|--------|
| GREEN | pass, 54898 of 54898 |
| Quality | pass: Aqua 11/11, JET clean, layering 2/2, coverage 99.0 percent |
| Size | pass: functions at most 25 lines, files at most 183 lines |
| Tests untouched | pass with the stated corrections: 3 assertions in `test/analysis/` (listed), 1 line in `test/experiments/ex_inputs.jl` (new layer) |
| Numerics | pass: no `@fastmath`; no accumulated logarithms in this module |

## Changes made in stage 04 (after the notes above)

The first run on the real store and the independent review found defects; they were fixed in `src/Analysis/` and recorded in `analysis-spec.md` ("Corrections found in stage 04") and in `analysis-verification.md`. In short: items with `breakdown_tau` are skipped by the finiteness check and listed in F7 (the real store has them; the fixture had none); the Y5 control without `K` no longer raises a `KeyError`; a manifest chunk with a status other than ok (failed included) refuses the run with a named cause; the ladder uses one common half point for the three precisions; the diagnostics flag is read on both L5 items; an S_D orbit without its long item or its directions, an S_L orbit without BigFloat items and a control with a positive exponent and no scope are refused; the input files are hashed again before the tables are written; the header names the module version and the commit; F2 has a pooled series over W; parameter ranges are checked; a failing figure drawer gives status `:figures_failed` and exit code 2; `FIGURE_DRAWER` replaced a generic function without methods (JET). After these changes: `src/Analysis` is 1231 lines in 14 files (largest file 201 lines, `pipeline.jl`), `test/analysis/` has 434 tests, the full `Pkg.test` gives 54916 pass in 6 min 43 s (Aqua and JET clean), coverage of `src/Analysis` is 98.3 percent (698 of 710 executable lines; uncovered: the invalid-TOML branch, the error constructor for an absent item, the threaded branch of `_pmap`, the `:figures_failed` branch and the success lines of the command).

Deviation 3 above (half point) now reads: the main and direction labels use the recorded step nearest N/2 of their own item; the ladder uses the step nearest N_B/2 that all three items recorded.

