# Implementation notes: analysis, Amendment D89 (stage 03)

Date of the run: 2026-10-29, Julia 1.13.2. Spec: `workspaces/module-build/stages/01-spec/output/analysis-spec.md`, "Amendment D89". Tests: `test/analysis/an_p2_*.jl`, plan `stages/02-tests/output/analysis-test-plan.md`. Code: `src/Analysis/p2_*.jl` (17 files, 2107 lines, none above 195 lines, no function of 50 lines or more), `pipeline.jl`, `cli.jl`, `Analysis.jl` (changed), `bench/analysis.jl` (four entries added), `src/Analysis/CONTEXT.md` (file table). No test, no `data/`, no frozen document, no other layer was touched.

## Result

Folder suite, `julia --project -t 4 --check-bounds=yes` (after the owner rework of 2026-10-29, tests reworked at stage 02): 1134 tests, all pass (including the 434 P1 tests). Aqua and JET clean in a temporary environment, `test/layering.jl` and D89t pass. Coverage measured before the rework: P2 files 99 percent, module 98 percent. The full `Pkg.test()` was not run.

## A1 (closed)

The stage 02 rework changed the base fixture slug; E21 is strict (`p2_dispatch.jl`: any slug but the two known ones is refused naming the slug and both known ones) and nothing special-cases a slug.

## Readings that differ from the spec text or that the tests fix

1. R2 version (E33), owner decision of the stage 03 checkpoint: the R2 version recomputes the detection times (M_X, Q, R_X), the verdict series and the H6 verdicts (directions, Float64 and Double64) from the R2 verdicts, for the primary sets (`r2`) and the extension (`r2_extension`). The R1 version stays the primary. The switch `P2_R2_REDETECT` was removed.
2. A2, owner decision of the stage 03 checkpoint: every key other than the 11 machine keys needs `quote` and `text`; `text` (the value as the section writes it: `1e-4`, `10^4`, `95 percent`, `5th`) must occur in `quote`, and `quote` in the named section (`p2_params.jl`, `_p2_entry` and `_p2_quotes`). A missing, non-text or not-in-quote `text` refuses naming the key. The string form of `value` is not matched to anything. The four label keys carry a quote and text of the section "Reference class rho(o)" that cites the P1 rules (they hold no number; the value is checked against the P1 table).
3. `p2_leaving_budget` is all or nothing in the pipeline: if the regular range is missing at any budget of an indicator and K, the leaving budget is `nothing` (as `an_p2_pipeline.jl` computes it). The spec reads "for every later index".
4. The extension H4 on the design ensemble uses the whole T_s^design (not restricted to E8) with the label set name `T_s_design_8` and set rule `E8` (the test reference does the same).
5. `p2_grid_starts(prereg)` has no panel argument in the test call. First id of the grid: keyword `first_id`, else one above the largest id of the panel file given by keyword `panel` (the pipeline passes it), else the single `*-panel.toml` next to the preregistration (the fixture), else 884.
6. Items of series kinds (tg, pr, st) must be finite in every row and column (table E31 of D86: every column required finite); window kinds finite unless the second flag bit is set, `reldig` also unless `reldig_undefined`; an infinite value anywhere refuses. The check covers rows that no budget uses (the test `nan-series` needs this). It passed on the real store.
7. The recorded steps of an item must equal those of the first orbit of its chunk and match `steps_recorded` (a missing row is detected there; a missing needed step cannot be detected otherwise because `n(b)` is the first step not before `b`).
8. A control with a positive exponent and no `control_scope` is not refused in P2 (the spec E36 does not list it); the stored `control_deviation` is checked to 1e-12.
9. `--out` inside the repository `data/` refuses for study p2-indicators only (P1 refusals stay as before, I11).
10. The old P1 command had a defect: `--threads N` always raised (`something(..., throw(...))` evaluates the throw). Fixed in `cli.jl` (`_threads`); no P1 test or output depends on it.
11. `analyze` is now the dispatcher (`p2_dispatch.jl`); the old body is `analyze_p1` (`pipeline.jl`) with the shared preparation `_p1_prepare`, also used for the in-memory P1 labelling. A P2 table given with a P1 preregistration is refused with its unknown and missing keys named (only when it has keys of the P2 set, so `read_params` itself behaves as before).
12. Tables: an undefined number is the text `nothing` in the markdown and in the digest (the test digest rule `string(x)`).
13. Spearman (E30) is reported in F1 (extra column `spearman`); F2 adds the rows `t_det`; T1 carries two text lines (omissions, drops).

## Performance (not authoritative)

Dry run on the real stores (`data/p2-indicators`, `data/p1-le-benchmark`, read only, output in the session scratch area) with a drafted real parameter table (58 keys with quote and text, accepted by the reader with the quote check on): 81 s wall at 23 threads (82 s at 8 threads earlier), peak resident memory 1.8 GiB, the tables file 1.0 MB; byte-identical at 8 and 23 threads (earlier table). Most of the time is serial (reading 695 chunk files, the P1 labelling). Set sizes equal the drafting prototype of the spec: E_all 2974, E8 2236, W2 214, S2 448, T_s 30, T_s^design 96, regular 1159, gate 80, S_D* 109, S_L* 51. No finite-value, identity, flag-count or coverage refusal fired. The drafted table is not a stage 04 artifact (not compared line by line to the preregistration) and no result of the run is recorded here.
`bench/analysis.jl`: `p2_bootstrap` H1 on 214 orbits at 1e4 resamples 0.61 s; H2 on 448 orbits 1.5 s; `p2_unanimity` 28 pairs 214 orbits 87 microseconds; `p2_auc` 300 against 150 scores 15 microseconds. About 35 bootstrap labels per run, parallel over labels.

## Known limits and not checked

- E36 o (input changed while the run goes): the hashes are taken at the start and compared before the tables are written (P2 inputs and P1 reference inputs); a test hook was not added.
- Not done here (stage 04): the independent Python layer (O6), the real P2 parameter table, D89v to D89x, the PDF figures (P11, exit code 2).
- The bootstrap runs on records through `p2_stat` (readable, same code as the point statistics); a faster index-based form was not needed.
- The coverage was measured on a run where the old P1 tests fail on the slug, so the P1 files are slightly under-covered.
