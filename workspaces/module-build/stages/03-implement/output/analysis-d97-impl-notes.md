# Implementation notes: analysis, Amendment D97 (module-build stage 03)

Julia 1.13.2, 4 threads. Spec, tests, `schema/`, `_config/`, `data/`, `verify/` and every other module untouched. No dependency added; layering test passes. This file is read by stage 04.

## Result

| Run | Pass | Fail | Error | Broken | Time |
|-----|------|------|-------|--------|------|
| RED (stage 02, analysis subset: D95 and D97 files) | 515 | 48 | 19 | 0 | |
| GREEN, same subset (driver `d97/driver.jl`) | 761 | 0 | 0 | 0 | 2 min 28 s |
| GREEN, whole package suite (`Pkg.test()`, 4 threads) | 1 404 780 | 0 | 0 | 1 | 22 min 54 s |

The broken result is the existing one (weighted-birkhoff). Quality 12 of 12 (Aqua, JET), layering 2 of 2, experiments 9 264, analysis 2 638. The golden tests of p1, p2 and `calibrate` (I47) pass inside the analysis line. Subset per file: fixture 41, shape 52, stats 327, params 49, log 45, refusals 5, dispatch 22, pipeline 111, d97 109.

Logs (scratchpad `d97/`): `d97-green.log`, `d97-pkg-test.log`.

## Change (all in `src/Analysis/`, line counts of `git diff`)

| File | +/- | Change |
|------|-----|--------|
| `p3t_store.jl` | 11/9 | E74: `_p3t_chunk` reads `k<i>:<class>:r<round>` (new field `class`), refuses a ref chunk of another form naming it; `p3t_cell_of(o)`; `p3t_cells` forms the expected cells from (`k_index`, `class`, `round`). The unused `p3t_cell_ids` is removed. |
| `p3t_inputs.jl` | 11/3 | E76: `_p3t_check_status` (the three forms, naming the status); E81: `P3T_ARMS` has the six arms. The per-chunk status check against the files is unchanged and gives the set equality. |
| `p3t_data.jl` | 19/8 | E75: `_p3t_missing_items` asks tg at N_c; `p3t_verdict_tables` (own copy of the loop of `p3_verdicts`) reads SALI, FLI, MEGNO at `horizon = N_c`, also for the Double64 verdicts of D; E77: failed cells of ind and ind_pa1 are cells of E74, matched with the cell of the orbit from the panel. |
| `p3t_run.jl` | 12/7 | E78: `p3t_measured_costs` reads tg, fq1, bk1 from `arms.ind.specs`, pa1 from `arms.ind_pa1.specs`; `_p3t_cost_l1` is 1e9 sum(thread_seconds) / sum(steps) over the ok L1 chunk entries of ref. |
| `p3t_tables.jl` | 25/20 | E79: `_p3t_chk_items` finds the chk items by content and refuses a horizon other than N_B; T3 without `labels agree`, with the omission text and the count sentence. E80: T2 rows and counts over the rounds up to the closing round; T11 orbits of a cell by `p3t_cell_of`. |
| `p3t_json.jl` | 4/4 | E81, P114, P115, P120 texts; the tg deviation is removed; the machine file lists `vcat(P3T_READINGS, P3T_DEVIATIONS)` under `readings`. |

Largest changed function: `_p3t_t3` (24 lines). All files stay under 400 lines (largest `p3t_stats.jl`, 216). The files keep their CRLF working-tree endings. No em dash, no `@fastmath`. `p2_*` and `p3_*` files are untouched, so p1, p2 and `calibrate` are unchanged (I47).

## Deviations from the spec (readings)

1. P120 in the machine file. The machine file has a fixed key order (test D95n) and no `deviations` key; the deviation text is appended to `readings` there (the test needs it in both files). The tables file keeps the two lists apart.
2. E75 is done inside `p3t_*` (a horizon keyword for SALI, FLI, MEGNO), not by changing `_p2_find`, which p2 uses with N_1.
3. E77 d64 and d64_pa1: a failed chunk leaves the orbits without Double64 items, so they leave D by the existing rule (P104); no extra count was added beyond the size of D in T2.
4. E77 chk: D_B is the set of orbits that carry a chk Double64 item; a failed chk chunk leaves a missing cell (`n.a.`) in T3 and refuses nothing (the test accepts both readings, ambiguity 1 of the test plan).
5. T2 sentence: "Orbits of the rounds up to the closing round: N (the panel file holds M); rounds R" (not pinned by the tests).

## Known limits

- Statistics with T_s > 0 on a real store, early close, class S, drift flags, a paused arm, kill and resume: not reached (spec "Out of scope").
- `p3t_verdict_tables` duplicates the six lines of `p3_verdicts` so that D93 stays byte-identical.

## Not done here

Coverage of the changed functions was not measured (the D97 tests exercise the changed branches by design; this was not checked by a coverage tool). No benchmark change. O28 and O29 (real fixture campaigns, independent reader) are stage 04. Nothing committed.
