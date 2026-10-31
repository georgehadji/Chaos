# Implementation notes: analysis, Amendment D93 (stage 03)

Date of the run: 2026-10-31, Julia 1.13.2. Spec: `analysis-spec.md`, "Amendment D93" (E46 to E58, I30 to I37, O15 to O19, criteria D93a to D93r, owner decisions of the stage 01 checkpoint). Tests: `test/analysis/an_p3_*.jl` (11 files from stage 02 and the new `an_p3_gate.jl`), plan `stages/02-tests/output/analysis-test-plan.md` (subsection "Rework at stage 03"). Code: 17 new files `src/Analysis/p3_*.jl` (1689 lines, largest 194 lines, largest function 26 lines), included from `Analysis.jl`. Changed: `cli.jl` (command `calibrate`, usage texts), `p2_dispatch.jl` (one line: `analyze` refuses the slug `p3-abstention`), `Analysis.jl` (includes, export of `calibrate`), `bin/chaoslab.jl` (routes `calibrate`), `bench/analysis.jl` (five entries), `src/Analysis/CONTEXT.md` (file table), `verify/analysis-p3.jl` (new, the moved cases, D38). The reference implementation of the tests (`an_p3_ref.jl`) was not read when the code was written and nothing was copied from it. Nothing under `data/` was read or written; the package suite never reads `data/`. Text by the implementing subagent, saved to this file by the coordinator (the subagent could not write report files).

## Result

| Suite | Pass | Fail | Error | Broken |
|-------|------|------|-------|--------|
| whole package, `Pkg.test()` (one thread, 23 min 14 s of tests, 24 min 03 s wall with precompilation) | 1401609 | 0 | 0 | 1 (weighted-birkhoff, existing) |
| quality (Aqua, JET), layering | 12, 2 | 0 | 0 | 0 |
| analysis (6 min 05 s) | 1877 | 0 | 0 | 0 |

Of the 1877 assertions of the analysis suite, 1287 are the earlier ones (D89 and before: 3 min 24 s, measured with the D93 files left out) and 590 are D93 (2 min 41 s): fixture self-test 49, rule 87, sets 44, select 58, PP 120, parameter table 96 (444 malformed files), dispatch 37, refusals 15, pipeline 81, gate 3. Every other module suite is unchanged. `verify/analysis-p3.jl` (`P3_FULL_RUN = 1`, 4 threads, `--project=verify`) runs the same files with all 54 refusal cases: 590 assertions, 0 fail, 0 error.

The first full run (before the last edits) failed one JET assertion with two possible errors in `p3_tables.jl` (`join` of an untyped table inferred as a union with `nothing`). Fixed with `::String` annotations; the second full run (the numbers above) is clean.

## Deviations from the spec

None. Where the spec leaves a form open the tests pin it (plan, "Readings pinned by the tests", readings 1 to 9) and the code follows the tests. The readings below are not in the spec text; none changes a threshold.

## Readings implemented

1. The pinned readings of the plan: a pair stream unit is written `a,b` (`pp-hb`: point index from 0, run index from 1; `pp-hc`: n_f and run; `optimism`: b from 1). An integer in 1 to m is `x mod m + 1` for the first output `x <= 2^64 - 1 - (2^64 mod m)`. The event "false-regular with probability p" is `u < p`. The nine uniforms of an orbit are drawn in the order rule, the seven comparators of I_7, decided. The type-I check uses the same streams as the main simulation at that point. K = 0 orbits are the panel orbits of K index `gate_k_index`. The family of C7 is theta* with `q_C` and `q_R` each moved by one step. Optimism is the mean over the valid resamples of the difference of the full-Cal and the resample false-regular shares, with `min_set` not applied inside a resample.
2. The p2 L1 item (Float64, N_1, `directions[1]`) of every non-control orbit of P and G is a required item (D89 coverage rule): E49 d takes lambda_ref from it for every orbit. A missing p2 L5 item of a grid orbit, a missing recorded step N_1/10, N_1/2 or N_1 of an L5 item, and a missing P1 L5 or L1 item (or its step N_1) are relabelling data missing: counted by source, K index and item, the orbit leaves Cal with reason "missing", no refusal. The test `item-missing-grid` pins the refusal for the p2 L1 item; E48 h with E49 e can be read the other way.
3. The verdict items of E51 (tg, fq1, pa1, bk1 at the primary horizons) are required for the orbits of Cal and of the gate set only (E48 h). A gate orbit outside Cal with a missing item refuses (P68; `an_p3_gate.jl`). An item listed in `not_applicable` that such an orbit lacks gives undecided verdicts at the ten budgets (the spec says nothing).
4. Only the chunks that hold an item of E49 or E51 are read from the p2 store, without the controls. The checks of D89 E36 run on these items only, so a NaN in an item that is not read changes nothing. The flags of an orbit are the union over the items read (p2 items, and for a panel orbit the P1 L5 and L1 items of C0).
5. P1 inputs: the checks of E48 f and the identity of the stored `ic` and `params` of the items read. The P1 items are not checked for finite values. The in-memory base labelling of the change counts (P46) tolerates an orbit whose items lack a recorded step (no old label; old class `missing`).
6. The p2 table is validated against the p2 preregistration by D89 E22 inside the reader of the p3 table (E47 iv). The p3 table is read with the P1 table for the label keys, as `analyze` does.
7. The input hashes of E48 j are those of the nine input files; the chunk files are not hashed again.
8. The single comparison of the type-I check of X is the conjunction of the comparisons of X alone, without the coverage condition and the other comparators; `single >= joint`.
9. PP for HB with no valid optimism resample: optimism taken as 0 (record `optimism: null` and the skipped count). PP for HC with no failure or no non-failure on Cal (plan A5): AUC undefined, PP for HC not run, `pp_hc` says "AUC not defined" with the two counts. A `:cap` target gives `n_T = y_t r_max`, `n_W = y_w r_max`; `power_cap` is the power of the last allowed point of the path.
10. The defect check of E54 (a feasible element with a capped non-LE comparator) is a refusal after the work, naming the element; F1 makes it unreachable.
11. `calibrate` without `out` writes into the folder of the preregistration (`workspaces/study/stages/01-design/output/`). `analyze` given the p3 preregistration refuses naming `p3-abstention` and `calibrate`.
12. Record: the items of E58 plus counts by K index, the first failed item of each element (`theta_space.first_fail`), shares by source, and the K > 0 values of theta* (P48). D93 accumulates no logarithm, so compensated summation has nothing to act on.

## Performance (not authoritative)

`bench/analysis.jl`, five entries (BenchmarkTools, `--project=verify`, one thread, 20 samples, minimum and median):
- Rule over 2088 elements and 3000 orbits (C3 to C6): 45.5 ms, 47.3 ms.
- One optimism resample: 45.7 ms, 47.5 ms, 21 MiB.
- One PP-HB path point (285 + 950 orbits, 100 runs): 3.36 ms, 3.57 ms.
- One PP-HC point (n_f = 65, 50 runs): 3.97 ms, 4.08 ms.

The McNemar critical-value table follows from two recurrences (30 ms for m up to 1235, equal to the direct binomial sums for m up to 400 at five levels; the first version took 200 ms per call).

Estimate for the real calibration: about 25 s serial (optimism 9 s, PP-HB with type-I check 11 s, PP-HC 12 s), 2 to 3 s at 23 threads; wall time dominated by reading the stores, estimated 1 to 2 minutes, memory below 3 GB, against the P60 bounds (30 minutes, 8 GB). Not measured on the real stores (stage 04).

## Quality

- Coverage (`--code-coverage=user`, analysis suite): 953 of 953 executable lines of the 17 files `p3_*.jl` (line coverage). The 103 `.cov` files of this run were deleted; 42 older ones left alone.
- Aqua clean; JET clean (0.12.3 in `Pkg.test()`, 0.13.1 in a temporary environment).
- Layering: `test/layering.jl` and D93q pass; no new package.
- Sizes: largest function 26 lines (`cli` 29), largest file 194 lines. No `@fastmath`.
- Exact arithmetic (I32): counts, shares, J, Theta_tol, the ceiling rule, the order of C6, the McNemar critical values and the comparisons with `rho_hb` and `c_min` are exact integers and rationals; binary64 only for the draws, the optimism mean, the quantile function and the Hanley and McNeil variance.

## Known limits

- The record layout below the top level is a proposal (plan A7).
- Sets use binary64 comparisons of stored lambda values (as D89 E29).
- A missing p2 L1 item of a grid orbit refuses (reading 2); a missing P1 item of a panel orbit does not.
- A dirty working tree is flagged "not a freeze run" in the tables file, not refused (P59).
- Byte identity across machines not claimed for PP-HC quantities (P58).
- E48 rule j is checked by hashes at the start and before writing; no test hook.
- Stage 04: O15 on the real stores, O19, D93r, time and memory at 23 threads, O16 at full size.
- The analysis suite (6 min 05 s) is over the D35 ceiling of 3 min per slug; the D89 part alone was already over (3 min 24 s).

## Rework list (tests, recorded in the plan)

1. New `test/analysis/an_p3_gate.jl` (owner decision 3, P68).
2. D38 move (owner decision 2): `an_p3_helpers.jl` gains `P3_FULL` and `P3_CASE_SAMPLE`, `an_p3_refusals.jl` one line; the package suite runs 30 of 54 refusal cases, `verify/analysis-p3.jl` all 54 (saved 17 s).
3. No assertion changed or weakened. Owner decision 4: all 54 refusal cases refuse for their own cause, named in the message (messages in the plan).

## Stage 03 audit

| Check | Pass condition | Result |
|-------|----------------|--------|
| GREEN | All tests pass, including every earlier module | pass (1401609 pass, 0 fail, 0 error, 1 broken existing) |
| Quality | Aqua and JET clean; layering; coverage at least 80% | pass (953 of 953 lines) |
| Size | Functions under 50 lines, files under 400 | pass (26, 194) |
| Tests untouched | `git diff` of `/test/<module>/` since 02 is empty | fail by the letter, recorded rework only (owner decisions 2 and 3); no assertion changed |
| Numerics | No `@fastmath`; compensated summation for accumulated logs | pass (no log accumulated) |
