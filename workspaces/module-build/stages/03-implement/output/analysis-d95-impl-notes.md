# Implementation notes: analysis, Amendment D95 (stage 03)

Date of the run: 2026-11-02 to 2026-11-03, Julia 1.13.2. Spec: `analysis-spec.md`, "Amendment D95" (E59 to E73, I38 to I46, O20 to O24, criteria D95a to D95q, owner decisions P87 to P115). Tests: `test/analysis/an_p3t_*.jl` (9 files of stage 02), plan `stages/02-tests/output/analysis-test-plan.md`, section "Amendment D95 (stage 02)". Code: 12 new files `src/Analysis/p3t_*.jl` (1442 lines, largest `p3t_stats.jl` 216 lines, every function under 45 lines), included from `Analysis.jl`. Changed: `Analysis.jl` (includes, export of `count_round`), `cli.jl` (command `count-round`, the p3 form of `analyze`), `p2_dispatch.jl` (the D93 refusal of `analyze` for `p3-abstention` lifted, E59), `bin/chaoslab.jl` (routes `count-round`), `bench/analysis.jl` (4 entries), `src/Analysis/CONTEXT.md` (file table, 78 lines). The reference of the tests (`an_p3t_ref.jl`) was not read when the code was written. Nothing under `data/` was read or written. Text by the implementing subagent (Sonnet), saved to this file by the coordinator.

## Result

| Suite | Pass | Fail | Error | Broken |
|-------|------|------|-------|--------|
| whole package, `Pkg.test()` (31 min 05 s), before the rework of `an_p3t_log.jl:39` | 1404597 | 1 | 0 | 1 (weighted-birkhoff, existing) |
| quality (Aqua, JET), layering | 12, 2 | 0 | 0 | 0 |
| analysis (8 min 41 s) | 2475 | 1 | 0 | 0 |

D95: 598 assertions pass (558 in the 33 theme testsets, 94.8 s; 40 in the 4 fixture testsets); the reworked D93a testset 5 pass. The analysis suite grows by 2 min 46 s (the plan estimated 8 to 12 min); no case is moved to `verify/`. The experiments suite took 12 min 54 s against 9 min 13 s at RED with no change to its code (machine load).

The one failure was the fixture constant of `an_p3t_log.jl:39` (`length(p3t_ref_hashes(S, 3)) == 22`): the toy campaign has 12 ref cells with L1 and L5 each, 24 ref chunks (plan: "24 ref chunks", 24 + 20 + 4 + 2 = 50 files). The assertion calls no module code. Changed to 24 as a recorded test rework (owner decision 2026-11-03); after the rework the analysis suite (`--project=verify -t 4`) gives 2476 pass, 0 fail, 0 error in 7 min 48 s. The whole package was not rerun after this one-line change; the other suites are untouched since the run above.

## Recorded test reworks

1. `test/analysis/an_p3_dispatch.jl`, testset "D93a ...": renamed "D93a analyze with the preregistration of this study and no round log is refused (exit code 1) naming the round log (D95 E59 lifts the refusal naming calibrate)". The API and command line calls now pass the p3 inputs with an absent round log, an absent test-side table and an absent record, and expect a refusal naming `round log`, exit code 1, no output folder and no log file. The "options of the P2 pipeline" assertion and the unknown-slug assertion for `p3-other` are unchanged. (Owner decision at the D95 stage 02 checkpoint.)
2. `test/analysis/an_p3t_log.jl:39`: 22 to 24 (above).

## Deviations from the spec

None known. The spec does not fix the order of the rules of E61; reading 1 fixes it. The tg item is read at N_1 (plan reading 4, E65 says N_c); the tables file lists this under "Deviations" until the study amendment P69 fixes the run block.

## Readings implemented

1. Order of the checks of `analyze`: `--out` and `--round-log` against `data/`, slug, hash rule, round log, 71-key table, test-side table, p2 inputs, panel, store scan, log verification (E63), manifest, closure completeness, identity, items. The round log comes before the tables so that an absent log is named first (D93a rework); the D95c cases pass with this order.
2. `count_round` treats a log whose last entry says `close` as closed and only verifies it. For an open log: no closure chunk, entries 1 to k verified, r = k + 1, no ref chunk above r, then round r complete. No ref chunk of round r at all gives `:nothing` (with an empty log: `counts` `nothing`, round 0); some but not all complete refuses naming the cells.
3. A log is split into entries at lines that start with `round=`; text before the first such line is entry 1. Each entry is parsed and compared byte for byte with its canonical rendering; the message names the entry and the first differing line.
4. Ref cells of round r are the (round, K index) pairs of the panel, each with an L1 and an L5 chunk. A failed chunk counts for the digest and excludes the panel orbits of its cell.
5. Lambda values become exact rationals through their shortest decimal text (`p3_rat`), as the thresholds.
6. HA and HC bounds are computed also below the floors (so `bounds` and `dropped` are reported); the decision is then "not decidable".
7. HB `margin` is `FR_X != 0 && FR_rule <= rho_hb * FR_X` (the D93 form); `failing` lists T_s first, then W2; `defects` lists capped non-LE comparators.
8. HC bootstrap on integer rank counts: Delta_AUC is an exact rational with denominator 2 n_f n_n. The set is the D orbits sorted by id.
9. Power at achieved sizes: model from `rule` and `comp` of the record; T_s comparators are the voters, W2 comparators I_7 (the `ts` and `w2` lists of the record model are not read). HC power uses the actual non-failure count; skipped if a class is empty.
10. Sensitivity (iii): measured cost from `arms.ref.specs` or `arms.ind.specs`; a missing value gives "not decidable" naming the spec. Sensitivity (i) uses the round-down budgets for LE too.
11. An orbit missing a verdict item and not in a chunk that failed twice refuses (E61 k). A failed ind chunk is found by its cell `k<K index>`; its orbits drop with "no verdict data". D_B is the set of orbits with the chk Double64 L5 item, whose `in_D` must be set. D is E_test with `in_D` and the Double64 items of V*.
12. T3: L5 value at the last recorded step not above N_B; N_B and the ladder bound from the P1 table; ladder test `|f64 - d64|` or `|d64 - bf|` above the bound; label agreement by `p3_relabel` on each precision's L5 with the ref L1.
13. Columns of T2 and T11 are the implementer's choice (the spec lists items, not columns).

## Performance (not authoritative)

`bench/analysis.jl`, `--project=verify -t 4`, median of 3: `p3t_ha` 7430 orbits, 100 replicates 6.5 ms (about 0.65 s for 1e4); `p3t_hc` 1309 orbits, 100 replicates 1.5 ms; `p3t_percentile` 1e4 exact values 0.52 ms; `p3t_counts` 7430 orbits 7.7 ms (file reads dominate a real `count_round`, not measured). A toy `analyze` run takes 1 to 3 s after compilation.

## Quality

- Coverage (`--code-coverage=user`, the D95 test files): 907 of 908 lines of the 12 new files (`p3t_tables.jl` 132 of 133, the rest 100%). The 70 `.cov` files of this run were deleted.
- Aqua and JET clean; layering and D95p pass; no new package; no `@fastmath`.
- Exact arithmetic (I42): counts, shares, p-values, bootstrap statistics, percentiles and threshold comparisons are rationals.

## Known limits

- No real test store exists: arm names, item horizons and manifest fields rest on the pinned readings (P115).
- Not run: the (V) parts (O20 independent Python, O23 alterations of the real D94 log, O24, 23 threads, time and memory at study size).
- Thread independence is by construction (one stream per replicate or run, integer sums); tested at 1, 2 and 4 threads.
- A failed write of the log is tested with a read-only file only.
