# Implementation notes: experiments, Amendment D94 (module-build stage 03)

Julia 1.13.2, 4 threads. Nothing committed. Tests (`test/`), spec, `schema/`, `_config/`, `CHANGELOG.md`, `src/Experiments/CONTEXT.md`, the decision log, `data/` and `verify/` untouched. No dependency added to `Project.toml`. This file is read by stage 04.

## Result

| Run | Pass | Fail | Error | Broken | Time |
|-----|------|------|-------|--------|------|
| RED (stage 02, D94 files) | 238 | 52 | 28 | 0 | 41 s |
| GREEN, D94 script (`run_d94.jl`: helpers, `ex_subsets.jl`, the five D94 files) | 2 457 | 0 | 0 | 0 | 51 s (68 s wall with start-up, 4 threads) |
| of which the D94 files alone | 2 390 | 0 | 0 | 0 | about 37 s of testset time inside the package suite (sum of the 28 D94 testsets), budget 120 s |
| GREEN, folder `test/experiments/` inside `Pkg.test()` | 9 243 | 1 | 0 | 0 | 9 min 30 s |
| GREEN, whole package suite (`Pkg.test(coverage = true)`, 4 threads) | 1 403 998 | 1 | 0 | 1 | 30 min 00 s |

The one failure is a conflict between two tests, not a defect of the code (section "Tests believed wrong"). The one broken result is the one of the earlier run (weighted-birkhoff), not in experiments. Quality testset 12 of 12 (Aqua and JET clean), layering 2 of 2, analysis 1 877 pass, diagnostics 245 pass. Every earlier experiments test passes except the one named below (it fails only because of the option the amendment must add).

Per file (D94 script, pass): grammar 187 (D94a, D94b, D94k, D94u and the guard), streams 127, panel 1 866, plan 155, run 55. A first JET run found three possible errors in the new code (regex captures typed `Union{Nothing, ...}` in `_ratio` and `_parse_entry`); fixed with `string(...)` and type assertions; JET then reports 0 (checked in a scratch environment with the package developed, and again in the final package run).

Logs (raw): scratchpad `d94-run1.log` to `d94-run4.log` (the D94 script while fixing), `pkg-test-run1.log` (first package run, one JET and one D86x failure), `pkg-test.log` (final package run), `d94-cov-run.log` (D94 script with coverage).

## Files

New in `src/Experiments/` (5 files, all under 400 lines, all functions under 50):

| File | Lines | Holds |
|------|-------|-------|
| `inputs_rounds.jl` | 120 | the new panel input keys (R22): source, points, `r_max`, `gate_rounds`, `quota`, `d_fraction`, `db_size`, `newton_grid_side`, `controls`; exact rationals |
| `rounds_streams.jl` | 106 | the four stream seeds, `random_screen_points`, `d_size`, `d_indices`, `db_allocation`, `db_subset`, `draw_pool_rounds` (E32 to E38) |
| `rounds_panel.jl` | 149 | the builder with the random source and the draw in rounds (R23, R24): round cells, orbit records, draw checks, pool sidecar, manifest keys |
| `round_log.jl` | 97 | read, parse and check of the round log (E40, R26 item 1): `read_round_log`, `round_state`, `check_round_log` |
| `run_rounds.jl` | 123 | the plan with rounds (E39, E41, R25): passes per phase and round, D and D_B at closure, the ref arm digest, plan text lines, manifest record, the request check (`arms`, `up_to`) |

Changed (lines added, removed in the diff): `Experiments.jl` (includes, one export line), `cells.jl` (`_draw_R` and `_draw_X` take `rng` and `n` keywords, defaults are the old behaviour), `inputs.jl` (new optional keys, `n_c` optional with `quota`, one range check), `panel.jl` (`_manifest` takes the cell rows, `_assemble` branches on the source, the pool sidecar is written after the panel), `panel_tables.jl` (`controls = "none"`), `run_block.jl` (`release`, subsets `D` and `D_B`), `run_plan.jl` (345 lines; `rounds` field, `round_log` and `check` keywords, passes), `run_campaign.jl` (`RunHooks.read`), `run_session.jl` (`round_log`, `arms`, `up_to`), `run_resume.jl`, `run_store.jl` (`onread`), `run_groups.jl`, `run_manifest.jl`, `run_cli.jl` (`--round-log`). `bench/experiments.jl` gains five entries (below). No `.md` file touched.

## Design in one page

- Inputs: the base reader first, then `_check_rounds`. `quota` replaces `n_c` (both given is refused naming both); without `quota` the table is derived from `n_c` (T and R at K = 0, the five classes elsewhere). Fractions are `Rational{BigInt}` from integers, decimal text, floats (their shortest decimal text) or `a/b`; `d_size`, the allocation and every count use integer arithmetic (I36).
- Builder: `screen_source = "lattice"` takes the old path untouched (golden bytes equal). The random source runs: Newton search on `newton_grid_side` (default `grid_side`, so the lattice path is unchanged), points of every K from the stream, `_screen_flags` on the points (the base function, index-write parallel distribution), `draw_pool_rounds` per K, R and X draws from the round streams with the base `_draw_R`, `_draw_X`, numbering by the base `_number_orbits!` (the order of E36 is the order of the cell list), draw checks (R24, I29), manifest keys, pool sidecar. Screening values are dropped at once (I34).
- Plan: a block with no `release` key is planned exactly as before (the log is ignored). A block with release keys on a panel of at least two rounds is planned in passes: for each phase (each_round arms, then closure arms when closed) and each round 1 to R_plan, one pass of the base item builder over that round's orbits. Items, chunks, ordinals, keys and hashes of the each_round phase are therefore a function of the round alone (I33). Cell names carry `:r<round>`.
- Log: parsed before the plan (R_plan from the entries), checked after the chunks exist (the digest is computed from the plan's chunk hashes, reading A of the pinned reading 6). The log is only read.
- Resume: files beyond the plan (ordinal above the chunk count) are stale and are not opened; the backend table check ignores run types the earlier session did not have, for plans with rounds only.

## Readings implemented (the 13 pinned readings of the test plan)

1. Recorded pool size: the candidate count at the moment of the draw (`draw_pool_rounds` returns it; the manifest copies it).
2. One manifest record per (round, K index, class), undefined and not drawn included (`round_cells`).
3. `quota` and `d_fraction` as tables by K index text and class letter; `d_fraction` accepts decimal text, a float and `a/b`; both refused without `r_max` where R22 says so (`d_fraction`, `db_size`).
4. A seed above the supported range is refused naming `seed`: the check is `v <= typemax(Int)` in `_int` (so a seed of 2^63 to 2^64 is refused too, the seed being an `Int` in the whole module; before this it was an `InexactError`).
5. Log text exactly as pinned (lines, count order free in the parser, chain by SHA-256 of the entry text, first entry the SHA-256 of the empty text); entry N is the 1-based position; text before the first `round=` line is "entry 1".
6. Ref arm: arms with `release = "each_round"`; digest over the chunks of the plan that have `:r<round>` in the cell name (controls chunks are not in it) and a round up to r; compared with the log at plan time. The comparison with the `chunk_hash` attribute of the files in the store is not made (the plan has no store).
7. `D` and `D_B` need `in_D` in the panel: refused at plan time naming `in_D` whatever the log says (also for a block without release keys).
8. `check = true` without a log: `R_plan` 0, not closed, no item, the text says "no round released (dry plan)"; the block and `in_D` are still validated.
9. Two blocks with keys in another order give equal plans but for the `log_sha256` line (the block text is canonical).
10. Chunk order: each_round pass for round 1, 2, ..., then closure passes round ascending (the closure order inside the phase is round ascending, not block order; P84 and E41 differ and the tests do not decide).
11. The canonical read-back of D94o leaves out `plan_sha256`: nothing to do in the code (the attribute is still written).
12. K = 0 and K > 0 use the same streams, indexed by K index.
13. Nothing in the builder needs `in_SD` and `in_SL`; the builder still sets them (the base rank rule runs on the round cell list, controls keep `in_SD`).

## Deviations from the spec (each with its reason)

None that change a number or a byte of a rule. Choices where the spec is silent or the tests ask for more:

1. `quota` with `screen_source = "lattice"` is refused (R22 item 2 says the lattice "ignores" it, which would give an empty panel with `n_c` absent). The other new keys are accepted with the lattice and ignored by it, except `controls` and `newton_grid_side`, which act on both sources.
2. An arm without a release key in a block that has release keys is planned with the each_round phase (the spec says only that it is "not subject to the log"). A panel with rounds and a block without any release key is planned as before (rounds ignored, one cell per K and class).
3. Controls and h4 starts (no round) go in round 1 of every phase, cell `controls`, and are not in the ref digest.
4. `RunHooks.read` is also called for the read-back of the temporary file of the writer (`:datasets`, then `:attrs`, path `...chunk-NNNNN.h5.tmp`), through a new `onread` keyword of `write_chunk` (default does nothing): D94n needs reads in the first round, when the store is still empty.
5. `classify_store` calls the read hook for every file it opens and treats a file whose ordinal is beyond the plan as stale without opening it (E43, I32: a closure file is never read before closure).
6. `check_backends!` ignores a (precision, method) the previous session's table does not hold, for plans with rounds only; otherwise the closure run after round-by-round runs refused (the old table had L1 only).
7. `run_campaign` gains `round_log`, `arms`, `up_to`. `arms` and `up_to` are validated (a closure arm without a closed log, a round above R_plan, an unknown arm name are `ArgumentError`s before anything is written) and then select the chunks that run. The entry point has only `--round-log`.
8. The manifest arm entry of an arm with no chunk in a plan with rounds has status `not released` (it was `complete` for an empty chunk list).
9. `_manifest` takes the cell rows instead of the cells; `totals.drawn` is `length(orbs) - length(controls)` (equal to the old sum).
10. `run_pilot` is not changed: a block with release keys and no log is refused there too (the pilot does not plan with rounds).

## Names against "API assumed by the D94 tests"

All adopted as written: `random_screen_points`, `screen_stream_seed`, `round_stream_seed`, `d_stream_seed`, `db_stream_seed`, `draw_pool_rounds`, `d_size`, `d_indices`, `db_allocation`, `db_subset`, `plan_run(...; round_log, check)`, `plan.rounds` `(R_plan, closed, log_sha256, ...)`, `run_campaign(...; round_log, arms, up_to)`, `RunHooks(read)`, `--round-log`, manifest keys, pool sidecar, chunk cell `k<i>:<class>:r<round>`, manifest `rounds` (`R_plan`, `closed`, `log_entries`, `log_sha256`, `D_size`, `DB_size`, `DB_allocation`, `release`), arm key `release` and subsets `D`, `D_B`. Extra names in the module (exported): `read_round_log`, `check_round_log`.

## Tests believed wrong

| Test | Reason |
|------|--------|
| `test/experiments/ex_d86_guards.jl:137` (D86x) | asserts that the usage text holds exactly the 16 options of D86 and names no other. D94u (`ex_d94_grammar.jl`, "usage text gains exactly one option") and the spec (R26 item 3, D94u) require `--round-log` in the usage text, once. The two assertions cannot both hold for any implementation. The precedent is D58, whose stage 02 added `--prior-record` to this list. Stage 02 of D94 should add `"--round-log"` to the list at line 138. Everything else in D86x passes. Owner decision at this stage: added (recorded rework, plan subsection "Rework at stage 03"); D86x then passes (8 of 8). |

## Performance

Measured on the host of the runs (4 threads), small sizes only; the study size runs are stage 04.

| Item | Time |
|------|------|
| `random_screen_points`, M = 262 144 (one K) | 71 to 84 ms (the 256-bit products dominate) |
| `draw_pool_rounds`, M = 262 144, r_max 4, sizes 250, 250, 25 | 3 ms |
| `d_indices` (n = 250) and `db_subset` (|D| = 1 300, b = 100) | 0.01 ms |
| `plan_run` of a hand panel of 3 000 orbits in 3 rounds, empty log (R_plan 1) | 29 ms |
| D94 testsets (28) in the package suite | 37 s in total; panel builds about 0.2 to 1.7 s each; the plan and run testsets 4 to 8 s |

`bench/experiments.jl` gains: `random_screen_points M = 262144`, `draw_pool_rounds M = 262144, r_max 4, quota 250`, `d_indices n = 250 and db_subset |D| = 1300`, and `plan_run in rounds, 3000 orbits, R_plan 1` (with its fixture functions). BenchmarkTools is not installed in this environment (the bench file could not be run as written); the numbers above come from the same calls timed with `@elapsed` (minimum of 3), and the bench file was run with a stand-in for the macro to check that its fixtures build and the entries run.

Memory: the screening of the random source holds the points of all K (2 M binary64 each, 4.2 MB at M = 262 144, 6 K: 25 MB) and the class bits (1 byte per point); the base chunking of the screening values is unchanged. The draw copies the candidate list per round cell (at most M entries).

## Quality

- Coverage (`--code-coverage=user`, D94 script and the package run give the same figures): `inputs_rounds.jl` 92 of 92 executable lines, `rounds_streams.jl` 50 of 50, `rounds_panel.jl` 92 of 92, `round_log.jl` 57 of 57, `run_rounds.jl` 72 of 75 (new files 363 of 366, 99.2 %); the three lines not run are the body of `_selected` (the `arms` and `up_to` selection of chunks, no test asks for a valid selection; the refusals around it are tested). Changed lines in the old files: 59 of 59 run. All of `src/Experiments/` in the D94 script alone: 98.4 % of the executable lines of the 44 files that loaded.
- Aqua and JET clean (quality 12 of 12 in the final package run); layering test passes (2 of 2), and D94u (no Analysis import in Experiments, no Experiments import below).
- Sizes: new files 97 to 149 lines; largest changed file `run_plan.jl` 345 lines (was 337), `run_block.jl` 352; every function under 50 lines (checked by script on all changed and new files). No `@fastmath`. No em dash in the new code or comments.
- Numerics: no accumulated logs in the new code (no float comparison decides a count, I36); the screening sum is the base compensated one.

## Known limits

- If an each_round arm and a closure arm contain an identical item (same method, precision, horizon, tau, direction, mode and orbit), the item is made once, in the each_round pass, with both arm names; its chunk is then a ref chunk and is planned before closure. The fixtures keep the arms apart by method. A study that wants the same method in both phases needs a distinct tau, horizon or direction (or a decision at stage 01).
- The digest of the ref arm excludes chunks without `:r<round>` in the cell name (controls, h4); a controls chunk is not covered by the guard.
- The log guard is against accident (P82): a hash chain and a digest of the plan's chunk hashes; it does not read the chunk files in the store.
- Seeds above 2^63 - 1 are refused (the module keeps the seed as an `Int`).
- `run_pilot` does not accept a block with release keys without a log.
- A K = 0 cell with a quota above 0 in a class undefined at K = 0 is refused; at K = 0 an absent or 0 quota for C, S, X is accepted.

## Stage 03 audit

| Check | Result |
|-------|--------|
| GREEN | D94 files 2 390 pass, 0 fail, 0 error; whole suite 1 403 998 pass, 1 fail, 1 broken. The one failure is `ex_d86_guards.jl:137`, which contradicts D94u (table above); the broken result is not in experiments. After the owner-approved rework (`--round-log` added to the D86x list) D86x passes 8 of 8; the whole suite was not rerun after this one-line test change. |
| Quality | Aqua and JET clean; layering 2 of 2; coverage of the new files 99.2 %, of the whole module in the D94 script 98.4 % (at least 80 %) |
| Size | functions under 50 lines, files under 400 (largest 352) |
| Tests untouched | fail by the letter, recorded rework only: one line of `ex_d86_guards.jl` (D86x option list gains `--round-log`, owner decision); no D94 test changed, no assertion weakened |
| Numerics | no `@fastmath`; no accumulated log added; exact integer and rational arithmetic for sizes, fractions and allocations |
