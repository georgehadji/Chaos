# Implementation notes: experiments, Amendment D86 (module-build stage 03)

Julia 1.13.2, 4 threads, `--check-bounds=yes`. Nothing committed. Tests, spec, `schema/`, `_config/` untouched. No dependency added to `Project.toml`. This file is read by stage 04.

## Result

Folder `test/experiments/` (outer testset, `runtests.jl` included as in the RED log):

| Run | Pass | Fail | Error | Total | Time |
|-----|------|------|-------|-------|------|
| RED (stage 02) | 4 901 | 13 | 32 | 4 946 | 3 min 59 s |
| GREEN (this stage) | 6 763 | 29 | 18 | 6 810 | 10 min 55 s |
| of which the old tests (D53 to D59, panel) | 4 861 | 0 | 0 | 4 861 | about 3 min 25 s |
| of which D86 | 1 902 | 29 | 18 | 1 949 | 7 min 20 s |

The old tests are unchanged (4 861 pass, as in the RED log). Every one of the 47 non-passing D86 results is a defect of a test or an assumption that cannot hold (table "Tests believed wrong", 17 entries); for each of them I checked with a scratch copy of the test (in the scratchpad, never in the repo) that the rest of the testset passes once the one defect is removed: D86i, D86j, D86k, D86l, D86m, D86s and D86u are then all green. No D86 criterion marked (S) fails because of the code as far as I can tell; the list is for the owner to confirm or to send back to stage 02.

Full suite (`Pkg.test()`, 4 threads): first run with the code before the JET fix below, 1 400 074 pass, 30 fail (29 experiments + 1 JET), 18 error, 1 broken (not in experiments), 23 min 38 s; all other modules pass. Second run (final code): 1 400 075 pass, 29 fail (all in experiments, the test defects below), 18 error (same), 1 broken (not in experiments), 19 min 40 s; quality 12 of 12 and layering 2 of 2 pass; experiments 6 763 / 29 / 18 in 11 min 20 s.

Per D86 testset (pass / fail / error): a 30, 5, 20, 5 / 0 / 0; b 8; c 169 and 5; d 206 and 20; e 777; f 7 / 4 / 0; g 3 and 34 / 16; h 8 and 4; i 1 / 0 / 1 and 3 / 1 / 0; j 9 / 2 / 12; k 2 / 1 / 1; l 0 / 0 / 1; m 0 / 0 / 1; n 19 and 17; o 5; p 3 and 5 / 3; q 218; r 17 / 1; s 1 / 0 / 1; t 144; u 103 / 0 / 1; v 29 / 1; w 17; x 8.

Red items that cannot be green in this stage by design: D86y, D86z (stage 04 and 05), the (V) and (G) parts of D86s, D86w, D86p, and the [full] part of D86b. They are not in the package suite.

## Files

New in `src/Experiments/` (7 files, all under 400 lines, all functions under 50 lines):

| File | Lines | Holds |
|------|-------|-------|
| `run_indicator.jl` | about 190 | grammar of `[[run.indicator]]` and `indicators`, validation in two passes (R11), `spec_text`, `spec_hash`, `ind_run_length` |
| `run_indicator_chain.jl` | about 190 | chains under the one-replacer rule, observables, reducer tuples, `orbit_bytes`, `group_capacity`, plan-time checks (R13) |
| `run_indicator_record.jl` | about 330 | the CPU loops of the six kinds, the tables of E24 and E25, frame checks, `hot_loop_bytes(spec)` |
| `run_indicator_device.jl` | about 210 | the two kernels (series kinds, window kinds), `ind_device_records` |
| `run_indicator_store.jl` | about 120 | `indicators/<id>/` datasets, provenance groups, pooled counts, root attributes |
| `run_indicator_health.jl` | about 40 | E31 finiteness, `complete` and `failed` defects |
| `run_indicator_report.jl` | about 70 | manifest table, flag counts, omitted orbits, per-spec cost lines |

Changed: `Experiments.jl` (imports `..Indicators`, includes, two exports `orbit_bytes` and `group_capacity`), `run_block.jl`, `run_plan.jl`, `run_backend.jl`, `run_record.jl`, `run_store.jl`, `run_health.jl`, `run_campaign.jl`, `run_groups.jl`, `run_device.jl`, `run_manifest.jl`, `run_session.jl`, `run_pilot.jl`, `run_resume.jl`, `CONTEXT.md` (66 lines). `bench/experiments.jl` gains six entries for the indicator recorders.

## Design in one page

- Block: two passes inside `run_block`, before anything else. Pass 1 validates each table alone (`_rb_indicator`), pass 2 validates specs against arms (`_ind_cross`: used by an arm, budget per precision of its arms, interval against tau, mode rules of E26, horizon equal tau times run length, one tau, no BigFloat). `control_methods` refuses a spec id. The canonical text of the block lists `run.indicator[i].*`; a block without specs has no `indicator` key, so its text and hash are unchanged (I18, D86b golden values pass).
- Items: `_ITEM_FIELDS` ends with `kind` and `spec`, printed only for indicator items; `direction` is printed only when present. The arm products enumerate methods, then specs; a spec without a direction is not swept over the directions. Applicability (Y5 has no torus: no `exp_i_q`, `cos`, `sin`) is decided at plan time and listed in `plan.not_applicable` and in the plan text. The checks of R13 item 1 run per (spec, orbit, horizon, tau, mode) and name the spec and the orbit.
- Recording: one function per kind calls the verified reducers (never a copy): the tangent, pair and stretching kinds advance between grid steps with `Tangent._orbit_body` and read the finals; the window kinds step one handoff at a time. The stretching workspace is one `UInt32` vector per orbit. The pair monitor reads orbit A (a method of `Tangent.reducer_step` for `RunMonitor{<:PairMap}`).
- Device: the existing machinery (`_launches!`, `_by_system`, `_todev`) with two new kernels; the grid rows are written through the same `_put!` methods as the host (`HostRow` and `DevRow`), so the two paths cannot diverge; the window tables are read on the host from the carried state. Stretching and frequency return the sixth fallback reason.
- Store: `OrbitRecord` gains `indicators::Dict{String, NamedTuple}`; `lyap` is optional on read; `write_chunk` writes the `provenance` groups (extracted from the root attribute dictionary) and `pooled/<id>` and compares them after the write.

## Deviations from the spec (each with its reason)

1. Window kinds step one handoff at a time and reset the deviation column to the constant `(1, 0)` after each handoff, for `phase_advance` and `birkhoff` as well (R14 item 4 says "the grid loop of item 1"). Reason: the unnormalised column of the one-column stream grows like exp(lambda n) and overflows Float64 in about 700 steps on a chaotic orbit, so every chunk of the real windows (8 000 handoffs) would fail with "deviation matrix". The reducers never read the column and the orbit does not depend on it, so the tables are bitwise those of the stand-alone reducer (D86l, D86m, D86q pass). The same loop is in the kernel.
2. The E1 buffer is per orbit (2 W_max complex values, filled by the orbit pass, transformed at once on the same worker, then dropped) and not per launch group (R14 item 3, R17 item 2). The stored bits are the same (sequential transform per orbit, I24); the peak memory is `threads x 32 W_max` bytes, below the group bound that E28 charges to `bytes(orbit)`. `orbit_bytes` still charges B as E28 says, so the capacity is the specified one.
3. `[run.pilot].horizon` and the rule "the pilot horizon is a multiple of every tau" ignore an arm that holds window specs only (R20 item 1: the pilot runs those at their own run length). Needed by D86m (K = 31, H = 1: run length 32 against the pilot horizon 40 of the fixture).
4. `launches_per_chunk` of a CPU chunk of an indicator spec is 1, whatever the capacity (for the methods it stays `ceil(orbits / capacity)`). Reason: D86n and D86q compare every root attribute except `chunk_size`, and I24 wants the file independent of the capacity; a CPU recording pass is one launch.
5. A launch capacity that does not fit one orbit of a spec (R17 item 3) creates the empty store directory before raising `ArgumentError` (D86n reads `readdir(store)` after the refusal). No file is written.
6. `first_full_step` in the stretching scalars is 0.0 when the handoff count never reaches `min_per_cell x bins` in the run (no grid step to name).
7. `window kinds`: `members` and `chain` in the file are the single name of the kind (`phase_advance` and so on), `provenance/<kind>` holds the reducer's keys; for `frequency` the group holds the settings of the transform. For stretching and Birkhoff the provenance, the range in the provenance and `observable_names` come from the first item of the chunk (the range of each orbit is in `scalars@range`; the names differ between a standard orbit and the torus-free control only for the `lifted` observable).
8. Stretching count matrices: the per-orbit counts of the whole chunk stay in memory (`pool`, 8 (bins + 2) bytes per orbit) until the chunk is written, to sum the pooled counts of E29; only the orbits of the named subset store theirs. These bytes are not in E28 (it charges one workspace per slot).
9. `_same(::OrbitRecord, ...)` compares the attributes with `isequal` (a failed pair stores `Dx_r = NaN`).
10. Frame preparation runs in `_open_session` right after the context is made, for all items of the plan, so a dependent frame refuses the start before any file or manifest (R13 item 2). `prime_frames!` keeps its other callers; the check runs once per (frame, spec).
11. The cost cap for a group of spec items reads the budget of the spec (`budget_line_ns(...; spec)`); the manifest arm entry gains `specs` only for arms with indicators.

## Names against "API assumed by the D86 tests"

All adopted as written, none renamed: `EX.spec_text`, `EX.spec_hash`, `plan.not_applicable`, `plan.steps.specs`, item fields `kind` and `spec` (`nothing` for the methods), `RunHooks(frame, series)`, `resolve_run_backend(...; kind, spec_sha256)`, `EX.orbit_bytes`, `EX.group_capacity`, `OrbitRecord.indicators`, `rerun_check` rows with `gated`, pilot rows with `horizon`, `bytes_per_orbit`, `capacity`, `flag_shares`, `orbit_pass_seconds`, `transform_seconds`. Additions: `RunPlan.omitted_cells` (orbit and cell of each omission, used for the per-chunk `not_applicable`), `chunk_entry(...; plan)`, `RunCtx.specs/hooks/bytes`, `hot_loop_bytes(spec; T)`, `block_spec`, `ind_*` helpers.

## Tests believed wrong or unreachable (file:line, reason)

Each was checked as stated; the cited expression cannot be true for any implementation of the spec, or contradicts another test.

| Test | Reason |
|------|--------|
| `ex_d86_chains.jl:142` | `PlanarObservable(1, 3)` with `d = 2` is valid: the index limit is `2d = 4` (the wrapped state has 4 entries). The refusal needs index 5 or 3 with a map of dimension 1. |
| `ex_d86_chains.jl:152-154` | D86f(d): at the lifted start 1e9 the first coordinate of B equals A (one ulp 1.2e-7), but the second coordinate separates by 1e-8 |u2|; the pair is refused only if |u2| < 1/2. With the seeded xi_0 of orbit 1 (`initial_deviations(run_seed, 1)`, u = (-0.544, 0.839)) Dx_r = 8.4e-9 lies inside [Dx/2, 2 Dx]. The test depends on the convention of A8. |
| `ex_d86_chains.jl:246` | MEGNO closed form of M10(a): the product runs over `i = 0 .. n-1` of `i^2 + 1` (megno-spec M10), the test sums `j = 1 .. n`. The recorder equals the stand-alone `MEGNO()` bit for bit (D86g, first testset); the reducer equals the spec form (n = 1 gives ln 2). 16 assertions. |
| `ex_d86_windows.jl:57` | `t[findfirst(...)]` is `t[nothing]` for the tau = 4 items of the block (25 handoffs against `min_per_cell x bins` = 40). Stops D86i before its checks. |
| `ex_d86_windows.jl:107` | "edge condition": with the derived range `hi = L (1 + 8 eps)` the stretching numbers of Y5 sit at the edge by construction, `u = 21/2 (1 + L/hi)` is 21 minus 1e-15; the following counts and D comparisons use the wrong cell. |
| `ex_d86_windows.jl:149, 150` | `big(...)` is called, but `big` is assigned a block at line 162 of the same testset, so it is a local there (UndefVarError). 12 errors. |
| `ex_d86_windows.jl:171` | D86j(c): the first call over the two control orbits (constant series) allocates about 340 KB of one-time cache inside the transform of Indicators (not a leak: later calls grow by 41 KB < 103 KB). The assertion passes when the warm-up includes a control. |
| `ex_d86_windows.jl:208` | D86j(f): the loop also visits the attributes `final@columns`, `final@rows`, `flags@bits` (vectors of strings), which the test itself reads in (b). |
| `ex_d86_windows.jl:268, 288` | D86k: `f0` has `mean_removal = false`; the fixed point (0, 0) gives the constant series 1, which is not degenerate (the pure function agrees, the test's own comparison passes). `degenerate` is set for `fm` and `fp` (checked: 2 per controls chunk). |
| `ex_d86_windows.jl:310, 328` | D86l: Y5 in wrapped mode has a 2-component state, `LiftedAngleObservable(1, 3)` is out of bounds; `plan.panel.orbits[o]` for the h4 starts 16 and 17 is out of bounds. |
| `ex_d86_windows.jl:342` | D86l: lifted K = 0 map, orbit 6 (p0 = 6.1): error 16.54 eps against the bound 16 eps (the lifted coordinate reaches 400, each increment is rounded at 4.5e-14). The rest passes. |
| `ex_d86_windows.jl:378` | D86m: `plan.panel.orbits[o]` for h4 ids (BoundsError). With this line removed, D86m is green (5 assertions, including K = 31 and the K = 0 map). |
| `ex_d86_run.jl:207, 212` | D86p: the canonical hash includes `backend`, `backend_fallback_reason`, `device`, `gpu_chunks`, which differ by design between a device run, a fallback run and a CPU run. All datasets are bit-identical (checked file by file). |
| `ex_d86_run.jl:318` | D86r: same, plus `attempt` 2 and `first_attempt_reasons` of the rerun. |
| `ex_d86_store.jl:57` | D86s reads `scalars` for every window kind, including `frequency`, for which E25 says "none"; the same test asserts `!haskey(lay, scalars)` for `frequency`, and D86j(f) forbids it. |
| `ex_d86_store.jl:168` | `all(a -> haskey(a, "specs") ..., man["arms"])` iterates the pairs of a `Dict` (MethodError); needs `values(...)`. |
| `ex_d86_store.jl:228` | D86v: the arm `a` has one run type and its four chunks form one launch group, so the pause cannot leave a chunk pending: `paused == ["a"]` holds and the cost cap acts (reason recorded), the status is `complete`. The old tests (D53l) rely on later groups. |

## Assumptions of stage 02

A8 true for D86h (the stored xi0 is a unit vector, the direction source equals the frame rule, the stand-alone pair run reproduces the columns) but D86f(d) depends on it (see above). A9 true. A10 true (SALI and FLI closed forms pass with the stated bounds). A11 true (16 eps for windows 32 and 64). A13 true (`undefined window` is the bit name, refined bits are zero). A14 true (`OrbitRecord.indicators`, mutated by D86r). A15 true (`RunHooks(frame, series)`, the clock rebuild). A18 true (the device batch is 8 with the 392 bytes of the chain, 9 with the allowance).

## Performance

Single thread, N = 1e5 steps (window kinds: their run length), standard map, wrapped, Float64, `bench/experiments.jl` functions timed by hand (BenchmarkTools is not in the project environment):

| Recorder | ns per orbit-step |
|----------|-------------------|
| tangent FLI, MEGNO, SALI (17 columns) | 142 |
| tangent SALI | 43 |
| stretching, 1000 bins | 51 before the follow-up, 47 after |
| frequency orbit pass (W = 1000, 4000, run length 8 000), transform not included | 31 |
| phase advance (1000, 4000) | 74 |
| Birkhoff K = 4 (1000, 4000) | 157 |

Hot loop allocation per kind after warm-up (`hot_loop_bytes(spec)`, Float64 and Double64): 0 bytes for every kind. Before the follow-up stretching allocated 4 096 bytes per orbit (20 bins, 23 grid steps) because `reducer_final(r, st, ws)` copies the count matrix at each grid step; the grid row now computes D with `spectral_distance` on views of the two count columns of the workspace (the same function and operands, so the bits are unchanged: D86i, D86q and the chain/device comparisons give the same results as before). Stretching at 1000 bins and N = 1e5: 1 316 503 bytes allocated per orbit before, 336 after (loop closure of the timing harness only); 51.0 to 46.9 ns per orbit-step. The other timings above are unchanged within noise (tangent 140.5, SALI 43.2, frequency pass 30.6, phase advance 73.2; Birkhoff K = 4 measured 157 before and 128 after with no code change, so that difference is machine noise); the few hundred bytes the harness reports for the other kinds are the timing closure, not the loops (`hot_loop_bytes` is 0).

Suite time: the D86 testsets add 410 s to the folder run (target 180 s, proposed in the criteria preamble). 275 s of it are two testsets whose cost is in the test: D86t 180 s (31 campaigns of 71 chunks, each about 2 s, plus a canonical read-back of every file for each) and D86d independent hashing 96 s (the test starts `sha256sum` 6 000 times to hash the spec of each seeded item). The rest is 135 s. A campaign of 71 chunks takes 2.0 s warm and 34 s the first time (compilation).

## Known limits

- The device path is verified on the KernelAbstractions CPU backend only; the GPU tag (`verify/gpu.jl`, `verify/experiments-device.jl`) must be rerun on the local GPU at stage 04 for the kernels `_ind_series_kernel!` and `_ind_window_kernel!` and the rows of `_put!`. Float64 only; the Double64 fallback on a real GPU is untested here.
- The per-orbit seconds of the transform are measured in the worker; `transform_seconds` is thread time.
- Reducer provenance of stretching and Birkhoff in the file is that of the first item of the chunk.
- A frequency series hook (`RunHooks.series`) runs inside the parallel region; it must be thread safe.
- The window kinds are recorded with one deviation column that is never read; the stored `diagnostics/det` and `symp` come from the monitor on the orbit.
- `prime_frames!` over the whole plan at the start costs one pass over the items; at 2 million items it is a few seconds (not measured).
- The plan-time checks (R13) are cached per (spec, orbit, horizon, tau, mode); the real P2 plan size was not planned here.

## Gates

| Gate | Result |
|------|--------|
| Layering test (`test/layering.jl`) | pass (2 assertions) |
| Aqua | pass |
| JET | the first run found three possible errors in my device code (a kernel variable shared by two closures, an untyped `Any` argument in the row writer); fixed (two small functions, argument types); rerun in a separate environment: pass (12 assertions of `test/quality.jl`) |
| Experiments folder | see Result |
| Function and file size | longest function below 50 lines, longest new file `run_indicator_record.jl`; `CONTEXT.md` 66 lines |
| No `@fastmath`, no em dash | checked in `src/Experiments` and `bench/experiments.jl` |
| Coverage | not measured here |
| Tests untouched | `git status` shows no change under `test/` |

## Final verification (orchestrator, 2026-10-28)

- Folder `test/experiments` on the final code and the corrected tests: 6830 pass, 0 fail, 0 error, 5 min 56 s (4 threads, `--check-bounds=yes`).
- Full `Pkg.test()`: 1,400,142 pass, 0 fail, 0 error, 1 broken (outside experiments, as before), 15 min.
- Coverage of `src/Experiments/` by the folder suite (`--code-coverage=user`): 96.1 percent of 3440 lines; the new `run_indicator*.jl` files 94 to 100 percent.
- Tests untouched since stage 02: the tests changed once after the stage 02 commit, in a recorded correction round (17 defects, see the test plan, "Correction round"); stage 03 itself changed no test.
- Open for the owner: the 11 deviations above (candidates for a D87 record) and the lifted-mode frequency bound of the correction round (16 eps plus max|q| eps).

