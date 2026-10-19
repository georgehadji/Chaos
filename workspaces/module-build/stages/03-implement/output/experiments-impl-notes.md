# Implementation notes: experiments (module-build stage 03)

Branch module/experiments, base ceb5227. Not committed. Julia 1.13.2 (1.10 is not installed on this machine), 2 threads for `Pkg.test`, 4 threads for the suite-alone runs.

## Result

`test/experiments/`: 2380 pass, 1 fail, 0 error (2381). The one failure is a test that contradicts another test of the same suite (see "Test believed wrong"); nothing in `src/` can satisfy both. Full suite (`Pkg.test(coverage=true)`, 2 threads): 52001 pass, 1 fail (the same assertion), 0 error, 5 min 13 s. Aqua 11/11 and JET 1/1 (quality testset 12/12), layering 2/2, every earlier module green (diagnostics 245, lyap-benettin 15009, lyap-gram-schmidt 18938, lyap-householder-qr 575, lyap-svd-stabilized 12915, lyap-two-trajectory 493, precision 150, standard-map 474, tangent-core 808).

## Files

- `src/Experiments/` (17 files, 1263 lines, longest file 149, longest function 27 lines): `Experiments.jl` (imports, exports, includes), `parallel.jl` (`_par`, index-keyed work distribution), `geometry.jl` (grid, torus distance), `inputs.jl` (E0 parsing and validation), `prereg.jl` (E18 hash, panel block), `periodic.jl` (residual, Newton, trace classes), `orbits.jl` (primitive period, representative, dedup, separation), `search.jl` (`periodic_set`, `periodic_distance`), `screening.jl`, `predicates.jl`, `draw.jl` (streams, seeds, placements, Y5), `subsets.jl`, `provenance.jl`, `cells.jl` (screening flags, per-cell draws), `panel_tables.jl` (tables, controls, ids, subsets, manifest pieces), `panel.jl` (`build_panel`), `portrait.jl`.
- `src/Chaoslab.jl`: include of `Experiments/Experiments.jl` after Lyapunov, `using .Experiments`, Experiments added to the re-export loop.
- `Project.toml`: `SHA` in `[deps]` (uuid ea8e919c-243c-51af-8825-aaa63cd721ce) and `SHA = "0.7, 1"` in `[compat]` (stdlib 0.7.0 on Julia 1.10, 1.0.0 on 1.13). `Manifest-v1.13.toml`: resolved (Chaoslab deps gain SHA, `julia_version` 1.13.2, `project_hash`). `Manifest.toml` (Julia 1.10) is not touched: 1.10 is not installed here, SHA is already in it as a stdlib, only its `project_hash` is stale. Run `Pkg.resolve()` on 1.10 once.
- `bench/experiments.jl`: four entries.
- No lower module changed. `test/experiments/` untouched.

## Design

- Imports: Foundation, Systems, Tangent, Lyapunov (for `safe_norm`), plus StaticArrays, LinearAlgebra, StableRNGs, SHA, TOML. Reused: `StandardMap`, `forward`, `tangent`, `wrap`, `LinearSymplecticMap`, `hyperbolic_block`, `CompensatedSum` (`add`, `total`), `safe_norm`. `orbit_seed` calls `Tangent._mix` (tangent-core's seed rule, unexported); exporting it from Tangent would remove the private access.
- Prescribed screening direction (owner decision P8): tangent-core already accepts `W0` in `propagate`, so no lower change is needed. Screening does not call `propagate`: with k = 1 the Householder reflection gives R_11 = ||W||, so `screen_point!` carries the vector itself (`safe_norm`, compensated sum of `ln R_11`, block exponents at multiples of w) with Systems `forward` and `tangent`. This keeps the hot loop allocation free (AC26) and independent of the propagate machinery; it is checked against the verified Benettin estimator (AC25, 200 cells, all within the S3 bound) and against the closed form at K = 0 (AC12).
- BigFloat precision: `build_panel` and `periodic_set` run inside `setprecision(BigFloat, bits)`; spawned tasks inherit the scoped precision on Julia 1.11 and later. On 1.10 `setprecision` sets a global, which is also safe here because nothing changes it while the tasks run.
- Work distribution: `_par` hands out item indices from an atomic counter to at most `ntasks` tasks (default `Threads.nthreads()`); results are stored by index, so the panel does not depend on the task count or the chunk size (AC19 passes for ntasks 1, 2, 4 and chunk 16, 4096). A first version shared one local variable between the spawned closures and gave a task-count dependent panel (found by AC19); the worker is now a separate function.
- Newton (`_newton_seed`): windings fixed from the seed, closed-form step of E3, accept at residual <= rho, polish towards rho^2 and stop at stagnation (residual not decreasing) or the cap. Abandon reasons counted in the manifest: singular (|tr J - 2| <= 4 eps), diverged or not finite (lifted coordinate beyond 1e4, proposed bound), cap. `q_max` above 12 is rejected (proposed supported range).
- Dedup (E9) and the separation rule use Float64 torus distances of the wrapped points (delta_dup of order 1e-8 is far above the Float64 resolution of the coordinates); the orbits, traces and windings are computed at 256 bits.

## Deviations from the spec

1. The E9 separation refusal (smallest distance between distinct orbit points below 10 delta_dup) is raised by `collect_orbits` (default) and by `build_panel`, but not by `periodic_set`, which records `min_separation` and returns. AC5 requires `periodic_set` at K = 2 pi (binary64) with q_max = 1 to return the two fixed points; at that K the residual test accepts a spurious pair near the saddle-node of the accelerator modes (see "Finding for the owner").
2. The periodic table of the panel file has rows for K = 0: one per reduced fraction a/b with b <= q_max, class parabolic, trace 2, exponent 0, windings (a, 0), b points (0, p), (p, p), ... (E7). AC24 requires at least one row per K. They are not controls and not Newton orbits.
3. The Y5 control carries `k_index = 0` (AC18 reads `k_index` of every control). It has no `K_text` or `K`.
4. Manifest and sidecar extras beyond the test plan: `counts` (periodic orbits per K, period and class), sidecar keys `backend_requested`, `threads`, `ntasks`, `julia_version`. `qualifying` is the candidate count after the S exclusion for C, the number of elliptic periodic points for R, 0 for X. A short T, C or S cell has reason "fewer qualifying grid cells than n_c".
5. `place_R` also returns the unrounded 256-bit sum (`xb`, `pb`) so that the decimal texts of R starts are the pre-rounding values (spec Algorithm 9). The orbit seed is the `Tangent._mix` key, written as a decimal string.
6. The default chunk is min(grid, free memory / 4 / bytes per cell) as the spec says; since the loop keeps one buffer per task, chunking only bounds the scheduling unit.

## Test believed wrong

`test/experiments/ex_panel.jl:57`, testset "panel structure and invariants": `@test allunique(cs) && isdisjoint(cs, w.(dr))` where `cs` are the wrapped starts of all controls, ignoring K. AC18 (same file, line 159) requires the fixed-point controls `[(2, 0.0, 0.0), (2, pi, 0.0), (3, 0.0, 0.0), (3, pi, 0.0)]`, and the Y5 control starts at `(0.0, 0.0)` (AC18, spec E19). The wrapped starts (0, 0) therefore appear three times and `allunique` is false in any conforming build (checked: `unique=false`, `disjoint=true`). Spec I3 says a control appears once as an orbit (distinct K, distinct system), not that starts are distinct across K. The assertion should compare (K index, system, start) or drop `allunique`. Not edited.

## Finding for the owner (study impact, not a module defect; superseded by D51, see the D51 section)

The preregistered K list contains "2pi", the binary64 value 6.283185307179586, which lies 2.4e-16 below the true 2 pi, the saddle-node of the accelerator-mode fixed points (sin x = 2 pi / K). No such orbit exists at that K, but the acceptance rule (residual <= rho) admits points around x = pi / 2 where the residual is about K delta^2 / 2 + 2.4e-16. Measured with the preregistered rho = 1e-12, delta_dup = 1e-8, q_max = 6 (N_g = 32): 312 orbits and a minimum separation of 1.08e-8, below 10 delta_dup = 1e-7, so `build_panel` refuses (E9) for the frozen inputs. The other K values give separations from 0.0008 to 0.12. Options for an amendment: lower delta_dup (the 10x margin then needs a smaller separation), accept the pair as one orbit, or handle trace-degenerate solutions explicitly (the trace there is 2 within about 1e-8, inside the parabolic band for rho = 1e-12 only when the band is widened). The frozen panel cannot be built until this is decided.

## Performance (Julia 1.13.2, this machine, other processes idle)

Bench (`bench/experiments.jl`, run from a scratch environment with BenchmarkTools, 4 threads):

| Entry | Time | Allocation |
|-------|------|------------|
| `screen_point!`, 1e5 steps, K = 5 | 4.24 ms (42.4 ns per orbit-step) | 0 B |
| `screen_point`, 2000 steps, K = 1.5 | 38.3 us | 272 B |
| Newton solve, q = 2, K = 5, 256 bit | 83.5 us | 56 kB |
| `build_panel`, 3 K, 8 x 8 grid, N_scr 2000, q_max 2 | 57.3 ms | 31 MB |

Suite time: experiments 26 to 30 s inside the full run (compile included), budget 3 min per slug (D35); whole `Pkg.test` 5 min 13 s on 2 threads, budget 30 min. Coverage of `src/Experiments`: 641 of 648 executable lines, 98.9 percent (`Pkg.test(coverage=true)`, `.cov` files counted then deleted); uncovered: four lines of `inputs.jl` (parse-failure branches), one in `orbits.jl`, one in `prereg.jl`, the `git` failure branch.

Full p1-le-benchmark panel (6 K, 256 x 256 grid, N_scr = 1e5, q_max = 6), estimated on 23 threads, assuming ideal scaling of the parallel parts:

- Screening: 6 x 65536 x 1e5 = 3.93e10 orbit-steps x 42.4 ns = 28 thread-minutes, about 75 s on 23 threads (the preregistration's Budget table used 33.09 ns per step; this loop is 28 percent slower per step). The class predicates add a distance scan over the periodic points only for cells that pass the cheap tests.
- Newton: measured 0.31 ms per seed solve including all periods (K = 0.5 to 10, N_g = 32: 1.8 to 2.3 s per 6144 solves on one thread); 5 x 65536 x 6 = 1.97e6 solves, about 10 thread-minutes, well below the spec's one thread-hour estimate. Measured at K = 10, N_g = 128: 28.4 s on one task, 9.8 s on four (parallel part 24.8 s, sequential dedup about 3.6 s). Extrapolated to N_g = 256: about 100 s parallel plus 15 to 30 s sequential at K = 10, less at lower K.
- Total: about 4 to 6 minutes wall for the whole build (extrapolated, not run: the frozen inputs are refused at K = 2 pi, see above).

## Known limits

- CPU only by design (P6); `backend` is validated and recorded in the sidecar as cpu with an empty reason. No GPU path.
- `numerics-defaults.toml` is found through `pkgdir(Chaoslab)`: the package must be used from its source tree (development install), as in this repository.
- `collect_orbits` (reduction to the primitive period and dedup) is sequential. Its cost grows with accepted solutions x kept points; fine up to the frozen size, a spatial hash would be the next step.
- `q_max` is limited to 12 (proposed range). `min_separation` is Float64; it assumes delta_dup well above 1e-15.
- Starts of displaced R and X orbits are stored unwrapped, as formed (plan A11); the portrait content wraps them.
- Newton seeds are the grid centres only (P9); an orbit whose basin misses all centres is not found (AC7 passes with N_g = 32 for the period-2 families at K up to 10).
- 1.10: not run here. `Manifest.toml` needs a `Pkg.resolve()` on 1.10 (project hash).

## Audit

| Check | Result | Evidence |
|-------|--------|----------|
| GREEN, including every earlier module | FAIL (1 test) | 52001 pass, 1 fail, 0 error; the failing assertion contradicts AC18 of the same suite (see "Test believed wrong"); all earlier modules pass |
| Aqua | pass | 11/11 checks in the quality testset |
| JET | pass | 1/1: `JET.test_package` over all own modules, no report (4 reports on `pkgdir`/regex `Nothing` unions fixed during the work) |
| Layering | pass | `layering.jl` 2/2 and AC28 2/2; imports are Foundation, Systems, Tangent, Lyapunov and declared dependencies only |
| Coverage of the module at least 80 percent | pass | 98.9 percent (641/648) |
| Functions under 50 lines | pass | longest 27 lines (`screen_point!`) |
| Files under 400 lines | pass | longest 149 lines (`inputs.jl`) |
| Tests untouched | pass | `git diff ceb5227 -- test/experiments` is empty |
| No `@fastmath` | pass | no match in `src/Experiments` |
| Accumulated logs use compensated summation | pass | `CompensatedSum` for S(n) in `screen_point!`; the only other logarithm is a single `log` per orbit exponent |
| Suite time inside D35 | pass | 26 to 30 s against 3 min; whole `Pkg.test` 5 min 13 s against 30 min |
| No em or en dashes in new files | pass | checked on `src/Experiments/*.jl`, `bench/experiments.jl`, this file |

## Amendment D51 (conditioning test, replaces the D50 merge; 2026-10-17)

D50 was withdrawn: its merge code, the `merged` field, the merged class, the manifest `merged` table, the portrait symbol and the `sn_trace` input are removed (nothing of them remains in `src/`). The coordinator's I3 test fix stays. Tests were changed first (stage 02 role, test plan section "Amendment D51"), then code.

Test diff against ceb5227 (uncommitted, includes the I3 fix): `ex_inputs.jl` (AC1: one case, `sn_trace` is an unknown key), `ex_panel.jl` (I3 key by K index; testset "D51 build at K = 2 pi (binary64)"), `ex_periodic.jl` (testset "D51 conditioning test"). `ex_helpers.jl` is unchanged against the stage 02 commit.

Code: `src/Experiments/newton.jl` (new): `error_bound(J, F) = ||F||_2 ||J - I||_2 / |2 - t|` (infinite at `t = 2`; `||J - I||_2` from the closed form of the largest singular value of a 2 x 2 matrix), `well_conditioned`, and `_newton_seed(K, x, p, q, rho, nmax, delta_dup)` moved there from `periodic.jl`. It evaluates the bound at the final iterate in the working 256-bit precision and returns code 4 (residual accepted, bound above `delta_dup / 10`) with `|2 - t|`; a singular step after reaching the residual also lands here (infinite bound). The exact fixed points bypass the Newton path (E5) and are kept. `periodic_set` counts code 4 in `stats.ill_conditioned` and per period in `PeriodicSet.ill` (count, smallest `|2 - t|`); the manifest has `ill_conditioned` rows (`k_index, q, count, min_abs_2_minus_t`) and `newton.ill_conditioned`. `search.jl`, `inputs.jl`, `orbits.jl`, `panel.jl`, `panel_tables.jl`, `portrait.jl` lost the D50 code. `bench/experiments.jl` passes `delta_dup` to `_newton_seed`.

Results: experiments suite 2597/2597 (27 s); full `Pkg.test` 52218 pass, 0 fail, 0 error, 4 min 55 s; Aqua, JET, layering clean; coverage 98.9 percent (656/663 lines).

### Real-panel build (repository preregistration with header hash 650b3558..., Julia 1.13.2, 23 threads, compile included)

- Wall time 406 s; file in the scratchpad `panel-try/panel.toml`. It builds with the frozen file: no refusal. Smallest separation between distinct orbit points 3.1e-4 (the refusal threshold is 1e-7).
- Newton: 1 966 080 seeds, 1 359 242 accepted, 73 963 ill-conditioned, 499 944 diverged, 32 931 cap, 0 singular.
- Ill-conditioned, all at K = 2 pi (K index 5): period 1: 44 032 (smallest |2 - t| 3.4e-8); period 2: 21 693 (1.5e-15); 3: 3 475 (6.1e-8); 4: 3 475 (6.2e-15); 5: 724 (8.0e-7); 6: 564 (1.4e-14). No other K has any rejection. The higher periods are the same ghost cloud seen through iterates of the period-1 and period-2 bifurcations.
- Periodic rows per K index (0, 0.5, 0.9716, pi/2, 2 pi, 10): 12, 24, 24, 32, 1056, 4135. K = 2 pi is not comparable with the lower K (it has many genuine orbits up to period 6, like K = 10), and all rows are elliptic, hyperbolic or parabolic. The earlier D50 build had 1086 rows at K = 2 pi with 28 ghost-class rows.
- Cells (n_c = 40; drawn / qualifying):

| K | T | R | C | S | X |
|---|---|---|---|---|---|
| 0 | 40 / 53248 | 40 | undefined | undefined | undefined |
| 0.5 | 40 / 63165 | 40 / 49 pts | 40 / 794 | 25 short / 25 | 40 |
| 0.9716 | 40 / 34295 | 40 / 37 | 40 / 22669 | 40 / 618 | 40 |
| pi/2 | 40 / 16462 | 40 / 31 | 40 / 41707 | 40 / 333 | 40 |
| 2 pi | 40 / 138 | 0 empty | 40 / 64596 | 5 short / 5 | 40 |
| 10 | 0 empty | 0 empty | 40 / 65536 | 0 empty | 40 |

  Class T is populated at every K except K = 10 ("no qualifying grid cell"). Empty: K = 2 pi R and K = 10 R ("no elliptic periodic point"), K = 10 T and K = 10 S ("no qualifying grid cell"). Short (reason "fewer qualifying grid cells than n_c"): S at K = 0.5 and K = 2 pi. R at K = 2 pi is now empty: the four elliptic points it had under D50 were ghosts.
- Orbits: 883 (870 drawn, 13 controls). Controls: ids 871 to 880 the fixed points of K index 2 to 6 (K = 0.5: hyperbolic 0.6931 and elliptic 0; 0.9716: 0.9496, elliptic; pi/2: 1.1831, elliptic; 2 pi: both hyperbolic 2.0993 and 1.3951; 10: 2.4779 and 2.0634); ids 881 and 882 Newton orbits at K = 10, period 1 (accelerator modes), hyperbolic, exponents 2.2697 and 1.7229; id 883 Y5 (exponent 0.96242). No ghost control.
- S_D 172 members, S_L 92, no cap applied.

### Audit rerun

| Check | Result |
|-------|--------|
| GREEN | pass: 52218 pass, 0 fail, 0 error |
| Aqua, JET, layering | pass (quality 12/12, layering 2/2) |
| Coverage | pass, 98.9 percent |
| Functions under 50 lines, files under 400 | pass (checked below) |
| Tests | changed on purpose, see the test diff above |
| `@fastmath`, dashes, compensated sums | pass |

## Amendment D52 (stage 04 review findings, 2026-10-17)

Tests first (test plan, section "Amendment D52"; new testsets in `ex_inputs.jl`, `ex_periodic.jl`, `ex_panel.jl`, no existing assertion changed), then code. Not run here by instruction: the study panel and `verify/experiments.jl`.

| Finding | Fix |
|---------|-----|
| M1 | `_snap` in `orbits.jl`: a stored orbit point coordinate whose binary64 value is not below 2 pi becomes 0 (identity key rule of E20); iteration of the orbit uses the unsnapped values. Controls are ordered by the 256-bit key `_key` of the orbit, not by stored binary64 values (`_controls` sorts orbit records first) |
| M2 | `collect_orbits`: a duplicate farther than `delta_dup / 5` from a kept point refuses (ArgumentError naming delta_dup, K and the period); closer ones merge as before (also when `check = false`) |
| M3 | manifest keys `julia_version` and `stablerngs_version` (the sidecar kept its `julia_version`); identity is promised for one machine and one Julia version (spec D52) |
| M4 | `_make_orbit` checks the return residual at every point of the orbit and returns `nothing` (orbit not accepted) above rho; the stored residual is the maximum over the points |
| L1 | `place_R0` returns the 256-bit `pb`; the decimal text of a K = 0 R start is `2 pi a / b` at 256 bits |
| L2, L3 | `K >= 0`; `N_scr` a multiple of 10 and `N_scr / 10` a multiple of `tau` (messages name the keys) |
| L4 | `_HASH_CELL` keeps the padding of the cell (trimmed text replaced by SELF); the real preregistration has single blanks, so its hash is unchanged |
| L5 | a solution reduced to a smaller primitive period is re-tested against the D51 bound at that period; rejected ones are counted in `ill_conditioned` (per period) |
| L6 | StableRNGs version in the manifest (`stablerngs_version`) |
| L7 | `_repo_path`: relative with `/` inside the repository, else `<outside repository>/<file name>` |
| L8 | `_search_k`: a refusal is rethrown with the prefix "K index i (K = ...)" |

Portraits (`verify/experiments.jl`): hyperbolic periodic points are no longer drawn (only elliptic circles and parabolic diamonds, plus controls by their start marker); each panel subtitle states the number of hyperbolic periodic points not drawn; the legend lost the hyperbolic entry. The file parses; it was not run.

Results: experiments suite 2632/2632 (29 s); full `Pkg.test` 52253 pass, 0 fail, 0 error, 5 min 09 s; Aqua, JET and layering clean; coverage of `src/Experiments` 98.2 percent (694/707). Longest function 27 lines, longest file 163 lines (`orbits.jl`); no `@fastmath`; no dashes.

Not fixed / caveats: the `delta_dup / 5` duplicate rule is untested on the full study inputs (the panel was not rebuilt); if a legitimate duplicate sits slightly farther than the linearized bound the build refuses with a message naming K index and period. M2's coordinate snap moves a point by at most about 4e-16 (its binary64 distance to 2 pi), well inside rho. `Manifest.toml` for Julia 1.10 still needs `Pkg.resolve()`.

## Amendment D53 (study runner and run store, 2026-10-18 and 2026-10-19)

Branch module/experiments-runner at 0082ef8, not committed. Julia 1.13.2 (1.10 is not installed here), 4 threads for the suite alone, 2 threads for `Pkg.test`, 23 threads for the real-panel smoke run. Tests untouched: `git diff 0082ef8 -- test/` is empty.

### Result

`test/experiments/` on 4 threads: 4279 pass, 7 fail, 3 error (4289); the 2632 earlier assertions all pass (AC28 included, see "Dependencies"). New part: 1647 pass, 7 fail, 3 error. The 10 non-passing results come from four defects of the tests themselves (next section), none from `src/`. With four minimal edits to a scratch copy of the tests (not the repository) the whole new part passes: 1668 of 1668. Full suite (`Pkg.test()`, 2 threads): 53900 pass, 7 fail, 3 error (53910), 4 min 07 s; the 10 non-passing results are the four test defects below, every earlier module passes. quality 12 of 12 (Aqua 11 checks, JET 1), JET clean (24 reports on the first run, all `Nothing` unions of regex matches and HDF5 `getindex`, fixed), layering 2 of 2.

### Tests believed wrong (stop at the point, report)

| Test | Where | Reason | Proposed fix |
|------|-------|--------|--------------|
| D53e, D53f (2 errors) | `ex_run_helpers.jl:411-414` (`ex_sys_x`), used by `ex_iter_states` (`ex_chunks.jl:80`) and `ex_direct` (`ex_chunks.jl:110`) | The fixture arm `wrap` (S_D, wrapped) contains the Y5 control, a `LinearSymplecticMap` (D53b requires that item: `ex_ref_items` takes it). The helper builds `torus_state(x)`, a 4-state, for it; `forward` and `tangent` of the linear map take 2-vectors only (`Systems/linear_symplectic.jl:30-32`), so the helper raises a MethodError whatever the runner stores. A linear map has no torus. The runner stores the Y5 control in wrapped mode in plane coordinates with zero winding counts. | `ex_sys_x`: `item.coordinate_mode == "wrapped" && item.system != "linear"` (and expect zero winding for linear items in D53e). |
| D53l wall cap (2 fails) | `ex_health.jl:151` and `:154` | The injected clock gives 150 s to a group of 1500 orbit-steps, 5e5 times the Float64 budget line of the fixture (200 ns), so the first-chunk cost cap of R5 item 5 pauses arm main after group 1 (and tau after group 3). That is exactly the rule the next testset (`D53l first-chunk cost cap`, 2.01 against 1.99) requires, so groups 2 and 4 are skipped and `ran` is 1 to 4, 9 to 12 and 17, not 1 to 12. The two testsets cannot both hold on the same fixture. | Give the wall-cap fixture a large `first_chunk_cost_factor` (1e9 in the scratch copy), or a clock cost proportional to the steps. |
| D53p input changed (1 error) | `ex_manifest.jl:117-118` and `:126` | `String(orig)` takes the memory of `orig::Vector{UInt8}` and leaves it empty, so `write(path, orig)` restores an empty file and the resumed run fails in `check_prereg` ("found 0" header rows). The 4 assertions before it pass. | `String(copy(orig))`. |
| D53q pilot (5 fails) | `ex_pilot.jl:38` and `:47` (4 times) | The injected clock starts at 1.7e9 s, where binary64 spacing is 2.4e-7 s; the injected group durations are 1.2e-5 to 4.5e-4 s, so any wall time formed as a difference of two readings has a relative error of 1e-4 to 1e-2. Measured 1.50005e-7 against 1.5e-7; the check is `rtol = 1e-9`. No implementation can pass. | `ex_clock(...; start = 0.0)` (passes in the scratch copy). |

### Files

- `src/Experiments/` gains 14 files, 2348 lines (longest file 309 lines `run_record.jl`, longest function 30 lines `run_group!`): `run_block.jl` (A4 block validation, canonical text), `run_grid.jl` (recording grid, deviation frame, exact closed forms, `code_version`), `run_plan.jl` (items, chunks, hashes, steps), `run_backend.jl` (R4, launch groups), `run_record.jl` (recording of one orbit, non-raising monitor, `hot_loop_bytes`), `run_store.jl` (atomic HDF5 chunk files, JSON), `run_health.jl`, `run_campaign.jl` (hooks, start checks), `run_manifest.jl`, `run_resume.jl`, `run_groups.jl`, `run_session.jl` (`run_campaign`), `run_pilot.jl` (`run_pilot`, `rerun_check`), `run_cli.jl` (`cli`). `Experiments.jl`: imports Diagnostics, exports 30 names, includes.
- `bin/chaoslab.jl` (7 lines, 4 of them code): `Chaoslab.Experiments.cli(ARGS)` and exit.
- `Project.toml`: `HDF5` 0.17 (uuid f67ccb44), `JSON` 1 (682c06a0), `Dates` 1 (stdlib) in `[deps]` and `[compat]`. `Manifest-v1.13.toml` resolved (HDF5 0.17.4, JSON 1.11.0, and their jlls). `Manifest.toml` (Julia 1.10) is not resolved here (no 1.10); `verify/Manifest.toml` needs a `Pkg.resolve()` for the same two packages (Chaoslab is a path source there).
- `schema/run-store.md`: one root attribute added, `setup_check` (needed so that the manifest can be rebuilt from the files, R8 item 2).
- `bench/experiments.jl`: 6 entries added.
- `src/Experiments/CONTEXT.md`: files and determinism.
- No lower module changed. P18 needed nothing: `propagate` takes `W0` (the first deviation vector) and the per-orbit non-raising report is built in Experiments from the pure residual functions of Diagnostics.

### Dependencies

- HDF5.jl is required by R6. JSON.jl 1.x (JuliaIO, maintained, one version line since 2025) over JSON3: JSON3 brings StructTypes and a typed-object API this module does not use; the module needs `JSON.json` of nested dictionaries and `JSON.parse(...; dicttype = Dict{String,Any})`, two calls. A hand-written writer was the third option and would need its own tests. `Dates` only for UTC strings.
- The old test AC28 (`ex_inputs.jl:148`) splits an import line at commas, so `import ..Tangent: a, b` counts `b` as a package. The new code therefore defines the monitor's methods as `Tangent.reducer_step(...)` and has no multi-name import line.

### Design

- Recording: the verified estimator reducers (`Benettin`, `QRSpectrum(alg)`, `SVDSpectrum`) run on the shared per-orbit body `Tangent._orbit_body` in segments between grid steps, with the reducer chain `(estimator, RunMonitor)` and the carried state; the value stored at a grid step is `reducer_final(...).lambda` of the same state sequence as a direct `propagate` to that step, which is why D53f holds with 0 difference (all nine methods, three precisions, wrapped mode). L2 does the same with `Lyapunov._shadow_body` and a carried end record. The two internals are reached by qualified name; exporting them would be the clean next step.
- Diagnostics (P18, P20): `RunMonitor` keeps the running maxima of the det and symplecticity residuals (`eps_det`, `eps_symp` of Diagnostics) and the largest `|M|`; orthogonality, `G_max` come from the estimator state; sum and pairing are evaluated at the grid steps. Thresholds in force are those of Diagnostics (`det_threshold`, `symp_threshold`, `sum_threshold`, `orthogonality_threshold`) at the largest `|M|` and `G` of the item; the first exceedance is the first grid step above them. Nothing raises: a health failure of the stream is recorded as a message, the remaining rows are NaN and the chunk fails (R5 item 1a).
- Plan: items in the order of R1 item 4 (precision, first arm, orbit); a chunk is one (run type, cell); the chunk hash adds the values that decide the status (the three health multipliers, `fail_fraction`, `control_bound`, `drift_bound`, `setup_check_tol`, `setup_check_steps`, `control_methods`, h1 flag). Closed-form K and extra starts are converted to 256-bit decimal texts (`_decimal_text`); decimals are read exactly in directions (`"0.3 + pi/5"` uses the decimal 0.3, not the binary64).
- Launch groups: consecutive chunks of one run type up to the capacity; `order = :reverse` reverses groups and the chunks inside them.
- Store: file written under `chunk-NNNNN.h5.tmp` with `complete` false, orbit groups, closed, read back and compared bit for bit, `complete` set, renamed. A replaced file goes to `aside/`, stale files to `stale/`, temporary files are listed in the manifest and removed. The manifest is written after every chunk (atomic, `.partial` then rename) and holds the session list; an initial manifest is written at start.
- Timing: the group clock events are read after the input and cap checks; the default clock is epoch time at the resolution of `time_ns` (`time()` ticks in milliseconds on this system). Before its first group a run type is compiled outside the timed region (`warm_up!`, at most 1000 steps), and a first group of an unmeasured arm of at most 50 000 orbit-steps on the real clock is timed five times and the shortest taken. Tasks are spawned only for 4000 orbit-steps each (below that the group runs on one task; results do not depend on it, I12).

### Deviations from the spec

1. GPU: the recorder state is not isbits, so every chunk runs on the CPU. The backend table of R4 is resolved and recorded as specified; a pair that resolves to the device is executed on the CPU and the chunk records the existing reason "reducer is CPU-only". No device path exists in this build (known limit; see the estimate below).
2. Sum and pairing running maxima are taken at the recorded grid steps, not at every handoff; for d = k = 2 pairing equals sum, and both use the sum threshold (the diagnostics spec leaves the pairing allowance to the caller).
3. The Y5 control in wrapped mode runs in its plane coordinates and stores zero winding counts (a linear map has no torus; see the table above).
4. A chunk that fails twice marks every orbit of the chunk `excluded`, not only the orbits named in the reasons. A control over its bound (R5 item 1c) writes the chunk as failed at attempt 1 without a rerun (a code defect, not a transient fault) and is not counted among the chunks that failed twice.
5. `setup_check` is a root attribute of the file (schema addition above); the manifest value is per chunk: `ok`, `failed` (one orbit differs), `no reference` (no Float64 item of the plan with the same orbit, method, tau, direction and mode whose grid holds the steps), `not applicable`. The Float64 reference is computed by direct iteration of the map, so the check does not depend on the order in which chunks run.
6. `code_version` hashes the working-tree files of `src/` and the lock in use (`Manifest-v<major>.<minor>.toml`, else `Manifest.toml`), CRLF read as LF; the spec says "committed content" and the tests change files without committing (test plan, choice 5).
7. The budget line of an arm in the manifest and in the cost cap is that of its first precision (all arms of the real block have one), CPU basis (the only basis that runs). `rerun_attempts` above 1 is read as 1.
8. A pilot writes the report and a store of the same layout under `<out>/store`; its parallel efficiency is the summed per-item time over wall time times workers, measured on the monotonic clock whatever the injected clock says.
9. Defaults of the command line `run`: with no `--store` the store is `<repo>/data/<slug>` (`default_store`), with no `--manifest` the manifest is `workspaces/study/stages/03-run/output/<slug>-manifest.json`.

### Performance (Julia 1.13.2, this machine, 4 threads, idle)

Recording loop of one orbit, steady state, `RunMonitor` included, ns per orbit-step (Float64 at N = 2e5, Double64 at 2e4, BigFloat 256 bits at 2e3; minimum of 3 and 2 runs):

| Method | Float64 | Double64 | BigFloat |
|--------|---------|----------|----------|
| L1 | 51.8 | 2432 | 20227 |
| L2 | 36.0 | 1548 | 13677 |
| L3 | 83.4 | 3114 | 36839 |
| L4 | 84.8 | 3123 | 36696 |
| L4b | 91.0 | 3324 | 40721 |
| L5 | 90.9 | 3156 | 43113 |
| L5b | 104.6 | 3156 | 40851 |
| L6 | 179.6 | n/a | n/a |
| L9 | 93.3 | 3333 | 46464 |
| mean | 90.6 | 2898 | 34823 |

Bench (`bench/experiments.jl`, BenchmarkTools, 4 threads, N = 1e5, Float64, minimum of 5): L1 4.69 ms (46.9 ns per orbit-step), L2 3.28 ms (32.8), L5 9.06 ms (90.6), L6 16.9 ms (169), L9 8.79 ms (87.9), all 0 B; one small campaign (panel of the bench fixture, one S_L orbit per cell, L5, 1e4 steps, files and manifest written) 210 ms, 7.3 MiB.

Against `[run.budget]` of Amendment A4 (rate, line with safety, first-chunk cap = 2 times the line): Float64 33.09, 66.18, 132.4; Double64 1553.8, 4661, 9323; BigFloat 18200, 54600, 109200. The measured means are 2.7, 1.9 and 1.9 times the rates and below every cap; only L6 at Float64 (179.6) is above the Float64 cap, and no arm of the real block starts with L6. The rates of the preregistration come from the screening kernel (42.4 ns per step measured in the panel notes), the runner also keeps the monitor (one `tangent` evaluation and the det and symplecticity residuals at every handoff) and the grid records. Allocation: 0 bytes in the recording loop of all nine methods after warm-up (D53v). Writing a chunk file with read-back costs about 1.4 ms per orbit (60 orbits: 85 ms); a manifest write about 9 ms at 1387 chunks.

Real plan (preregistration with Amendment A4, real panel of 883 orbits): 15 741 unique items, 1387 chunks, 3.067e11 orbit-steps after the merge (3.170e11 before, saving 1.032e10), plan time 4.9 s. Real-panel smoke run with horizons scaled down (3200 to 9600 steps) and the caps relaxed, 23 threads: 1100 chunks in 85 s of which 50 s writing and 14 s manifest, i.e. about 2 minutes for the 1387 chunks.

Campaign estimate with the measured rates (unique items of the real plan): Float64 2.933e11 steps, 7.2 thread-hours; Double64 1.308e10 steps, 11.3 thread-hours; BigFloat 2.76e8 steps, 2.8 thread-hours; total 21.3 thread-hours (the budget lines of A4 without safety give 9.8). On 23 threads: 0.93 h with ideal scaling, 1.03 h at 90 percent efficiency, well inside the wall cap of 12 h. With the CPU alone that is the figure; a device path for the Float64 methods other than L6 (2.656e11 steps) would take about 0.18 h at the recorded 2.468 ns per step and leave Double64, BigFloat and L6 (15.5 thread-hours, 0.67 h on 23 threads), about 0.7 h in all, but it does not exist yet.

Suite time (4 threads): new part 101.6 s (budget 120 s, spec), the whole slug 2 min 14 s (D35 ceiling 3 min); under `Pkg.test` with 2 threads the slug takes 2 min 07 s and the whole suite 4 min 07 s (budget 30 min).

### Coverage

`Pkg.test(coverage=true)`, 2 threads: `src/Experiments` 2047 of 2101 executable lines, 97.4 percent (new files `run_*`: 1353 of 1394, 97.1 percent); counted from the `.cov` files (merged over the two processes of the entry-point test), then deleted. Caution for that mode: instrumented code is several times slower, the fixtures' budget lines are 100 to 200 ns, and the first-chunk cost cap then pauses arms in the real-clock tests (134 failures in that run, 37 of them show the status `paused arms` directly and the rest follow from incomplete stores); the totals above come from the plain run. The same risk exists on a slow CI machine; it is a property of the fixture budgets, not of the module.

### Known limits

- No device path (deviation 1); the GPU fallback machinery of tangent-core is not used by the recorder.
- The cost cap on the real clock is noisy for groups of a few thousand orbit-steps and is only informative from groups of about 1e5 orbit-steps; the injected clock path is exact.
- `rerun_check` reruns on the CPU only (no device chunk exists in this build).
- Manifest rewrite is O(chunks) per chunk (9 ms at 1387, 12 s in all); a journal would be needed beyond about 1e5 chunks.
- The panel file is hashed (6 MB) before every launch group; a stat signature would save about 40 ms per group.
- 1.10: not run; `Manifest.toml` and `verify/Manifest.toml` need `Pkg.resolve()`.

### Audit (D53)

| Check | Result | Evidence |
|-------|--------|----------|
| GREEN, including every earlier module | FAIL (10 results, 4 test defects) | see "Tests believed wrong"; scratch copy with four edits: new part 1668 of 1668; all earlier modules pass |
| Quality | pass | Aqua 11 of 11, JET no report, layering 2 of 2 and AC28 2 of 2; coverage 97.4 percent (97.1 percent for the new files) |
| Functions under 50 lines, files under 400 | pass | longest function 30 lines (`run_group!`), longest file 309 lines (scratch lint over `src/Experiments/run_*.jl`) |
| Tests untouched | pass | `git diff 0082ef8 -- test/` empty |
| No `@fastmath`; compensated logs | pass | no match; the sums are the `CompensatedSum` of the estimators, the runner adds none |
| No em or en dashes in new files | pass | checked on `src/Experiments/run_*.jl`, `bin/chaoslab.jl`, `bench/experiments.jl`, this section |
| Spec numbers | pass | grid 201 and 221 points; plan 52, 66 and 223 items; study caps plan 1100 x 9 x 1e7 |
