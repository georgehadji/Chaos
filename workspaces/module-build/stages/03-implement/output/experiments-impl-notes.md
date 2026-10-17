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
