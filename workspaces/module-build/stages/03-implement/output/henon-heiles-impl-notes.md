# Implementation notes: henon-heiles (catalog row F3; spec D116)

Code (new): `src/Systems/henon_heiles.jl` (100 lines), `src/Integrators/section.jl` (79), `src/Integrators/escape.jl` (49), `bench/henon-heiles.jl` (42). Edited: `src/Systems/flow_models.jl` (the `HenonHeilesFlow` docstring and comments: F3, no "T2 test-only"), `src/Systems/Systems.jl` and `src/Integrators/Integrators.jl` (exports, includes), `src/Systems/CONTEXT.md` (79 lines), `src/Integrators/CONTEXT.md` (37), `docs/src/api/systems.md` and `integrators.md` (one sentence each; the names appear through `@autodocs`). Nothing under `test/`, `_config/` or `data/` was touched. Tests: `test/henon-heiles/`, plan `stages/02-tests/output/henon-heiles-test-plan.md` (the API proposal is implemented as written).

## Status

| Run | Result |
|-----|--------|
| New suite, first run on the code (driver `red_hh.jl`, 1 thread, Julia 1.13.2) | 6899 pass, 1 fail, 0 error, 0 broken, 3 min 45 s. The failure: AC17 `escape_sector` allocated 160 bytes (Double64: `mod` and `T(pi)` allocate) |
| New suite after the fix (`escape_sector` by dot products, no `mod`, no `pi`) | 6900 pass, 0 fail, 0 error, 0 broken, 2 min 5 s (three runs with coverage restricted to one new file each) |
| Aqua, JET (`test/quality.jl`), layering (`test/layering.jl`) and the integrator-catalog suite (`test/integrator-catalog/runtests.jl`), one process after the fix | 7334 pass, 0 fail, 0 error, 2 broken (the two extension skips of the catalog suite), 8 min 24 s |
| Coverage, `--code-coverage=user`, new suite alone | 109 of 110 executable lines of the three new files, 99.1 percent (`henon_heiles.jl` 43/43, `section.jl` 45/46: the unreachable safeguard return after 200 iterations of `_solve_substep`, `escape.jl` 21/21); `--code-coverage=user` over the whole suite did not finish in 75 minutes, so each file was measured by `--code-coverage=@<file>`; the `.cov` files were deleted |
| Whole package suite (`Pkg.test()`, main session, 2026-11-11) | 1426407 pass, 0 fail, 0 error, 3 broken (the two catalog extension skips, one weighted-birkhoff), 59 min 30 s; henon-heiles 6900 |
| Not run | `verify/gpu.jl`, `HH_FULL_RUN` sizes (stage 04), `--check-bounds=yes` |

## Design

- Systems (`henon_heiles.jl`). Generic in the element type; `critical_points`, `escape_energy` and `d3_elements` are built in `T` (`sqrt(T(3))/2`, `T(1)/6`; correctly rounded). `section_interval` bisects `u(q2) = E` to the resolution of `T` (same formula as Equation 7), closed form `(-1/2, 1)` at `E = 1/6`. `energy_start` forms `p1 = sqrt(2 (E - p2^2/2 - u(q2)))` with the operation order of Equation 8. `random_energy_start` is the rejection loop of Algorithm 2 over `energy_start`. `energy_grid` visits cell centres, `i` outer, `j` inner.
- One root solver (`_solve_substep`, `section.jl`) serves `locate_crossing` (`q1`, derivative `p1`) and `locate_escape` (`r - R`, derivative `q.p/r`): Newton from `tau/2`, a step outside the current bracket is replaced by the midpoint, the bracket is kept from the signs. It accepts at `abs(g) <= 4 eps` (`4 eps R` for the radius), so the tests' bound of 16 eps holds with margin; a collapsed bracket returns its upper end. The state returned is always `advance(m, flow, z, s)` of the returned `s` (AC7b asserts equality). Newton needs 2 to 4 evaluations, bisection only near grazing.
- `section_map`: the grid states are the plain `advance` orbit (not re-based at a crossing), crossing `i` has `t = k tau + s`. With `jacobian`, the frame `Phi` is advanced across each grid step by `tangent_advance`; at a crossing the leg frame is `Phi` through the sub-step `s`, the 4 x 4 return map `DP = (I - f e1' / p1) Phi` (Equation 12) is applied to the two energy-surface displacements `(0, 1, -g2/p1, 0)` and `(0, 0, -p2/p1, 1)` of the previous section point, and the new leg starts from `tangent_advance(zc, I, tau - s)`.
- Escape: `escape_orbit` declares escape at the first step with `r > R_esc`, located with `locate_escape`; `escape_sector` takes the largest dot product of `q` with the unit exit vectors (the nearest centre angle, Equation 14), with no trigonometry.

## Deviations from the spec and the test plan

1. `energy_start` also returns `:outside` for `abs(q2) > 1`. The spec text is only the energy condition; for `E <= 1/6` a point with `q2 > 1` satisfies it (the escape channel beyond `u = E`) but is not in the origin component, and Algorithm 2 requires the origin component. One rule (energy and `abs(q2) <= 1`) covers `E <= 1/6` (origin component) and `E > 1/6` (`r <= 1`).
2. `section_interval` and `energy_grid` throw `ArgumentError` unless `0 < E <= 1/6` (the origin component of the section line is unbounded above it).
3. `section_map` with `jacobian` and a start off the section: `J[1]` is a placeholder identity (there is no section point 0); the legs of the later crossings are correct. A grazing bracket is flagged and not recorded; the leg continues to the next crossing.
4. `section_map` `max_crossings`: the start counts as a crossing; the status is `:buffer_full` as soon as `length(t) == max_crossings`.
5. `escape_sector` uses dot products with the unit exit vectors instead of the polar angle: the same labels (the largest cosine is the nearest angle), allocation free for Double64, bitwise reflection equivariance. Ties go to the lower label.
6. The uniform draws of `random_energy_start` are binary64 (`rand(rng)`) converted to the element type, so a Double64 start has 53 random bits.
7. `escape_orbit` returns `t = 0` and `:escaped` for a start with `r > R_esc`; `R_esc` is converted to the element type.

No compensated summation is needed: the module accumulates no logarithm (times are `T(k) tau + s`, not running sums). No `@fastmath`.

## Performance (`bench/henon-heiles.jl`, benchmark environment, 1 thread, medians; indicative: the coverage run of the suite was on the machine, BenchmarkTools resolution 100 ns)

| Entry | Time | Allocations |
|-------|------|-------------|
| KahanLi8 step / tangent step, Float64 | 0.7 / 1.3 us | 0 |
| KahanLi8 step, Double64 | 2.8 us | 0 |
| `locate_crossing`, `locate_escape` | 2.1 us, 3.5 us | 0 |
| `escape_sector` | 0.1 us | 0 |
| `section_map`, T = 1000 (2 10^4 steps), plain / with Jacobian | 14.5 ms / 40.4 ms | 16 / 23 (the record vectors) |
| `escape_orbit`, E = 1/4, up to T = 200 | 0.15 ms | 0 |
| `random_energy_start`, 10^4 draws | 0.24 ms | 0 |
| `energy_grid`, 100 x 100 | 0.30 ms | 16 (the result) |

The section map costs about 0.73 us per step, the step itself 0.7 us: the crossing solves add little; the Jacobian triples the cost (one tangent step per grid step).

## Known limits

- The crossing is a first-passage test on the grid: two crossings inside one step are missed; the tests use `tau = 1/20`.
- The section Jacobian is that of the numerical return map with the sub-step (not exactly symplectic with the grid map, spec Equation 11); the determinant defect is what AC9 measures.
- `section_map` stores crossings in `Vector`s (host only). The GPU kernel with a bounded record buffer is stage 04 and later (`verify/gpu.jl` not run).
- The portrait (classification at `T_f`) is not a library function (spec P2, test plan U9); the tests build it from `energy_grid` and the existing reducers.
- The root tolerance is `4 eps` absolute in `q1`; for Double64 `effective_eps` is `2^-100`.

## Audit

| Check | Result |
|-------|--------|
| GREEN | new suite: see Status; integrator-catalog suite 7334 pass (with quality and layering), 0 fail; whole package suite 1426407 pass, 0 fail, 3 broken |
| Quality | Aqua and JET clean, layering passes (the combined run above); coverage of the new lines: 109 of 110 executable lines of the three new files, 99.1 percent (`henon_heiles.jl` 43/43, `section.jl` 45/46: the unreachable safeguard return after 200 iterations of `_solve_substep`, `escape.jl` 21/21); `--code-coverage=user` over the whole suite did not finish in 75 minutes, so each file was measured by `--code-coverage=@<file>`; the `.cov` files were deleted |
| Size | longest functions: `section_map` 28 lines, `_solve_substep` 15, `escape_orbit` 17, `random_energy_start` 10; files 100, 79, 49 lines; `src/Systems/CONTEXT.md` 79 lines, `src/Integrators/CONTEXT.md` 37 |
| Tests untouched | `git diff -- test/` empty |
| Numerics | no `@fastmath`; no accumulated logarithm; crossing and escape tolerances are multiples of `effective_eps` |
| Em dashes | none in the new or edited files |
