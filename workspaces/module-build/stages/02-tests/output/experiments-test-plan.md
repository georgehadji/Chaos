# Test plan: experiments

Spec: `workspaces/module-build/stages/01-spec/output/experiments-spec.md` (29 acceptance criteria, owner decisions P1 to P9 of 2026-10-17), with the preregistration `p1-le-benchmark-preregistration.md`, Amendment A1 (SELF hash rule, machine-readable panel block). Tests: `test/experiments/runtests.jl` including `ex_helpers.jl`, `ex_inputs.jl` (API, AC1, AC2, AC28), `ex_geometry.jl` (AC3, AC4), `ex_periodic.jl` (AC5 to AC9), `ex_screening.jl` (AC11, AC12, AC25, AC26), `ex_predicates.jl` (AC10), `ex_draw.jl` (AC14 to AC18), `ex_subsets.jl` (AC21), `ex_panel.jl` (build-level AC12, 13, 15 to 18, 21 to 24, 29, invariants I1 to I6, I10), `ex_determinism.jl` (AC19, AC20). Helpers are prefixed `ex_` / `EX_`; nothing in `src/`; no top-level call into `Chaoslab.Experiments` (`EX.name` resolves the module when called, so a missing module errors the testsets only).

Conventions. Fixture inputs are the suite's own (`ex_block_default`): K in {0, 1.5, 5}, grid side 8 (build) or 16/32 (unit), N_scr = 2000, w = 100, n_c = 3, q_max = 2, rho = 1e-10, delta_dup = 1e-8, synthetic thresholds. No study value is copied into a test; the real preregistration is read by one testset (AC2a) for the hash rule and the block. Fixture preregistrations are generated with the true SELF hash computed by the platform `sha256sum` (fallback `certutil`, `shasum`), independent of the module (O3). Seeds: `EX_BASE = 20261012 + 27000`, `ex_seed(ac, sub) = EX_BASE + 50 ac + sub` (sub < 50); run seed in builds is the default `seed` of `numerics-defaults.toml`. Every random draw uses StableRNGs with those seeds. Regression values: none from memory; analytic oracles are evaluated at 1024 bits (family A, B, C of E6, fixed points, K = 0 screening, Y5), O1 uses the benchmarks.md row "Island centre and unstable 4/9 orbit point" (ContopoulosHarsoula2008 Sec. 2 p.4), K_c = 0.971635406 the row "Golden-mean KAM torus breakup". ChaosTools comparison is stage 04 (see mapping).

## New API (all in `Chaoslab.Experiments`, exported through `Chaoslab`)

Any 2-indexable (tuple, SVector) is accepted where a point is passed. Names introduced by the tests, checked for collisions against the existing exports (none collide; `Experiments` has no existing names).

| Name | Contract tested |
|------|-----------------|
| `panel_inputs(dict)` | validates the E0 inputs (string keys, as in the TOML block), returns an object with the key names as fields (`K::Vector{Float64}`, `n_c`, `grid_side`, `N_scr`, `alpha`, `mu`, `w`, `theta_T`, `rho`, ...); `ArgumentError` naming the key. Optional keys `tau`, `seed`, `bigfloat_bits` (defaults from `numerics-defaults.toml`) |
| `prereg_hash(text)`, `check_prereg(path)`, `read_prereg(path)` | E18; `check_prereg` returns `(header_sha256, computed_sha256, amendments)` and throws `ArgumentError` on any mismatch, missing or duplicated row, non-hex value; `read_prereg` adds `.inputs` |
| `build_panel(prereg_path; out, chunk=nothing, ntasks=nothing, backend=:auto)` | writes the panel TOML file and `out * ".provenance.toml"` (keys `backend`, `backend_fallback_reason`); nothing written on a bad input |
| `portrait_content(panel_or_path; n_plot)` | one NamedTuple per K: `k_index, K, starts, iterates, periodic, notes` |
| `grid_centres(N_g)` | `Vector` of centres, index = id + 1, id = j N + i |
| `torus_delta(s)`, `torus_distance(u, v)` | E8 |
| `periodic_set(K; q_max, N_g, rho, newton_max, delta_dup, bits=256)` | `(orbits, lines, stats, min_separation)`; orbit fields `q, m, j, points (SVector{2,BigFloat}, representative first), trace, class (:elliptic/:parabolic/:hyperbolic), exponent, residual, grid_id`; `lines` = (a, b) pairs at K = 0 |
| `collect_orbits(K, sols; rho, delta_dup, bits=256)` | E9 on solutions `(q, z, id)` (id -1 = exact fixed point); returns the orbit vector; `ArgumentError` naming `delta_dup` when distinct points are closer than 10 delta_dup |
| `classify_trace(t, rho)`, `orbit_exponent(t, q)`, `accept_residual(res, rho)`, `newton_correction(J, F)` | E4, E3 acceptance, Newton step (`nothing` at tr J = 2) |
| `periodic_distance(z, pset; elliptic_only=false)` | E8 distance to the periodic set (lines at K = 0); parabolic counts as elliptic |
| `screen_buffer(N, w)`, `screen_point!(buf, K, x, p, N, w, alpha, tau)`, `screen_point(K, x, p; N_scr, w, alpha, tau)`, `screen_bytes(N, w)` | L5 screening; fields `lambda, blocks, S_tenth, S_N` |
| `block_median`, `stuck_length(blocks, f_S, w)`, `pred_T(lam, S10, SN, d; theta_T, g_T, r_T)`, `pred_C(lam, d; theta_C, r_C)`, `pred_S(lam, blocks; theta_S, w, f_S, L_S)` | E11 predicates |
| `cell_stream_seed(run_seed, i, class)`, `cell_rng(...)`, `orbit_seed(run_seed, id)` | E15, E16 (class as String or Symbol) |
| `draw_candidates(rng, cands, n_c)`, `place_X(rng, lo, hi)`, `place_R(rng, pts, delta_R)`, `place_R0(rng, q_max)` | E14, E12, E13 (`place_*` return NamedTuples `x, p, r, phi` / `index, phi, x, p` / `a, b, x, p`) |
| `y5_control(mu)` | `(sys::LinearSymplecticMap, x0, exponent)` |
| `select_subsets(members; s_D, cap_D, s_L, cap_L)` | E17 on `(id, cell, rank, control)`; returns `SD, SL, removed_D, removed_L` |
| `git_state(dir)`, `package_version()` | `(sha, dirty)`; version string of `Project.toml` |

Panel file layout fixed by the tests (the spec leaves layout to the build): `[manifest]` with `format_version, package_version, git_sha, git_dirty, stream_rule, draw_order, K_binary64, min_separation`, sub-tables `prereg (path, header_sha256, computed_sha256, amendments)`, `inputs.<key> (value, source)` including `seed`, `bigfloat_bits`, `tau` with source "default", `newton`, `subsets (capped_D, capped_L, removed_D, removed_L)`, `totals (orbits, drawn, controls)`, arrays of tables `cells (k_index, class, qualifying, drawn, status, reason)`, `periodic (id, k_index, q, m, j, class, trace, exponent, residual, points = [[x_dec, p_dec], ...])`, array `controls`; then one `[[orbit]]` per orbit with `id, class, control, control_kind (fixed/newton/y5), system, x0, p0, x0_dec, p0_dec, run_seed, orbit_seed, in_SD, in_SL` plus `k_index, K_text, K` (not for Y5), `cell_seed, draw_index` (drawn only), `a, b` (R at K = 0), `source_orbit, source_point, phi` (R at K > 0), `r, phi` (X), `period, wind_m, wind_j, stability, exponent` (controls), `mu` (Y5). 64-bit seeds are decimal strings. Starts are stored as formed (unwrapped sum for displaced R and X starts).

## Criterion to test

| AC | Test (file: testset) | Sample | Pass |
|----|---------------------|--------|------|
| 1 | inputs: "AC1 input validation" | 35 malformed dictionaries plus 3 malformed preregistration files through `build_panel` | each `ArgumentError` names the key; no panel or sidecar file |
| 2 | inputs: "AC2 preregistration hash" | (a) real frozen file (header = recomputed = independent, block parsed); (b) 20 seeded one-byte changes (3 forced: slug row, first and last hash digit); (c) no row, two rows, 63 digits, non-hex; (d) CRLF copy; (e) module value = `sha256sum` value | exact |
| 3 | geometry: "AC3 grid" | N_g in {1, 2, 3, 5, 32, 256} | ids, order, 0 ulp against 1024-bit value, neighbour distance within 2 ulp of 2 pi |
| 4 | geometry: "AC4 torus distance" | known pair, 2000 symmetric pairs, 10^4 triples, K = 0 lines, `torus_delta` | see ambiguity A4 for the first tolerance |
| 5 | periodic: "AC5 fixed points" | 8 K values, `q_max = 1`, N_g = 16 | class, trace to 64 eps_BF, exponent to 1e-60, listed once |
| 6 | periodic: "AC6 parabolic band ..." | synthetic traces and their negatives, `orbit_exponent`, trace and class at every point of every K = 5 orbit | exact classes |
| 7 | periodic: "AC7 period-2 families" | K in {0.5, 1, 1.5, 2, 2.5, 3, 5, 7, 10}, `q_max = 2`, N_g = 32 | families A, B, C vs 1024-bit closed forms, traces to 1e-16, O1 within 1e-4 unit, exponents |
| 8 | periodic: "AC8 no spurious periodic point" | all orbits at K = 1.5 and 5 recomputed at 1024 bits; synthetic acceptance and Newton cases | residual <= rho (1 + 2^-200), windings, primitive period, det bound at 256 bits |
| 9 | periodic: "AC9 deduplication" | (a) to (e) on synthetic solutions | as stated |
| 10 | predicates: "AC10 predicates" | 46 assertions, each threshold at the value and one step either side, median odd and even, e_med <= 0, parabolic as elliptic | hand computed |
| 11 | screening: "AC11 windowed estimate" | 51 seeded grid cells over K in {0, 1.5, 5} | blocks N/w, mean = lambda to 1e-12 |
| 12 | screening: "AC12 K = 0 ..." ; panel: "AC12 K = 0 cells" | all 64 cells; thresholds either side of the closed value in a K = 0 build | closed form, T cell full or empty (the [full] part with study inputs: stage 04) |
| 13 | panel: "AC13 K = 0 resonant class" | R starts of a K = 0 build | gcd 1, b <= q_max, p0 nearest, return within b ulp; C, S, X undefined |
| 14 | draw: "AC14 X starts" | 10^4 draws | log10 r range, KS <= 1.95 / 100 for log r and phi |
| 15 | draw: "AC15 R placement"; panel: "AC15 R starts at K > 0" | 20 + 30 seeded placements; build starts rebuilt from (point, phi) bit for bit; K = 10, q_max = 1 build for the empty cell | as stated |
| 16 | draw: "AC16 draw"; panel: "AC16 draw in the build" | reference partial Fisher-Yates, 2000 seeds with M = 1000, n_c = 10 (5 sd), short and empty, S/C disjoint, starts are screening points | as stated |
| 17 | draw: "AC17 streams and seeds"; panel: "AC17 seeds in the file" | 500 distinct stream seeds; E15 re-derived with the platform hash; per-orbit seed keys `initial_deviations` | see ambiguity A5 |
| 18 | draw: "AC18 Y5 control"; panel: "AC18 controls" | Y5 eigenvalues and exponent; order of controls, cut at n_ctrl_max (2, 4, 7), Newton controls ordered and checked at 1024 bits | as stated |
| 19 | determinism: "AC19 determinism" | rerun, `ntasks` 1, 2, 4, `chunk` 16 and 4096, forbidden keys | byte-identical, SHA-256 equal |
| 20 | determinism: "AC20 backend" | cpu, gpu and auto with a fake device present and absent | byte-identical panel; sidecar backend cpu, empty reason |
| 21 | subsets: "AC21 subsets (synthetic)"; panel: "AC21 subsets in the build" | the 202/200 population, ladder cap, short cells; capped build against the reference procedure | removed ids as stated |
| 22 | panel: "AC22 counts and file" | base build | counts, fields, round trip of decimals, inputs with source |
| 23 | panel: "AC23 seed fields and recorded reasons" | base build | seeds present; every cell populated or with a reason |
| 24 | panel: "AC24 portrait content" | `portrait_content` with `n_plot = 20` | markers, symbols, ranges, notes, no screening names |
| 25 | screening: "AC25 screening against L1" | 200 seeded cells, K in {0, 1.5, 5}, same orbit, W0 = (cos alpha, sin alpha) | within (B_S + B_S) / N_scr (S3 as amended); the K = 10 line vs L5 at 10 N_scr is printed, not asserted |
| 26 | screening: "AC26 zero allocation" | `@allocated` after warm-up; `screen_bytes` | 0 bytes; 48 + 16 + 8 N/w |
| 27 | stage 04 only | ChaosTools 3.5.4 needs its own environment (D27) | unit proxy: AC7 and AC5 give the analytic orbits and classes the oracle is compared with |
| 28 | inputs: "AC28 layering of src/Experiments" | scans `src/Experiments/*.jl` imports (only Foundation, Systems, Integrators, Tangent, Lyapunov, Diagnostics and declared `Project.toml` deps, so `SHA` must be declared) and that nothing else imports Experiments; the existing `layering.jl` checks the layer order | pass |
| 29 | panel: "AC29 provenance" | scratch git repo (head, dirty after edit); manifest git sha and dirty flag against `git`; package version against `Project.toml` | exact |

Invariants: I1 to I6 and I10 in "panel structure and invariants"; I7 is AC8; I8 is AC19; I9 is AC17 (cell stream independent of other cells); I4 is AC19 (forbidden keys) and AC24 (no screening names).

## Deviations and ambiguities (owner decisions needed)

A1. AC2(b) says "including inside the header row". Bytes of the header value's parenthetical are replaced by SELF before hashing, so changing them cannot be detected by the rule. The test forces changes in the slug row and in the digits, and excludes the parenthetical.
A2. AC12 "8 eps_53 relative" for the K = 0 closed form is below the rounding of the sum of 2000 logarithms (about N eps / S). The test uses an absolute 8 eps_53 on lambda_scr, the scale of the Householder suite (S(n) to 8 n eps).
A3. AC9(c) "at distance 10 delta_dup kept" sits on the refusal edge ("below 10 delta_dup"). The test uses 10.5 delta_dup (kept) and 5 delta_dup (refused).
A4. AC4 first bullet "4 eps_53 relative" is not reachable because the Float64 value of 2 pi - 0.1 is already off by 4.4e-16 (error 1e-15 measured on a reference). The test uses 4 eps_53 times 2 pi absolute.
A5. AC17 asks for a stored list of 20 uniforms. A list of generator output is not a benchmark row, so the test re-derives the stream from the E15 rule with the platform hash and compares the first 20 uniforms. Stage 03 may record a list in the 04 verify report.
A6. AC19 asks for runs at 1, 4 and all threads. The package suite varies `ntasks` (concurrency cap of `build_panel`) in one process; real `julia -t` runs go to stage 04.
A7. The manifest records the preregistration path, so panels from different paths differ; the determinism tests reuse one path. The spec says the file is a function of the preregistration bytes only; the owner may prefer a relative or no path.
A8. E13 does not fix the order of the reduced fractions (K = 0) or of the elliptic points (K > 0). Tests assume (b, a) ascending and the order of the supplied list.
A9. AC1 "q_max beyond the supported range": the range is not stated; the test uses 10^6. AC1 "a cap below the number of controls" is tested as cap_D < n_ctrl_max + 1.
A10. AC7 completeness depends on the seeds: the Newton search from a 32 x 32 grid must find every family C root at K = 5, 7, 10 with q_max = 2. If it misses one, N_g in the test is raised, not the tolerance.
A11. Starts of displaced R and X orbits are tested as stored unwrapped; the module may wrap instead if the owner prefers (the test then compares the wrapped value).
A12. The portrait test defines the first plotted point as the start itself.
A13. AC25 uses the S3 bound with the n^2 term as the spec states (loose at N_scr = 2000, K = 5: about 1e-7 per exponent).

## Time

Measured against throwaway stubs of the pure-function parts (scratchpad, not in the repo): unit parts about 30 s including compilation; AC7 about 12 s with a slow stub Newton search; the 9 builds of the panel tests are the unknown (fixture: 64 grid cells, 2 non-zero K, q_max = 2; about 12 builds including determinism and backend variants). Estimate for the real module: under 3 minutes (D35). Cut order if over: AC7 K list to {0.5, 2, 5, 10}, chunk 4096 build, `ntasks` 2.

## Not in the package suite (stage 04, `verify/`)

AC2(a) with the study inputs is present here (hash), the [full] runs of AC12 (T cell at K = 0), AC21 (full population) and the full-panel counts; AC27 (ChaosTools orders 1 to 3, sub-grid 64 x 64, class by automatic differentiation); O7 (separate-process re-evaluation); thread-count runs (`julia -t 1, 4, N`); the portrait rendering (verify script, CairoMakie).

## RED log

Command (as the earlier plans, suite alone; the module does not exist yet): `julia --project -t 4 -e 'using Test, Chaoslab; @testset "experiments" begin include("test/experiments/runtests.jl") end'` from the repository root. Result: 110 reported results, 0 Pass, 27 Fail, 83 Error, exit 1, 5 s. All 83 errors are `UndefVarError: Experiments not defined in Chaoslab`; the 27 Fail are 24 `@test_throws ArgumentError` that received that `UndefVarError`, `isdir(src/Experiments)` (AC28) and the two export checks. No parse error, no error inside a helper. (Counts are per assertion reached; most testsets stop at their first call into the module, so the full count appears at stage 03.) Static size: 51 testset blocks (plus 35 generated in AC1), 333 assertion lines.

Helper and test-logic checks done against throwaway stubs (not kept): the pure functions (grid, distance, draw, placement, predicates, subsets, streams, Y5, screening), a stub Newton search and `collect_orbits` pass AC3 to AC12, AC14 to AC18 (except build tests), AC21 (synthetic), AC25, AC26 with 1400 passing assertions after fixing three test errors found that way (a Newton check with det J != 1, the AC4 tolerance, the Y5 start direction). The build-level tests (`ex_panel.jl`, `ex_determinism.jl`, AC1, AC2, AC28) have not been run against any implementation.

Full suite: `julia --project -t 2 -e 'using Pkg; Pkg.test()'`: 49731 results: 49621 Pass, 27 Fail, 83 Error; the 27 Fail and 83 Error are exactly those of the experiments suite (110 results), so every other suite (quality, layering, diagnostics, lyap-*, precision, standard-map, tangent-core) passes with 0 Fail and 0 Error. Exit 1 from the experiments RED only. Log: scratchpad `red-full.log`.

## Checkpoint (2026-10-17)

Owner approved the tests and the gap choices above. The manifest records the preregistration path relative to the repository root, so the same file gives the same panel wherever the checkout lives (spec I8).
