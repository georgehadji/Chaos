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

## Correction after stage 03 (2026-10-17)

`ex_panel.jl` I3 check: control starts are now keyed by K index before the uniqueness test. Unkeyed, it contradicted AC18, which requires the fixed point (0, 0) as a control at every K > 0 (and the Y5 control starts at (0, 0) too), so no conforming build could pass both.

## Amendment D51 (2026-10-17)

Spec Amendment D51 (decision D51) supersedes D50; preregistration Amendment A3 removes `sn_trace` from the panel block. The D50 tests are removed.

| Criterion | Test | Change |
|-----------|------|--------|
| D51b | ex_periodic.jl "D51 conditioning test" | new. `error_bound(J, F)` equals `||F||_2 ||J - I||_2 / |2 - t|` within 1e-12 relative on 200 seeded random 2 x 2 matrices of determinant 1; `tr J = 2` gives an infinite bound |
| D51a | same testset | new. Same residual 1e-12: `J` with `2 - t = -1e-6` has a bound above delta_dup / 10 and is not `well_conditioned`; `2 - t = -1` is; at 256 bits `_newton_seed` at K = 2 pi (binary64) from (pi/2 + 1e-7, 0) returns code 4 (ill-conditioned, small `|2 - t|`) and at K = 5, period 2, near the family A orbit returns code 0 |
| D51c | ex_panel.jl "D51 build at K = 2 pi (binary64)" | new small-grid proxy (K = ["2pi"], q_max 1, 16 x 16 grid, rho 1e-10): the build completes, `manifest.ill_conditioned` has a period-1 row with count > 0 and smallest `|2 - t|` < 1e-5, `newton.ill_conditioned` covers it, no accepted period-1 point lies within 1e-5 of (pi/2, 0) or (3 pi/2, 0), `min_separation` >= 10 delta_dup. Full study inputs: stage 04 |
| D51d | same testset | new. No `merged` table, no `sn_trace` input, periodic and control classes are elliptic, parabolic or hyperbolic only, no `SN_CLASS` or `MergedPair` in the module; the base build lists only ill-conditioned rows with count > 0 |
| input validation | ex_inputs.jl AC1 | the four D50 `sn_trace` cases are replaced by one: `sn_trace` is now an unknown key |
| fixtures | ex_helpers.jl `ex_block_default` | `sn_trace` removed again (net change against the stage 02 commit: none) |
| D50 tests | ex_periodic.jl, ex_panel.jl | the D50 merge testset and the D50 build proxy are deleted |

New API: `error_bound(J, F)`, `well_conditioned(J, F, delta_dup)` (exported); `_newton_seed(K, x, p, q, rho, nmax, delta_dup)` returns `(code, x, p, den)`; manifest `ill_conditioned` (`k_index, q, count, min_abs_2_minus_t`) and `newton.ill_conditioned`. The I3 correction of the earlier section stays.

## Amendment D52 (2026-10-17)

Stage 04 review findings (decision D52). Tests added, none of the existing assertions changed:

| Criterion | Test | Content |
|-----------|------|---------|
| D52a | ex_periodic.jl "D52 coordinates, duplicates, every point" | a 256-bit coordinate 2 pi - 1e-78 (and +1e-78) in a solution: stored coordinates all below the binary64 2 pi, the first one exactly 0 for the minus sign, the same orbit order for both signs |
| D52a | ex_panel.jl "D52 panel: ..." | in the base build no control or T, C, S start and no periodic point reaches the binary64 2 pi; no orbit coordinate equals it |
| D52b | ex_periodic.jl (same testset) | duplicate at 0.5 delta_dup refuses naming delta_dup and "period 1"; at delta_dup / 20 it is merged silently |
| D52c | same | period-2 orbit at K = 10 (found by a scan) whose first point meets rho and the second does not: not accepted; a genuine orbit still is |
| D52d | ex_inputs.jl "D52 input checks" | K negative, N_scr not a multiple of 10, N_scr / 10 not a multiple of tau: each rejected naming the key; K = 0 stays valid |
| D52d | ex_inputs.jl "D52 SELF rule keeps the padding..." | a header cell with extra blanks: the independent hash keeps them, the module agrees, differs from the fixed-padding hash, and a padding change refuses |
| D52e | ex_panel.jl "D52 panel: ..." | K = 0 R decimal text equals 2 pi a / b at 256 bits; manifest `stablerngs_version` and `julia_version`; a path outside the repository is stored as `<outside repository>/name` (no drive letter, no `\`), one inside as a relative path (`_repo_path`); a refusal names the K index (delta_dup = 0.5 gives "K index 2") |

L5 (conditioning at the primitive period) is covered by construction in `collect_orbits` (only reduced solutions are re-tested; the D51 tests are unchanged). Internal names used by the tests: `_repo_path`.

## Amendment D53 (2026-10-18)

Second scope of the module: the study runner and the run store (spec Amendment D53 with the owner decisions P10 to P25; preregistration Amendment A4). Tests are added; no existing assertion or helper is changed. Stage 02 for the amendment, RED.

### Files

New files in `test/experiments/`, registered at the end of `runtests.jl` in this order:

| File | Holds |
|------|-------|
| `ex_run_helpers.jl` | fixtures (panel, run blocks `:small`, `:full`, `:alt`, `:noref`, `:c20`, `:c40`, `:big`), TOML emitter with a chosen key order, hooks (`ExKill`, fake clock, fault injection), canonical read-back bytes, reference code written from the spec (recording grid at 1024 bits, Cartesian product with merge, backend table, direct estimator calls, c factor at 256 bits) |
| `ex_run_inputs.jl` | API export check, D53a (65 malformed blocks, 4 file-level refusals and a zero-item plan), D53r (real A4 block, source scan) |
| `ex_plan.jl` | D53b, D53c (hashes, 1000 seeded items, code version), D53d, D53g (planning part), D53h |
| `ex_chunks.jl` | deviation frame, D53e, D53f, D53g and D53i (three full-fixture campaigns), D53r (second block), D53s, D53t, D53v |
| `ex_store.jl` | D53m, D53n |
| `ex_resume.jl` | D53o |
| `ex_health.jl` | D53j, D53k, D53l |
| `ex_manifest.jl` | D53p, D53u, D53x (proxy) |
| `ex_pilot.jl` | D53q |

48 testsets (44 top-level at run time, D53a expands to 65 nested cases), about 380 `@test` expressions in the source before loops.

### Conventions

- Seeds: `ex_s53(letter, sub) = 20261012 + 28000 + 50 (letter index) + sub`, letter = the criterion letter; every random draw uses `StableRNG`.
- Fixtures are generated by the helpers, never copied from the study. The panel file is written by `ex_write_panel` (not by `build_panel`), so the suite does not pay for a panel build. Orbits: 12 drawn in cells (k1,T), (k2,C), (k3,S) at K in {0.5, 1.5, 2.5} (ids 1 to 12), three controls (13 elliptic fixed point at K = 1.5, 14 hyperbolic fixed point at K = 2.5, 15 Y5 with mu = (3 + sqrt 5) / 2), S_D = the first two of each cell plus the controls (9 orbits), S_L = orbit 1, two extra starts from the block (ids 16 and 17).
- Run blocks (dictionaries of the `[run]` table, written to a fixture preregistration by `ex_emit_run`, which fixes the textual key order):

| Block | Arms | Items | Chunks | Use |
|-------|------|-------|--------|-----|
| `:small` | main (P, L1, L5, h1_check), tau (S_D, L4, tau 1 and 2), h4f (L5), h4d (Double64 L5, horizon 40) | 52 | 18 (6 launch groups) | health, caps, resume, store, manifest, pilot |
| `:full` | all nine methods, merges (dup, dirs, epsw, tauw), wrapped, h4 Float64 and Double64, BigFloat on orbit 1; horizons 100 and 200 | 223 | 76 | series, direct call, thread independence, Y5, rerun |
| `:alt` | differs from `:small` in every value of R1 to R5 (tau 2, eps 1e-6, direction 0.3 + pi/5, other grid, fractions, bounds, backend, budget, k19, pilot) | 66 | 22 | D53r |
| `:noref` | main (L1) and Double64 on the extra starts without a Float64 counterpart | 17 | 5 | setup check "no reference" |
| `:c40`, `:c20` | 9 methods x 4 cells + 4 extra-start chunks; 5 methods x 4 cells | 143, 75 | 40, 20 | stop rule |
| `:big` | one cell of 40 orbits (panel `:big`, no control) | 40 | 1 | fail_fraction boundary |

- Control bound of the fixtures is 0.05 (study value 1e-8): at horizons of 100 to 200 steps the bias ln(c) / n of the K = 2.5 hyperbolic fixed point at direction pi/4 is 6e-4, so a bound of 1e-8 would stop every fixture campaign (this is P17 at fixture scale, not a defect). Injected control faults add 1.0.
- Launch capacity of all comparisons is 64 orbits unless a test sets it (the `chunk_size` attribute then does not depend on the free memory of the host); thread override and dirty-tree override are passed by the `ex_campaign` wrapper because the repository under test is dirty and the suite runs on 4 threads while a fixture block may name another k19 value.
- Reference code (independent of the module): the recording grid of R2 item 1 at 1024 bits (`ex_ref_grid`), the Cartesian product and merge of R1 items 2 and 3 (`ex_ref_items`), the backend table of R4 (`ex_expect_backend`), the direct estimator call of O9 (`ex_direct`, lower-layer API only: `propagate` with `W0`, `propagate_shadow_orbit`), direct iteration of the map (`ex_iter_states`), the exact finite-time exponent of a linear map at 256 bits and the factor c (`ex_lin_lambda`, `ex_c_factor`), the platform SHA-256 over many texts in one process (`ex_sha256_many`). All were run against throwaway stubs in the scratchpad (grid counts 201 and 221, fixture counts, TOML round trip of the emitter in both key orders, direct calls for all methods and precisions, Y5 error 0 against bound 2e-15).

### New dependencies

- `HDF5` (package, `[deps]` and `[compat]` of `Project.toml`): required by the store. The tests never `using HDF5`: they read and write only through `read_chunk`, `write_chunk` and `chunk_layout`, so the RED failure is a missing Chaoslab function. D53u asserts that `HDF5` is a declared dependency and that every package imported under `src/Experiments/` is declared.
- A JSON reader and writer for the manifest and the pilot report: a package (JSON3, JSON) or a small writer in the module, stage 03 chooses; the tests go through `read_json`. The hash library is `SHA` (already a dependency); `Dates` (UTC strings) must be declared if used.
- No test dependency is added: the suite uses `Test`, `StableRNGs`, `TOML`, `StaticArrays`, `LinearAlgebra`, the platform `sha256sum` (fallback `certutil`) and `git`.

### New API (all in `Chaoslab.Experiments`, exported through `Chaoslab`; checked for collisions against all exports of the lower layers: none)

| Name | Contract tested |
|------|-----------------|
| `run_block(dict)`, `read_run_block(prereg_path)` | validate the `[run]` table (string keys, as parsed from the TOML) and return an object; `ArgumentError` naming the key, before any computation. `read_run_block` first applies the E18 check, reads LF and CRLF files, takes the `run` table of the TOML block that has one |
| `block_text(rb)`, `block_hash(rb)` | canonical text of the parsed block (independent of the textual key order of the file) and its SHA-256 (hex) |
| `plan_run(prereg, panel; code_version=nothing, defaults=nothing, ntasks=nothing)` | the plan: `.items`, `.chunks`, `.counts (items, chunks)`, `.steps`, `.hash`. Errors: header hash, panel hash outside the accepted list, block errors, zero items. `code_version` defaults to `code_version()`; `defaults` is the path of `numerics-defaults.toml` (default: the repository file) |
| item (NamedTuple in `plan.items`) | `id, orbit, system, param, x0, p0, run_seed, orbit_seed, method, precision, horizon, tau, eps, direction, coordinate_mode, vectors, grid_sha256, code_version, hash_version, hash, arms, design_class, cell`. `precision` is the label with bits (`Float64`, `Double64`, `BigFloat(256)`); `x0, p0, param, direction` are texts; `eps` is `nothing` except for L2, else the text `string(Float64)`; `system` is `standard` or `linear` (`param` is K or mu); `arms` lists arm names in block order |
| `item_text(item)`, `item_hash(item)` | LF text of the 18 hashed fields, one `name=value` line each in the order `orbit, system, param, x0, p0, run_seed, orbit_seed, method, precision, horizon, tau, eps (L2 only), direction, coordinate_mode, vectors, grid_sha256, code_version, hash_version`; hash = SHA-256 of that text. Extra fields of the item are ignored |
| chunk (NamedTuple in `plan.chunks`) | `ordinal, key, hash, cell, run_type (method, precision, horizon, tau, eps, direction, coordinate_mode), items (item ids), orbits (orbit ids), arms`; cell is `k<k_index>:<class>`, `controls` or `h4` |
| `chunk_text(plan, ordinal)` | LF text whose SHA-256 is the chunk hash: the key, the sorted item hashes, the values that decide the status (the three health multipliers of the defaults, `fail_fraction`, `control_bound`, `drift_bound`, `setup_check_tol`) |
| `plan_text(plan)` | canonical bytes of the plan (counts, chunk table with hashes, panel SHA-256, block SHA-256, defaults, code version); no path, host, time or preregistration file hash; `plan.hash` is its SHA-256 |
| `plan.steps` | `arms[name] = (items_before, steps_before, steps_after)`, `precisions[label] = (steps_before, steps_after)`, `total_before`, `total_after`, `saving`; "after" assigns a merged item to its first arm |
| `record_grid(N, tau; points_per_decade, first_steps, extra_fractions)`, `grid_hash(grid)` | R2 item 1; hash is a 64-digit hex SHA-256 of the canonical grid text |
| `deviation_frame(T, direction_text, k)` | `SMatrix{2,k,T}`: first column (cos a, sin a) with a evaluated at 256 bits, second column the counterclockwise quarter turn `(-W[2,1], W[1,1])` exactly |
| `code_version(root=<package root>)` | content identifier (hex) of the files under `src/` and the lock `Manifest.toml`; not of other files |
| `launch_capacity(free_bytes; f, M, reducer, columns, grid_steps)` | `floor(f free / (M + reducer + 8 columns grid_steps))` |
| `launch_groups(plan; capacity)` | consecutive chunks of one run type up to the capacity: `(chunks, orbits, launches)`; a chunk larger than the capacity gives `launches > 1` |
| `resolve_run_backend(request, T, method; threads, record)` | `(backend, device, reason, note)`; `record` is `nothing` or `(threads, gpu_faster)`; the device is `Foundation.DEVICE[]` |
| `cpu_only_method(id)` | true for L6 (library Householder path) and false for the other ids tested |
| `backend_table(plan; threads, record)` | `Dict{(precision label, method) => (request, backend, device, reason, note)}` |
| `RunHooks(; clock, kill, mutate, observe)` | test hooks, see below |
| `run_campaign(prereg, panel; store, manifest, threads, thread_override=false, allow_dirty=false, repo=nothing, code_version=nothing, defaults=nothing, capacity=nothing, free_bytes=nothing, order=:plan, stale=:stop, retry_failed=false, cap_override=String[], device_record=nothing, prior_record=nothing, hooks=RunHooks())` | returns `(status, reasons, ran, skipped, failed, paused, stale, manifest, plan)`; refusals at start throw `ArgumentError` (nothing written); stops and pauses return a status |
| `read_chunk(path)`, `write_chunk(path, attrs, orbits)`, `chunk_layout(path)` | store reader and atomic writer, see "Store contract" |
| `read_json(path)` | parsed JSON (objects as `Dict{String,Any}`, integers as `Int`, other numbers as `Float64`) |
| `run_pilot(prereg, panel; out, threads, thread_override, allow_dirty, repo, code_version, capacity, device_record, prior_record, hooks)` | returns `(report, path, starts)`; writes `<out>/pilot-report.json` and the pilot store `<out>/store/chunk-NNNNN.h5` |
| `rerun_check(prereg, panel, store; code_version, threads, thread_override, allow_dirty, sample=nothing)` | CPU reference rerun of `[run.verify].rerun_sample` items per arm; vector of `(item, arm, backend, bitwise, max_abs_diff)` |
| `exit_status(status)`, `cli(args; io=stdout)` | `0` only for `complete`; `cli` returns the exit code of `panel`, `run`, `pilot` |
| `default_store(prereg; root=<repo>/data)` | `<root>/<study slug>` |
| `hot_loop_bytes(method; T=Float64)` | bytes allocated by the recording loop of one orbit after warm-up |

### Contracts fixed by the tests

- Status strings (campaign): `complete`, `complete with failed chunks`, `stopped: control bound`, `stopped: failed chunks`, `stopped: wall cap`, `stopped: preregistration changed`, `stopped: panel changed`, `paused arms`. Arm status in the manifest: `complete`, `paused: cost cap`. Chunk status in the manifest: `ok`, `failed`, `pending`, `paused`. Setup check in the manifest chunk entry: `ok`, `failed`, `no reference`, `not applicable`.
- Reasons (substrings): `non-finite` with `orbit <id>`, `diagnostics exceeded`, `setup check`, `control bound` with `orbit <id>` and the method id. Flags: `diagnostics exceeded`, `drift`, `excluded`. Exceedance and drift are strict (`>`); the failure fraction rule is a count strictly above `fail_fraction x n`.
- Exit codes: 0 complete, 1 refusal or usage error (message on `io`), 2 complete with failed chunks, 3 stopped, 4 paused arms. `cli`: `panel <prereg> --out <file>`; `run <prereg> <panel> --store <dir> --manifest <file>`; `pilot <prereg> <panel> --out <dir>`; options `--threads N`, `--thread-override`, `--allow-dirty`, `--repo <dir>`, `--code-version <text>`, `--capacity N`, `--cap-override a,b`, `--stale move`, `--retry-failed`. `bin/chaoslab.jl` is at most 12 code lines, calls `Chaoslab.Experiments.cli(ARGS)` and exits with its value.
- Directory names inside a store: `aside/` (a file replaced by a rerun or `retry_failed`), `stale/` (files moved by `stale = :move`); temporary files end in `.tmp` and sit in the store directory; they are listed in `manifest["removed_temporaries"]` at the next start and removed.
- Order of processing: the plan order; `order = :reverse` processes the chunks in reverse plan order; a launch group never mixes run types.

### Test hooks (D53l, D53n, D53o, injected faults)

`RunHooks(; clock=nothing, kill=nothing, mutate=nothing, observe=nothing)`; every field optional; the module calls a hook only through this struct.

| Hook | Call | Contract |
|------|------|----------|
| `clock(event, info) -> Float64` | events `:session_start`, `:group_start`, `:group_end`, `:session_end` | the only source of time (seconds since the Unix epoch, default `time()`); a group's wall seconds are `clock(:group_end) - clock(:group_start)`; UTC strings of the manifest come from it; for group events `info = (group, chunks, arms, orbits, steps, backend)` with `steps` the orbit-steps of the group and `arms` the union of the arm names of its chunks; the wall cap is checked before each group from the manifest's sessions plus the current session |
| `kill(point, info)` | `:group_start` with `info = (group, chunks)` before the input checks and the cap checks of that group; `:before_data`, `:during_data`, `:after_data`, `:before_rename`, `:after_rename`, `:chunk_done` with `info = (ordinal, path, chunks_done)` | called from the writer, which is the main task; temporary file created before `:before_data`, half of the orbit groups written at `:during_data`, closed and read back at `:after_data`, `complete` set at `:before_rename`; a thrown exception propagates unchanged (no catch, no cleanup, no wrapping); side effects are allowed (D53p changes an input file at `:group_start`) |
| `mutate(info, orbits)` | after a chunk attempt is computed, before the health checks | `info = (ordinal, attempt, backend, chunk, items, arms)` with `backend` the symbol of the attempt (`:cpu` on the rerun); `orbits` is the vector of orbit records below; arrays and dictionaries are mutated in place |
| `observe(info, attrs, orbits)` | after checks and flags, just before the write | `info = (ordinal, attempt)`; the file read back equals what was observed, bit for bit |

Orbit record (also returned by `read_chunk` and accepted by `write_chunk`): `(id, attrs, t, state, winding, lyap, diagnostics, thresholds, first_exceedance, cost)` with `t::Vector{Int64}`, `state::Matrix{Float64}` (n x 2), `winding::Union{Nothing,Matrix{Int64}}` (wrapped mode only), `lyap::Dict{String,Matrix{Float64}}` keyed by method, `diagnostics::Dict{String,Vector{Float64}}` running maxima keyed `det, symp, orth, sum, pair`, `thresholds::Dict{String,Float64}`, `first_exceedance::Dict{String,Int}` (0 = none), `cost::Dict{String,Vector{Float64}}` = `[steps, vectors, wall_seconds, thread_seconds, resident_bytes]` per method. Diagnostics present: L1 `det, symp`; QR family (L3, L4, L4b, L5, L5b, L6) `det, symp, orth, sum, pair`; L9 `det, symp, sum, pair`; L2 none. `resident_bytes` is per orbit and does not depend on the launch capacity (I12).

### Store contract (D53m)

`read_chunk(path)` returns `(attrs::Dict{String,Any}, orbits)`; string lists are `Vector{String}` (empty when none); `complete` and `pilot` are `Bool`. Root attributes: every name of the table "Root attributes" of `schema/run-store.md` except `python_version`; `perturbation` only for L2. `chunk_layout(path)` returns `Dict{path => (eltype = "Int64", size = (n,))}` read from the HDF5 file: `t` Int64 (n), `state` Float64 (n, 2), `winding` Int64 (n, 2) only in wrapped mode, `lyap/<method>` Float64 (n, k), `diagnostics/<name>` Float64 (n), `cost/<method>` Float64 (5), nothing under `/grid`. File names `chunk-NNNNN.h5` with the zero-padded ordinal of the plan. Orbit attributes: `system, params, ic, class_label, seed, item_hash, design_class, flags, steps_recorded`; `seed` is the per-orbit seed (text or integer).

### Manifest contract (D53p)

Top-level keys: `study, schema_version, plan_sha256, plan_counts (items, chunks, steps_before, steps_after), status, reasons, sessions (git_sha, clean, start_utc, end_utc, threads_used), git (sha, clean), versions (language, packages_sha256, hdf5, module), host (cpu, cores), threads (used, available), backend_table (precision, method, request, backend, reason, note), device, device_record, prior_record, prereg (sha256_start, sha256_end, header), panel (sha256_start, sha256_end, build_hash), block (text, sha256), defaults, launch (free_bytes, f, capacity), chunks, arms, wall (total_seconds, cap_hours), overrides, removed_temporaries`. `device_record` and `prior_record` are objects `(threads, cpu_ns_per_step, gpu_ns_per_step, gpu_faster, date)` or the text `absent`. Chunk entries: `ordinal, key, hash, arms, status, attempt, backend, fallback_reason, orbits, steps, wall_seconds, first_attempt_reasons, flags (diagnostics exceeded, drift, excluded), setup_check`. Arm entries: `steps, wall_seconds, thread_seconds, cost_ns_per_orbit_step, budget_ns, status`. Overrides: `dirty (paths, diff_sha256)`, `threads (used, k19)` only when the thread count differs from k19, `cap` (arm names).

### Criterion to test

| Criterion | Test (file: testset) | Sample | Pass |
|-----------|---------------------|--------|------|
| API | ex_run_inputs: "API exported through Chaoslab (D53)" | 30 names | exported by the submodule and by `Chaoslab` |
| D53a | ex_run_inputs: "D53a block validation" | 65 malformed blocks: missing sections and keys (16), unknown keys (4), unknown ids and requests (7), empty lists (4), duplicate arm name, horizon not a multiple of tau, N / 10 not a multiple of tau, extra fractions (4), grid integers (4), fail and stop fractions (4), safety, cost factor, wall cap, control_bound, drift_bound, setup_check_tol, k19 threads, rerun_sample, setup_check_steps (2), duplicate h4 label, pilot (2), accepted list (2), h1_check type, eps, tau, missing budget | each `ArgumentError` names the key |
| D53a | ex_run_inputs: "D53a refusals before any computation, nothing written" | preregistration hash not equal to its header, panel hash outside the accepted list, thread count differing from k19 without override, malformed block in a file, zero-item plan | each refused, no store directory, no manifest |
| D53b | ex_plan: "D53b plan on fixture small / alt / full" | hand counts 52 and 66, reference product with merge for three blocks | item and chunk counts, ids 1..M, key set and arm lists, order (precision, first arm, orbit), steps per arm before and after, saving, chunk membership, plan bytes identical for ntasks 1, 4 and all and for the reversed textual key order |
| D53b | ex_plan: "D53b merged items list every arm" | `:full` | items of four, two and three arms; chunk keys and hashes unchanged when arms are renamed |
| D53b | ex_plan: "D53b study block on a synthetic panel of the study caps" | the real A4 block, 1100 / 200 / 100 orbits, plan only | item count, steps per arm equal the reference, main arm 1100 x 9 x 10^7 |
| D53c | ex_plan: "D53c hashes" | alt fixture | hash equals the text hash, 25 items by the platform tool, 17 field cases plus eps, exclusions (backend, k19, arm names, ntasks), code version, chunk hash with four status values, chunk text |
| D53c | ex_plan: "D53c independent hash on 1000 seeded items" | 1000 random items | module hash = platform `sha256sum` of the module text (one process) |
| D53c | ex_plan: "D53c code version" | scratch tree | changes with `src/` and the lock only |
| D53d | ex_plan: "D53d recording grid" | N in {1e3 .. 1e7}, tau in {1, 2, 4, 8, 16} | 0 differences against 1024 bits, 201 and 221 points, distinct grid hashes |
| D53e | ex_chunks: "D53e series content" | `:full` campaign, 223 items | grid, shapes, diagnostics sets, monotone maxima, thresholds, finite, dtypes, lifted state 0 ulp, wrapped state and winding exact, cost record, attributes |
| D53f | ex_chunks: "D53f direct call (O9)" | 5 seeded grid steps per item, all items of the run with capacity 64, a third of the items of the runs with capacities 1 and 7 | lambda equal to the direct estimator call, 0 difference |
| D53g | ex_plan: "D53g chunk ids, keys and launch capacity" | 3 free-memory values, f = 1/2 and 1/4 | exact formula, halving, groups partition the plan, ids and hashes unchanged |
| D53g | ex_chunks: "D53g stored datasets bit-identical ..." | capacity 1 and reverse order; capacity 7 | canonical datasets equal the base campaign |
| D53h | ex_plan: "D53h backend resolution (R4)" | 3 precisions x 3 requests x 4 device states x 4 records x 2 reducers | equals the written table, no tenth reason, thread-count note, plan table |
| D53i | ex_chunks: "D53g ...; D53i threads" | 1, 2 and 4 threads, all threads on the small block | datasets bit-identical, thread count in the manifest |
| D53j | ex_health: "non-finite values", "diagnostics fraction", "exact-exponent controls", "setup check and drift" | NaN, Inf; 40 orbits at fail_fraction 0.05, 0.1, 0.25 (count at the limit and one above); controls in h1_check and another arm, exponent zero; setup check at 2 and 0.5 setup_check_tol, no reference; drift at bound and bound (1 + 1e-9) | as in the criterion |
| D53k | ex_health: "D53k failure and rerun" | transient and persistent fault, `:c40` with 2 and 3 failed chunks, `:c20` with 1 and 2 | attempt 2 on CPU, reasons, stop at the strict count, manifest statuses sum to the plan (I17) |
| D53l | ex_health: "D53l wall cap ...", "D53l first-chunk cost cap ..." | fake clock 150 s per group against 450 s; costs of 2.01 and 1.99 times the line | stop before the fourth group, resume adds the sessions; pause at 2.01, continue at 1.99, override resumes and is recorded |
| D53m | ex_store: "D53m store layout", "D53m write and read back bit for bit" | `:small` and `:alt` stores, observe hook, special values | every attribute and dataset of R10 with dtype, file names, schema document, read-back equals what was written |
| D53n | ex_store: "D53n atomic write killed at 20 seeded points" | 20 (chunk, point) pairs, each point four times | no incomplete final file, temporaries listed and removed, chunk recomputed, store equals the reference |
| D53o | ex_resume: (a), (b), (c), (d) | finished store; stale by code version and by fail_fraction; 10 seeded k with another thread count; failed chunk | 0 chunks and 0 files touched; stale refused and with the option moved to `stale/`; 0 differences and equal SHA-256; failed final, retried on request, old file in `aside/` |
| D53p | ex_manifest: three testsets | scratch repository, fake records, changed inputs | every field of R7, clean flag, refusals on tracked change and untracked file, override recorded, input change stops with a resumable store |
| D53q | ex_pilot: "D53q pilot" | `:small`, fake clock 150 ns per orbit-step | isolation, marked files, starts at least 1e-6 from panel starts and deterministic, report fields, projection equals steps x measured cost within 1e-9 relative, later study run bit-identical, dirty override recorded |
| D53r | ex_run_inputs: "D53r ..."; ex_chunks: "D53r the runner with the second fixture block" | real A4 block; two blocks differ in every value; `:alt` campaign (grid, direct call, diagnostics, 1 and 3 launch capacity) | the real block is accepted; no distinctive number or string of the study block in the new module files |
| D53s | ex_chunks: "D53s Y5 and exact exponents (O11)" | Y5 and the K = 2.5 hyperbolic fixed point, 8 and 7 methods, every grid step of horizon 200 | `abs(lambda - exact) <= B_S / n` (form of D35) with the exponent from 256 bits; asymptotic form ln mu + ln(c) / n for n >= 20 |
| D53t | ex_chunks: "D53t CPU reference rerun (O12)" | rerun sample 3 per arm, 2 and 1 threads | sample size, same sample for both thread counts, CPU items bitwise |
| D53u | ex_manifest: "D53u layering and declared dependencies", "D53u entry point" | Project.toml, imports, script inspection, `cli`, one subprocess | HDF5 declared and all imports declared; script is at most 12 lines without numerics; exit codes; the script with no arguments exits nonzero |
| D53v | ex_chunks: "D53v time and memory of the run itself" | nine methods; the largest chunk file | zero bytes after warm-up; writer allocation below series bytes + 64 MB |
| D53x | ex_manifest: "D53x isolation ..." | default store of the fixture | named `<data>/ex-fixture`, empty after the suite |

### Stage 04 mapping and unit proxies

| Criterion or part | Why not in the package suite | Unit proxy here |
|-------------------|------------------------------|-----------------|
| D53b [full]: Budget table rows and totals with the real panel | needs the real panel file | plan only with the real A4 block and a synthetic panel of the caps; counts equal the reference and 1100 x 9 x 10^7 |
| D53f: all nine methods at Double64 and BigFloat | suite time (D38 policy) | Double64 and BigFloat carry L5 only; L6 is Float64 only by construction; the full product runs in `verify/` with the same assertions |
| D53f chunk sizes 1, 3 and all | suite time | capacities 1, 7 and all on the full fixture, 3 on the second block |
| D53h with a real device | needs the GPU ([local]) | simulated devices; the record resolution on the local device in D53w |
| D53i at the process level (`julia -t 1/4/N`) | needs separate processes | task counts 1, 2, 4, all inside one process |
| D53m: the independent reader (O13) | Python or HDF5 tools outside the package | read-back through the module, dtypes through `chunk_layout` |
| D53q factor 2 against a measured run of the real arms | timing noise | injected clock: projection equals steps x measured cost within 1e-9 |
| D53s on the real arms at 1e7 | long run | the same form at horizon 200 |
| D53t on the device | [local] | CPU items only |
| D53u: script run as a subprocess for each command | compile time | `cli` calls plus one subprocess with no arguments |
| D53w reduced campaign, D53x on the real run | stage 04 by the spec | none (D53x: the default store is empty) |
| O14 independent plan and hash script in `verify/` | stage 04 | independent hashing here (platform tool) |

### Ambiguities and choices made by the tests

1. Run-type nesting (R1 item 4): the spec orders items by precision, arm and "run type in block order" without naming the nesting of horizon, mode, direction, tau, eps and method. The tests check precision-major, then the first arm, then the orbit order inside a run type, and the exact set of items.
2. Plan hash and the preregistration: I11 lists the preregistration bytes among the inputs, D53b requires equal plan bytes for another textual order of the tables. The tests require the plan text to hold the panel SHA-256 and the block SHA-256 but not the preregistration file hash (the manifest holds it).
3. Setup check span: "the first setup_check_steps grid steps" is read as grid entries with n <= setup_check_steps.
4. K19 gate: applied to the request `auto` only; an explicit `gpu` request on a functional device resolves to the device without a record.
5. Code version: a content identifier of the working-tree files under `src/` and the lock (the spec says "committed content"); the tests change files without committing.
6. Fixture size: three controls instead of two (an elliptic fixed point with exponent zero and a hyperbolic fixed point at K = 2.5 are both needed); Double64 and BigFloat arms carry L5 only.
7. Control bound 0.05 in the fixtures (see Conventions).
8. Exit codes other than 0 are fixed by the tests (2, 3, 4); the spec only requires nonzero.
9. Stop rule granularity: the failed-chunk stop is taken at a launch-group boundary; the tests put the failing chunks in different groups and bound the last chunk run by the end of the group that holds the third failure.
10. D53r source scan: the files of the first scope (panel builder) are excluded from the scan, because they hold the literal 1e9 for another purpose; numbers with at least four significant digits and values of at least 1e5 of the real block are searched as numeric literals in the other files of the module.
11. The thread override flag is required by `run_pilot` and `rerun_check` too, because they resolve the thread count by the same rule.
12. Existing suite: the AC2 testset "frozen p1-le-benchmark" (`ex_inputs.jl` line 82) errors in this working tree before any D53 test is reached: the working copy of the preregistration has CRLF line endings (index LF, `core.autocrlf = true`), and that test matches ```` ```toml\n ```` on the raw text. With an LF copy the same testset passes (35 passes in AC2 instead of 31 and one error, so the baseline 2632 is restored). Not changed here; `read_run_block` must read both line endings, as `check_prereg` does.
13. Suite time: not measurable before the module exists; the estimate is about 60 s of campaigns (three full-fixture campaigns of 76 chunks, about 120 small-fixture campaigns of at most 18 chunks) plus compile time of the nine methods at Float64 and of Double64 and BigFloat L5 already paid by the other suites; stage 03 measures it against the 120 s ceiling and moves cases to `verify/` (D38) if it is exceeded.

### RED log

Command: `julia --project -t 4 -e 'using Test, Chaoslab; @testset "experiments" begin include("test/experiments/runtests.jl") end'` (Julia 1.13.2, 4 threads, 1 min 23 s including the existing suites). Result: `experiments | Pass 2628 | Fail 12 | Error 61 | Total 2701`.

- Existing suites: 2628 passes and 1 error. The error is the pre-existing one of ambiguity 12 (CRLF in the working tree; unrelated to D53); with an LF copy of the preregistration the same testset passes and the count is the 2632 of the baseline. No other existing testset changed.
- New: 44 testsets (48 in the source, D53a counted once), Pass 0, Fail 12, Error 60.
- Errors (60): 59 `UndefVarError` in `Chaoslab.Experiments`, namely `run_campaign` (14), `plan_run` (14), `hot_loop_bytes` (9), `launch_capacity` (7), `RunHooks` (5), `read_chunk` (2) and once each `run_block`, `read_run_block`, `item_text`, `record_grid`, `deviation_frame`, `resolve_run_backend`, `code_version`, `default_store`; 1 `SystemError` opening the missing `bin/chaoslab.jl`.
- Fails (12): 2 export checks (names not defined, ex_run_inputs "API exported"), 8 `ex_rejects` expressions that are false because the call raises `UndefVarError` and not `ArgumentError` (D53a refusals 7, D53p clean tree 1), 1 `HDF5` not in `[deps]` (D53u layering), 1 `bin/chaoslab.jl` missing (D53u entry point).
- No parse errors, no errors inside helpers: every fixture, writer and hook builder ran (the fixtures and files are generated before the first missing name).
- Not run against an implementation: the code behind the first missing name of each testset (all campaign-level testsets). Checked in the scratchpad against throwaway stubs (not in the repository): `ex_ref_grid` (201 and 221 points, 77 assertions of D53d), the product with merge (counts 52, 66, 223, 17, 143, 75, 40, hand counts agree), the TOML emitter in both key orders (round trip equals the dictionaries), `ex_c_factor` and `ex_lin_lambda` against the direct `propagate` (Y5 error 0, fixed point within the bound), `ex_direct` and `ex_iter_states` for all nine methods at Float64, Double64 and BigFloat in lifted and wrapped mode, the deviation-frame testset (36 assertions), the code-version, launch-capacity and backend-matrix testsets against stubs built on `Foundation.resolve_backend` (matrix of 288 cases agrees with the written table), `ex_sha256_many` (1000 texts in one process).

### Correction during D53 stage 02 (2026-10-18)

`ex_inputs.jl` AC2 "frozen p1-le-benchmark": the real preregistration text is read with CRLF normalised to LF before the regex and hash helpers. A Windows checkout with `core.autocrlf` writes CRLF, and the block regex then found nothing; the module hash rule already normalises (E18).

### Corrections during D53 stage 03 (2026-10-19)

Four test defects found by the implementation, fixed in the tests (no src change):

- `ex_run_helpers.jl` `ex_sys_x` and `ex_chunks.jl` D53e: the Y5 control (linear map) has no torus; in wrapped mode it runs in plane coordinates with zero winding counts.
- `ex_health.jl` D53l: the wall-cap fixture sets `first_chunk_cost_factor` to 1e9 so the first-chunk cost cap (its own testset) does not pause the arm first.
- `ex_manifest.jl` D53p: the restore writes back a copy of the original bytes (`String(orig)` had emptied the buffer).
- `ex_pilot.jl` D53q: the injected clock starts at 0, since binary64 spacing at 1.7e9 s exceeds the injected durations and made rtol 1e-9 unreachable.

Result: experiments suite 4300 of 4300.

## D53 device path (return of the runner to stage 03, 2026-10-19)

The owner sent the D53 runner back to stage 03 to add the device path for the Float64 kernel methods. Stage 02 role first: new file `test/experiments/ex_device.jl` (registered in `runtests.jl`, last include); no existing assertion changed. The package suite has no GPU, so the device is simulated by `ExHookDev`, a device type whose backend is the KernelAbstractions CPU backend (the kernel CUDA would run executes in the test), with a hook before every launch (`Tangent._before_launch`), a one-wave record value, a free-memory figure and a compile flag. The real device is stage 04 (`verify/experiments-device.jl`).

New interface under test: `Chaoslab.Experiments.device_records(ctx, items; hook = true) -> (recs, prov)`; `recs === nothing` when the group falls back, and `prov` holds `backend`, `device`, `reason`, `gpu_chunks`, `chunk_size`, `launch_steps`, `launches`. Campaigns and pilots on a present device record these per chunk (`backend`, `device`, `backend_fallback_reason`, `gpu_chunks`, `chunk_size`, `steps_per_launch`, `launches_per_chunk`).

| Testset | Assertions | What it checks |
|---------|-----------|----------------|
| methods with a kernel path | 10 | the recorder state is isbits exactly for the methods that resolve to the device (L1, L3, L4, L4b, L5, L5b, L9; L2 through `ShadowEnd`); L6 (library QR) is the only CPU-only method |
| records equal the CPU loop on the KA CPU backend | 4 | every run type of the full fixture at Float64 for the eight kernel methods (lifted and wrapped, tau 1, 2, 4, directions, eps, the Y5 control in the group): the device records equal `compute_records` bit for bit (t, state, winding, lambda, diagnostics, thresholds, first exceedance, resident bytes); one launch of the full horizon |
| relaunch boundaries not aligned with grid steps | 6 | forced launch lengths 1, 7 and 13 (tau 1) and 28 (tau 1, 2, 4): records identical to the CPU loop, `launches == cld(N, L)`, the hook runs once per launch and batch; a length that is not a multiple of tau is refused (`ArgumentError`) |
| chunk size from free device memory | 3 | batches from the free memory (`free = 8000` gives batches smaller than the group, one batch per system type with plenty), records unchanged |
| a failing orbit is recorded, not thrown | 8 | a state corrupted between launches (launch 2: failure in a segment that ends on a grid step; launch 16: failure in the tail of a launch that holds no grid step, recorded at the next grid step): NaN rows from that grid row, finite and bit-equal rows before it, the message names the state and the step, the other orbits untouched |
| fallback reasons of a group | 5 | out of memory (fourth reason), compile error (third), lost device (ninth, the device is dropped), launch over the budget (eighth, `ns = 1e9`), an unclassified exception is rethrown |
| mixed backend table, records and provenance of a campaign | 236 | full fixture, device present, record at the run's thread count: Float64 rows of eight methods `gpu`, L6 `cpu` with "reducer is CPU-only", Double64 and BigFloat `cpu` with an empty reason; every chunk file equal bit for bit to the CPU campaign (`ex_store_canon`); per chunk `backend`, `device`, reason, `gpu_chunks`, `steps_per_launch`, `launches_per_chunk`, `chunk_size`; `rerun_check` on the device store (O12): bit for bit on the KA CPU backend |
| relaunch inside a campaign | 19 | forced launch 14 on the small fixture: store equal to the CPU campaign, `steps_per_launch == 14`, `launches_per_chunk == cld(N, 14)` |
| fallback reasons per pair in a campaign | 9 | out of memory in the first launch only (that group CPU with the fourth reason, later groups return to the device), lost device (ninth reason for every Float64 chunk), compile error, watchdog; the CPU records equal the reference store |
| K19 record thread count against the run thread count | 5 | record at 8 threads, run at 4: `cpu`, seventh reason, note names the mismatch; `gpu_faster` false: seventh reason; record at the run's thread count: `gpu`; an explicit `gpu` request needs no record |
| failed chunk reruns on the CPU, cost cap uses the device line | 6 | a fault on the device attempt: attempts logged `(1, :gpu)`, `(2, :cpu)`, file attempt 2 with `backend = cpu` and the fifth reason; 100 ns per orbit-step passes the CPU line and pauses the arm on the device line (`ns_gpu x safety`) |
| pilot on the device | 4 | cost rows `gpu` for Float64 and `cpu` for Double64, one worker for the device (thread seconds equal wall seconds), chunk files with `backend` `gpu` |

Choices: (1) bit equality on the KA CPU backend, the form the tangent-core tests use there; across the real device only backend_agreement holds (stage 04). (2) The simulated device omits hooks at the warm-up launch (`hook = false`), so a fault injected by a test hits a real group. (3) The Y5 control is a linear map and runs in its own batch (a kernel takes one system type), so `gpu_chunks` counts batches per system type. (4) A device fault falls back for that group only; a lost device keeps every later group on the CPU.

### RED log (before the code)

Command: `julia --project -t 4 -e 'using Test, Chaoslab; @testset "experiments" begin include("test/experiments/runtests.jl") end'` (Julia 1.13.2, 4 threads, 2 min 36 s). Result: `experiments | Pass 4349 | Fail 240 | Error 6 | Total 4595`. Existing suites: 4300 of 4300 pass. New part: 12 testsets, 49 passes (the isbits testset, which needs no new code, and the campaign assertions that hold on the CPU path), 240 fails (every chunk recorded as `cpu` where `gpu` is required, launch counts and reasons absent) and 6 errors: `UndefVarError: device_records not defined in Chaoslab.Experiments` in 5 testsets, and one test defect (a two-argument anonymous function over dictionary pairs in the backend-table check) fixed before the code was written. After the code: 4615 of 4615 (315 new assertions).

## Amendment D54 (per-method record on the recording chain, 2026-10-19)

New file `test/experiments/ex_record.jl` (registered in `runtests.jl`, last include); no existing assertion changed. Interface: a record NamedTuple with a `recorder` field (method id to `(threads, gpu_faster, ...)`) is read per method by `resolve_run_backend` and `backend_table`; `Experiments.gpu_record(root)` returns the field from the file.

| Testset | Assertions | Criterion |
|---------|-----------|-----------|
| resolution reads the per-method record | 8 | D54a, D54c: device only for the methods with an entry at the run's threads and `gpu_faster` true; top-level flag ignored; absent method, other thread count (note) and not faster give the seventh reason; L6 the sixth; explicit requests; Double64 and BigFloat; empty table |
| matrix of the resolution with per-method records | 1 (a list of mismatches that must be empty, 3 precisions x 3 requests x 4 device states x 4 record states x 3 methods) | D54b |
| gpu_record reads the layout of the file | 4 | D54d |
| campaign: backend table and chunks follow the per-method record | 21 | D54e |

### RED log

Before the code, `resolve_run_backend` ignored the field: `julia --project -t 4` on the quick runner (helpers plus ex_device.jl plus ex_record.jl): first testset Pass 3, Fail 5 (every case that reads a per-method entry or the absence of one); the remaining three testsets were not reached in that run because a failed top-level testset ends the runner. Running them against the unchanged code is not informative either: `gpu_record` returned no `recorder` field and the campaign used the single record. After the code: 34 of 34 (8, 1, 4, 21). D54f is the verify script.


## Amendment D53 review fixes (stage 03 rework after the stage 04 independent review, 2026-10-19)

Stage 02 role first: the findings H1, H2, M1 to M9 and the cheap LOW items of `experiments-runner-verification.md` ("Independent review") get tests before any code. New file `test/experiments/ex_review.jl` (12 testsets, 108 assertions), registered as the last include of `test/experiments/runtests.jl`. No existing assertion is changed. Two helper constants/calls in `ex_run_helpers.jl` change because the new behaviour makes them necessary (they are not assertions): `EX_ROOT_SKIP` gains `"session"` (the session index differs between a killed-and-resumed store and an uninterrupted one, like `created_utc`), and `ex_campaign` passes `manifest_reset = true` (D53o (a), (b) copy a finished store without its manifest and resume it; M4 now refuses that unless asked, and `ex_review.jl` calls `run_campaign` with `manifest_reset = false` to test the refusal).

New interface under test (all in `Chaoslab.Experiments`): `run_campaign` keywords `retry_control_stop`, `backend_override`, `manifest_reset`; root attributes `dirty`, `session`, `control_stop`; manifest keys `control_stops`, `sessions[i].index`, `overrides.retry_control_stop`, `overrides.backend_table`, `overrides.manifest_reset`; manifest backups `<manifest>.1` to `.3`; internal seams `_move_retry`, `_write_json(...; keep, mover)`, `_fsync`, `_attr_diff`, `_check_threads`; the group info gains `precision`; command-line flags `--retry-control-stop`, `--backend-override`, `--manifest-reset`.

| Finding | Test (ex_review.jl: testset) | Sample | Pass |
|---------|------------------------------|--------|------|
| H1 | "H1: a control-bound stop is not a result" | fixture small, control of chunk L5/controls bumped by 1.0 | file `control_stop = true`, status failed, no `excluded` flag; manifest `control_stops` record; plain resume refused and nothing touched; changed code version gives the stale refusal; `retry_control_stop = true` recomputes (file to `aside/`, override and `retried` recorded, status complete) |
| H2 | "H2: the backend table is compared at resume" | fake device with record, killed after 2 chunks, resumed without device / with the same device / without manifest | refusal naming the table, nothing written; `backend_override` goes on and records the changes; same table accepted; with no manifest the kept chunk written on `gpu` against a CPU table is refused |
| M1 | "M1: renames retry ..." | mover that throws EBUSY or EACCES twice, ENOENT, EBUSY always; four rewrites of a JSON file; two writes to one chunk name | retries 3 calls, ENOENT not retried, failure after `tries`; failed replace leaves the previous file whole; backups `.1` to `.3` newest first; taken chunk name moved to `.old1`, new content at the final name |
| M2 | "M2: root attributes in the read-back ..." | a string attribute with an embedded NUL character (HDF5 strings end at the NUL), `_attr_diff` cases, `_fsync` | write stops naming the attribute, no final file; fsync on a file and an error on a missing one |
| M3 | "M3: the session clock moves ..." | clock that ticks only on `:manifest_write`, killed after 2 chunks | `wall_seconds` of the manifest equals the number of manifest writes |
| M4 | "M4: chunk files without a manifest" | killed store with unreadable or missing manifest, with and without backups, first run on an empty store | refusal naming `manifest_reset`, nothing written; option goes on and is recorded; backup used (two sessions kept); no option needed on an empty store |
| M5 | "M5: code version covers the code and nothing else" | temporary root with `src`, `ext`, `Project.toml`, lock | `*.md`, `.cov`, `.txt`, `docs/` unchanged; CRLF content equal; `ext/` and `Project.toml` change it; no `ext/` defined |
| M6 | "M6: panel-derived values are part of the chunk hash" | control exponent changed by 1e-3, stability key removed (hyperbolic with positive exponent, elliptic with zero), panel file with one comment added | item hashes equal, chunk hash of the chunks holding the control differs; no-stability hyperbolic control is the same plan; elliptic one refused; kept files of another panel file are stale |
| M7 | "M7: the thread count must be the thread count of Julia" | `_check_threads` with k19 equal to `nthreads()` and to `nthreads() + 1`; `run_campaign` with `nthreads() + 1` | accepted when equal; refusal naming `Threads.nthreads()` and `thread_override`; override accepted; nothing written |
| M8 | "M8: a cost-cap pause is on disk ..." | clock at 2.01 times the line, killed before the first chunk file of the paused group | manifest says `paused: cost cap`; the resumed run does not start the arm |
| M9 | "M9: chunk files record the dirty tree and the session" | clean scratch repository, killed after 2 chunks, dirty tree, resume | files 1, 2: `dirty = false`, `session = 1`; the others `true`, `2`; sessions carry `index`; the schema table lists the three new attributes |
| LOW | "LOW: ..." | free device memory 0; arm with precisions Double64 and Float64; sha_start; stop fraction at the start of a resumed session; unknown `--cap-override` through the module and the CLI; L6 at Double64 | out-of-memory fallback (reason 4); pause by the Float64 group's line (200 ns); hashes equal the plan's; `stopped: failed chunks` at once; refusal naming the arm, exit code 1; plan refuses L6 at Double64 |

### RED log (review fixes)

Command: the new file against the tree at 05792a3 (`git archive HEAD` unpacked in the scratchpad, the helper files of 05792a3 and the helper definitions of `ex_chunks.jl`, `ex_store.jl`, `ex_manifest.jl`, `ex_device.jl`, `ex_review.jl` inside one outer testset; Julia 1.13.2, 4 threads). The copy of `ex_review.jl` used there drops the arguments `manifest_reset = false` (they restate the new default and the old `run_campaign` does not know the keyword; the remaining uses are the M4 and H2 calls that need the option). Result: `Pass 17, Fail 22, Error 30, Total 69`; a testset stops at its first error, so the totals count what ran. By finding (Pass/Fail/Error): H1 1/1/3 (no `control_stop` attribute, the plain resume is not refused, `retry_control_stop` is an unsupported keyword), H2 2/2/1 (no refusal, `backend_override` unsupported), M1 0/0/1 (`_move_retry` undefined), M2 2/2/3 (no attribute comparison, `_attr_diff` and `_fsync` undefined), M3 1/1/0 (the manifest wall time stops at the group boundary), M4 1/2/1 (no refusal, no backup), M5 2/7/0 (every file under `src/` enters the version, `ext/` and `Project.toml` do not), M6 2/3/0 (chunk hash blind to the exponent, a control without stability key accepted, a changed panel file kept), M7 2/2/0 (`nthreads()` never checked), M8 2/2/0 (pause not on disk), M9 2/0/20 (attributes missing), LOW 0/0/1 (`csize = 0` reaches `partition`). After the implementation: 108 of 108.


## Amendment D56 (preregistration A5 and runner support, 2026-10-19)

Stage 02 role first: `test/experiments/ex_a5.jl` (6 testsets, 48 assertions, registered as the last include of `runtests.jl`), then the code. Spec: `experiments-spec.md` Amendment D56, criteria D56a to D56g. Three existing places change because the real preregistration now has two run blocks: `ex_real_run()` (ex_run_inputs.jl) returns the A4 and A5 run tables merged by an independent merge written in the helper; D53r asserts the budget values of A5 (207.6, 6634.6) instead of 33.09 and 1553.8; D53b (ex_plan.jl) asserts the 9 methods as the sum of the arms `main` and `main_l2` (and `long`, `long_l2`) instead of one arm. No other assertion changes.

New interface under test: `merge_run`, `read_run_block` over several blocks; arm key `stop_exempt`; `[run.health]` key `control_scope`; work-item field `control_exact`; `breakdown_of`, `chunk_exempt`; orbit flag `breakdown_tau` with attributes `breakdown_step`, `breakdown_kind`; orbit attributes `control_scope`, `control_deviation`; manifest chunk entry keys `control_info` and the flag list entry `breakdown_tau`.

| Criterion | Test (ex_a5.jl: testset) | Sample | Pass |
|-----------|--------------------------|--------|------|
| D56a | "merge of the run tables of successive amendments" | fixture A4 block plus a partial A5 block (scope, one budget key, a replaced arm, a new arm) and three invalid A5 blocks | block hash equals the hand-merged table; arm order, replaced arm without its old `tau` list, `ns_gpu` kept; one-block file has no new keys; refusals name `control_scope`, `bogus`, `stop_exempt` |
| D56b | "the real preregistration: A5 over A4" | real file | 5 amendments, header hash equals computed; L2 only wrapped; main and long without L2; main_l2, long_l2; only tau_sweep exempt; scope; budget 207.6, 6634.6, 296311.3 with safety 2, 3, 3 and `ns_gpu` 2.468; panel hash of A4 accepted; equals the independent merge |
| D56c | "exactness of a control start is decided at plan time" | small panel with orbit 13 changed to a hyperbolic control at (pi, 0); a control at x0 = 0; a Newton control with winding (1, 1) at K = 10 | (0, 0) and Y5 true, (pi, 0) false, winding false, non-controls false; chunk hash changes with the flag |
| D56d | "control scope: the bound for exact starts, INFO for the others" | three campaigns on the small fixture | without the key the campaign stops at orbit 13; with `exact_start` it completes, orbit 13 `info` (deviation above 0.5), orbit 14 `bound` (below 0.05), manifest lists (13, L1) and (13, L5), also after a resume; a bump of orbit 14 stops |
| D56e | "orthonormalisation breakdown in an exempt arm" | one orbit at K = 10, L5, tau = 16, horizon 1600 | exempt: chunk ok, flag, step 48, kind, rows blank from the failure; not exempt: chunk failed; `breakdown_of` on three texts |
| D56f | "an exempt arm is outside the failed-chunk stop rule" | NaN injected in the L4 chunk of arm tau, stop fraction 0.01 | not exempt: `stopped: failed chunks`; exempt: `complete with failed chunks`, chunk still failed; a chunk shared with `main` is not exempt |

### RED log (D56)

Command: the new file against the tree at d787d25 (`git archive HEAD` in the scratchpad, with the new `ex_a5.jl` and `ex_run_inputs.jl`; the old tree has the A4-only preregistration and none of the code; Julia 1.13.2, 4 threads). Result: `Pass 4, Fail 6, Error 6, Total 16` (a testset stops at its first error): D56a 2/6/1, D56b 0/0/1 (no `main_l2` arm in the file), D56c 0/0/1 (`control_exact` missing), D56d 0/0/1 (key refused), D56e 0/0/1 (`stop_exempt` unknown key), D56f 2/0/1 (`chunk_exempt` undefined). After the implementation: 48 of 48. The control-exactness test first used the torus residual of the first return and marked control 881 exact; the plan-only check on the real panel showed it must be the lifted bit-for-bit return (the lifted state moves by the binary64 2 pi), which is the rule now in the spec and in the winding assertion of D56c.

## Amendment D57 (per-method control bound)

Spec: Amendment D57. New file `test/experiments/ex_a6.jl` (registered after `ex_a5.jl`); `ex_a5.jl` D56b reads the amendment count as at least 5 (A6 adds a sixth).

| ID | Property | Input | Expected |
|----|----------|-------|----------|
| D57a | block key read, merged and validated | single block, A4 plus a second block with the table, five bad tables | table kept; second block merges; a block without the key has none; refusals name `health.control_bound_method` and the method |
| D57b | method bound replaces `control_bound` for its method | exact control 14, deviation +0.01 injected in the L5 control chunk, table `L5 = 1e-3`; same without the table; the L1 chunk with the table; table `L5 = 2.0` with +1.0 | `stopped: control bound` naming method L5 and `above 0.001`; `complete`; `complete`; `complete` |
| D57c | table enters the chunk hash | plan with and without the table | hashes differ with the table, equal to the plan without the key when absent |
| D57d | real preregistration | `p1-le-benchmark-preregistration.md` | `check_prereg` accepts 6 amendments; `control_bound` 1e-8, `control_bound_method` `{L9 = 5e-8}`, scope `exact_start` |
| D57e | [local] real horizon | verify section CTRL | stage 04, not in the package suite |

## Amendment D58 (--prior-record)

Spec: Amendment D58. New file `test/experiments/ex_prior.jl` (registered after `ex_a6.jl`).

| ID | Property | Input | Expected |
|----|----------|-------|----------|
| D58a | record reader | fixture files in the format of `data/gpu-benchmark.toml`: the 4-thread record, a file with a per-method `recorder` table, a missing file, a text file, a file without the device, a device without Float64 | fields as in the file, empty recorder without a table, recorder entries with a table; four refusals naming the path or the device; `gpu_record` of a root without a file is `nothing` |
| D58b | option | `EX._cli_common` with and without the option and a given device name; `EX.cli(["pilot", ... "--prior-record", f])` in the test environment (no device) | `prior_record` key with threads 4 only when given; usage names the option; exit status 1 and a message naming `--prior-record` and the device |
| D58c | manifest | campaign with the record read from the fixture and `EX_REC` as device record; one without | `prior_record` has threads 4, cpu 8.254622999999999, date 2026-10-17 and `device_record` is `EX_REC`; without the option both stay as before |

## Amendment D59 (manifest backups)

Spec: Amendment D59. New file `test/experiments/ex_backup.jl` (registered after `ex_prior.jl`).

| ID | Property | Input | Expected |
|----|----------|-------|----------|
| D59a | retry | copier failing twice with EBUSY, EACCES, EPERM then copying; writer failing twice then writing | backups equal those of an undisturbed rotation; the file is written |
| D59b | best effort | copier always failing with EACCES; copier failing with ENOSPC | file replaced with the new content, older backups kept, one failure text per skipped copy; the ENOSPC error is raised |
| D59c | manifest | small campaign with the failing copier installed in `EX.BACKUP_COPIER`; one without | status `complete`, `backup_failures` non-empty; no such key without failures |
| D59d | no change | the tests of M1 in `ex_review.jl` | unchanged |


## Amendment D86 (stage 02, 2026-10-28)

Spec: `experiments-spec.md`, "Amendment D86" (E21 to E31, R11 to R21, I18 to I27, O16 to O22, criteria D86a to D86z, P26 to P39 approved). New files in `test/experiments/`, registered at the end of `runtests.jl` under the comment "Amendment D86": `ex_d86_helpers.jl` (fixtures, independent references), `ex_d86_grammar.jl` (D86a to D86d), `ex_d86_chains.jl` (D86e to D86h), `ex_d86_windows.jl` (D86i to D86m), `ex_d86_run.jl` (D86n to D86r), `ex_d86_store.jl` (D86s to D86v), `ex_d86_guards.jl` (D86w, D86x), `ex_d86_golden.jl` (golden values of D86b, generated, see below). One existing line changes: `allowed_internal` of "AC28 layering of src/Experiments" in `ex_inputs.jl` gains "Indicators" (D86x says Indicators is a declared import of Experiments; the change is harmless before the code and necessary after). Nothing under `src/`, `schema/`, `_config/` or the spec is touched.

Revision of the same day: the owner approved the corrections 1 to 7 appended to the spec ("Corrections at the D86 stage 02 checkpoint"); the tests and this section follow them (the items A1 to A7, A12, A16, A17, A19, A20 below are marked "resolved by correction N"). The corrections changed the fixtures (fraction 0.1 back in the grid), the D86a and D86d tests, and added two testset lines; the helper file now memoizes the reference runs and the campaigns are shared.

Fixtures. All values are the suite's own (no study value). The panel is `ex_panel_orbits(:std)`: 12 drawn orbits in 3 cells, controls 13 (elliptic fixed point), 14 (hyperbolic fixed point at (0, 0)), 15 (Y5), and two extra starts (16, 17). S_D = ids 1 2 5 6 9 10 13 14 15, S_L = id 1. The grid of the D86 blocks is `points_per_decade 5, first_steps 10, extra_fractions [0.1, 0.5]`; the fraction 0.1 is allowed for an arm with indicator specs only whatever the run length (correction 1), but still refused for an arm that names a catalog method, so the L1 twins of the window kinds run at a horizon with N / 10 a multiple of tau; the worked numbers of R17 use the preregistration grid (20, 100, [0.1, 0.5]). Window kinds: windows [16, 24, 40] for frequency and weighted Birkhoff (run length 80, with L1 twins at the horizon 80), [16, 24, 32] for phase advance (run length 65, no twin), [16, 32] in the all-kinds block (64, 65), tau 1 (tau 2 and 4 for the tangent chain with its L1 twins at the horizon 80); the K = 0 blocks use [32, 64] (frequency, phase advance) and [1000] (weighted Birkhoff, run length 2000; the one check that keeps the window 1000, correction 2). The panel has the 3 controls of the existing fixture (correction 6a). Spec dictionaries come from `d86_spec(kind; ...)` and every key is written (I27). Seeds: `d86_seed(letter, sub) = 20261012 + 29000 + 50 (letter - 'a' + 1) + sub` (sub < 50); testset D86<letter> draws only from `StableRNG(d86_seed(letter, k))` (D86d 1 items, D86f 1 directions, D86g 1 grid steps, D86h 1, D86i 1, D86k 1 K = 0 starts, D86m 1, D86t 1 kill points; helper starts `d86_seed('z', k)`). No reference value from memory: references are the stand-alone Indicators reducers and pure functions run on the identical orbit (`d86_stand_alone`, `d86_stretch_final`, `d86_pair_ref`, `d86_ref_freq`, `d86_ref_pa`, `d86_ref_bk`), a plain-array iteration of the map (`d86_points_T`), the independent enumeration of the items (`d86_ref_keys`), the independent bytes formula of E28 (`d86_ref_bytes`), independent text and hash of items, specs and chunk keys (`d86_item_text`, `d86_spec_text`, `d86_chunk_key`, hashed by the platform tool), BigFloat(256) closed forms, and hand counts derived in comments. The reference helpers were run in the scratchpad on the pre-amendment code against fabricated items: the independent bytes formula reproduces the worked numbers of R17 item 4 (18 864; 136 000 and, in Double64, B = 256 000; 16 896; 8 544; capacities 424 088, 58 823, 473 484, 936 329 at 16 GB and f = 1/2), and the independent enumeration of the all-kinds block gives 157 spec items with per-spec items and steps tg 38 / 2 756, st 36 / 2 628, prd 36 / 2 628, pr 18 / 1 476, fq 10 / 640, bk 10 / 640, pa 9 / 585, which are the literals asserted in D86c.

D86b golden values (code version before the amendment, commit 38eee93, fixture `ex_fix(:small)`, `EX_CODE`): `ex_d86_golden.jl` holds the block hash, the plan hash, 52 item hashes, 18 chunk keys and hashes, and the SHA-256 of the canonical read-back of each of the 18 chunk files (all datasets and root attributes except created_utc, cpu, threads, session, chunk_size, git_sha, dirty, julia_version, packages). It was produced by a scratchpad script: `include` of `ex_helpers.jl`, `ex_run_helpers.jl`, `ex_d86_helpers.jl`; `ex_fix(:small)`; `plan_run`; two campaigns (threads 4 / capacity 64 and threads 2 / capacity 7 / reverse order, which gave equal file hashes); then a print of the named tuple. The values are numbers of this tree and platform (they are a regression guard, not a published reference).

### Criterion to test

"@test" is the number of static `@test` statements of the testset; loops over orbits, cases and files multiply it (a green run reports more assertions than the RED log).

| Criterion | Tag | File: testset | Asserted | @test |
|-----------|-----|---------------|----------|-------|
| D86a | S | grammar: "D86a valid blocks read and give the canonical text" | the all-kinds block and one block per kind are read; canonical text lines and hash (platform tool); the arm `indicators` line; an arm without `methods` and an arm with `methods = []` (correction 7); key order of the file does not matter | 12 |
| D86a | S | grammar: "D86a the fraction 0.1: exempt for an indicator-only arm ..." | the horizon 8001 of phase advance with the windows 250, 1000, 4000 and the grid fraction 0.1 is accepted; an indicator-only arm with N / 10 not a multiple of tau is accepted; an arm with the catalog method L1 and the horizon 8001 is still refused naming `extra_fractions` (a guard that already passes), and so is an arm with a method and a spec together | 5 |
| D86a | S | grammar: "D86a malformed blocks are refused with the key path" | 69 malformed blocks (at least 2 per kind, covering every refusal class of R11 item 6, correction 6b; 14 common, 55 per kind or arm), each refused before computation with a message that holds the key path and the offending text; 6 file-level refusals through the campaign write 0 files and no manifest | 5 (+69 cases, +6 campaigns) |
| D86a | S | grammar: "D86a merge of two run blocks by id ..." | a replaced spec keeps its place and is replaced as a whole, a new spec is appended, the arm list merges (unknown keys and `control_methods` are in the malformed list) | 5 |
| D86b | S | grammar: "D86b absence of change ..." | block text and hash, plan hash and text hash, 52 item hashes, 18 chunk keys and hashes, 18 canonical file hashes equal the golden values; no `kind` in old items | 8 |
| D86b | [full] | not here | the real P1 plan hash and the 1363 chunk hashes: stage 04 | 0 |
| D86c | S | grammar: "D86c plan of the fixture block with all six kinds" | 159 items (157 spec items plus L1, L5), items per arm before the merge, direction multiplication (3 for tangent, stretching, pair with direction source; 1 for window kinds and seeded pair), omitted list `{fq: [15], bk: [15]}` and its presence in the plan text, merged item lists, order, per-spec step sums and totals (11 553 after, 12 553 before, saving 1 000), repeated planning for 1, 4 and all threads and reversed key order, chunk keys of the second form, chunk hashes | 27 |
| D86c | S | grammar: "D86c refusals ..." | window arm with two horizons or two tau refused; tau 2 with horizon tau L accepted; a subset in which the spec is undefined for every orbit is an error naming the spec | 5 |
| D86d | S | grammar: "D86d hashes ..." | one case per hashed field for each of 7 specs (more than 24 cases), 24 parameter edits each changing the item hash and the spec hash; budget, key order, arm name, backend request, thread count do not change them; old chunk key unchanged | 13 (+ about 150 in loops) |
| D86d | S | grammar: "D86d independent hashing (O19)" | 1000 items of each of 6 kinds built as NamedTuples from seeded values of every field of E27 (not plan items, correction 6c) and passed to `item_text` and `item_hash`; text equals the independent text and the hash equals `sha256sum`; spec hash and text from `EX.spec_*` equal the independent ones | 14 (+ 12 000 in loops) |
| D86e | S | chains: "D86e chains ..." | forms of all 13 chains, attributes `chain`, `members`, `spec`, `spec_sha256` of every file, a provenance group per member with no overwrite, package read-back equal; mutation check | 12 |
| D86f | S | chains: "D86f checks before the stream" | (a) dependent frame refused naming "initial frame" with 0 files; 10 000 seeded directions pass the frame checks; (b) unequal pair block frame fails its chunk naming the pair id; (c) plan-time refusals named by spec id and orbit id; run of the wrong length; Indicators-level index check; (d) separation outside [Dx/2, 2Dx] fails the chunk, the first rerun repeats it | 11 |
| D86g | S | chains: "D86g tangent family ..." | 7 tangent chains, lifted and wrapped, tau 1, 2, 4, Float64 and Double64: `t` equals the grid; columns equal E24; stored values at 5 seeded grid steps and the last equal the stand-alone member bitwise; orbit and monitor columns equal those of the L1 item; FLI nondecreasing; SALI in [0, sqrt 2 (1 + 2 eps)]; finite; no `lyap` | 3 (+ over 1 000 in loops) |
| D86g | S | chains: "D86g closed forms ..." (O18) | K = 0 shear, direction 0: FLI_1 and MEGNO column 1 exactly 0; FLI_B, SALI and MEGNO column 2 against their BigFloat closed forms within the bounds of the cited criteria | 5 (+ loops over the grid) |
| D86h | S | chains: "D86h pair ..." | 5 pair specs (RLI; RLI, MEGNO; MEGNO, RLI; direction source), both modes, 9 orbits: columns equal the stand-alone pair run bitwise at 6 steps; Dx_r; BigFloat(256) check (lifted); attributes `Dx`, `Dx_r`, `mode`, `xi0`, `frame_source`, `direction_rule`; xi0 unit and equal to the frame rule for the direction source; orbit equals the L1 item; RLI = 0 for Y5, lifted and wrapped bitwise equal; seeded and direction sources give different xi0 | 4 (+ loops, 90 items) |
| D86h | S | chains: "D86h pair on the K = 0 shear ..." | RLI and DL exactly 0 at every grid step, wrapped and lifted bitwise equal | 3 |
| D86i | S | windows: "D86i stretching: series, scalars, counts, ..." | D at 6 grid steps and the 18 orbit scalars equal the stand-alone final with workspace bitwise; range equals `stretching_range`; counts only for S_L, uint32 (bins+2) x 2, columns sum to the handoff count, equal stand-alone; pooled counts equal the exact integer sum over the chunk and are unchanged for reversed chunk order; an orbit recorded alone equals the same orbit in the group | 4 (+ loops) |
| D86i | S | windows: "D86i stretching closed forms ..." (O18) | Y5 in its eigenframe: single-cell counts (SN8b) and D = sqrt(2/h) within 4 eps; shear: column 1 in the central cell, D equals SN5 from the stored counts and lies in [0, sqrt(2/h)] | 4 |
| D86j | S | windows: "D86j frequency feed ..." (naff-fma criterion 16) | (a) manifest capacity equals floor(f free / bytes(orbit)) for 3 free values, B = 32 / 64 W_max, halving; (b) 96 bytes per window; (c) resident bytes after a chunk (full GC) at most the stored results plus 64 kB; (d) table equals `frequency_map` on the independent series, 2 observables, 3 windows, Float64 and Double64, flags decoded; (e) orbit points equal the L1 item; (f) no dataset holds the series | 11 |
| D86k | S | windows: "D86k frequency values ..." | (correction 6d: no Y5 rotation block; the Y5 origin start covers `undefined window` in D86l) K = 0 map, p_0 in {0.1, 1, 2, 3.3, 5, 6.1, pi, ...}: nu within 16 eps of p_0 / (2 pi) mod 1 (BigFloat), below_floor set; hyperbolic fixed point gives `degenerate` and an ok chunk counted in the manifest; island orbit: flags equal the pure function and the fixture shows mirror and near_zero; injected NaN sample fails the chunk naming orbit and sample | 4 |
| D86l | S | windows: "D86l phase advance ..." | table equals the stand-alone reducer bitwise (planar and lifted) and the pure functions on the independent series (planar); scalars `handoffs`, `near_pi`, `complete`; Y5 origin gives `undefined window`; K = 0 lifted observable gives p_0 / (2 pi) within 16 eps in both modes; wrong run length refused (D86f) | 3 |
| D86m | S | windows: "D86m weighted Birkhoff ..." | 31 accepted, 32 refused; K = 4, H = 3 table equals the pure functions bitwise; K = 31, H = 1 at 20 seeded starts; absdig and reldig bounds; constant observable gives delta 0 bitwise and below_floor; S = 0 gives zero and below_floor; K = 0 map (criterion 3, window 1000); cos / sin in lifted mode refused; observable names | 5 |
| D86n | S | run: "D86n sizes of the reducers and the worked numbers of R17" | Float64 sizes SALI 72, FLI 96, MEGNO 200, RLI 88, stretching 104 + 8 016, phase advance 96 + 64 per extra horizon, weighted Birkhoff 16 + 8 H (2 + 5 K), monitor 24, M 112 / 80 / 176; all four worked examples of R17 item 4 including (d) phase advance with three windows at the run length 8001: 141 grid steps, R = 248, bytes 8 224, capacity 972 762 | 16 |
| D86n | S | run: "D86n capacity equals the formula ..." | capacity equals the formula for the all-kinds plan at 3 free values and for each kind alone; chunk keys and hashes do not depend on the free memory; one orbit that does not fit refuses the start naming spec, bytes and free memory; stored data equal for capacities 1, 7, all and reversed order | 7 |
| D86o | S | run: "D86o backend matrix ..." | 7 specs x 3 precisions x 3 requests x 4 device states x 6 record states resolve as R16: sixth reason for stretching and frequency whatever the request, seventh reason for no record, other thread count, not faster, other spec hash (with the note), no tenth reason; backend table rows per (precision, spec id) | 5 |
| D86p | S | run: "D86p device path: records equal the CPU loop ..." | tangent, pair, phase_advance, birkhoff run types on the KA CPU backend: records bit-identical to the CPU loop; stretching and frequency never reach the device | 3 |
| D86p | S | run: "D86p injected faults ..." | out-of-memory and compile faults send groups to the CPU with reasons 4 and 3 and the data of a CPU run; the device batch size uses the real chain bytes (392 against the 256 allowance) | 6 |
| D86p | G | not here | local GPU, O21, GPU tag: stage 04 | 0 |
| D86q | S | run: "D86q determinism ..." | all-kinds block, Float64 and Double64: canonical read-back identical for 1, 2, 4 and all threads, capacities 1, 7, 64 and all, plan and reverse order, file by file; manifest thread count | 6 (+ files) |
| D86r | S | run: "D86r health: injected cases" | NaN in a tangent column fails; Inf in a window table fails naming orbit and column; NaN with its flag passes (degenerate, undefined window, reldig_undefined, zero) and is counted; without the flag it fails; complete false fails; status failed fails; all orbits below_floor passes; monitor exceedance 2 of 40 passes, 3 fails, flagged in both; CPU rerun after a device fault gives the data of a CPU run; the control rule does not read indicator items | 17 |
| D86s | S | store: "D86s store layout ..." | root attributes, dtypes, shapes and attribute lists of every dataset of every kind from the layout listing and the raw file; `lyap` absent; counts only for the subset; pooled group; pair orbit attributes; a version 2 old file reads; round trip | 5 (+ loops over files) |
| D86s | V | not here | independent reader of O20: stage 04 | 0 |
| D86t | S | store: "D86t atomic write killed at 20 seeded points ..." | 20 seeded kill points (5 write hooks) leave no incomplete final file and the resumed store equals the uninterrupted one; 10 seeded k with another thread count, canonical read-back and pooled counts equal; a changed spec parameter makes the old files stale; a final failed pair chunk is not retried | 13 (+ loops) |
| D86u | S | store: "D86u manifest ..." | `indicators` table (kind, spec, hash, chain, windows, run length, frame source, direction rule, observable names, budget, bytes terms), backend rows, arm cost lines, chunk `indicator_flags` recounted from the files, `not_applicable` lists, `transform_seconds`; a manifest without specs has the keys of R7 | 13 |
| D86v | S | store: "D86v pilot at the own run length ..." | pilot rows per spec with horizon (own run length for window kinds, pilot horizon for the others), bytes, capacity, flag shares, orbit pass and transform seconds for frequency; "not study data"; pilot files marked; cost cap on the spec line (2.01 pauses, 1.99 continues); CPU rerun bitwise; a simulated-device item is reported with its difference and not gated | 11 |
| D86w | S | guards: "D86w the suite runs with two fixture blocks ..." | generic checks (items equal the enumeration; tangent, stretching, pair, frequency, phase advance and birkhoff values equal the references) on the all-kinds block and on an alternative block whose specs differ in every value; no threshold, class or detection name in any new dataset or attribute; unknown keys refused | 6 (+ generic checks) |
| D86w | V | not here | `verify/` script that searches the source for the p2 values: stage 04 | 0 |
| D86x | S | guards: "D86x layering and entry point" | declared imports of Experiments include Indicators and nothing outside the layers and dependencies; nothing below imports Experiments; the usage text lists no new option; `run_cli.jl` holds no indicator or numerics name; every stop status has a nonzero exit status; a refused indicator block exits 1 naming the key | 8 |
| D86y | V, G, [full] | not here | stage 04 | 0 |
| D86z | stage 05 | not here | document checks at stage 05 | 0 |

### API assumed by the D86 tests

Stage 03 may adopt or rename these; the tests then change with them. Names with a `d86_` prefix are test helpers and are not API.

| Assumed name | Used in | Contract assumed |
|--------------|---------|------------------|
| `[[run.indicator]]` as `rb.raw["indicator"]` (list of normalized tables), arm key `indicators`, arm `methods` may be absent | a, c, all | as R11; the message of a refusal holds `indicator[i].<key>` or `arm[i].<key>` and the offending text |
| `EX.spec_text(rb, id)`, `EX.spec_hash(rb, id)` | d | canonical text of the spec without `budget`: lines `indicator.<key> = <value>`, keys sorted, lists in order, tables of a list as `key[i].name`, values as in `block_text`; the hash is its SHA-256 |
| plan item fields `kind`, `spec` (spec hash), `method` = spec id, `direction === nothing` for window kinds and seeded pairs, `eps === nothing`; `item_text` of E27 | c, d | `direction` printed only when not `nothing`, then `kind=` and `spec=` after `hash_version` |
| `plan.not_applicable::Dict{String, Vector{Int}}`; the plan text holds the list | c, w | omitted orbit ids per spec id |
| `plan.steps.specs[id]` with fields `items`, `steps` | c | unique items and horizon steps per spec |
| chunk key of a new kind | c, d | `method=<id> kind=<kind> spec=<16 hex> precision=... horizon=... tau=... eps=- direction=<d or -> mode=... cell=...` |
| orbit group layout | e to v | `orbits/<id>/indicators/<spec>/{series, final, flags, scalars, counts}`; attributes `columns` (series, final, scalars), `rows` (final), `bits` (flags); `pooled/<spec>/counts` with attribute `orbits`; orbit attributes `Dx`, `Dx_r`, `mode`, `xi0`; `provenance/<member>` as a group with the member's attributes; `scalars@range` for stretching |
| root attributes | e, s | `kind`, `method`, `spec`, `spec_sha256`, `chain` (strings such as `"FLI(replace=false)"`), `members`, `windows`, `run_length`, `frame_source`, `direction_rule` (`"u = xi_0"`), `separation`, `observable_names`; `schema_version` 2 |
| column names | e, g, h, i, s | `<member>.<field>`: `SALI.sali`, `SALI.L`, `FLI.fli`, `FLI.fli1`, `MEGNO.Y1`, `MEGNO.Ybar1`, `MEGNO.Ybar20_1`, `MEGNO.Y2` ..., `RLI.DL`, `RLI.RLI`, `stretching.D`; window tables as E25 (`nu1` ... `D`; `avg1` ... `dig`; scalars `handoffs`, `near_pi`, `complete` / `handoffs`, `complete`) |
| flag bit order | j to m, r | the order of E25 (`below_floor`, `degenerate`, `failed`, `tie1`, `tie2`, `near_zero1`, `near_zero2`, `mirror1`, `mirror2`, `refined1`, `refined2`; `below_floor`, `undefined window`; `below_floor`, `zero`, `reldig_undefined`); bit i is 2^(i-1), listed in `flags@bits` |
| stretching scalars | i | 18 values in `scalars` (handoffs; sum, comp, sum of squares, comp, mean, std, below, above for both columns; the first grid step with min_per_cell x bins handoffs); the test looks the values up, it does not assume their order |
| in-memory orbit record | r, t, v | `OrbitRecord.indicators::Dict{String, NamedTuple}` with mutable `series`, `final`, `flags`, `scalars`, `counts`; the `mutate` hook may edit them |
| `RunHooks(frame = (kind, info, W) -> W)`, kinds `:frame` and `:pair_frame` | f | a test seam applied to the frames before the frame checks |
| `RunHooks(series = (info, orbit, buffer) -> ...)` | k | called after the orbit pass of an orbit filled its series buffer (complex vector of 2 W_max entries), before the transforms |
| `EX.resolve_run_backend(request, T, id; threads, record, kind, spec_sha256)`; records `recorder[id] = (threads, gpu_faster, spec_sha256)` | o | CPU with the sixth reason for `stretching` and `frequency`; another `spec_sha256` gives the CPU, the seventh reason and a note containing "another chain" |
| `run_campaign(...; free_bytes, capacity = nothing)` reading the capacity of E28 (existing keywords); refusal as `ArgumentError` naming spec id, bytes(orbit) and the free memory | j, n | manifest `launch.capacity` |
| device batch size | p | `device_records(ctx, items).prov.chunk_size` from the real reducer bytes of the chain |
| manifest | u | `indicators[id]` with `kind`, `spec`, `spec_sha256`, `chain`, `windows`, `run_length`, `frame_source`, `direction_rule`, `observable_names`, `budget`, `bytes_per_orbit = (M, R, S, B, total)`; chunk entries `indicator_flags`, `not_applicable`, `transform_seconds`; `arms[i].specs[id] = (measured_ns_per_orbit_step, budget_line_ns)` |
| pilot report | v | cost rows for spec ids with `horizon`, `bytes_per_orbit`, `capacity`, `flag_shares`; `orbit_pass_seconds` and `transform_seconds` for frequency |
| `rerun_check` rows | v | field `gated::Bool` (false for device items) |
| existing internals used as in `ex_device.jl` | p, j, v | `EX.make_ctx`, `EX.compute_records`, `EX.device_records`, `ex_hdev`, `ExDev`, `ex_rec_same` |

### Mutation and negative cases

- D86a: 69 malformed blocks (the 45 required covering every class, at least two per kind): missing, malformed, duplicate and catalog-equal ids; unknown and missing kind; missing budget, budget without a precision of the arms, safety below 1, non-positive cost; keys named `threshold`, `class`, `detect`; per kind the refusal classes of R11 item 6 (members empty, unknown, repeated, SALI not last, key of another kind; separation zero, negative, NaN, missing; unknown frame source; pair member SALI; range with lo >= hi or one entry; bins; interval; interval not a multiple of tau; min_per_cell; counts not a subset; windows empty, not increasing, below the minimum, not integers; order, pad, mirror_ratio, floor_multiplier, observable, positions; exp_i_q and cos in lifted mode; birkhoff_order 0 and 5; 0 and 32 observables, unknown observable kind, frequency 0 or not an integer, coordinate 3); arm with neither methods nor indicators, with an empty list, with an unknown id; a spec used by no arm; window arm with two tau or another horizon; phase advance at the wrong run length; frequency window below 8 (order + 1); BigFloat; an id in `control_methods`. The fraction 0.1: accepted for an indicator-only arm at the horizon 8001, refused for an arm with a catalog method or with a method and a spec. A scratchpad run against the pre-amendment code confirmed that none of the 68 first cases passes for the wrong reason (all fail); the 69th was added afterwards and fails in the same way.
- D86c, D86d: window arm refused for horizon or two tau; every hashed field and parameter changed; the things that must not change the hash.
- D86e: sabotage chains computed in the test (SALI replacing and first; two replacing members) must change the values of the later members, so that the comparison of D86g catches them; the correct chain must not.
- D86f: dependent frame, unequal pair block frame, injected non-finite start, derived range for another system, Dx_r outside [Dx/2, 2 Dx], run of the wrong length.
- D86r: NaN, Inf, flagged and unflagged NaN, `complete` false, status failed, all-below_floor chunk, 2 and 3 exceeding orbits of 40, device fault, control rule.
- D86t: kills at the five write points, kills after k chunks, stale files, final failed pair chunk.
- D86v: 2.01 and 1.99 times the spec line.

### Not covered here

- D86b [full] (real P1 plan and 1363 chunk hashes), D86s (V) (O20), D86w (V) (script in `verify/`), D86p (G), D86y, D86z: stages 04 and 05, as the task states.
- D86f (c): "interval not a multiple of tau" is refused earlier, at block validation (D86a); "observable index beyond the dimension" cannot be written in the grammar (coordinates 1 and 2 only, the state has two coordinates), so only the Indicators-level `check_run` is asserted (it passes before the code); the run of the wrong length is tested through the horizon refusal.
- D86k: the Y5 rotation block of naff-fma criterion 7(a): the Y5 of the panel is the hyperbolic block (mu), not a rotation, and the panel format cannot express a rotation block (correction 6d, former A7). The Y5 origin start is used for the `undefined window` and `degenerate` bits instead. The elliptic fixed-point control is not exact in binary64 (sin(pi) is 1.2e-16 and p drifts), so only the hyperbolic one (exact) is used for constant series.
- D86h: Double64 pairs (the stored xi0 is two binary64, so the Double64 start cannot be rebuilt bitwise from the file); the seeded xi0 convention is not recomputed (A8).
- D86g: the Y5 closed forms of SALI, FLI and MEGNO and the second-column counts of the shear are not tested through the recorder; only the K = 0 shear forms of O18 and the Y5 and shear forms of stretching and RLI.
- D86i: the shear counts of the second column (closed form SN8c) are not recomputed; D is checked from the stored counts (SN5 algebra) and its bounds.
- D86n: the two launches (58 823 and 41 177) of worked example (b) are not asserted; example (d) is asserted for both weighted Birkhoff (8 544, 936 329) and phase advance (8 224, 972 762).
- D86p: only the four device kinds, Float64; the display-watchdog fault is not repeated for the new kinds.
- D86q: "alone and inside a group" is covered by capacity 1 against all orbits and, for stretching, by D86i.
- Suite time: the added time of at most 180 s on 4 threads cannot be measured at RED (the D86 testsets stop at their first error); measure at stage 03.
- tau in {2, 4} for the window kinds: only tau 1 (and tau 2 in the alternative block of D86w).

### RED log

Command (repository root, Julia 1.13.2, `-t 4 --check-bounds=yes`): the experiments folder only; `runtests.jl` is `include`d (with `@__DIR__` replaced by the path) inside one outer testset. "Before" is the same file cut at the line "Amendment D86" (the old D53 to D59 tests; the only change in them is the list of AC28); "after" is the full file. The scratchpad runner is `runner.jl` (reads `runtests.jl`, cuts or not, `include_string`). Final run, after the corrections of the checkpoint:

| Run | Pass | Fail | Error | Broken | Total | Time |
|-----|------|------|-------|--------|-------|------|
| before (old tests only) | 4 861 | 0 | 0 | 0 | 4 861 | 3 min 25 s |
| after (old tests plus D86) | 4 901 | 13 | 32 | 0 | 4 946 | 3 min 59 s (3 min 22 s to 4 min 12 s in the runs; the spread is machine noise) |

The first RED run (before the corrections) gave 4 899 / 12 / 30. The 40 passes of D86 are all guards that must hold before and after the code: D86b 8 (absence of change), D86a malformed 13 (the counts of cases and kinds, 6 "nothing written" checks of the file-level refusals), D86a fraction 0.1 2 (the grid holds the fraction; an arm with the catalog method L1 at the horizon 8001 is refused, which the existing code already does), D86n sizes 9 (the verified reducer sizes), D86f 2 (10 000 seeded frames pass the frame checks; no file written), D86x 6 (the declared imports are legal, nothing imports upward, no new option, no numerics in the entry point, stop statuses nonzero). Every old testset has the same pass count before and after, 4 861 in all, with 0 failures and 0 errors.

First error class of each D86 testset (all are missing functionality; 31 of the 32 errors are `ArgumentError: run block run.indicator: unknown key` from the block validation, one is a `MethodError` for the missing keyword):

| Testset | Outcome | First error |
|---------|---------|-------------|
| D86a valid blocks | Error | `ArgumentError`, unknown key `run.indicator` |
| D86a the fraction 0.1 | 2 Pass, 1 Fail, 2 Error | the accepted blocks error with the unknown key; the refusal of the mixed arm fails on the message |
| D86a malformed blocks | 7 Fail (13 Pass) | the list of cases not refused for the right reason holds all 69; 6 file-level refusals fail on the message |
| D86a merge | Error | `ArgumentError`, unknown key |
| D86b | 8 Pass | guard, passes before the code by design |
| D86c plan | Error | `ArgumentError`, unknown key |
| D86c refusals | 2 Fail, 1 Error | the message does not name `arm[i].horizons` and `tau`; then unknown key |
| D86d hashes; D86d independent hashing | Error; Error | unknown key |
| D86e | Error | unknown key |
| D86f | 1 Fail, 1 Error (2 Pass) | no refusal naming "initial frame" (the seam is absent); `MethodError: no method matching RunHooks(; clock, frame)` |
| D86g (2), D86h (2) | Error x 4 | unknown key |
| D86i (2), D86j, D86k, D86l, D86m | Error x 6 | unknown key |
| D86n sizes | 1 Error (9 Pass) | unknown key (the sizes of the verified reducers hold) |
| D86n capacity, D86o, D86p (2), D86q, D86r | Error x 6 | unknown key |
| D86s, D86t, D86u, D86v, D86w | Error x 5 | unknown key |
| D86x | 2 Fail (6 Pass) | `imports_indicators` is false; the refusal of the entry point reads "unknown key" instead of naming `indicator[1].members` |

No parse error, no error inside a helper, no `UndefVarError` of a test name. Defects of the tests themselves found and fixed earlier: a campaign that failed inside an `@test` followed by a read of a store that does not exist (reads are now guarded by `d86_files` and `d86_nofiles`, and campaigns are called outside `@test`); `pair_start` takes the lifted start (found by running the reference helper); the reference of the lifted phase advance uses the stand-alone reducer instead of a difference of rebuilt lifted values (different rounding); the golden file excludes `git_sha` and `dirty`. Defects caused by the return of the fraction 0.1: the L1 twins of the window kinds and of the tau 2 and 4 chain had to move to horizons with N / 10 a multiple of tau (80), the phase advance twins were dropped (run length 65 has no twin), and the D86d edit of `order` became 0 (the window 16 is below 8 (order + 1) for order 2).

Suite time. The D86 testsets sum to about 6 s in the RED runs (they stop at their first error), so the time of the green suite and the gain from the memoization cannot be measured here; measure at stage 03. The memoization (below) removes repeated reference runs and repeated campaigns, so the green time is expected to be lower than the same suite without it, but no number is claimed.

### Test economy (revision)

No assertion and no seed changed. In `ex_d86_helpers.jl`: the stand-alone reducer runs (`d86_stand_alone`, `d86_stretch_final`, `d86_pair_ref`), the pure-function references (`d86_ref_freq`, `d86_ref_pa`, `d86_ref_bk`) and the independent series (`d86_complex_series`) are memoized by the key (kind, orbit data, precision, mode, tau, direction, member or spec text, step); the raw read of a chunk file (`d86_raw`) is cached by (path, modification time, size); plans are cached by fixture tag (`d86_plan`). The all-kinds indicator-only campaign is run once (`d86_all_ind_run`, threads 4, capacity all) and shared by D86q, D86t, D86v and D86w; the chain and window campaigns are shared by D86e to D86h and D86j to D86m; D86p reads the plans of those campaigns instead of planning again.

### Resolution of the ambiguities

Resolved by the corrections 1 to 7 of the spec ("Corrections at the D86 stage 02 checkpoint"): the tests and comments follow them. The rest are test assumptions for stage 03 to confirm or change.

| Item | Status | Text |
|------|--------|------|
| A1 | resolved by correction 6a | 3 controls (elliptic fixed point, hyperbolic fixed point, Y5); the existing fixture is used and the hand counts use ids 13, 14, 15 |
| A2 | resolved by correction 1 | the 0.1 rule applies only to an arm with a catalog method; fixtures use `extra_fractions = [0.1, 0.5]`; phase advance at the run length 8001 has 141 grid steps, bytes 8 224, capacity 972 762, asserted |
| A3 | resolved by correction 7 | an arm with `indicators` may omit `methods` or give `methods = []`; both are tested |
| A4 | resolved by correction 7 | a window is at least 8 x (order + 1); a malformed block with order 2 and the window 16 is in the list |
| A5 | resolved by correction 3 | root path `indicator`, one line per key, keys sorted, 1-based list index, no `budget` |
| A6 | resolved by correction 4 | `provenance/<member>` is an HDF5 group, one attribute per key |
| A7 | resolved by correction 6d | no Y5 rotation block; exact fixed points and the K = 0 map; the Y5 origin start covers `undefined window` and `degenerate` |
| A8 | test assumption | the seed convention of the seeded xi_0 is not stated; the tests read the stored `xi0`, check that it is a unit vector, that the direction source equals the frame rule, and that the stand-alone pair run started from it gives the stored columns |
| A9 | test assumption | D86h: Dx_r is a rounded norm; one eps relative to the BigFloat norm in lifted mode, `pair_start` bitwise in both modes |
| A10 | test assumption | O18 / D86g: B_SALI(n) of sali A15 is taken as 8 eps Lambda_n + SALI_cf * 2 B_S + 4 eps with B_S(n) = 8 eps (n + 3) (G_first = 3 for the shear) and B_meas = the SALI floor 4 eps; the stage 04 B_meas is not known here |
| A11 | test assumption | D86k: the bound 16 eps of naff-fma criterion 7(b) is stated for N in {1000, 4096}; the fixture windows are 32 and 64 |
| A12 | resolved by correction 2 | the K = 0 weighted Birkhoff check uses the window 1000, birkhoff_order 1, run length 2000 |
| A13 | test assumption | E25: the bits `refined1` / `refined2` ("false only") are not compared; `undefined window` is assumed as the exact bit name |
| A14 | test assumption | the in-memory layout of an indicator result is not fixed; see the API table; D86r edits `final`, `flags`, `scalars`, `series` in place |
| A15 | test assumption | D86f (a), (b): the seam is `RunHooks(frame = ...)`; `ex_campaign` rebuilds hooks when `clock` is missing, so every D86 test with hooks also passes a clock |
| A16 | resolved by correction 6b | at least 45 blocks covering every refusal class, several per class allowed; 69 are used |
| A17 | resolved by correction 6c | D86d hashes 1000 items per kind built from seeded values of every field with `item_text` and `item_hash`; `item_text` accepts any value with the properties of E27, so no plan item is needed |
| A18 | test assumption | D86p capacity: free memory chosen so that 9 orbits fit with the old 256-byte allowance and 8 with the 392 bytes of FLI, MEGNO, SALI plus the monitor |
| A19 | resolved by correction 5 | the attribute `range` of `scalars`: two binary64, lo then hi (the test checks the element type) |
| A20 | resolved by correction 4 | the key `indicator` is not merged; three members show three values, each read from its own group |

New points found in this revision (test assumptions):

- B1. Correction 1 names "an arm with indicator specs only" as exempt; an arm that holds methods and specs is read as not exempt (the arm "names at least one catalog method"), and the test for the mixed arm assumes the refusal.
- B2. Because the 0.1 rule still holds for the L1 twins, a window kind with an odd run length (phase advance, 2 W + 1) has no L1 twin: the orbit points of the phase advance items are not compared with an L1 item (I21 for that kind rests on the stand-alone reducer).
- B3. The grid target round(0.1 N) is not a handoff when tau does not divide it; the tests' independent grid applies the rule of D53 R2 item 1 (the first handoff not before the target), as correction 1 states.


### Correction round (2026-10-28)

Source: `experiments-d86-impl-notes.md`, "Tests believed wrong or unreachable" (17 items), with "Deviations from the spec" read for context. The folder suite gave 6 763 pass, 29 fail, 18 error; the D86 part 1 902 / 29 / 18. Each item below was checked against the spec text and against the code (by running the test, a scratch script, or the Indicators source). In all 17 items the test was wrong or the fixture could not test the criterion, and none of them shows a deviation of the code from the spec, so none is left red and nothing goes back to the implementer. Two points of the corrections where the fix needed a basis beyond the criterion text are marked "basis".

| # | Test (old line) | Old claim and finding | Spec line | Fix |
|---|-----------------|-----------------------|-----------|-----|
| 1 | chains 142 | `check_run(WeightedPhaseAdvance(.., PlanarObservable(1, 3)), 1, 33; d = 2)` must throw. Code: the limit is `1 <= i <= 2d` (`naff_reducer.jl`: "observable indices must lie in 1:2d"), so 3 is valid for d = 2. The claim of the implementer holds. | D86f(c): "an observable index beyond the dimension is refused at plan time"; the state of a wrapped map has 2d entries | the refusal is asserted for index 5, and the acceptance of index 3 |
| 2 | chains 152-154 | D86f(d) depends on the seeded xi_0: with u = (-0.544, 0.839) the second coordinate separates by 8.4e-9, inside [Dx/2, 2 Dx]. The claim holds (A8). | D86f(d): "A pair with Dx_r outside [Dx / 2, 2 Dx] (injected start) fails its chunk with kind separation and the pair id; the first rerun repeats it (R18 item 5)" | the start has both lifted coordinates at 1e9 (ulp 1.2e-7 > 2 Dx = 2e-8): B rounds back to A in both coordinates whatever u is (Dx_r = 0); orbit id, "separation", attempt 2 and the repeated first-attempt reason are still asserted |
| 3 | chains 246 | MEGNO closed form used `prod j = 1..n`. Wrong. | megno-spec M10(a): "Y(n) = ln(n^2 + 1) - (1/n) ln prod_{i=0..n-1} (i^2 + 1)" | product over i = 0..n-1 |
| 4 | windows 57 | `t[findfirst(...)]` is `t[nothing]` for the tau = 4 items (25 handoffs against min_per_cell x bins = 40). The spec names the entry only when it exists. | E25: "the first grid step at which the handoff count reaches min_per_cell x bins" (silent when it never does; the code stores 0.0, deviation 6 of the notes) | that entry is compared only when the run reaches the count; the length 18 is still asserted |
| 5 | windows 107 | The "edge condition" fails for Y5 with the derived range: the range is the exact bound padded by 8 eps, so z = ln mu sits at its edge by construction. Two further defects found when this was fixed: the test used the default 20-bin spec for the reference counts, and the eigenframe (direction pi/4) cannot give single cells: the column along the contracting direction is rounded off it and the rounding grows like mu^k, so after about 20 handoffs the recorded column is the expanding one (the stand-alone reducer does the same). | SN8b: counts "equal those of the closed form ... tolerance 0 under the edge condition (otherwise within the certificate)", D "equal to the closed form within 4 eps" | given range [-2, 2] (allowed by correction 5), direction 0.3, the reference counts and D from the closed form z_k = ln(norm(A^k v) / norm(A^(k-1) v)) in BigFloat(256) with the margin to the cell edges checked (> 1e-9), the spec of the plan used for the reference |
| 6 | windows 149, 150 | `big(...)` used before the local `big` of the same testset. Test defect. | none | local renamed |
| 7 | windows 171 | The assertion needs a warm-up that includes the control orbits: the constant series allocate about 340 KB of one-time cache in the transform. The claim holds (growth of later calls is 41 KB < 103 KB). The conjunct "the buffers would exceed 64 kB" was false for the 8 items of the fixture (64 000 bytes). | D86j(c): "After a chunk the resident host bytes exceed those before it by at most the stored results plus 64 kB (the buffer is released; measured by the allocation counters after a full garbage collection)" | one full pass over the same items as warm-up, then the measured pass; the fixture has 40 orbits (320 kB of buffers if they were kept) |
| 8 | windows 208 | The loop over the keys of the file visited attributes (`final@columns` and others are vectors of strings). Test defect. | D86j(f): "No dataset of the chunk file holds the series" | attributes (keys with `@`) are skipped |
| 9 | windows 268, 288 | `f0` has `mean_removal = false`, so the constant series at the fixed point (0, 0) is not degenerate. The claim holds. | naff-fma 11(b): "A window is degenerate when its samples, after the optional removal of the window mean, are all equal, and either mean removal is on or the common value is zero ... A window of equal nonzero samples with mean removal off is a line at frequency 0: status ok" | `degenerate` is asserted for `fm` (mean removal on), a zero frequency with `near_zero1` and no `degenerate` for `f0`; the manifest count is read from the `fm` controls chunk |
| 10 | windows 310, 328 | Y5 in wrapped mode has a state of two components, so `LiftedAngleObservable(1, 3)` is out of range; `plan.panel.orbits[o]` fails for the extra starts 16 and 17. Test defects. | E26: lifted "is built from the state coordinate and its winding count in wrapped mode, or the lifted coordinate in lifted mode"; the torus-free control has no count | `LiftedAngleObservable(1)` for a linear system; `d86_pinfo` returns `nothing` for the extra starts |
| 11 | windows 342 | Lifted K = 0 map, orbit 6 (p_0 = 6.1): error 16.54 eps against 16 eps. The recorder equals the stand-alone reducer bitwise (the main D86l comparison passes for the lifted observable in both modes), so the recorder is right. Basis: naff-fma 7(b) gates 16 eps for the "wrapped state with winding counts ... lifted post-kick momentum rebuilt from the winding counts"; in the lifted mode the coordinate itself is rounded at each step and 7(c) says its error is "about |q| eps, not gated". | naff-fma 7(b), 7(c); D86l: "within the bound of naff-fma criterion 8" (criterion 8 gates the algebra, 7(b) the K = 0 answer) | wrapped mode 16 eps; lifted mode 16 eps + max abs(q) eps (q read from the stored state, up to 400); this is a widening with the basis of 7(c) only, stated here for the owner |
| 12 | windows 378 | `plan.panel.orbits[o]` for the extra starts. Test defect. | none | the unused line is removed |
| 13 | run 207, 212 | The canonical hash included `backend`, `backend_fallback_reason`, `device`, `gpu_chunks`, which differ by design. The claim holds; all datasets are bit-identical. | D86p: "records bit-identical to the CPU loop"; I24 (datasets and attributes independent of backend and capacity) | the comparison is strict on every dataset and on every attribute except `backend`, `backend_fallback_reason`, `device`, `gpu_chunks`, `attempt`, `first_attempt_reasons` (`D86_BACKEND_ATTRS`); I first added `steps_per_launch` and `launches_per_chunk` and removed them again: the strict comparison passes without them. Mutation check kept: a changed dataset value, and a changed attribute outside the list (`horizon`), change the hash; a changed `backend` does not |
| 14 | run 318 | same, plus the rerun attributes `attempt` 2 and `first_attempt_reasons` | D86r: "the CPU rerun after an injected device fault gives the data of a CPU run" | the same list; the rerun is asserted separately (attempt 2, backend cpu, reason 5, first-attempt reason recorded) |
| 15 | store 57 | `scalars` was read for every window kind; E25 gives none for `frequency`. The claim holds. | E25, table row `frequency`: "Orbit scalars: none (planes = 1)" | `scalars` is read for phase advance and Birkhoff only; absence is asserted for frequency |
| 16 | store 168 | iterating the pairs of a `Dict`. Test defect. | R19 | `values(...)` |
| 17 | store 228 | The arm `a` had one run type, so its four chunks formed one launch group and the pause could not leave a chunk unrun. The claim holds. The first rewrite also asked for the manifest status "pending"; the code (and D53l) records "paused". | D86v: "a first group at 2.01 times the line pauses the arm, at 1.99 continues, as D53l"; D53l: the later chunks are paused and not run | the arm has two run types (horizons 100 and 80); asserted: status "paused arms", later chunks "paused" and not in `ran`, arm status "paused: cost cap", the other arm complete, the first group ok; at 1.99 every chunk ok |

Speed (no assertion weakened, no seed changed). `ex_d86_grammar.jl`: the spec hash was computed with one external `sha256sum` per item (6 000 processes); it is now computed once per kind and memoized, and the 1 000 items per kind are hashed by one external process per kind (`ex_sha256_many`), the independent implementation staying external. `ex_d86_store.jl` D86t: the 20 kill points, the 10 resumes and the stale and failed-pair checks run on a fixture of the same six kinds, both precisions, the direction sweep and the pooled counts, on the orbit 1 (S_L) and the extra starts (about 23 chunks instead of 71); the uninterrupted reference is one campaign. The seeds are unchanged; because the number of chunks M changed, the same seeds now pick other kill points.

Results, D86 files only against the current `src/` (the loader evaluates the helper definitions of the old files and skips their testsets; each D86 file is one testset of an outer testset; Julia 1.13.2, `-t 4`, `--check-bounds=yes`):

| File | Before (pass / fail / error, s) | After (pass / fail / error, s) |
|------|--------------------------------|--------------------------------|
| `ex_d86_grammar.jl` | 468 / 0 / 0, 114.1 | 468 / 0 / 0, 33.0 |
| `ex_d86_chains.jl` | 833 / 20 / 0, 24.0 | 854 / 0 / 0, 21.1 |
| `ex_d86_windows.jl` | 15 / 4 / 16, 13.2 | 44 / 0 / 0, 17.1 |
| `ex_d86_run.jl` | 284 / 4 / 0, 40.2 | 292 / 0 / 0, 39.7 |
| `ex_d86_store.jl` | 277 / 1 / 2, 110.4 | 286 / 0 / 0, 43.9 |
| `ex_d86_guards.jl` | 25 / 0 / 0, 3.7 | 25 / 0 / 0, 3.3 |
| D86 total | 1 902 / 29 / 18, 305 | 1 969 / 0 / 0, 158 (2 min 36 s in the last run) |

The "before" counts reproduce the implementer's numbers (1 902 / 29 / 18). The first compilation of each campaign is inside these times. The windows file has fewer assertions than its code paths because most of its checks are collected in `bad` lists (one `@test isempty(bad)` per loop); the testsets now run to the end, which they did not before (the errors stopped them). Left red: none. The old tests were not rerun in this round (the coordinator runs the full suite).

New points (test assumptions): C1. Y5 with a symmetric hyperbolic block and the direction pi/4 cannot be used for SN8b in binary64 (item 5); the closed form is tested for a generic direction. C2. The lifted-mode bound of D86l (item 11) is wider than 16 eps by the scale of 7(c); if the owner wants the criterion 16 eps in the lifted mode, that cannot hold for a lifted coordinate of 400.
