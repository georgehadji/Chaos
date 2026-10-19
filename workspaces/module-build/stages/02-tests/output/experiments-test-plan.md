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
