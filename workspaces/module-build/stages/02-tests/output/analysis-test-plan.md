# Test plan: analysis

Spec: `workspaces/module-build/stages/01-spec/output/analysis-spec.md` (approved 2026-10-20 with P1 to P8 as proposed; decision D60). Suite: `test/analysis/` (`runtests.jl` includes one file per group). Stage 02, RED: the module `Chaoslab.Analysis` does not exist yet, so every test that calls it errors with `UndefVarError: Analysis not defined in Chaoslab`.

Seeds. Every random draw uses a StableRNGs seed from `an_seed(ac, sub) = 20261012 + 28000 + 50 ac + sub` (the repository seed of `_config/numerics-defaults.toml`, offset for this slug). Bootstrap streams are keyed by the label of the statistic and `seed` (E12); the tests that need a seed pass `AN_SEED = 20261012`.

No memory values. Regression numbers come from the preregistration (through the parameter table E0), `benchmarks.md` rows (SalaArtuso2015 0.298 and 1.172; Shevchenko2014 Eq. (25) 0.1333 K, Eq. (26) 0.469 (K - 1.037)^(1/2) and ln(K/2) + 1/K^2) or from closed forms (the integrable map, preregistration A1 item 2). Every other expected value is computed in the test from the inputs it builds.

## API under test

Names are chosen here (the spec has no Julia names, Pattern 10). `p` is the parameter table, a NamedTuple with the 37 keys of E0 (below). All names live in `Chaoslab.Analysis`; the tests reach them through `AN.<name>`.

| Name | Call | Returns |
|------|------|---------|
| `read_params` | `read_params(path)` | NamedTuple of the 37 keys; `ArgumentError` naming the key for a missing, unknown or mistyped key and the path for a missing file |
| `h1_orbit` | `h1_orbit(lam::Dict{String,Float64}, p)`, keys L1, L3, L4, L4b, L5, L5b, L6, L9, L2 (values at N_1) | `(delta_qb, diff_l9, diff_l2, pass_a, pass_b, pass_c)`; Delta_QB over L1 L3 L4 L4b L5 L5b L6 only; bounds inclusive |
| `decide_h1` | `decide_h1(orbits, p)`, orbits `(id, pass_a, pass_b, pass_c)` | `(verdict, n, fractions=(a,b,c), passes=(a,b,c))`, verdict `:supported`, `:falsified` or `:not_decidable` (n below `min_set`) |
| `control_deviations` | `control_deviations(lam, exact)` | `Dict{String,Float64}` of absolute deviations |
| `check_controls` | `check_controls(rows, p)`, rows `(id, exponent, scope, lam::Dict, stored::Dict or nothing)` | `(ok, defects=[(id, method, deviation, bound)], info=[(id, deviations)], mismatches=[ids])`; scope "bound" is gated per method (`control_bound`, L9 `control_bound_L9`), "info" never gates, exponent 0 is skipped; a stored deviation differing from the recomputed one by more than 1e-12 is a mismatch |
| `orbit_sets` | `orbit_sets(orbits, p)`, orbits `(id, K, class, control, exponent, flags, instability)` | `(E_ch, W, S, S_nonpositive, dropped=[(id, reason)])`; reasons "drift" and "excluded" |
| `relative_spread` | `relative_spread(lam::Dict{String,Vector}, lam_ref)` | `Vector` R(n) over the eight methods of V (L2 is ignored if present) |
| `decide_h2` | `decide_h2(R, p; seed)` | `(verdict, median, lo, hi, n)`, verdict `:supported`, `:not_supported`, `:not_decidable` |
| `convergence_time` | `convergence_time(steps, lam, lam_ref, band)` | `(n, defined)`; undefined gives `n` = last step; `lam_ref <= 0` is undefined |
| `orbit_convergence` | `orbit_convergence(steps, lam::Dict, lam_ref, band)` | `(n, defined)` maximum over V, undefined if any method is |
| `decide_h3` | `decide_h3(nconv_W, nconv_S, p; seed)` | `(verdict, ratio, lo, hi, n_w, n_s)`; the interval end is compared strictly (`>`), the ratio with `>=` |
| `stream_seed` | `stream_seed(label, seed)` | `UInt64`: first 8 bytes, big-endian, of SHA-256 of `chaoslab-analysis-1|<seed>|<label>` |
| `bootstrap_rng` | `bootstrap_rng(label, seed)` | `StableRNG` seeded with `stream_seed` |
| `bootstrap_interval` | `bootstrap_interval(values, label; seed, resamples, level)` | `(point, lo, hi)`, median, percentile type 7 |
| `bootstrap_ratio` | `bootstrap_ratio(w, s, label_w, label_s; seed, resamples, level)` | `(point, lo, hi)` for median(w) / median(s) with independent resamples |
| `decide_h4` | `decide_h4((halfpi=(L1,L5), twopi=(L1,L5)), p)` | `(verdict, halfpi=(pass,L1,L5), twopi=(...))`, verdict `:supported` or `:did_not_reproduce` |
| `shevchenko_fit`, `chirikov_line` | `shevchenko_fit(K)`, `chirikov_line(K)` | the fit of the K range or `nothing` (between 1 and 1.1); `log(K/2)` |
| `growth_ratio` | `growth_ratio(lam_N, lam_N10)` | `10 lam_N / lam_N10`, `nothing` if `lam_N10 <= 0` |
| `base_label` | `base_label(d, p)` | `(instability, structure)` before any test |
| `orbit_label` | `orbit_label(d, p)`, d `(lam_N, lam_N10, rho1, rho2, periodic, residual, estimator_fail, diag_exceeded, excluded, ladder, direction)` | `(instability, structure, forced::Vector{Symbol}, raw, ladder_same)`; forced symbols `:ladder :estimator :direction :diagnostics :excluded` in that order |
| `decide_h5` | `decide_h5(orbits, p)`, orbits `(id, in_SL, ladder_same, instability, structure, excluded)` | `(verdict, h5a, h5b, h5c)`, each sub-result `(fraction, n, pass::Union{Bool,Nothing}, verdict)`; excluded orbits are left out |
| `read_store` | `read_store(manifest, store)` | `Vector` of items `(orbit, method, precision, horizon, tau, eps, direction, mode, arms, t, lam, state, winding, flags, attrs, diag, cost)`; `lam` is the first exponent; `eps` is `nothing` except L2 |
| `find_item` | `find_item(items; orbit, method, precision, horizon, tau, eps, direction, mode)` | the item found by content, or `nothing` |
| `analyze` | `analyze(; prereg, manifest, panel, params, store, out, figures=false, threads)` | `(status, tables, decisions, sets, labels, controls, cells, digests, deviations, readings, header, defects)`; status `:ok`, `:figures_skipped` or `:control_defect` |
| `cli` | `cli(args; io)` | exit code `0`, `1` (refusal or usage), `2` (tables written, figures skipped) |

`analyze` result fields used by the tests. `decisions`: NamedTuple `H1 .. H5` in this order, or `nothing` after a control defect. `H3` also carries `undefined_share`; `H4` carries `crosscheck` (one `(float64, double64)` row per K, pi/2 first). `sets`: as `orbit_sets`. `labels`: `Dict{Int, ...}` of every panel orbit with `instability`, `structure`, `forced`, `raw`, `ladder_same`, `applicable` (the tests that could be applied: `:ladder` for S_L, `:direction` for S_D, `:estimator` and `:diagnostics` for all). `controls`: rows `(id, scope, deviations)`. `cells`: rows `(k_index, class, orbits, dropped, labels::Dict{Tuple,Int})` with a cell "controls". `header`: `prereg_header_sha256`, `panel_sha256` (plain), `manifest_sha256`, `params_sha256`. `digests`: `Dict("F1" .. "F7" => 64 hex)`. Refusals are `ArgumentError`; a control defect is a status, not an exception (a defect report `<study>-defect.md` is the only output).

Parameter keys (E0, flat TOML): `N_1 N_long N_B n_H2 h1a_bound h1a_fraction h1b_bound h1b_fraction h1c_bound h1c_fraction h1c_eps control_bound control_bound_L9 h2_median h2_interval bootstrap_resamples bootstrap_level h3_band h3_ratio h3_interval h4_ref_halfpi h4_tol_halfpi h4_ref_2pi h4_tol_2pi h4_T h4_T_cross h5a_min h5b_max h5c_max lambda_threshold growth_ratio rotation_tol return_tol ladder_tol drift_bound min_set K_c`.

Readings the tests pin where the spec leaves a choice (stage 03 follows them; none changes a frozen threshold):
- Periodic points are the panel controls of kind `fixed` or `newton` (the spec says "control kind periodic"); the return residual is the torus distance between the stored point at the recorded step q = `period` and the start. Kind `y5` is not periodic.
- `control_deviation` and `control_scope` are read per item (per method) from the orbit group of the `main`-arm items of the controls; the recomputed value is |lambda_m(N_1) - exponent|.
- Flags: "drift" and "excluded" count for an orbit if any of its items has them; "diagnostics exceeded" is read on the primary L5 item of arm `main`. The drift flag of an item must equal the recomputation from its stored points (some |p| above `drift_bound`; a wrapped item adds `2 pi` times the winding count, the rule of the runner).
- A start of the h4 arms is recognised by an id above the panel's largest id; its K comes from the orbit attribute `params`, matched to pi/2 and 2 pi within 1e-9.
- Item identity is by content: the attributes `method`, `precision`, `horizon`, `tau`, `perturbation` (L2), `direction`, `coordinate_mode` of the chunk file and the orbit id.
- The tables file is `<study>-tables.md` and the defect report `<study>-defect.md` (study slug from the preregistration header), with sections `## Header`, `## T1` to `## T7` and `## F1 data` to `## F7 data` (each T and F heading carries the word "confirmatory"; each F section has a line `digest: <64 hex>`), `## Readings` (one `- ` line per reading), `## Deviations`, `## Exploratory items` (the word "none"). No clock, no absolute path.
- H5 verdict: `:supported` if all three sub-rules are decidable and pass; `:not_supported` if any decidable one fails; otherwise `:not_decidable`.
- Excluded orbits stay in the cell counts and T6 (I7) but are left out of every fraction.

## Scenarios (the closed forms of the fixture)

The synthetic study (`an_helpers.jl`) has 20 panel orbits and 2 h4 starts. Series: `lambda_m(n) = lref (1 + c_m / n^q)` for the chaotic classes and controls, the integrable closed form `ln sqrt((cos a + n sin a)^2 + sin^2 a) / n` for class T and the elliptic control. Horizons are scaled (N_1 = 1e4, N_long = 1e5, N_B = 1e3, n_H2 = 100); thresholds are the frozen values, except `lambda_threshold = 1e-2` (the scaled regular floor is 9e-4, above the frozen 1e-4) and `min_set = 3`. Items go through the Experiments store writer, one chunk per method and parameter set (55 chunks), with the arms of the real plan (identical items merged).

| Scenario | Construction | Expected |
|----------|--------------|----------|
| default | c(L9) = 0.02, c(L2) = 0.001, others 0 | E_ch = {1 to 4, 11 to 14, 17, 18, 20}; W = {1 to 8}; S = {11 to 14}; H1 supported (11 of 11); R(100) = 2e-4 on W: H2 not supported; n_conv = 1 on W and S: ratio 1, H3 not supported; H4 supported; H5 supported (n = 4, 20, 5); labels positive and irregular (chaotic orbits and Y5), none detected and quasiperiodic (T), positive and periodic (17), none detected and periodic (19), positive and irregular for the info control 18 (return residual 2e-6) |
| weak | c(L9) = 3e5 / n^3 on orbits 1 to 8 | R(100) = 0.3 on W: H2 supported; n_conv(W) = 200, n_conv(S) = 1: ratio 200, H3 supported; H1 still supported (Delta at N_1 about 1e-7) |
| force | orbit 2 c(L9) = 100 (fails H1b); orbit 9 BigFloat +0.05; orbit 5 BigFloat +0.01; orbit 7 L5 direction "0.3" constant 1e-5; orbit 13 flag "diagnostics exceeded" on L5 | forced: 2 estimator, 9 and 5 ladder, 7 direction, 13 diagnostics; W = {1, 3, 4, 6, 8}; H1 falsified (b: 10 of 11); H5a = 3/4 (raw labels), H5 not supported; H5b = 0.25, H5c = 0 on 4 |
| flags | orbit 3 flag "drift" (lifted items, |p| = 2e6), orbit 4 "excluded", orbit 5 "diagnostics exceeded" | dropped {3 drift, 4 excluded}; E_ch loses 3 and 4; W = {1, 2, 6, 7, 8}; H1 n = 9; H5b n = 19; cells still sum to 20 |
| defect | control 17, method L9, +6e-8 at N_1 | status `:control_defect`, defect (17, L9), only the defect report is written |

## Criterion to test

| AC | Spec | Test file, testset | Notes |
|----|------|--------------------|-------|
| 1 | E17, E0 | `an_pipeline.jl` "AC1 refusals" (17 malformed inputs: status, chunk status, preregistration text, manifest hash, panel, missing chunk, pilot, incomplete, chunk schema, manifest schema, chunk hash, extra file, NaN, unflagged drift, wrong control deviation, missing and unknown parameter key; each refused with the named cause and no output directory); `an_params.jl` "AC1 parameter table" (37 missing keys, unknown key, wrong type, missing file) | the `stale/` area is ignored (same testset) |
| 2 | E1, E2, E17 | `an_pipeline.jl` "AC2 reader" (all 55 chunks, every orbit: t, exponent, point, flags, arms, attributes, diagnostics names, control attributes; winding of a wrapped item; item found by content) | n differs from k throughout |
| 3 | Oracle O1 | not in the RED suite: real store, stage 04 (verify script) | needs the 475 MB store |
| 4 | E4 to E6 | `an_h1.jl` "AC4 H1 metrics" | bounds at the value and one float below |
| 5 | E8 | `an_h1.jl` "AC5 H1 decision" | 98 / 98 / 90 of 100, 19 against 20 |
| 6 | E9 | `an_h1.jl` "AC6 exact-exponent controls"; `an_pipeline.jl` "AC6 control defect stops the analysis" | L9 own bound, info, exponent 0, stored deviation |
| 7 | E7, E10 | `an_sets.jl` "AC7 orbit sets", "AC7 finite-time spread", "AC7 H2 decision"; `an_pipeline.jl` "AC7 and AC12 flags" | K_c inclusive, S and X need K > 0, drift and excluded listed |
| 8 | E11 | `an_sets.jl` "AC8 convergence time" | analytic first step 501 for c = 20, late excursion, undefined, lam_ref <= 0, maximum over V |
| 9 | E12, O2 | `an_bootstrap.jl` (five testsets: stream, determinism, exact enumeration of 6^6 resamples against the module at 4 standard errors, coverage 900 to 1000 of 1000, level and ratio) | independent SHA-256 of the stream label |
| 10 | E11, E12 | `an_sets.jl` "AC10 H3 ratio" | ratio inclusive, interval end strict, independence of the W and S resamples |
| 11 | E13 | `an_h4.jl` (decision, comparators); `an_pipeline.jl` default scenario (H4, cross-check rows) | fit branches by K range, "no fit" between 1 and 1.1 |
| 12 | E14 | `an_labels.jl` (growth ratio, instability axis, structure axis, disagreement tests); `an_pipeline.jl` "AC12 and AC13 disagreement tests end to end" | about 40 label assertions, each test alone and together |
| 13 | E15 | `an_h5.jl` (H5a, H5b and H5c, verdict); `an_pipeline.jl` force scenario | H5a from raw labels, not 1.0 by construction |
| 14 | I1, I5 | `an_pipeline.jl` "AC14 a parameter changes only what depends on it" (six perturbations), "AC14 weak scenario", default scenario key order | decisions in the order H1 to H5 |
| 15 | E19, I7, I8 | `an_pipeline.jl` default scenario (sections, confirmatory marks, digests, readings at least 8, no exploratory item, cells sum to 20, labels once, header hashes, no clock or path) | |
| 16 | I3 | `an_pipeline.jl` "AC16 determinism" (threads 1 and 4, another output place; bytes equal, digests equal) | the run at 23 threads is stage 04 |
| 17 | I2 | `an_pipeline.jl` "AC17 read-only" | after a run and after a refusal |
| 18 | E19, Compute | `an_pipeline.jl` "AC18 figures" (status, digests, no PDF without the extension); `an_cli.jl` "AC20" (extension wiring in `Project.toml`, file `ext/ChaoslabCairoMakieExt.jl`) | the seven PDF files and their signature need CairoMakie: stage 04 |
| 19 | O3 | not in the RED suite: real store, independent recomputation, stage 04 | |
| 20 | Layering | `an_cli.jl` "AC20" (imports of `src/Analysis`, no plotting package, layer list, weak dependency and compat entry, verify environment) | `test/layering.jl` gained the layer in this stage (D60) |
| 21 | Compute | not in the RED suite: wall time and memory of the real store, stage 04 | |
| 22 | Command | `an_cli.jl` "AC22 command" (exit 2 with figures skipped, help, 5 usage refusals, refusal with cause, entry point routing) | exit code 0 needs CairoMakie: stage 04 |

Coverage of the spec: every acceptance criterion maps to at least one test, except 3, 19 and 21 (real store, stage 04 by the spec's own wording) and the figure files and exit code 0 of 18 and 22 (they need the plotting package, which the test environment lacks; stage 04 runs them under `verify/`).

## RED run

Command: `julia --project -e 'using Test; @testset "analysis" verbose = true begin include("test/analysis/runtests.jl") end'` (2026-10-20, Julia 1.13.2).

| Result | Count | Cause |
|--------|-------|-------|
| Pass | 20 | 16 fixture self-test checks (`an_fixture.jl`: the fixture is valid without the module); 1 own-enumeration check (6^6 resamples counted); 1 case count (17 refusal cases); 2 preconditions already true (the layer list of `test/layering.jl` has the layer, edited in this stage; `verify/Project.toml` has CairoMakie) |
| Error | 80 | all 80: `UndefVarError: Analysis not defined in Chaoslab` |
| Fail | 17 | checks on artifacts that do not exist: `src/Analysis`, `ext/ChaoslabCairoMakieExt.jl`, the `CairoMakie` weak dependency and compat entry, the `analyze` route in `bin/chaoslab.jl`, and the results of `cli` calls that did not run |
| Total | 117 | |

No error comes from a helper, a parse error or a broken fixture (checked: the only distinct exception text among the 80 errors is the `UndefVarError`). Test time before the module exists: 20 s (the fixture builds 5 synthetic studies of 55 chunks, 1.5 s each, plus 8.7 s for the first one with compilation).

The 17 fails by testset: AC22 command 9 (cli call results, entry point), AC20 layering 8 (`src/Analysis`, files, imports, plotting, weak dependency, compat, extension file).

## Corrections found in stage 03

Three assertions of this stage were wrong; the rest of the suite is unchanged. Each was found when the code first ran against the tests and is a defect of the test, not of the spec.

| File, testset | Defect | Correction |
|---------------|--------|------------|
| `an_sets.jl` "AC8 convergence time" | the case `lam = lref (1 + 2 / n)` with band 0.05 puts the boundary `2 / 40 = 0.05` on a grid step; in binary64 `abs(lam - lref) <= 0.05 lref` is false there (0.015000000000000013 against 0.015), so the expected step 40 depends on rounding | constant 2.2: inside the band from k = 44, expected grid step 50 |
| `an_sets.jl` "AC8 convergence time" | the same boundary with band 0.2 at k = 10 | constant 2.2: inside from k = 11, expected grid step 13 |
| `an_bootstrap.jl` "AC9 determinism of the interval" | two assertions compared the end points of one interval for two labels and for two seeds; the end points of a bootstrap median of 40 values are discrete and coincided for the fixed data | five labels and five seeds are compared and must not all give the same end points; the first resampled medians of the streams were checked to differ |

The pinned reading on drift (above) was widened from "lifted items" to all items; the test of the flags scenario is unchanged.

## Tests added in stage 04

Two files were added after the approved suite and the approved assertions were not changed. They are regressions of defects found on the real store and by the independent review.

| File, testset | Covers |
|---------------|--------|
| `an_breakdown.jl` | an item of the tau sweep with the flag `breakdown_tau` and blank rows (t = 0, NaN) does not refuse the run; F7 counts it and the tables file lists its step and kind; a NaN in an unflagged item is still refused |
| `an_review.jl` "failed chunk" | a manifest chunk with status failed is refused, naming chunk and status |
| `an_review.jl` "ladder half point" | a spike of x at a step that only the BigFloat item records does not change the ladder labels (one common half point) |
| `an_review.jl` "control scope" | a control with a positive exponent and no `control_scope` is refused |
| `an_review.jl` "parameter ranges" | level outside (0, 1), `N_1`, `min_set` below 1 and a NaN value are refused naming the key; a manifest without a chunk list is refused |
| `an_review.jl` "header" | module version and commit in the header; pooled series over W in F2 |

The fixture helper changed with them: the Y5 control has no `K` and no `K_text` (as in the real panel), a breakdown item is written with t = 0 and NaN rows, and a hook `state3(id, st, steps)` can change the stored points by step. Suite: 434 tests in `test/analysis/`.

