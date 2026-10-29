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


## Amendment D89 (stage 02)

Spec: `workspaces/module-build/stages/01-spec/output/analysis-spec.md`, section "Amendment D89 (2026-10-29)", approved with P9 to P34 as proposed (decision D89). Date of this stage: 2026-10-29, Julia 1.13.2. The module `Chaoslab.Analysis` has none of the D89 functionality yet: the names below do not exist, the P1 path is unchanged. Suite files (all registered LAST in `test/analysis/runtests.jl`, comment "D89 stage 02"): `an_p2_ref.jl`, `an_p2_helpers.jl`, `an_p2_study.jl`, `an_p2_build.jl`, `an_p2_refstore.jl` (definitions only), `an_p2_fixture.jl` (self-test), and the theme files `an_p2_dispatch.jl`, `an_p2_params.jl`, `an_p2_verdicts.jl`, `an_p2_detection.jl`, `an_p2_agreement.jl`, `an_p2_sets.jl`, `an_p2_stats.jl`, `an_p2_bootstrap.jl`, `an_p2_refusals.jl`, `an_p2_pipeline.jl`, `an_p2_variant.jl`, `an_p2_reported.jl`, `an_p2_output.jl`. No existing test file, no `src/`, no `data/`, no frozen document was touched.

Seeds. Every random draw uses a StableRNGs seed `p2_seed(letter, sub) = 20261012 + 31000 + 50 n + sub` (n = position of the criterion letter, sub < 50): `p2_seed('i', 1)` (random verdict tables), `('i', 2)` (Spearman), `('j', 1)` (AUC), `('m', 1)` (statistics), `('o', 2)` (coverage), `('o', 3)` (bootstrap against the reference), `('o', 4)` (exact enumeration). The bootstrap streams are keyed by `seed = 20261012` and the label `p2|<item>|<set>|<set rule>|<verdict rule>|<variant>`; four stored constants pin the stream (platform independence): label `p2|H1|W2|E_all|R1|all` has stream seed 3444927369723232618 and first outputs `0x986079ec7c98ba24, 0x1086a453d2baa0ac, 0x2c374cbb5b651d9f, 0x8a250f9f138a33db`; label `p2|H3|S2|E_all|R2|minus-E1` has stream seed 17828611372690436285 and first outputs `0x0943a72e539ace21, 0xe51db9067f4b005b, 0x8ccebb39745cc993, 0x346ec6584d67d5f2` (computed with StableRNGs 1.x from the rule E12, SHA-256 of the label text; the independent SHA-256 of the suite, the platform tool `sha256sum` or `certutil`, gives the same seed).

No memory values. Every expected number is (a) an analytic or hand value with its arithmetic in the test or in the section "Hand-computed expected values" below, (b) a number of the preregistration or of the spec text (E26 budget steps, the pilot table of `p2-indicators-pilot-evidence.md`, section 3, 192 verdicts), (c) a value computed by the independent reference implementation of `an_p2_ref.jl` and `an_p2_refstore.jl` (written from the amendment text, never calling the module, validated against the hand values by the self-test), or (d) the golden value of the P1 path taken from the module before the amendment. No real `data/` store is read; the only file of `data/` that the suite touches is the listing of the tree (relative path, size, modification time) for the read-only check.

### API proposed by the tests

The spec has no Julia names (Pattern 10); the names below are chosen here, in the style of the existing `Chaoslab.Analysis` (`decide_h1`, `bootstrap_interval`, ...: lower case with underscores, NamedTuple results, `nothing` for an undefined quantity, `ArgumentError` for a refusal), prefixed `p2_` so that the P1 API stays untouched. All live in `Chaoslab.Analysis`; the tests call `AN.<name>`. `p` is the parameter table (58 keys, below). Indicator symbols: `:LE :SALI :FLI :MEGNO :E1 :E2 :WBA :RLI` (`:STR` only for scores); verdicts `:chaotic :regular :undecided` (`:na` for an orbit without the item); rules `:R1 :R2`; hypotheses `:H1 .. :H6`.

| Name | Call | Returns |
|------|------|---------|
| `P2_PARAM_KEYS` | constant | tuple of the 58 keys of E22 |
| `p2_read_params` | `p2_read_params(path; prereg, p1_params = nothing, defaults = nothing)` | NamedTuple of the 58 keys (integer keys `Int`, reals `Float64`, `budgets`, `windows` `Vector{Int}`, `directions`, `gate_classes`, `grid_k_values` `Vector{String}`); `ArgumentError` naming the key for every defect of E22 (missing, unknown, repeated, mistyped, out of range, machine key different from the run blocks, quote or section not found, `text` missing, not a text or not in the quote (A2), `budgets`/`windows`/`N_c`/`b_primary` inconsistent, `bootstrap_seed` different from the `seed` of `defaults` (default: the repository `_config/numerics-defaults.toml`), a label key different from the key of the P1 table `p1_params`); the path for a missing file |
| `p2_budget_step` | `p2_budget_step(t, b)` | the first recorded step of `t` not before `b`; `ArgumentError` naming "budget" if there is none |
| `p2_nearest_step` | `p2_nearest_step(t, x)` | the recorded step nearest to `x`, a tie going to the larger step |
| `p2_window_row` | `p2_window_row(windows, b)` | position of `W = b / 2` in `windows`; `ArgumentError` naming "window" for an odd `b` or `W` not in `windows` |
| `p2_growth`, `p2_ratio` | `p2_growth(X, n, value)`, `p2_ratio(g, gp)` | `g_LE = n value`, `g_SALI = -value` (column `SALI.L`), `g_FLI = value`, `g_MEGNO = value`; `g / gp`, or `nothing` unless both are finite and `gp > 0` |
| `p2_verdict` | `p2_verdict(X, rule, v, p)` | the verdict by table R1 or R2. `v` per indicator: LE `(n, lambda, r)`; SALI `(r, sali)`; FLI `(n, r, fli)`; MEGNO `(r, ybar)`; RLI `(rli,)`; E1, E2 `(D, flags)` with the stored flag word; WBA `(dig::Vector, flags::Vector)` of the four rows of the window; `r` may be `nothing`. The bit of the j-th flag name has the value 2^j, so degenerate, undefined window and zero are the bit of value 2. `ArgumentError` for `:STR` and for an unknown rule |
| `p2_score` | `p2_score(X, v)` | the chaos score of E31 (`nothing` if unavailable; `-Inf` kept for RLI = 0); `X = :STR` takes `(D,)` |
| `p2_detection` | `p2_detection(V, rho, p)` | `(t_det, outcome, censored)` from the ten verdicts and `rho` in `:C :R`; `outcome` in `:correct :wrong :undetected`; `t_det` `nothing` if undefined |
| `p2_median` | `p2_median(v)` | middle value, mean of the two middle values for an even count |
| `p2_contingency` | `p2_contingency(vx, vy)` | `(cc, cr, rc, rr, n, A, kappa)`; `A` and `kappa` `nothing` if `n = 0`; kappa 0 if `p_e = 1` |
| `p2_unanimity` | `p2_unanimity(V::Dict{Symbol, Vector}, I)` | `(U, used, dropped, A_min, A::Dict{(X, Y)}, kappa::Dict, kappa_mean)`; orbits decided for every indicator of `I` enter `U`; pairs in the order of `I` |
| `p2_spearman` | `p2_spearman(x, y)` | Pearson correlation of the midranks, `nothing` for fewer than 3 values or a constant rank vector |
| `p2_auc` | `p2_auc(c, r; min_set = 0)` | the rank-sum AUC (ties one half), `nothing` if either class has fewer than `min_set` scores |
| `p2_regular_range`, `p2_leaving_budget`, `p2_share_below` | `p2_regular_range(scores, p)`, `p2_leaving_budget(scores, upper, p)`, `p2_share_below(scores, m)` | `(lo, hi, median, n)` from the `range_lower_percentile`-th and `range_upper_percentile`-th percentile (type 7), `nothing` below `min_set`; the budget `b_k` of the smallest k from which `scores[j] > upper[j]` for all j up to 9 (scores `Vector{Union{Nothing, Float64}}`), `nothing` if none; the share of scores below `m` |
| `p2_sets` | `p2_sets(orbits, p)`, orbit `(id, source, kindex, K, class, control, flags, rho, lambda_ref, lambda_ns, in_SD, in_SL)`, `source` `:panel` or `:grid` | `(E_all, E8, W2, S2, T_s, T_s_design, regular, gate, S_D, S_L, overlap, dropped)`; sorted id vectors; `dropped` is a vector of `(id, reason)` with reasons `"control"`, `"drift"`, `"excluded"`, `"unresolved"` (precedence in that order); `S_D`, `S_L` are S_D* and S_L*; `ArgumentError` naming "gate" if the gate set size differs from `gate_size` |
| `p2_decidable` | `p2_decidable(n, p)` | `n >= p.min_set` |
| `p2_rho`, `p2_grid_label`, `p2_rotation_numbers`, `p2_lift` | `p2_rho(instability, structure)`; `p2_grid_label((lam5, lam5_tenth, lam1, rho1, rho2), p)`; `p2_rotation_numbers(x0, xh, xN, h, N)`; `p2_lift(x, winding)` | `:C`, `:R`, `:U`; `(instability, structure, rho)` with the label texts of P1; `(rho1, rho2)`; the lifted coordinate `x + 2 pi winding` |
| `p2_stat` | `p2_stat(h, sets::Tuple, p; I = I_7, kp)`, a record is `(id, V::Dict{Symbol, 10 verdicts}, T::Dict{Symbol, censored time})` | H1 `(n, used, dropped, U, A_min, A, kappa_mean)`; H2 `(n, M, Q)`; H3 (sets `(W, S)`) `(n_w, n_s, M_w, M_s, R)`; H4 `(n, phi, decided, dphi)` (LE left out); H5 `(n, F, decided, F_med)`; `nothing` for an undefined quantity |
| `p2_h6_stats` | `p2_h6_stats(dir, prec, p)`; `dir` rows `(id, v::Dict{X => 5 verdicts at k = 9})`, `prec` rows `(id, f64::Dict, d64::Dict)` | `(dir = (n, delta, dropped), prec = (n, delta, dropped))`, dictionaries by indicator; undecided orbits dropped and counted |
| `p2_flips` | `p2_flips(V)` | the flip count Phi_X |
| `p2_decide` | `p2_decide(h, s, p)`, `s` the statistics and the lower end of the interval (`U_lo`, `Q_lo`, `R_lo`, `dphi_lo`, `F_med_lo`; H6 `(dir, prec)`) | NamedTuple with `verdict` in `:supported :not_supported :not_decidable` and the inputs echoed; collections by indicator may be vectors or dictionaries |
| `p2_robust` | `p2_robust(primary, extension)` | `:ROBUST`, `:NOT_ROBUST`, `:NOT_COMPARABLE` |
| `p2_label`, `p2_resample`, `p2_bootstrap` | `p2_label(item, set, set_rule, rule, variant)`; `p2_resample(label, m, B; seed)`; `p2_bootstrap(h, sets, labels, p; I, seed = p.bootstrap_seed)` | the label text; the `m x B` matrix of indices (column r = replicate r, drawn with `rand(rng, 1:m)` in order from the stream of the label); `(interval, dropped, defined)`, `interval` a dictionary by scalar (`:U`, `:Q`, `:R_<X>`, `:dphi`, `:F_med`) of `(lo, hi)` or `nothing` when more than 5 percent of the replicates are undefined. `labels` is a string, two strings (W then S) for H3 |
| `p2_read_store` | `p2_read_store(manifest, store)` | one NamedTuple per orbit record of every chunk: `orbit, id, kind, precision, horizon, tau, direction (nothing if none), mode, arms, chunk, t, state, winding, flags, attrs, lam, series, columns, final, rows, fcolumns, flagwords, bits, scalars, scolumns, range` (Julia shapes; `lam` the first exponent or `nothing`) |
| `p2_find_item`, `p2_flag`, `p2_grid_starts` | `p2_find_item(items; orbit, id, precision, horizon, tau, direction)`; `p2_flag(item, row, name)`; `p2_grid_starts(prereg)` | the item found by content or `nothing`; whether bit j of the word of that row is set; the starts of the `[[run.h4]]` tables `(id, label, kindex, i, j, K, x0, p0)` in table order, ids from 21 in the fixture (884 in the real study: the first id after the panel) |
| `analyze` | `analyze(; prereg, manifest, panel, params, store, out, figures = false, threads, p1_prereg, p1_manifest, p1_params, p1_store)` | for a P2 preregistration the result below; the four `p1_*` keywords are an `ArgumentError` for a P1 preregistration |
| `cli` | `cli(args; io)` | exit code 0, 1 (refusal, usage), 2 (tables written, figures skipped); the options `--p1-store`, `--p1-preregistration`, `--p1-manifest`, `--p1-parameters` are a usage error (exit 1, message naming `p2-indicators`) for a P1 preregistration |

Unchanged and used by the tests: `read_params` (37 keys), `analyze`/`cli` for `p1-le-benchmark`, `stream_seed`, `bootstrap_rng`.

Result of `analyze` for study `p2-indicators` (the fields the tests read; `H` below is the NamedTuple of the six decisions in the order H1 to H6):

| Field | Content |
|-------|---------|
| `status`, `study`, `tables` | `:ok` or `:figures_skipped` (with `figures = true`, no drawer); `"p2-indicators"`; path `<out>/p2-indicators-tables.md` (the only file written) |
| `header` | `prereg_header_sha256, manifest_sha256, panel_sha256, params_sha256, p1_prereg_header_sha256, p1_manifest_sha256, p1_params_sha256, plan_sha256, chunks, items, seed` |
| `rho`, `lambda_ref`, `lambda_ns`, `p1_labels` | dictionaries by orbit id (P and G; the labels of the panel orbits are `(instability, structure)` of the P1 pipeline) |
| `grid` | rows `(id, label, kindex, i, j, K, x0, p0, instability, structure, rho, ref_weak)` |
| `sets` | as `p2_sets` |
| `verdicts`, `scores`, `steps`, `detection` | dictionaries keyed `(id, X, rule)`, `(id, X)` (also `:STR`), `(id, X)` (the step n(b_k), or the window W, per budget index), `(id, X)` for the orbits of E_all |
| `decisions`, `extension`, `r2`, `r2_extension`, `robust` | `H` on E_all with I_7; `H` on E8 with I_8 (sets W2_8, ...); `H` with verdicts of R2 (E33: detection times, M_X, Q, R_X and the H6 verdicts of the direction and Double64 items all with R2); the same for the extension; `(H1 ... H6)` of ROBUST flags. `H.H1 = (verdict, n, used, dropped, U, U_lo, U_hi, A_min, A, kappa_mean, U_S2)`, `H.H2 = (verdict, n, M, Q, Q_lo, Q_hi)`, `H.H3 = (verdict, n_w, n_s, M_w, M_s, R, R_lo, passing)`, `H.H4 = (verdict, n, phi, decided, dphi, dphi_lo, design = (same fields for T_s_design))`, `H.H5 = (verdict, n, F, decided, F_med, F_med_lo)`, `H.H6 = (verdict, dir, prec)` |
| `gate`, `second` | rows `(X, decided, chaotic, undecided, na, share, violated)`; `Dict(X => (decisions, extension))` for each violating indicator, with X removed from the indicator set |
| `auc`, `range`, `leaving`, `below_median` | `(X, k, K index, source)` => `(n_C, n_R, auc)` (`source` `"pooled"`, `"panel"`, `"grid"`; `auc` `nothing` below `min_set`); `(X, K index, k)` => regular range or `nothing`; `(orbit, X)` => leaving budget; `(X, K index)` => share |
| `sensitivity`, `packages` | `(pair, k)` => `(n, A, kappa, median_diff, max_diff)` for the pairs `:raw_vs_fq1 :qp_vs_fq1 :planar_vs_pa1 :cos_vs_wba`; `(set, (X, Y), k)` => dictionary of joint verdicts to counts |
| `flag_counts`, `omissions`, `controls_table`, `k0line`, `cost` | spec id => flag => count; spec id => omitted orbit ids; rows `(id, scope, exponent, lambda, deviation, stored_deviation, verdicts9)`; window => `(fq1, pa1)` largest circle distance; spec id => `(measured_ns, budget_line_ns, bytes_per_orbit, transform_seconds)` |
| `items`, `digests`, `readings`, `deviations` | `T1 .. T17`, `F1 .. F6` each `(columns, rows)` (columns pinned where the tests read them: T1 `source, K index, orbits, C, R, U, flagged`; T2 `indicator, k, budget, set, source, K index, chaotic, regular, undecided, not applicable, step`; F1 `set, budget, X, Y, cc, cr, rc, rr, N, A, kappa`; F2 `set, indicator, kind, value, count`; F3 `set, indicator, k, decided, fraction`; F4 `indicator, k, K index, source, n_C, n_R, AUC`; F5 `id, label, K index, i, j, rho` and `<X> b_P`, `<X> N_c`; F6 `class, step, count, p5, p25, p50, p75, p95`); `Dict("F1" .. "F6" => hex)`; vectors of text |

The tables file `<out>/p2-indicators-tables.md` has the sections `## Header` (the hashes of the preregistration header, manifest, panel, parameter table and of the P1 preregistration, manifest and parameter table, the plan hash, the seed, the module version, the git commit, the 58 parameters), `## Decisions`, `## T1 ...` to `## T17 ...` and `## F1 data ...` to `## F6 data ...` (each heading carries the word "confirmatory"; each data section a line `digest: <64 hex>`), `## Readings`, `## Deviations`, `## Exploratory items` (the word "none"). The digest is the SHA-256 of the data table written as markdown (`| a | b |` lines, header line and separator line first, numbers as `round(x; sigdigits = 12)`, lines joined by LF, no final LF), the P1 rule. The second versions of a violating gate are the sections `## T3 without E1 ...` to `## T8 without E1 ...` and `## T13 without E1 ...`.

### Readings pinned by the tests (where the spec leaves a choice; none changes a threshold)

1. Stream labels for the extension: the resampled set is named with the suffix 8 (`W2_8`) and the set rule is `E8`: `p2|H1|W2_8|E8|R1|all`; H4 on the design ensemble uses the item `H4design` and the set `T_s_design`; a second version uses the variant `minus-<X>` (P20, P21).
2. `dropped` of the sets: an orbit is listed once; precedence control, drift or excluded, unresolved (I16, "every orbit removed from a set appears in a drop list with its reason").
3. A statistic that is undefined (empty denominator, no available score, interval undefined) makes its decision `:not_decidable` (P27, P33); the statistics are still printed.
4. A2 (owner decision of the D89 stage 03 checkpoint, replaces the string(value) reading of the first version of this plan): every key other than the 11 machine keys carries `section`, `quote` and `text`; `text` is the value as the section writes it (`1e-4`, `0.20`, `-4`, `10^3`, `95 percent`, `5th`), must occur in `quote`, and `quote` in the section; the string form of `value` is not required to occur anywhere and the value is not matched to its text. The fixture writes the section style (47 keys, many different from Julia's `string(value)`); the style `:julia` (text = string(value)) is tested as an accepted variant.
5. T_s_design, S_D*, S_L*, E8 are intersected with E_all; the gate set is not (P26). The extension uses the sets of E_all restricted to E8 with the suffix 8; H6 (a) on the extension is the primary H6 (a) restricted to E8, H6 (b) uses the seven indicators with a Double64 item (P23).
6. The gate share has the decided orbits at k = 9 in the denominator; an indicator with no decided orbit has no share (`nothing`); `min_set` is not applied to it (P26).
7. The digest rule ("as E19", E37): the base rule above.
8. Flag counts per spec are the sums over the chunks of the manifest's `indicator_flags`; the module refuses a chunk whose recomputed counts differ.
9. The percentile keys `range_lower_percentile` and `range_upper_percentile` are in percent (5 and 95).
10. `ref_weak` is true for every grid orbit; the Y5 control has K index 0 and no K (P28); the controls count in T1 by their own K index.

### Synthetic study "p2-indicators" (fixture, scaled)

Built with the real chunk writer (`Experiments.write_chunk`, schema version 2 with the D86 extension, 33 chunk files, 466 orbit records) in `an_p2_study.jl` and `an_p2_build.jl`; the fixture preregistration, the parameter table and the run blocks are checked by the self-test with the real run-block reader (`Experiments.read_run_block`).

- Scale. Budgets b_k = 40 x 2^k (k = 0 to 9), windows 20 x 2^k, `N_c` = 20480 = 2 x 10240, `b_primary` = 2560 (k = 6), `persist_divisor` 4 (so b_7 = 5120 = N_c / 4 is the last stable start), `censor_factor` 2 (t~ in {40, ..., 5120, 40960}), `N_1` = 80000, `N_dir` = 40000, `n_s` = 800, `min_set` = 5, `gate_size` 2 (class T at K index 1, K = 0.5 in this panel), `bootstrap_resamples` 1000. Every threshold of R1, R2, the hypotheses, the sets and the labels is the frozen study value (the label instability floor is the P1 fixture value 0.01).
- Recorded grids (E26). GA (L1 and tangent items, horizon N_1): 1 to 10, then n_k = 44 x 2^k = 44, 88, 176, 352, 704, 1408, 2816, 5632, 11264, 22528 (the first step not before b_k) and 4, 9, 18, 35, 70, 141, 282, 563, 1126, 2253 (the step nearest to n_k / 10), plus 800, 8000, 40000, 80000. GAD = GA up to 40000 (direction and Double64 items). GB (pair, stretching, E1 and WBA items, horizon N_c): GA up to 5632, then 2048, 10240, 20480; n(b_k) = 44, ..., 5632, 10240, 20480. GB1 = GB and 20481 (E2).
- Panel and grid. The 20 orbits of the P1 fixture panel (16 orbits, 4 controls, Y5 = 20; S_D = {1, 2, 5, 7, 9, 11, 12, 15} and the controls, S_L = {1, 5, 9, 11}) and a grid of 24 starts (ids 21 to 44, K texts `0.5` and `6.283185307179586`, a 3 x 4 lattice, labels `grid-k<K>-i<ii>-j<jj>`). The P1 labels come from a P1 reference study (slug `p1-le-benchmark`, 55 chunks, the P1 fixture with orbit 2 failing the estimator test and orbit 13 the diagnostics test): rho C for 1, 3, 4, 5, 6, 7, 8, 11, 12, 14 and the control 18 and Y5; R for 9, 10, 15, 16; U for 2, 13 and the periodic controls 17, 19. Grid classes by the reduced rule: C 21 to 26 and 33 to 40, R 27 to 31 and 41 to 43, U 32 (the second half rotation differs by 8e-7 > 1e-7) and 44 (|lambda_L1 - lambda_L5| = 1e-6 > 1e-8).
- Designed verdict sequences. Every orbit and indicator has a designed sequence of ten verdicts (strings of `c r u`, built from an onset index per profile: S (strong), S2, W (weak), W2, T (sticky), N (neither), R (regular), plus overrides for the regular orbits), and the stored numbers realise it by the rules of E26 and E27 (growth ratios r = n g(n') / n' ..., RLI 1e-3 or 1e-12, D = -1 or -9, dig minimum 3 to 5 or 11 and more, flag words with every kind of bit, `degenerate` / `undefined window` / `zero` for the undecided rows). The self-test reads the chunk files back with the real reader and applies the independent rule R1: 334 of 334 sequences agree.
- Sets (hand): lambda_ref and lambda_ns are placed at the boundaries (S2: 1.0 and prevfloat(1.0); W2: 0.1 and nextfloat(0.1); T_s: lambda_ns = 0.1 against 0.2 x 0.5 and nextfloat(0.1)): E_all = 36 orbits, W2 = {1, 3, 5, 7, 8, 21, 22, 23, 24, 26, 39} (11), S2 = {11, 12, 33 to 38} (8), T_s = {5, 23, 24, 39} (4, below `min_set`: H4 not decidable), T_s_design = {5, 6}, regular = 12, E8 = 26, gate = {9, 10}, S_D* = {1, 5, 7, 9, 11, 12, 15} (7), S_L* = {1, 5, 9, 11} (4: H6 (b) not decidable).
- Scenarios. `default`: no flagged orbit, gate not violated, H1 SUPPORTED, H2 SUPPORTED, H3 SUPPORTED (6 of 7), H4 NOT DECIDABLE, H5 NOT SUPPORTED, H6 SUPPORTED (part a); ROBUST for H1, H2, H3, H5, NOT COMPARABLE for H4 and H6. `var`: orbit 8 flagged drift (|p| = 2e6 on every item), orbit 25 sticky (T_s = 5, H4 SUPPORTED), gate orbit 10 chaotic for E1 at k = 9 (share 0.5 > 0.10: E1 violates, second versions without E1), RLI of the orbits 1, 21, 22, 26 chaotic only from k = 8 (the extension H5 is SUPPORTED, the primary NOT SUPPORTED: NOT ROBUST).
- Double64 and direction items differ from the primary at k = 9 on purpose: Double64 orbit 5 E2 and orbit 11 LE (delta_prec = 0.25); direction: orbit 7 SALI in direction 3 and orbit 12 LE in direction 5 (delta_dir = 1/7), orbit 15 FLI undecided in direction 2 (dropped and counted).

### Hand-computed expected values

The arithmetic is in the comments of the tests; the self-test (`an_p2_fixture.jl`) checks the reference implementation against these values before the module is involved.

| Quantity | Value and arithmetic | Test |
|----------|----------------------|------|
| n(b_k) on GA, GB | b_k = 40 x 2^k; the first grid step not before b_k: 44, 88, 176, 352, 704, 1408, 2816, 5632, 11264, 22528 (GA), the same up to 5632 then 10240, 20480 (GB); n' = nearest to n / 10 = 4.4, 8.8, 17.6, 35.2, 70.4, 140.8, 281.6, 563.2, 1126.4, 2252.8 gives 4, 9, 18, 35, 70, 141, 282, 563, 1126, 2253 | `an_p2_verdicts.jl` D89e, self-test |
| n(b), n' on the real grids | 200 x 2^k on the grids of `Experiments.record_grid`: 200, 447, 891, 1778, 3548, 7079, 14125, 28184, 56234, 112202 (horizons 1e7, 1e6) and ..., 51200, 102400 (102400); n' = 20, 45, ..., 2818, 5623, 11220 and 5012, 10240 (the numbers of E26 and P15) | D89e |
| kappa | [[20, 5], [10, 15]]: p_o = 35/50 = 0.7, p_e = (25 x 30 + 25 x 20) / 2500 = 0.5, kappa = 0.2 / 0.5 = 0.4; [[25, 25], [25, 25]]: 0; one verdict only (p_e = 1): 0; a = d = 0, b = c: -1 | `an_p2_agreement.jl` |
| U, A_XY, kappa of the 3-indicator table | LE c c r c c, SALI c c r u r, FLI c r r r r: used 4, split 2, U = 0.5; A = 3/4, 2/5, 3/4; kappa 0.5, 0.08/0.68, 0.5; mean 0.37254901960784315 | `an_p2_agreement.jl` |
| Spearman with ties | x = (1, 2, 2, 3), y = (1, 3, 2, 4): midranks (1, 2.5, 2.5, 4), (1, 3, 2, 4): Pearson = 4.5 / sqrt(4.5 x 5) = sqrt(0.9) | same |
| AUC | separated 1; identical score sets 0.5; C = {1, 2, 3} against R = {2}: (0 + 0.5 + 1) / 3 = 0.5; {-Inf, 1} against {-Inf}: (0.5 + 1) / 2 = 0.75 | same |
| Regular range | type 7 on 1..20: h = 19 x 0.05 + 1 = 1.95 and 19 x 0.95 + 1 = 19.05; median 10.5 | same |
| Detection (b_7 = 5120 = N_c / 4) | `rrrrrrrccc` stable from k = 7: t = 5120; `rrrrrrrrcc` undetected (censored 40960); `cccccccrcc` undetected; `uucccccccc` t = 160; `cccrrrrrrr` for rho R: 320; `rrrrrrrrrr` for rho C: wrong, 40960; the number of sequences per outcome of the 3^10 for rho C: correct 3^7, wrong 3^7, undetected 3^10 - 2 x 3^7 | `an_p2_detection.jl` |
| F_X, F_med | four orbits, (b_P, k = 9) pairs: F = 1/2, 0, 1/3, 0, 1, 0, 1/2, median (4th of 7) 1/3; with RLI 1/4: (1/4 + 1/3) / 2 | `an_p2_stats.jl` |
| phi_X, Delta_phi | five orbits: 3/4, 0, 1, 2/5, 2/5, 4/5; Delta_phi = 1 - 0 = 1 | same |
| M_X, Q, R_X | S: 80, 40, 160, 160, 2560, 160, 1280, Q = 2560 / 40 = 64; W (even count): 480, 80, 320, 240, 25600, 1280, 40960; R = 6, 2, 2, 1.5, 10, 8, 32 | same |
| delta_dir, delta_prec | four orbits: LE 1/4, SALI 0 (one orbit dropped), FLI 2/4, MEGNO 0; five orbits: LE 2/5, FLI 1/4 (one dropped) | same |
| Sets of the default study | E_all = 36, W2 = {1, 3, 5, 7, 8, 21, 22, 23, 24, 26, 39} (11), S2 = {11, 12, 33 to 38} (8), T_s = {5, 23, 24, 39}, T_s_design = {5, 6}, regular 12, E8 26, gate {9, 10}, S_D* 7, S_L* 4; boundaries: S2 at 1.0 and prevfloat(1.0), W2 at 0.1 and nextfloat(0.1), T_s: 0.1 <= 0.2 x 0.5 (equal) and nextfloat(0.1) | `an_p2_fixture.jl`, `an_p2_sets.jl`, `an_p2_pipeline.jl` |
| T1 | rho counts (orbits, C, R, U): panel K index 1 (13, 8, 2, 3), K index 2 (6, 3, 2, 1), K index 0 (1, 1, 0, 0), grid K 1 (12, 6, 5, 1), K 2 (12, 8, 3, 1); sums 44, 26, 12, 6 | `an_p2_pipeline.jl` |
| H1 .. H6 of the default study | U = 1.0 (interval [1, 1]), A_min = 0 (FLI against E1); Q = 16, interval [16, 64]; R_X = 64, 16, 16, 32, 8, 16, 128, passing 6 (E1 8 < 10); H4 not decidable (4 < 5); F_med = 0; delta_dir(LE) = delta_dir(SALI) = 1/7, delta_prec(LE) = delta_prec(E2) = 1/4, part (b) not decidable (4 < 5) | `an_p2_pipeline.jl` |
| AUC counts | K index 1: 13 chaotic, 7 regular orbits of E_all (RLI 8 and 6); K index 2: 11 and 5 (RLI 9 and 3: not computed); by source only the grid at K 1 (6 and 5) | same |
| T2 rows | LE at b_P on W2: 7 chaotic, 4 regular, step 2816; RLI on W2: 8 chaotic, 3 not applicable; E1 on W2: 11 regular, window 1280 | `an_p2_output.jl` |
| F1 row | W2 at b_P, (LE, E1): cr = 7, rr = 4, A = 4/11, kappa = 0 | same |
| T11 | fq1_raw: D + 0.5 on every row, so |dD| = 0.5; fq1_qp: D - 0.25 on odd windows; pa_planar: D + 0.1 on odd windows; k = 3: four regular, A = 1, kappa 0; k = 5: kappa 1 | `an_p2_reported.jl` |
| T15 | W2 at b_P, (MEGNO, SALI): (c, c) 7, (r, r) 4 | same |
| T16 | K = 0 line: orbit 9 has nu1 = p0 / (2 pi) + k x 1e-10, orbit 10 twice that (pa1: 3e-10 per window index): largest per window 2 x 9 x 1e-10 = 1.8e-9 (fq1), 5.4e-9 (pa1) at the last window, 0 at the first | same |
| Pilot | 24 starts, 192 verdicts of the pilot evidence, section 3 | `an_p2_verdicts.jl` D89y |
| R2 version (E33) | LE chaotic iff n lambda >= 3 ln n; the stored n lambda at n_k is r x 0.2 x n'_k, so R2 needs r >= 15 ln(n_k) / n'_k = 14.19, 7.46, 4.31, 2.51, 1.41, 0.77 for k = 0 .. 5; a chaotic row (r >= 8.0) is chaotic from k = 1, a regular row (r 1.2 to 1.45) from k = 5 (k = 4 when id mod 5 = 4). S2: every LE onset is 1, M_LE(S2) = 80 (R1: 40); SALI: ids 12, 33, 36 (id mod 3 = 0) hold sali = 1e-9 on regular rows, chaotic under R2, giving t~ = 40 six times and 80 twice: median 40. W2 (11 orbits): LE onsets 5 (nine orbits) and 4 (orbits 24, 39): t~ = 1280 nine times and 640 twice, M_LE(W2) = 1280; R_LE = 1280 / 80 = 16 (R1: 64). H6: the regular LE row of orbit 12 in direction 5 is chaotic under R2 (r >= 1.2 > 0.067 at k = 9), so delta_dir(LE) = 0 (R1: 1/7); the SALI flip of orbit 7 stays (sali 0.1): 1/7; the undecided FLI row of orbit 15 holds fli = 3, regular under R2: nothing dropped; Double64: delta_prec(LE) = 0, delta_prec(E2) = 1/4 | `an_p2_pipeline.jl` |

### Criterion to test

Tags of the spec: (S) package suite, (V) verification environment with the real stores. Every (S) criterion has at least one test; the (V) criteria and the (V) halves are for the stage 04 scripts (`verify/`).

| Criterion | Test file, testset | Notes |
|-----------|--------------------|-------|
| D89a (S) | `an_p2_dispatch.jl` "D89a P1 path unchanged" (golden values taken from the module before the amendment: the SHA-256 of the tables file without the two lines "module version" and "git commit", the seven figure digests, exit codes 2, 0 and 1, two refusal messages, the 37 keys), "D89a dispatch on the study slug" (unknown slug, header slug against manifest `study` for both studies, the two parameter tables exchanged, the four `--p1-*` options and keywords with the P1 preregistration, `--out` inside `data/`) | the first testset is a guard: it passes before and after. (V): the real P1 inputs, stage 04 |
| D89b (S) | `an_p2_params.jl`: valid table; each of the 58 keys missing and mistyped (116 files); for each of the 47 non-machine keys `text` missing, of the wrong type, and not in the quote (141 files, A2); unknown and repeated key; 12 out-of-range values; empty `budgets`; four inconsistencies of `budgets`, `b_primary`; the 11 machine keys altered; a quote (with a text that occurs in it) not in the section, a section not found, a quote changed alone (the text is no longer in it); value formats: the section style (1e-4, 0.20, 10^3, 95 percent, 5th) and the `string(value)` style are both read, a value changed alone is read as given; the seed; the label keys (they carry quote and text) | the `text` rule replaces the value-text rule of the first version |
| D89c (S) | `an_p2_refusals.jl` "D89c more than 40 malformed inputs" (54 cases: slug absent, several, invalid; header hash row absent, repeated, wrong; each manifest hash; chunk hash attributes; schema, status, chunk status failed and pending, git, session, control stops, chunk and item counts; plan and panel hashes (manifest, chunk, file, build); missing, pilot, incomplete, wrong-version, wrong-hash and extra chunk; omission not listed, listed with an item, item missing; four identity mismatches; missing step and window row; non-finite without flag (frequency, series), `complete` 0; drift unflagged and spurious; missing table key; wrong `control_deviation`; wrong flag count; five P1 reference failures); "D89c things that are not refusals" (overrides listed, `stale/` ignored) | rule o (an input file changed while the run goes on) needs a seam in the module and is not tested; the case `N_c` not twice the largest window cannot be isolated (the run-block reader refuses first) and is covered by the machine key `N_c` altered |
| D89d (S) | `an_p2_refusals.jl` "D89d reader and identity" (466 items against the real reader, shapes, flag bit rule, grid starts) and the refusal cases of E24 and E25 | the real HDF5 extents are reversed in the file: not observable through the Julia API, covered by the shapes and by stage 04 O9 |
| D89e (S) | `an_p2_verdicts.jl` "D89e budget rule" | the two real grids built by `Experiments.record_grid` |
| D89f (S) | `an_p2_verdicts.jl` "D89f growth ratio and verdict R1", "D89 chaos scores"; `an_p2_pipeline.jl` (flags in the store do not change a verdict; RLI absent outside S_L and G) | |
| D89g (S) | `an_p2_verdicts.jl` "D89g verdict R2" | |
| D89h (S) | `an_p2_detection.jl` (118098 sequences against an independent characterisation; hand cases; medians) | the golden file of the spec is replaced by the independent suffix characterisation of the test (the Python implementation O6 is stage 04) |
| D89i (S) | `an_p2_agreement.jl` (kappa hand tables, unanimity hand table, 1000 random tables, Spearman with ties) | |
| D89j (S) | `an_p2_agreement.jl` (AUC hand, random, min_set 19 and 20, negation; regular range; leaving budget; share below), `an_p2_pipeline.jl` AUC, range, leaving budget end to end | |
| D89k (S) | `an_p2_sets.jl`, `an_p2_pipeline.jl` "the sets", `an_p2_variant.jl` (drift orbit) | |
| D89l (S) | `an_p2_sets.jl` "D89l", `an_p2_pipeline.jl` (panel rho against the base pipeline's labels; the grid orbits) | |
| D89m (S) | `an_p2_stats.jl` (hand values and 200 random record sets) | |
| D89n (S) | `an_p2_stats.jl` "D89n decisions", `an_p2_pipeline.jl` (order, extension, ROBUST, the R2 version and its extension end to end, gate), `an_p2_variant.jl` (NOT ROBUST, second versions) | "a later decision never changes an earlier output" is checked by the order of the fields and by the independence of the statistics, not by two runs. E33 R2 version (owner decision): the independent reference recomputes the detection times, M_X, Q, R_X and the H6 verdicts with R2 (`Det2`, `p2_ref_recs(rule)`, `p2_ref_h6(rule)`); hand checks in the section "Hand-computed expected values" |
| D89o (S) | `an_p2_bootstrap.jl` | |
| D89p (S) | `an_p2_output.jl` (T2, F1 to F6, header, sections, digests, readings, deviations), `an_p2_pipeline.jl` (T1) | |
| D89q (S) | `an_p2_output.jl` "D89r read-only ... D89q" (1, 4 and 3 threads, two folders, regeneration after deletion) | 23 threads: stage 04 |
| D89r (S) | `an_p2_output.jl` (inputs, P1 reference, store and the data tree listing unchanged after a run), `an_p2_refusals.jl` (every refused case leaves the bytes of all inputs unchanged), `an_p2_dispatch.jl` (`--out` in `data/`) | (V): the real stores, stage 04 |
| D89s (S) | `an_p2_output.jl` "D89s command" | |
| D89t (S) | `an_p2_dispatch.jl` "D89t layering and dependencies" | the layering test `test/layering.jl` already lists Analysis |
| D89u (S) | `an_p2_output.jl` "D89u" and the status in the determinism testset | the PDF files need the extension (P11): not in the suite |
| D89v (V) | stage 04 | prefix consistency, flag counts and E24 on the 695 real chunks, the K = 0 line |
| D89w (V) | stage 04 | independent recomputation (O6, O10) on the real stores |
| D89x (V) | stage 04 | the real-run report, wall time and memory |
| D89y (S, V) | `an_p2_verdicts.jl` "D89y" (192 verdicts from the 24 tabulated rows) | (V) half only if a pilot store is made available (P30) |

(V)-only: D89v, D89w, D89x, and the (V) parts of D89a, D89r, D89u, D89y.

### What the tests cannot cover

- Rule o of E36 (an input file changed while the run goes on): the run is a single call; a test would need a hook inside `analyze`. Left to stage 04 or to a seam the implementer may add (`analyze(...; hook)` is not part of the proposed API).
- 23 threads (D89q), peak memory, wall time, real stores, the figure files (P11), the independent Python layer: stage 04.
- The actual extents of the HDF5 datasets (the Julia API shows the Julia shapes): stage 04.
- The wording of the messages beyond the words the tests require (the cause, the key, the orbit, the file).
- That a statistic computed on 695 chunks equals one computed on 33: scale effects (the 8 GB and 15 minute limits).
- Every flag name in the store: the chunks carry below_floor, degenerate, tie1, mirror1, refined1 and refined2 of the frequency table, below_floor and undefined window of the phase advance, below_floor, zero and reldig_undefined of the Birkhoff table; the eleven bit values of the frequency table are exercised as words by the unit tests of `p2_verdict`.

### Spec ambiguities and decisions for the owner or stage 03 (none edits the spec)

- A1. The slug dispatch (E21) refuses any slug other than the two known ones, but the base fixture of the existing suite (`an_helpers.jl`, `an_study`) writes the slug `an-fixture`, the manifests with `study` `an-fixture` and expects the file `an-fixture-tables.md` (in `an_pipeline.jl`, `an_cli.jl`, `an_breakdown.jl`, `an_review.jl`). After the amendment the tests that build a study with `an_study` are refused unless the fixture slug becomes `p1-le-benchmark` (and the file names with it). D89a says the golden values are those of "the P1 fixture run of the base suite". This stage did not edit the existing files (instruction), so the golden guard uses its own P1 fixture with the slug `p1-le-benchmark` (`p2_p1_study`, a copy of `an_study` with a slug argument). Stage 03 has to change `an_helpers.jl` (one slug and the expected names) or the owner has to state another rule.
- A2 (resolved by the owner at the D89 stage 03 checkpoint): the value text is a required field `text` of each non-machine key, written as the section writes it and found in the quote; the tests follow this rule (reading 4).
- A3. The label grammar for the extension sets (E34): the set named with the suffix 8 together with the set rule E8 is a reading (above, item 1); P20 leaves the wording open.
- A4. The digest rule (E19, E37) is stated as "numbers written with 12 significant digits, LF line ends": the tests use the base rule (markdown table lines joined by LF, no final LF, `round(x; sigdigits = 12)`).
- A5. E22 gives `range_lower_percentile` and `range_upper_percentile` without a unit; the tests use percent (5 and 95, as "5th and 95th percentile").
- A6. E32 H1/H4/H5 statistics that are undefined (empty denominator) make the decision NOT DECIDABLE (P27, P33); E33's gate "decided" counts: the tests pin the share as chaotic over decided with undecided and not applicable reported separately.
- A7. The in-memory P1 labelling (E23 a) must give the same labels as the P1 `analyze` on the same inputs; the tests compare `p1_labels` with the labels of the P1 `analyze` result. A control defect in the P1 reference run is expected as a refusal naming "control" (the base `analyze` returns a status instead).
- A8. The pure functions take records and parameter tuples as documented above; if stage 03 prefers other shapes the tests change, the contract (what is computed) does not.

### RED log (2026-10-29)

Command, from the repository root: `julia --project -t 4 --check-bounds=yes -e 'using Test, Chaoslab; @testset "an" begin include("test/analysis/runtests.jl") end'` (Julia 1.13.2, 4 threads). Wall time 1 min 47 s (the 434 old tests take 68 s; the new files take about 40 s before the module exists, of which 25 s are the two synthetic studies and their P1 reference; the studies are cached for the whole run).

| Result | Count | Cause |
|--------|-------|-------|
| Pass, old tests | 434 | the P1 analysis suite, unchanged (it passed 434 of 434 before this stage, same count after) |
| Pass, new tests | 110 | 66 fixture self-test checks (`an_p2_fixture.jl`: the study, the preregistration, the run blocks, the reference implementation, none calls the module); 12 guards of the P1 path (golden values taken before the amendment); 4 layering checks and 5 dispatch checks that hold already; 5 of the stream test (existing `stream_seed` and `bootstrap_rng`, rule E12 reused) and 1 length check of the enumeration, 15 usage-error checks of `cli` (the P1 command already refuses unknown and missing options); 1 count of the refusal cases and 1 of the machine keys |
| Fail, new tests | 11 | dispatch 6: 4 `cli` message checks (the message of a `--p1-*` option does not name `p2-indicators`), 2 `an_rejects` checks (a manifest `study` that differs from the slug is not refused; a P2 table given with the P1 preregistration is refused without naming both keys); output 5: 2 exit codes (`cli` for a P2 preregistration returns 1, not 2, because the options are unknown), 1 help text without the `--p1-*` options, 2 reads of tables that were not written |
| Error, new tests | 206 | 183 `UndefVarError: p2_<name> not defined in Chaoslab.Analysis` and 23 `MethodError: no method matching analyze(; ..., p1_prereg::String, ...)` |
| Broken | 0 | |
| Total | 761 | |

First line of each distinct error class:

- `UndefVarError: \`p2_read_params\` not defined in \`Chaoslab.Analysis\`` (and the same for `p2_auc`, `p2_bootstrap`, `p2_budget_step`, `p2_contingency`, `p2_decide`, `p2_detection`, `p2_grid_label`, `p2_growth`, `p2_h6_stats`, `p2_label`, `p2_median`, `p2_nearest_step`, `p2_ratio`, `p2_read_store`, `p2_regular_range`, `p2_resample`, `p2_rho`, `p2_rotation_numbers`, `p2_score`, `p2_sets`, `p2_spearman`, `p2_stat`, `p2_unanimity`, `p2_verdict`, `p2_window_row`; and `P2_PARAM_KEYS` in `an_p2_dispatch.jl`): 183 errors.
- `MethodError: no method matching analyze(; prereg::String, manifest::String, panel::String, params::String, store::String, out::String, p1_prereg::String, ...)`: 23 errors (the keyword arguments of the P2 study).
- `Test Failed` (11): expressions `occursin("p2-indicators", msg) && ...`, `AN.cli(p2_cli_args(st; ...)) == 2`, `an_rejects(...)` and `all(o -> occursin(o, h), ("--p1-store", ...))`, all about behaviour that the amendment adds (the P1 command is unchanged).

Per file (pass / fail / error): `an_p2_fixture.jl` 66 / 0 / 0; `an_p2_dispatch.jl` 21 / 6 / 6; `an_p2_params.jl` 1 / 0 / 51; `an_p2_verdicts.jl` 0 / 0 / 42; `an_p2_detection.jl` 0 / 0 / 11; `an_p2_agreement.jl` 0 / 0 / 19; `an_p2_sets.jl` 0 / 0 / 27; `an_p2_stats.jl` 0 / 0 / 24; `an_p2_bootstrap.jl` 6 / 0 / 7; `an_p2_refusals.jl` 1 / 0 / 3; `an_p2_pipeline.jl` 0 / 0 / 6; `an_p2_variant.jl` 0 / 0 / 2; `an_p2_reported.jl` 0 / 0 / 3; `an_p2_output.jl` 15 / 5 / 5. A testset stops at its first error, so the counts of errors are counts of testsets and test statements reached, not of assertions; the full suite has about 680 assertions once the module exists.

Statement for each file: every failure and error of the new files is of the missing-functionality kind (a name of the proposed API that is not defined, or an `analyze` call with the keywords of the P2 study, or a check of behaviour that the amendment adds). No parse error, no error inside a helper or a fixture (the first RED run found three defects of this kind, all corrected: a helper defined in a later file, an unsorted multidimensional comprehension, a file read after a failed command; the validation below found more, all corrected, among them a weak hash of the file bytes in the read-only checks, a raw regex substitution and a flag value that the library had not yet written to the disk; the last two RED runs are clean). The passing new tests are the self-test, the guards of the unchanged P1 path, the checks of behaviour that the P1 command already has, and counts of cases.

How the tests themselves were validated (they cannot be run against the module): a throwaway implementation of the proposed API, assembled from the reference implementation (`p2ref_*`) and kept outside the repository (the session scratch area), was loaded into `Chaoslab.Analysis`; with it all pure-function files (322 assertions) and the pipeline-level files (261 of 276; the 15 that remain are the dispatch refusals and `cli`, which the throwaway code does not implement) pass, the 54 refusal mutations all change the files they name and are restored, and the reader and parameter-table tests pass. This validates the plumbing, the hand values and the field names of the contract, not the module.

Time. Once the module exists the D89 part of the suite is expected to take 4 to 6 minutes (the target of 3 minutes is exceeded): two P2 studies of 33 chunk files, each with its P1 reference of 55 (about 35 s), about ten full runs of the P2 `analyze` that each re-read the 88 chunk files and run the P1 labelling (default 1, determinism 3, command 2, overrides 1, variant 1, accepted run after the refusals 1, plus the P1 run of the label comparison), 54 copies of the study for the refusal cases, and 1000 resamples per bootstrap label. The run count can be reduced at stage 03 if the owner prefers; no assertion depends on it.

Checkpoint (stage 02): present this section and the RED run to the owner; the owner decides whether to approve the tests or add cases, and the owner is asked for the decision on ambiguity A1 (the slug of the base fixture) before stage 03.

### Owner decisions at the D89 stage 02 checkpoint (2026-10-29)

- A1: stage 03 changes the base fixture slug of the P1 tests (`an-fixture`) to one of the two known slugs. E21 stays strict: any other slug is refused.
- A2: a numeric `value` is matched to its quote by Julia's `string(value)`, as the tests pin it.

### Rework after the stage 03 checkpoint (2026-10-29, base commit a3ca6dd)

Two defects of the tests were corrected (the owner decisions are in the spec, "Owner decisions at the D89 stage 03 checkpoint"); `src/`, `data/` and the frozen preregistration are untouched.

1. A2, the `text` field. `an_p2_helpers.jl`: `p2_texts(p; style)` gives the text of each key (section style or Julia's `string(value)`); `p2_quotes`, `p2_write_params` and `p2_prereg_text` carry it (`style` keyword; `edit` with `nothing` removes a field). `an_p2_fixture.jl` checks that every non-machine key has a text that occurs in its quote and that at least 20 texts differ from `string(value)`. `an_p2_params.jl` has the new cases (criterion D89b). No other test built a quote from `string(value)`.
2. E33, the R2 version. `an_p2_refstore.jl`: `p2_ref_all` also returns `Det2` (the detection of every orbit and indicator from the R2 verdicts), `p2_ref_recs(...; rule = :R2)` takes the censored times from `Det2`, `p2_ref_h6(...; rule)` takes the direction, Float64 and Double64 verdicts with the rule. The R1 expectations are unchanged. `an_p2_pipeline.jl`, testset "extension ..., the R2 version": the R2 decisions and the R2 extension (`res.r2_extension`, a field the implementation already has) are compared with the reference, and the hand values of the section "Hand-computed expected values" (M_LE(S2), M_SALI(S2), M_LE(W2), R_LE, delta_dir, delta_prec under R2) are asserted as literals. The hand computation was done from the design; the reference agrees with it.

Runs against the implementation at a3ca6dd (`julia --project -t 4 --check-bounds=yes -e 'using Test, Chaoslab; @testset "an" begin include("test/analysis/runtests.jl") end'`):

| Run | Pass | Fail | Error | Meaning |
|-----|------|------|-------|---------|
| Official (tables with the section-style texts) | 864 | 25 | 20 | RED as expected. The implementation still requires `string(value)` in the quote, so it refuses the parameter table of every fixture study: all tests that need an analysed study stop at their first call (`ArgumentError: ... key ratio_threshold has the value 5.0, whose text does not occur in its quote`), and the loops of `an_p2_params.jl` get the message of `ratio_threshold` instead of the key under test. The other fails are the checks of the `cli` output file that follow, and one dispatch check, all consequences of the same refusal |
| Verification (same files, tables written with `text` = `string(value)`, which the old rule accepts) | 1121 | 11 | 2 | every test that does not concern A2 or the R2 version passes: the 434 P1 tests, the fixture, dispatch, verdicts, detection, agreement, sets, stats, bootstrap, refusals, pipeline (R1, extension, gate, AUC), variant, reported and output tests. The 13 non-passing: 5 asserts of the R2 version (the implementation constant `P2_R2_REDETECT = false` keeps the R1 detection times and H6 verdicts), 2 self-test checks of the text variants and 6 checks of `an_p2_params.jl` that need the section-style texts (this run used the other style) |

The first run is the one to read for the A2 field; the second shows that nothing else is broken. After stage 03 implements the `text` rule and the R2 recomputation, the expected result is the full suite green.


### Review findings (stage 04) (2026-10-29, base commit 9ed8855)

The independent review of the D89 code found defects that the tests missed. The new tests are written from the spec (E22, E24, E26, E31, E36, T16), not from `src/Analysis/p2_*.jl`; `src/`, `data/`, the spec and the frozen preregistration are untouched. New file `test/analysis/an_p2_review.jl`, registered last in `runtests.jl`. Findings 8, 10 and 11 became readings and have no test.

Fixtures. Study "review": the default study with the grid K list swapped (`grid_k_values = ["0.75", "0.5"]`, grid K index 1 = 0.75, index 2 = 0.5) and design variant `:review`. The K values are 0.5, 0.75 and 2 pi. In the fixture the panel orbits have K index 1 = 0.5 and 2 = 2 pi and the grid orbits K index 1 = 0.75 and 2 = 0.5, so panel orbit 1 (index 1) and grid orbit 33 (index 2) share the value 0.5 with different indexes, and panel orbit 1 and grid orbit 21 share the index 1 with the values 0.5 and 0.75 (panel orbit 11 and grid orbit 33 share the index 2 with the values 2 pi and 0.5). Study "review-dg": design variant `:review` on the default grids (K index and K value agree), so that finding 5 and finding 6 are not mixed up with the pooling of finding 1. Variant `:review` makes FLI undecided for the regular orbits 41 and 42 at k = 0 to 2, FLI of the sticky orbit 39 chaotic from k = 3, and nu1 = NaN of the gate orbit 10 at window 320 (E1 and E2 undecided there). Reference code: `p2_ref_Kvals`, `p2_ref_ranges_K`, `p2_ref_auc_K`, `p2_ref_leaving_K`, `p2_ref_below_K` in `an_p2_refstore.jl` (pooling by the binary64 K value, a budget without a range counts as "not above", Y5 left out); `p2_ref_k0` filters and counts non-finite distances; `p2_del_attr` in `an_p2_build.jl`; the `texts` override of `p2_texts`, `p2_prereg_text` and `p2_pair`.

Contract proposed by the tests (where the spec fixes the semantics but not the key form): `res.auc[(X, k, K, source)]`, `res.range[(X, K, k)]` and `res.below_median[(X, K)]` are keyed by the K value (Float64), and the "K index" column of T2 to T10, T13, T14 and F4 holds the K value or "all" (the tests find the column by the name "K" or "K index"). `res.k0line[W]` gains `fq1_nonfinite` and `pa1_nonfinite` (counts of the distances left out). The older tests accept both forms of the key and label (the default study has K index and K value in the same order), so they pass before and after the fix.

| Finding | Spec | Test (testset of `an_p2_review.jl`) | What it pins | Why it is RED now |
|---------|------|-------------------------------------|--------------|-------------------|
| 1 K pooling by value | E31, E24 | "review 1: pooling by K value, not by K index" | AUC rows (n_C, n_R, AUC by source), regular range, leaving budget, share below the median, T9, F4, T14 of the swapped study equal the reference grouped by K value; hand counts for LE at k = 6: K = 0.5 pooled (15, 5), panel (7, 2), grid (8, 3); K = 0.75 (6, 5); K = 2 pi (3, 2) | keys and labels are K indexes: 18 fail (2 guards pass) |
| 1 (tables and decisions) | E31, E24 | "review 1: the by-K rows of T2, T3 to T8, T10 and T13 ..." | by-K rows of T2 equal the reference counts of the set orbits with that K value; hand counts of "orbits" in T3 to T7 (W2 {0.5: 6, 0.75: 5}, S2 {0.5: 6, 2 pi: 2}, T_s {0.5: 2, 0.75: 2}); K labels in T8, T10, T13; the pooled decisions H1 to H6 (R1, R2, extension) and the sets unchanged | labels are K indexes: 11 fail; the decisions and sets (6 tests) pass now and must keep passing. T8, T10 and T13 have no by-K rows in these studies, so their label test is a guard only |
| 2 infinite percentile | E31 | "review 2: type 7 percentile with infinite neighbours is never NaN" | 100 scores with 10 of them -Inf (5th percentile -Inf, 95th 85.05, median 40.5); one -Inf next to a finite value (lower end -Inf, upper 3.8); a 15th percentile with the neighbours -Inf and 1 is not NaN; all-finite value unchanged | -Inf minus -Inf gives NaN: 4 fail |
| 3 stripped attributes | E36 b, d, e, f | "review 3: a chunk with a stripped attribute is refused in P2" | each of `plan_sha256`, `prereg_sha256`, `panel_sha256`, `status`, `study` deleted from one chunk is refused with the attribute named and no output written; the P1 analysis of a chunk without `panel_sha256` and `prereg_sha256` stays accepted | none of the five is refused: 1 fail; the P1 guard passes |
| 4 text against value | E22, I14, A2 | "review 4: a numeric text must equal its value" | the text (and its quote) differs from the value for `persist_divisor` ("8"), `min_set` ("50"), `e1_d_threshold` ("-5"), `h2_q_min` (value 100, text "10"), the last budget, `bootstrap_level` ("0.5"), `bootstrap_resamples` ("1e4"): refused naming the key; the formats "0.2000", "1e3", "0.0001", "-4.0", "5" that equal the value are accepted; "95 percent", "5th" and "10^3" stay unchecked | none of the seven is refused: 1 fail; the accepted formats pass |
| 5 leaving budget | E31, P25 | "review 5: the leaving budget when the regular range does not exist ..." | unit: `p2_leaving_budget` with no range at k = 0 to 2 gives b_3 = 320, a gap at k = 5 gives b_6, a gap at k = 9 or no range at all gives nothing; end to end (study "review-dg"): the regular set at K = 2 pi is {15, 16, 41, 42, 43}, FLI has no range at k = 0 to 2 and five scores from k = 3, orbit 39 is above from k = 3: `res.leaving[(39, :FLI)] == 320` | the pipeline gives nothing: 1 fail; the unit tests (4) pass |
| 6 K = 0 line | T16 | "review 6: non-finite nu1 distances of the K = 0 line are filtered and counted" | window 320 of the gate has one NaN distance for fq1 and one for pa1: the largest distance is finite (4e-10 and 1.2e-9, equal to the reference), the counts `fq1_nonfinite`, `pa1_nonfinite` are 1 at that window and 0 at the others | the maximum is NaN and the count fields do not exist: 3 fail, 2 error |
| 7 unknown kind and bit names | E36 flags | "review 7: an unknown indicator kind and wrong bit names are named refusals" | kind "bogus" is an ArgumentError naming the kind; the name of the bit is verified at the positions failed (frequency table, position 3) and reldig_undefined (Birkhoff table, position 3), and also at degenerate, zero and undefined window (position 2, guards) | `KeyError: key "bogus" not found`; renamed `failed` and `reldig_undefined` are accepted (for reldig_undefined the manifest count is renamed too, so that the count recomputation does not refuse the chunk): cases kind, bit-failed, bit-reldig fail; the three position 2 cases pass |
| 9 control without deviation | E36 n | "review 9: a control item without control_deviation is refused" | controls 17 and 20 (positive exponent, bound scope) without the attribute are refused naming `control_deviation` | not refused: 1 fail |
| 12 growth ratio overflow | E26, E31 | "review 12: a growth ratio that overflows to Inf ..." | `p2_ratio(1e308, 1e-10)` is Inf, the verdict is chaotic and `p2_score` is nothing for LE, SALI, FLI, MEGNO; for -1e308 regular and no score; RLI = 0 keeps -Inf | the score is Inf: 5 fail; the verdicts pass |

What the tests cannot cover. Finding 12 is a unit test (no store fixture produces a ratio that overflows with the real indicators). The `^` texts of finding 4 are not compared (a form such as "10^4" is not a pure number). That a refusal names the attribute, key, kind or bit is a reading of "a named refusal" (E36); the words used are the attribute or key name, `control_deviation`, the kind, or the expected or renamed bit name. Finding 5 end to end uses the default grids (study "review-dg") so that it fails for its own reason and not for the pooling of finding 1.

RED log (`julia --project -t 4 --check-bounds=yes -e 'using Test, Chaoslab; @testset "an" begin include("test/analysis/runtests.jl") end'`, 3 min 52 s):

| Run | Pass | Fail | Error | Total |
|-----|------|------|-------|-------|
| Before (HEAD 9ed8855 plus the working tree, folder suite) | 1134 | 0 | 0 | 1134 |
| After (new file added, four old testsets made tolerant of the K key form) | 1165 | 46 | 2 | 1213 |

The 79 new tests: 31 pass (guards and unchanged behaviour listed in the table), 46 fail and 2 error, all in `an_p2_review.jl` and for the reasons of the table. All 1134 old tests pass; the testset "D89 default scenario: AUC, regular range, leaving budget" keeps its 24 tests. After stage 03 fixes the findings the expected result is the full suite green (1213).
