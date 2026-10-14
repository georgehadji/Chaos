# Test plan: lyap-gram-schmidt

Spec: `workspaces/module-build/stages/01-spec/output/lyap-gram-schmidt-spec.md` (D30, as amended by its Amendment D35 for AC10; the Double64 coverage note is there too). API: the Batch B contract `api-batch-b.md`, final under D34 (`QRSpectrum(CGS()/MGS()/RGS())`, `orthonormalize`, `window_G`, `window_Gw`, `orthogonality_threshold`, `safe_norm`, `Benettin`, `effective_eps`, the seven-field `NumericalHealthError`, `reducer_provenance`). Float32 is dropped (D24).

Files in `test/lyap-gram-schmidt/`: `runtests.jl` (outer testset, includes the rest; the outer testset keeps one failing criterion from aborting the others when the folder is included alone), `gs_helpers.jl` (shared with `test/lyap-svd-stabilized`, which includes it through a guard variable), `gs_api.jl`, `gs_ac01_05.jl`, `gs_ac06_09.jl`, `gs_ac10_14.jl`. Helpers are prefixed `gs_` and are functions only: no Lyapunov name is touched at include time, so a missing `Chaoslab.Lyapunov` errors the testsets and never the include.

Seeds: StableRNG / run seed = 20261012 + 20000 + 100 x criterion (`gs_seed(ac)`, range 20100 to 21400 for this slug; the SVD slug uses 30100 to 31200; no Batch A offset is reused). Starts: x, p uniform in [0, 2 pi) drawn in Float64 and widened (lyap-benettin S2, S7). Element types Float64 and Double64 where the spec says so; theta_sum = 16 x effective_eps(T) x (G_max / tau + M_max) (G13) is formed in the test: G_max from the estimator final (d = k = 2, where it is the exact condition ratio), or, for the 4-dimensional linear map A_S, from `cond(A_S^tau)` (the exact G of every window, API Q5); M_max is replayed from the stored orbit (largest tangent entry) or is the matrix maximum.

## Criterion to testset

| Criterion | Testset | Sample sizes | Tolerances |
|-----------|---------|--------------|------------|
| (API) orthonormalize | "API1 orthonormalize contract" | Float64, shapes 2x2 and 4x3, 4 well conditioned frames each, 3 variants; zero and NaN column (4x3); `@allocated` | W = Q R to 16 k eps, Q'Q = I to 8 k eps (RGS) or 4x that (CGS, MGS), exact zeros below the diagonal, R_ii > 0, same R for all variants (G3), status `QR_OK`; zero or NaN column gives `QR_PIVOT` with column 3; 0 bytes |
| (API) conditioning, thresholds | "API2 window conditioning and thresholds" | 12 (Float64) and 2 (Double64) upper triangular 2x2, one 3x3 | `window_Gw` = ‖R‖_F ‖R^-1‖_F and `window_G` = kappa (256-bit reference), `G_w >= G`, k = 1 gives 1, k = 3 `window_G >= cond`; `orthogonality_threshold` = tol eps_eff k g to 1e-12 relative (g = G_w for CGS, MGS); `effective_eps` |
| (API) safe_norm, orth_loss | "API3 safe_norm and orth_loss" | exponents -1000..1000, 5 vectors; one Float64 matrix | relative error <= 2 eps (L1 AC8); zero, NaN, Inf; `orth_loss == eps_orth` bitwise |
| AC1 cat map | "AC1 cat map" | Float64 with 3 variants, Double64 with RGS, tau in {1, 2, 4}, n in {40, 400}, 20 seeded frames (orbit ids), z_0 = 0; eigenbasis frame | lambda_1 vs G12a (256-bit), lambda_1 + lambda_2, eigenframe lambda_1 vs ln mu: each <= theta_sum; RGS/MGS/CGS orthogonality bound (AC6a, b) |
| AC2 K = 0 | "AC2 K = 0" | Float64, 3 variants, n in {1, 10, 1000, 1e4}, tau = 1, frame (e_2, e_1) | S_1 vs (1/2) ln(n^2 + 1) and S_1 + S_2 <= 8 n eps; lambda_1(1e4) <= ln(1e4)/1e4 + 1e-12 |
| AC3 bitwise structure | "AC3 bitwise structure" | (a) k = 2, K in {0.5, 2, 5, 10}, 100 orbits, n = 1200, tau in {1, 3}: sum, comp, final Q; 600 random d x 2 matrices ((Float64, d = 2), (Float64, d = 4), (Double64, d = 2); 200 each) through `orthonormalize` for Q and R. (b) A_S, k = 4, tau = 3, 200 windows | tolerance 0 (`isequal`); (b) Q differs in at least one entry (own CGS and MGS chains; also the stream with `orthogonality_tol = Inf`) |
| AC4 256-bit recomputation | "AC4 256-bit recomputation" | Float64: K in {0.5, 2, 5, 2 pi}, 20 orbits, tau in {1, 2, 4}, n = 1000, 3 variants. Double64 (RGS only): K in {2, 2 pi}, 4 orbits, tau in {1, 4} | lambda_1, lambda_2 vs the 256-bit RGS chain on the stored orbit points (same widened Q_0, tangent at 256 bits at the same points; one chain per orbit, tau independent, one logarithm of the pivot products) <= theta_sum |
| AC5 cross variant and tau | "AC5 cross variant and cross tau" | K in {0.5, 2, 5, 10, 2 pi}, 50 orbits, n = 1000, 9 (variant, tau) pairs (36 pairs); A_S, k = 4, n = 600, tau in {1, 2}, 10 frames | pair bound theta^A + theta^B for every i; a CGS run of A_S that raises the orthogonality error is reported (tau = 2 only) and not scored; MGS, RGS and CGS at tau = 1 must run |
| AC6 orthogonality | "AC6 orthogonality" | (a), (b): K in 5 values, tau in {1, 2, 4}, 20 orbits, n = 1000, k = 2 (also embedded in AC1, AC4, AC5, AC8); (c): A_S, k = 4, 200 windows, 20 frames, tau = 3 and tau = 1 | (a) `orth_max <= 8 k eps`, (b) `<= 8 k (G_max + 1) eps` (G_w <= G + 1); (c) tau = 3: CGS raises "orthogonality" in 20 of 20 (column 0, step multiple of tau), MGS and RGS in 0 of 20, loss ratios >= 100 and >= 100 (loss measured with `orthogonality_tol = Inf`); tau = 1: loss(CGS) > loss(MGS) > loss(RGS) in 20 of 20 |
| AC7 breakdown | "AC7 breakdown" | Float64 A^20 (3 variants) and Double64 A^40 (RGS), tau = 1, 200 seeded frames each; A^16 / A^32 no breakdown, 20 frames, n = 200; injected zero column (4x4 identity map, k = 4, column 3 and column 1; Float64, all variants), through the per-orbit body and through `reducer_update`; NaN in a 4x4 W | kind "orthonormalisation breakdown", `orbit`, `step` = 1, `column`, `backend`; no-breakdown runs satisfy the orthogonality bounds, CGS = MGS bitwise; NaN: inside the stream tangent-core's finiteness check fires first ("deviation matrix", step = tau), the estimator alone (`reducer_update`) raises the breakdown with threshold Inf and the column of the NaN |
| AC8 wrapped state | "AC8 wrapped state" | WrappedStandardMap, K in {2, 5}, tau in {1, 3}, n = 1200, 3 variants, 20 orbits; 10^6 steps, K = 5, 8 orbits, tau = 1, RGS | whole final record, x and W bitwise equal to a plain per-orbit loop driving `reducer_step` by hand, and S_i bitwise equal to an `orthonormalize` + `add` loop; 10^6 run: finite, `|lambda_1 + lambda_2| <= theta_sum` (M_max <= 1 + K), `@allocated propagate_orbit` = 0 |
| AC9 thread, chunk, allocation | "AC9 thread and chunk invariance, zero allocation" | K = 2 pi, n = 1000, k = 2, tau = 1, chunk in {64, 1000, default} with 4096 orbits for RGS; chunk in {64, default} with 512 orbits for CGS and MGS; A_S k = d = 4, 10^4 steps | `propagate` equals the serial map of `propagate_orbit` bitwise (x, W, finals); 0 bytes |
| AC10 Benettin | "AC10 agreement with Benettin" | K in {0.5, 2, 5, 2 pi}, 50 orbits, n = 1e4, tau in {1, 4}, k = 1 for the 3 variants and k = 2 on the first 10 orbits | bound 2 B_S(n)/n with B_S of `lyap-benettin` S3 as amended (`gs_BS`, K_s = abs(K) for the standard map; Amendment D35); the k = 1 run is also bitwise equal to `Benettin` (contract: shared `safe_norm` and compensated sum) |
| AC11 volume | "AC11 volume and sum identities" | A_S, k = 4, tau in {1, 2, 3}, n = 60, 10 frames, Float64, 3 variants | sum of the first p exponents vs (1/n) ln vol_p (exact rationals, exact Gram determinant, logarithm at 256 bits) <= p theta_sum, p = 1..4; total <= theta_sum; CGS at tau >= 2 may raise "orthogonality" (reported) |
| AC12 CPU vs GPU | stage 04 | | [local], tangent-core AC12 tolerances; fallback reasons verbatim |
| AC13 oracle package | stage 04 | | [local], ChaosTools full-spectrum routine; lambda_1 pathwise from lyap-benettin AC13 |
| AC14 provenance | "AC14 provenance" | 3 variants x 3 `orthogonality_tol` values for `reducer_provenance`; propagate with each variant (Float64); with a chained `ResidualMonitor` | exact equality of the NamedTuple and of the record fields; no `iteration_path`, `implementation` or `d` key |
| Stage 04 only | (no testset) | | Double64 parts of AC2, AC10, AC11 (not in the package suite, see "Shrunk" below); G8 sandwich at n = 12 and 512 bits (AC11), measured loss / (eps G_w) ratios and pass histograms, scratch constants of Tolerances and AC4/AC6/AC7 (RGS minimum eps G_w 1.74, 19 of 200 fifth passes), the O3 constant, recorded line R3, the vacuous-check band of P1 (eps^-1/2 < G_w < 1/(2 eps)) |

## Finding for the owner (library, not a test issue)

`exact(Double64, n)` for an `Integer` above 2^53 returns `Double64(n)`, which rounds through Float64 (`lo = 0`): the precision contract ("integers ... exact") does not hold for the A^40 entries of AC7. The tests avoid it with `gs_exactint`; stage 03 may want to fix `exact` or state the limit.

## Spec ambiguities resolved

- Orthogonality bound for CGS and MGS (AC6b): the estimator records `G_max` (the condition ratio G), not the per-handoff G_w; at d = k = 2 `G_w` lies in `[G, G + 1]` (G10), so the test uses `8 k (G_max + 1) eps`. RGS uses `8 k eps`.
- `orth_max` of the final record is the maximum over all handoffs (API: "recorded every handoff"), so one comparison covers "at every handoff". The same quantity with `orthogonality_tol = Inf` measures the loss in AC3(b) and AC6(c), as the spec allows.
- AC3 "Q and R bitwise": R is not part of the final record; the pure `orthonormalize` is compared on 600 seeded random matrices (CGS vs MGS, k = 2), and the stream is compared on S_i (sum and compensation) and final Q.
- AC4 reference: the 256-bit RGS chain depends on neither tau nor the variant (exact S_i are tau independent, G7), so one chain per orbit serves all 9 runs; the pivots are multiplied and one logarithm is taken at the end.
- AC5 and AC11, Y5 pair: the estimator's `G_max` is G_w for k >= 3 (API Q5, an upper bound), so theta_sum is formed from `cond_2(A_S^tau)`, which is exactly G of every window (R_j has the singular values of A_S^tau Q_{j-1}).
- AC7 Double64: the entries of A^40 exceed 2^53, and `Double64(::Int64)` rounds through Float64; `gs_exactint` builds hi + lo, otherwise the matrix is not the integer matrix and no breakdown occurs. If A^40 does not break down at GREEN, the spec's fallback power 44 is the change (conditional, not applied: the stub, a throwaway estimator written from the contract, breaks down in 200 of 200 frames at A^40 with eps_eff G_w >= 19). The threshold of the zero-pivot case is not asserted (Q11); the NaN case asserts `isinf(threshold)`.
- AC7 NaN: through `propagate` tangent-core fires first ("deviation matrix"); the spec asks to report which kind fires, the test asserts it.
- AC10: the spec states the bound and expects 0 deviation; the contract says the k = 1 run is bitwise equal, so both are asserted.
- AC14: `QRSpectrum(...).orthogonality_tol` is a Float64 (`8.0` default); the record compares with `===`.
- Sweeps use the per-orbit body (`propagate_orbit`) through `gs_serial`, which is bitwise the stream (AC9) and compiles one kernel specialisation less per variant and type; `propagate` itself is exercised in AC8, AC9, AC14. This keeps the compiled specialisations to about 25.
- Shrunk against the spec for the 3-minute budget (the suite time is almost entirely compilation of type specialisations, so the cuts remove element-type and variant combinations, not samples): Double64 runs only AC1 (cat map), AC3 (direct `orthonormalize`, d = 2), AC4 (4 orbits, RGS) and AC7 (A^40, A^32), all with RGS; AC2, AC5, AC6, AC8, AC9, AC10, AC11 and AC14 run Float64 only, so the Double64 parts of AC2, AC10, AC11 are not covered in the package suite and should be added to the stage 04 list; AC9 uses 512 orbits for CGS and MGS (RGS keeps 4096) and omits chunk 1000 for them; AC7 injected cases run Float64 only. Everything else follows the spec sample sizes.

## Runtime

Against a throwaway stub (a copy of `src/` plus a 250-line `Lyapunov` module written from the contract, kept in the scratchpad, never in the repo) the whole folder takes about 195 s of CPU (wall 190 s, 4 threads, a shared machine). About 80 % of that is compilation of type specialisations (per testset, CPU / compile seconds: API1 30 / 34, API2 7 / 8, API3 2 / 2, AC1 26 / 27, AC2 2 / 2, AC3 11 / 10, AC4 18 / 10, AC5 16 / 13, AC6 5 / 4, AC7 5 / 5, AC8 21 / 11, AC9 28 / 10, AC10 14 / 6, AC11 5 / 5, AC14 8 / 8); the stub allocates in every step (the real kernel will run the step loops faster) and the real compile cost is unknown until stage 03. By this stub measurement the folder is OVER the 3-minute budget (about 190 s wall, roughly 80 % compilation); the stub is not the real kernel, so compile and run cost both move at GREEN in an unknown direction: re-measure at GREEN and cut then. If it is over, the cheapest cuts are the API1 to API3 units, the Double64 parts of AC1 and AC4, chunk 1000 for RGS in AC9, and the 10^6-step run of AC8 (move to stage 04); only the last three touch a spec criterion.

Stub results (what the tests were validated against, not evidence about the real code): every package criterion passes against a contract-conforming implementation except the checks a stub that allocates cannot meet (zero allocation in API1, AC8, AC9) and the NaN-threshold check of AC7 (the stub does not set `Inf`). The helper references were checked in isolation as well: exact Gram determinant against a Float64 determinant at n = 1, the cat-map G12a formula against `ln |A^n q| / n`, the 256-bit chain against the stub.

## RED log

Re-run after the Amendment D35 edits (`gs_BS`, AC10 pointer): 0 passed, 0 failed, 81 errored in 15 testsets (57 s), the same names (`CGS` 13, `window_G` 29, `window_Gw` 28, `safe_norm` 8, `effective_eps` 1, `Benettin` 1, `Chaoslab.Lyapunov` 1); no error from a helper or the include.

Command (repo root, 4 threads, test env with `Chaoslab` developed): `julia -t 4 --project=<scratchpad>/m2/testenv -e "using Test, Chaoslab; include(raw\"test/lyap-gram-schmidt/runtests.jl\")"`. Log: `scratchpad/m2/red-lyap-gram-schmidt.log`. Result: 0 passed, 0 failed, 81 errored in 15 testsets (42 s). Every error is an `UndefVarError` on a name the contract adds (`CGS`, `window_G`, `window_Gw`, `safe_norm`, `effective_eps`, `Benettin`, `Chaoslab.Lyapunov`); no error comes from a helper or from the include.

| Testset | Errors | First error |
|---------|--------|-------------|
| API1 orthonormalize contract | 2 | UndefVarError `CGS` |
| API2 window conditioning and thresholds | 58 | UndefVarError `window_G` (29), `window_Gw` (28), `effective_eps` (1) |
| API3 safe_norm and orth_loss | 9 | UndefVarError `safe_norm` (8), `Chaoslab.Lyapunov` (1) |
| AC1 cat map | 1 | UndefVarError `CGS` |
| AC2 K = 0 | 1 | UndefVarError `CGS` |
| AC3 bitwise structure | 1 | UndefVarError `CGS` |
| AC4 256-bit recomputation | 1 | UndefVarError `CGS` (after the 256-bit orbit references, 9 s) |
| AC5 cross variant and cross tau | 1 | UndefVarError `CGS` |
| AC6 orthogonality | 1 | UndefVarError `CGS` |
| AC7 breakdown | 1 | UndefVarError `CGS` |
| AC8 wrapped state | 1 | UndefVarError `CGS` |
| AC9 thread and chunk invariance, zero allocation | 1 | UndefVarError `CGS` |
| AC10 agreement with Benettin | 1 | UndefVarError `Benettin` |
| AC11 volume and sum identities | 1 | UndefVarError `CGS` (after the exact volume references, 9 s) |
| AC14 provenance | 1 | UndefVarError `CGS` |

Testset count: 15 (3 API units, 12 criteria: AC1 to AC11 and AC14). Deferred to stage 04: AC12, AC13, the G8 sandwich of AC11, the scratch reproductions listed above.

