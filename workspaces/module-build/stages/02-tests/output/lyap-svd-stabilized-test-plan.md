# Test plan: lyap-svd-stabilized

Spec: `workspaces/module-build/stages/01-spec/output/lyap-svd-stabilized-spec.md` (D30, as amended by its Amendment D35 for AC3(c) and AC5; its definitions of theta_sum, the window factorisation, the panel and the identical orbit point to `lyap-gram-schmidt`). API: the Batch B contract `api-batch-b.md`, final under D34 (`SVDSpectrum`, `log_singular_values`, the reducer final `(n, sum, comp, lambda_qr, lambda, c, u, G_max, orth_max, count_max)`, `reducer_provenance`, `QRSpectrum(RGS())` for the bitwise comparison). Float32 is dropped (D24).

Files in `test/lyap-svd-stabilized/`: `runtests.jl` (outer testset; includes the shared helper file `test/lyap-gram-schmidt/gs_helpers.jl` once, through its guard variable, so the folder also runs alone), `svd_helpers.jl` (`sv_` helpers, the dimension-4 test map `sv_Boom4`, the exact guard message), `svd_ac01_05.jl`, `svd_ac06_12.jl`. No Lyapunov name is touched at include time.

Seeds: run seed and StableRNG = 20261012 + 30000 + 100 x criterion (`sv_seed(ac)`, range 30100 to 31200 for this slug; the GS slug uses 20100 to 21400; no Batch A offset is reused). Starts as in `lyap-benettin` S2 and S7. theta_sum = 16 x effective_eps(T) x (G_max / tau + M_max) (GS G13), G_max from the estimator final (exact for d = k = 2), M_max replayed from the stored orbit or the matrix maximum; the QR-family comparisons use the pair bound theta^A + theta^B.

## Criterion to testset

| Criterion | Testset | Sample sizes | Tolerances |
|-----------|---------|--------------|------------|
| (API) V6 closed form | "API V6 log_singular_values" | Float64: 4 values of c x 20 random (S_1 in [0, 15], S_2 in [-15, 0]); Double64: c = 1.27 x 4; S_1 = 1e5 | (ln sigma_1, ln sigma_2) vs the 2x2 singular values at 256 bits to 16 eps (1 + |S_1| + |S_2| + ln(1 + c^2)); finite and ordered for S_1 = 1e5, sum S_1 + S_2; sigma_1 = sigma_2 (S = 0, c = 0) is finite |
| AC1 symmetric map | "AC1 symmetric map" | cat map, T in {Float64, Double64}, tau in {1, 2, 4}, n in {40, 400}, 20 seeded frames (orbit ids), z_0 = 0 | lambda_1^SVD vs ln mu and lambda_2^SVD vs -ln mu <= theta_sum; lambda_1^QR vs G12a (256-bit) <= theta_sum; any two frames and tau differ by <= theta^A + theta^B; `lambda_1 >= lambda_2` |
| AC2 K = 0 | "AC2 K = 0" | T in {Float64, Double64}; frames (e_2, e_1), (e_1, e_2), one seeded; tau = 1; n in {1, 10, 100, 1000, 1e4} | lambda_1^SVD vs ln((n + sqrt(n^2 + 4))/2)/n (256-bit), lambda_2 + lambda_1 <= theta_sum; frame (e_2, e_1): SVD minus QR vs [ln sigma_1 - (1/2) ln(n^2 + 1)]/n <= theta_sum; lambda_1(1e4) <= ln(1e4)/1e4 + 2e-12 |
| AC3 identity and sandwich | "AC3 identity with Gram-Schmidt and sandwich" | Float64, K in {0.5, 2, 5, 2 pi}, tau in {1, 2, 4}, 50 orbits, n = 1000 (and the n = 500 prefix runs) | (a) sum, comp, lambda_QR bitwise equal to `QRSpectrum(RGS())` on the same orbit and frame; (b) V7: -theta_sum <= lambda_1^SVD - lambda_1^QR - ln(1 + c^2)/(2n) <= u/(2n) + theta_sum, `lambda_1^SVD >= lambda_1^QR`; (c) classes by S_1(500) over every orbit of the panel (Amendment D35): >= 25: `|c_1000 - c_500| <= 8 eps max(1, |c|)`; 5 to 25: `<= 0.06 max(1, |c_500|)`; < 5: reported, not scored (class counts and the worst relative change printed per K, tau) |
| AC4 256-bit recomputation | "AC4 256-bit recomputation" | Float64: K in 4 values, 20 orbits, tau in {1, 2, 4}, n in {12, 20, 100, 1000}. Double64: K in {2, 2 pi}, 4 orbits, tau in {1, 4} | lambda_1, lambda_2 vs the 256-bit tangent product at the stored points with det Y = 1 used analytically (sigma_1 = (s + sqrt(s^2 - 4))/2, s = sqrt(‖Y‖_F^2 + 2), sigma_2 = 1/sigma_1) <= theta_sum |
| AC5 negative control | "AC5 negative control and overflow" | K = 2 pi, one seeded orbit (the first of the seeded starts with ln sigma_1(20) >= 25, Amendment D35), tau = 1, n = 20 and 2000 | `svdvals` of the plain Float64 product at n = 20: \|ln sigma_2^naive - ln sigma_2^ref\| > 10; the plain product has a non-finite entry before n = 700; the estimator is finite and within theta_sum of the 256-bit value at n = 20 and 2000 |
| AC6 degenerate | "AC6 degenerate and regular cases" | rotation by 0.7311, T in {Float64, Double64}, tau in {1, 10}, n = 1e4, 4 frames | \|lambda_i\| <= theta_sum, nothing non-finite, no health error (K = 0 is in AC2) |
| AC7 pairing and sum | "AC7 pairing and sum" | cat map, tau = 1, n = 1e5, 4 frames, two types; K in 4 values, 10 orbits, n = 1e4 | \|lambda_1 + lambda_2\| <= theta_sum, equal to (S_1 + S_2)/n within 8 eps, `eps_pair == eps_sum` (Diagnostics), `lambda_1 >= lambda_2` |
| AC8 wrapped state | "AC8 wrapped state" | WrappedStandardMap, K in {2, 5}, tau in {1, 3}, n = 1200, 20 orbits; 10^6 steps, K = 5, 8 orbits, tau = 1 | whole final record (S_i, compensation, c, lambda) and x, W bitwise equal to a plain per-orbit loop driving `reducer_step`; 10^6 run: finite, `\|lambda_1 + lambda_2\| <= theta_sum` (M_max <= 1 + K), 0 bytes allocated after warm-up |
| AC9 thread, chunk | "AC9 thread and chunk invariance" | 4096 seeded orbits, K = 2 pi, n = 1000, tau = 1, chunk in {64, 1000, default} | `propagate` equals the serial map of `propagate_orbit` bitwise (x, W, finals incl. c); 0 bytes per-orbit body |
| AC10 CPU vs GPU | stage 04 | | [local], tangent-core AC12 tolerances plus \|c_cpu - c_gpu\| <= 1e-12 |
| AC11 oracle package | stage 04 | | [local], ChaosTools full-spectrum routine |
| AC12 guard, provenance | "AC12 scope guard and provenance" | (d, k) in {(4, 4), (4, 2), (4, 1)} on a map that throws if stepped, and (2, 1) on the cat map; provenance for three tol values (Float64), with a chained `ResidualMonitor` | `ArgumentError` with the exact message, thrown by `propagate` and by `reducer_init`; `reducer_provenance` and the run record equal `(estimator = "svd", orthonormalizer = "RGS", d = 2, orthogonality_tol)` plus the run fields; no `iteration_path` / `implementation` key |
| Stage 04 only | (no testset) | | Double64 parts of AC3, AC5, AC8, AC9, AC12; measured constant (1/(2n)) ln(1 + c_n^2) over the panel, c of the slowly converging orbits (class S_1(500) < 5 of AC3c, reported only), H1 to H3 of the held general-d algorithm (none scored in M2) |

## Spec ambiguities resolved

- AC3(c): resolved by `lyap-svd-stabilized` Amendment D35. The first draft scored the spec's 0.02 on four frames only, which narrowed the sample. The spec now defines the classes by S_1(500) (>= 25, 5 to 25 with 0.06, below 5 reported with the reason: the increments of c do not decay near elliptic points and accelerator modes), and the test scores every orbit of the panel by its class.
- AC5: resolved by Amendment D35 (the orbit is the first seeded start with ln sigma_1(20) >= 25, the control needs a strongly stretching orbit; the overflow statement is unchanged).
- AC1 "any two frames and any two tau": checked with the pair bound theta^A + theta^B of GS G13 (each value is already within theta of ln mu, so the stronger statement is the pair bound).
- AC3(a): the QR-family part compares `sum`, `comp`, `lambda_qr` of `SVDSpectrum` with `sum`, `comp`, `lambda` of `QRSpectrum(RGS())` (`isequal`, tolerance 0).
- AC12 "k != d" inside d = 4 and "d != 2" with k = 1: the guard message is asserted for `propagate` (before any work, proved by a 4-dimensional map whose `forward` and `tangent` throw) and for `reducer_init` directly. The d = 2, k = 1 case uses the cat map.
- AC4 Double64 and the effective eps: theta_sum uses `effective_eps(Double64)` = 2^-100 (GS G13, D34 Q6). The report states it.
- Sweeps use the per-orbit body (`propagate_orbit`, bitwise the stream by AC9) through `gs_serial` to keep the compiled kernel specialisations low; `propagate` is exercised in AC8, AC9, AC12.
- Shrunk against the spec for the budget (the suite time is almost entirely compilation): AC4 Double64 uses 4 orbits, K in {2, 2 pi}, tau in {1, 4}; the V6 unit runs Double64 on one value of c; AC12 runs the provenance record in Float64 only. The spec header says "Float64 and Double64" for every criterion, but AC3, AC5, AC8, AC9 and AC12 run Float64 only; their Double64 parts are not in the package suite and belong on the stage 04 list (spec Amendment D35, Coverage; owner question 6). AC1, AC2, AC4, AC6 and AC7 run both types.

## Runtime

Against the same throwaway stub (see the GS plan) the whole folder takes about 129 s of CPU (wall 126 s, 4 threads), mostly compilation (per testset, CPU / compile seconds: V6 7 / 10, AC1 22 / 21, AC2 9 / 8, AC3 6 / 4, AC4 14 / 5, AC5 2 / 2, AC6 8 / 8, AC7 8 / 1, AC8 30 / 19, AC9 19 / 4, AC12 4 / 4). The stub passes every package criterion except the two zero-allocation checks it cannot meet (AC8, AC9); the V6 closed form, the V9b values, the V7 sandwich on the seeded panel, the 256-bit references and the c convergence rule of AC3(c) (all 12 combinations of K and tau) were checked against it. Run in one process after the GS folder, part of the compilation is shared.

## RED log

Re-run after the Amendment D35 edits (AC3(c) classes, AC5 comment): 3 passed (the premise assertions of AC5), 0 failed, 11 errored in 11 testsets (21 s), every error an `UndefVarError` (`SVDSpectrum` 10, `log_singular_values` 1).

Command (repo root, 4 threads): `julia -t 4 --project=<scratchpad>/m2/testenv -e "using Test, Chaoslab; include(raw\"test/lyap-svd-stabilized/runtests.jl\")"`. Log: `scratchpad/m2/red-lyap-svd-stabilized.log`. Result: 3 passed, 0 failed, 11 errored in 11 testsets (34 s). The 3 passes are premise assertions of AC5 that name no Chaoslab function (the plain Float64 product loses sigma_2 and overflows); the testset still errors on `SVDSpectrum`. Every error is an `UndefVarError`.

| Testset | Errors | First error |
|---------|--------|-------------|
| API V6 log_singular_values | 1 | UndefVarError `log_singular_values` |
| AC1 symmetric map | 1 | UndefVarError `SVDSpectrum` |
| AC2 K = 0 | 1 | UndefVarError `SVDSpectrum` |
| AC3 identity with Gram-Schmidt and sandwich | 1 | UndefVarError `SVDSpectrum` |
| AC4 256-bit recomputation | 1 | UndefVarError `SVDSpectrum` |
| AC5 negative control and overflow | 1 (3 passes) | UndefVarError `SVDSpectrum` |
| AC6 degenerate and regular cases | 1 | UndefVarError `SVDSpectrum` |
| AC7 pairing and sum | 1 | UndefVarError `SVDSpectrum` |
| AC8 wrapped state | 1 | UndefVarError `SVDSpectrum` |
| AC9 thread and chunk invariance | 1 | UndefVarError `SVDSpectrum` |
| AC12 scope guard and provenance | 1 | UndefVarError `SVDSpectrum` |

Testset count: 11 (1 API unit, 10 criteria). Deferred to stage 04: AC10, AC11, the measured constant, H1 to H3.

