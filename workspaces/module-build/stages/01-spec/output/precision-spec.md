# Spec: precision

Owner checkpoint 2026-10-13 (decision D24): Float32 (binary32) is dropped from M2. Every Float32 item in this spec (element-type lists, the screening tier, the Float32 GPU gate and Float32 test cases) is inactive until revisited at M4; stage 02 writes no Float32 tests.

Library layer: Foundation (element-type handling). Draft for the 01-spec checkpoint, branch m2-p1-modules.

## Scope

- Catalog rows: `_config/catalogs/lyapunov.md` section E ("Precision variants, apply to every row"). No numbered row; section E applies to all rows.
- Decisions: D20 (Double64 is the extended-precision type; BigFloat only for short CPU oracle runs; owner decision, M2 kickoff), D17 (GPU first, CPU fallback; `_config/architecture.md` "Compute"), risks K15 (weak Float64 on the GTX 1050, compute capability 6.1) and K18 (published lambda(K) values are loose sanity checks only).
- Element types specified: Float64 (default), Double64 (double-double, about 32 significant digits), BigFloat (arbitrary precision, working precision in bits recorded), Float32 (screening tier, conditional, see Algorithm 6).
- Out of scope: the estimators (lyap-* specs), the tangent step and the backend resolver (tangent-core), the dataset dtype in the store (Experiments spec, see Open questions).
- Prerequisites: standard-map (the test system) and tangent-core (backend resolution, fallback record, tangent step). Both are M2 peers; this spec's tests run once they are verified.

## Equations

Notation: T is the element type; eps_T is the spacing of T at 1 as reported by the implementation (criteria use the reported value); u_T = eps_T / 2 is the unit roundoff for the round-to-nearest models. Nominal values: Float32 2^-23 = 1.19e-7; Float64 2^-52 = 2.22e-16; Double64 about 2^-104 = 4.93e-32; BigFloat with b bits 2^(1-b), so b = 256 gives 2^-255 = 1.73e-77.

1. Canonical map (`_config/research-question.md`, "Ground truth for variational methods (standard map)"): p' = p + K sin x, x' = x + p', lifted coordinates, wrap only for display. Tangent map DT = [[1 + K cos x, 1], [K cos x, 1]], det = 1.
2. Conversion of the McMillan form of SalaArtuso2015 to the canonical form (derived). The benchmarks row gives (x, y) -> (2x + K sin x - y, x). Set p_n = x_n - y_n. Then p_{n+1} = x_{n+1} - y_{n+1} = 2x_n + K sin x_n - y_n - x_n = p_n + K sin x_n, and x_{n+1} = x_n + p_n + K sin x_n = x_n + p_{n+1}. Start (x, y) = (1, 2) x 10^-3 becomes (x, p) = (10^-3, -10^-3). That row (K = pi/2: 0.298; K = 2 pi: 1.172) is used here only to size horizons (K18).
3. First-order round-off propagation (derived). Let z_n be the orbit in the reference (higher) precision and e_n the difference to the orbit computed in T from the same exactly representable data. Then e_{n+1} = DT(z_n) e_n + d_n, where d_n is the round-off of one step. One step has four rounded operations (sin, multiply by K, two additions), all in magnitudes at most M_n = max(1, max_{k<=n} |x_k|, |p_k|), so |d_n| <= (K + 3) M_n u_T. Hence |e_n| <= (K + 3) u_T M_n Gamma_n, with Gamma_n = sum_{k<n} || DT(z_{n-1}) ... DT(z_{k+1}) ||_2 (empty product = identity), computed on the reference orbit. For a chaotic orbit with finite-time exponent lambda_n, Gamma_n is of order exp(lambda_n n) (within a factor 1/(1 - exp(-lambda_n)), at most 2 for lambda_n >= 0.7). For K = 5 the prefactor is (K + 3) u = 4 eps.
4. Orbit agreement bound (criteria 4, 9): |e_n|_inf <= c eps_T M_n Gamma_n with c = 8 (twice the derived constant 4 at K = 5). Perturbative horizon n_H(T) = largest n with c eps_T M_n Gamma_n <= 10^-4, so nonlinear saturation is excluded. Sizing (loose, K18): Chirikov law lambda about ln(K/2) = 0.92 at K = 5 (benchmarks.md, "Standard map largest Lyapunov exponent, large K", Shevchenko2004), M about 20, giving n_H of about 24 (Float64) and about 63 (Double64). For Float32, 8 x 1.19e-7 x 20 = 1.9e-5 leaves room for exp(0.92 n) <= 5.3, so n_H is about 2: Float32 has no useful perturbative horizon and is not tested pathwise.
5. Consequence of a Float64-rounded parameter (derived). pi/2 and 2 pi rounded to Float64 have relative error 3.90e-17 (absolute 6.12e-17 and 2.45e-16; computed from 60-digit pi). An orbit of the map with that K differs from the true-K orbit by about delta exp(lambda n); it reaches order 1 at n* about ln(1/delta)/lambda. With lambda = 0.298 (K = pi/2): Float64 round-off alone gives about ln(1/1.1e-16)/0.298 = 123 steps; a Double64 run with K widened from Float64 gives 125 steps, i.e. no gain over plain Float64; K computed in Double64 gives about ln(1/2.5e-32)/0.298 = 244 steps. At K = 2 pi (lambda = 1.172) the widened-K run decouples after about 31 steps. Same effect for any parameter or initial condition given as a decimal literal or an expression in pi.
6. Accumulated tangent growth (criterion 6): renormalizing every step in this test, S_n = sum_{k<n} ln(||DT(z_k) u_k|| / ||u_k||) with u_k the unit tangent vector; this is ln||DT^n v_0|| - ln||v_0||. Its Double64-versus-reference difference is bounded in criterion 6.

## Algorithm

1. Element-type request: the run configuration names one precision; the resolved element type is recorded and never changed by a fallback. Values: `Float64`, `Double64`, `BigFloat(b)` (b >= 64 in MPFR sense, significand bits including the leading bit, rounding to nearest, b recorded as an integer inside the string), `Float32`.
2. Conversion rules.
   - Exact data (a Float64 number from a seeded generator, a dyadic number such as 5 or 0.5) widens exactly into T. Seeded initial conditions and tangent vectors are drawn once as Float64 and widened, so that runs at different precision start from the identical point (criteria 4, 6, 9).
   - Symbolic or decimal data (K = pi/2, K = 2 pi, 0.1, 1e-3, any rational) is built in T from its definition: pi from the T constant, decimals parsed from text, divisions done in T. Result: within 1 ulp of the true value (criterion 1). Rounding a Float64 value to T, or T to Float64 and back, is a defect.
   - Narrowing (Float64 for display, plot or store) is allowed only on output, rounds to nearest, and never feeds back into a computation.
3. Evaluation rules. sin, cos, log, exp and sqrt are evaluated in T (their argument and result both T). Every constant inside the step (the 1 in DT, 2, K) is T. Accumulators for log growth and the compensation term are T. No fast-math reassociation (`_config/coding-conventions.md` "Tests"). No operation silently promotes or demotes (criterion 12).
4. Compatibility and fallback (`_config/architecture.md` "Compute"; reasons are the strings listed there). Resolution order: element type first, then device.

   | Element type | CPU | GPU | Reason recorded when the run is on CPU instead of a requested or `auto` GPU |
   |---|---|---|---|
   | Float64 | yes (reference) | yes, slow on cc 6.1 (K15) | device reasons only |
   | Double64 | yes | attempted; must compile for the device | "kernel does not compile for the element type" if it does not, else device reasons |
   | BigFloat(b) | yes, CPU only | never attempted | "kernel does not compile for the element type" (always, with or without a GPU present) |
   | Float32 | yes | yes (screening tier) | device reasons only |

   Device reasons: "CUDA.jl not loadable", "no functional device", "device out of memory", "numerical health error raised only on the GPU". The `precision` attribute is the requested type in every case.
5. Recording per run (`schema/run-store.md`, root attributes `precision`, `backend`, `device`, `backend_fallback_reason`): `precision` is exactly `Float64`, `Double64`, `Float32` or `BigFloat(<b>)`; the config hash includes precision, so two runs differing only in precision never share a hash. BigFloat bits is the working precision of every operation of the run.
6. Float32 admission (screening tier). Argument from K15: on this GPU the Float64 rate is a small fraction of the Float32 rate (K15 states "weak Float64 throughput"; the ratio is to be measured, not assumed), so Float32 is the only element type that can make the GPU beat all CPU threads. Against it: with eps = 1.19e-7 an orbit decouples from its exact counterpart after about 2 steps at the 10^-4 level at K = 5 (Equations 4; saturation to order 1 takes about ln(1/1.2e-7)/0.92 = 17 steps), and a Float32 orbit lives on a finite lattice of about 2^46 points per unit cell (heuristic, derived: a random map on N states has cycles of length about sqrt(N), here about 10^7), so long orbits are not trustworthy. Rule: Float32 is admitted only for GPU ensembles of short orbits used to screen (grids, indicator maps), never as a reported estimate unless the study's preregistration says so; every Float32 store carries `precision = Float32`. It ships only if criterion 10 passes; otherwise it stays specified but not shipped.
7. Oracle runs: BigFloat at b = 256 by default (about 77 digits). A BigFloat run with n steps loses about lambda_max n / ln 2 bits; it is valid only if it agrees with the run at 2b bits to the claimed digits (criterion 5 is the standing check for the default).

## Invariants

1. Round trip: Float64 -> T -> Float64 is the identity for every finite Float64 (T in Double64, BigFloat(b >= 64), Float32 excluded: it narrows).
2. Element-type closure: every operation of the map step, tangent step and log accumulation returns T when fed T.
3. Precision in force: the result of a BigFloat(b) run depends on b only through rounding; a Float64 constant leaking in shows as an error of order 10^-16 Gamma_n in a BigFloat(256) run (criterion 5).
4. Fallback never changes the numerics request: the same configuration gives the same `precision` string on every backend.
5. A Float32 store never backs a reported estimate unless preregistered.
6. Symplectic invariants (det DT = 1 within eps_T scale) are owned by diagnostics, not here.

## Compute

- Data-parallel: ensembles of independent orbits (grids, sweeps) in Float32, Float64 and Double64 are GPU candidates. CPU only: BigFloat (heap-allocated, size grows with b, not a fixed-size inline value), single long orbits, oracle runs.
- Packages (names allowed here): DoubleFloats (Double64 is a fixed-size inline value, unlike BigFloat, so the same kernel runs on both backends: D20), MPFR through BigFloat, KernelAbstractions, CUDA. Whether the Double64 kernel compiles for cc 6.1 and whether CUDA fused multiply-add contraction alters the double-double error-free transformations are untested hypotheses; criteria 9 and 11 test them.
- Memory per orbit (derived count, maximal-exponent run): 6 scalars (state 2, tangent vector 2, compensated log sum 2) plus a step counter: Float32 24 B, Float64 48 B, Double64 96 B. BigFloat is not budgeted (about 6 heap objects whose size scales with b); chunk size from free memory as in `architecture.md`.
- Zero-allocation rule of `coding-conventions.md` ("Tests") applies to Float32, Float64, Double64; BigFloat is exempt (every operation allocates), stated here as the deliberate exception.
- CPU-GPU agreement: orbit level, same ICs, |difference| <= 2 c eps_T M_n Gamma_n over n_H (twice the single-run bound, both sides round); estimator level at the preregistered estimator tolerance, deferred to the lyap-* specs. Backends are never compared bitwise.
- Stage 05 setting: `_config/environment.md`, section "High precision", becomes "Backend for extended precision: Double64 (double-double, about 32 digits, same code on CPU and GPU); BigFloat (256 bits default) for short CPU oracle runs only" (replaces the placeholder).
- Proposed `numerics-defaults.toml` key: `bigfloat_bits = 256` (default working precision of oracle runs, Algorithm 7). `precision = "Float64"` stays; no other key.
- CI has no GPU: the kernel runs on the CPU backend; GPU comparison runs locally and goes into the 04-verify report.

## Oracle

BigFloat at b = 256 (MPFR), same map, same exactly widened data, as the reference for Float32, Float64, Double64; BigFloat at 512 as the reference for 256. Analytic: the map and tangent formulas of Equations 1-2. No external package is needed for the reference.

Placement (owner decision, M2 kickoff): the package test suite holds the small self-contained checks, criteria 1 to 8, 12, 13 (BigFloat is part of the language; at most 200 orbits x 100 steps and 10^4 function points, seconds). The separate verification environment (module-build stage 04) holds the heavy and external items: criteria 9 to 11, and the full estimator comparison across precisions (Benettin lambda_1 at Float64, Double64, BigFloat, with the ChaosTools maximal-exponent routine as cross-method check), which is deferred to the lyap-benettin stage 04. This spec's criteria 6 and 9 are pathwise only; long-orbit estimator agreement is a statistical statement not made here (K18: no published lambda(K) is a hard target).

## Acceptance criteria

Common data unless stated: K = 5 (exactly representable, so no conversion enters), 200 seeded ICs uniform in [0, 2 pi)^2 as Float64 (seed from `numerics-defaults.toml`), reference BigFloat(256). "Suite" = package tests, "Verify" = stage-04 environment.

1. (Suite) Parameter conversion. For K in {pi/2, 2 pi, 0.1, 1e-3} and T in {Float32, Float64, Double64, BigFloat(256)}: |K_T - K_ref(512 bits)| <= 0.5 ulp_T for Float32, Float64 (correctly rounded) and <= 2 eps_T K for Double64 and BigFloat. Fails if any value was built through Float64 (Double64 error would be 3.9e-17 K).
2. (Suite) Exact widening. 10^4 seeded Float64 values, including subnormals and values near 2 pi, round-trip Float64 -> T -> Float64 bitwise for T in {Double64, BigFloat(256)}, and equal the BigFloat value exactly.
3. (Suite) Transcendentals in T. For 10^4 Double64 points with all 106 mantissa bits random: sin and cos on [-2 pi, 2 pi] and on [-10^3, 10^3], log on log-uniform [10^-3, 10^3]. Pass: |f_D - f_ref(256)| <= 16 eps_D max(1, |f_ref|) for every point. The 16 ulp bound is a proposal; the library's real accuracy is not documented in the sources read (Open questions). A Float64-evaluated function gives errors near 10^-16 and fails by 16 orders. Same test at Float32 on CPU, bound 4 eps_32.
4. (Suite) Pathwise orbit agreement. T in {Float64, Double64} versus BigFloat(256) over n_H(T) (Equations 4): max over the 200 orbits and n <= n_H of |e_n|_inf / (eps_T M_n Gamma_n) <= 8. Floors so the test cannot pass vacuously: realised n_H >= 15 (Float64), >= 40 (Double64) for every orbit; orbits with a smaller n_H are reported and at most 5 percent may be excluded.
5. (Suite) Reference self-consistency, detects precision leakage. BigFloat(256) versus BigFloat(512), n = 100: |e_n|_inf <= 8 eps_256 M_n Gamma_n for all n, and |e_10|_inf <= 10^-60 (a Float64 constant inside the 256-bit run would give about 10^-12 at n = 10).
6. (Suite) Tangent growth. For the same 200 orbits, unit tangent vector from a seeded Float64 pair, renormalized every step in this test (independent of the tangent-core default): for n <= min(n_H(Double64), 40), |S_n(Double64) - S_n(ref)| <= 64 (1 + K) eps_D M_n Gamma_n. The constant 64 is a first-order proposal (Equations 6 has no derivation of the local stretch factor); stage 04 reports the realised ratio and the owner may tighten it. Full estimator comparison across precisions: lyap-benettin stage 04.
7. (Suite) Fallback matrix. For every (T, requested backend in {auto, gpu, cpu}) the resolved backend and reason equal the table in Algorithm 4 (CI: "CUDA.jl not loadable" or "no functional device" for GPU-capable types; BigFloat with an auto or gpu request: "kernel does not compile for the element type" regardless of device presence; any cpu request: empty reason); the stored `precision` equals the request in all 12 cases.
8. (Suite) Recording. For each T the `precision` string matches exactly one of the four forms; BigFloat bits parse back to the integer b. Store-side checks (a store whose `precision` mismatches the request raises an error; config hashes differ when only precision differs) belong to the Experiments spec (M3), not here.
9. (Verify, local GPU) Cross-backend orbit agreement. GPU versus CPU for Float64, Double64 (if it compiles) and Float32 (if admitted; n <= 2 only): same ICs, n <= n_H(T), max |difference| <= 2 c eps_T M_n Gamma_n with c = 8. If Double64 does not compile, the recorded outcome is the fallback reason of criterion 7, and this is a pass of the fallback, not of the kernel.
10. (Verify, local GPU) Float32 gate. Throughput of the Float32 GPU ensemble (10^5 orbits x 10^3 steps) on the GTX 1050 is >= 4x the Float64 CPU all-threads throughput on the same ensemble. Below 4x, Float32 is not shipped (the threshold is a proposal). The Float64 and Double64 GPU/CPU ratios are recorded as well (K15).
11. (Verify) Performance, reported not gated. Time ratio Double64/Float64 on CPU (1 thread and all threads) and on GPU where it compiles, plus Float32/Float64; same ensemble as criterion 10. Pass = the ratios, the hardware and the commit are in the 04-verify report; no threshold.
12. (Suite) Closure and allocation. For each T the map step, tangent step and log accumulation return T (no promotion or demotion); allocated bytes in the hot loop after warm-up are 0 for Float32, Float64, Double64 (BigFloat exempt).
13. (Suite) Conversion consequence is visible (Equations 5). K = 2 pi, 100 seeded ICs, Double64 orbit with K built in Double64 versus K widened from Float64: the first step at which max-norm separation exceeds 10^-8 is <= 60 for at least 80 percent of the ICs (derived estimate 15 to 30 steps; the percentage allows accelerator-mode islands). The measured median n* is reported.

## Sources

- `_config/catalogs/lyapunov.md`, section E; `_config/architecture.md`, "Decision" and "Compute"; `_config/environment.md`, "High precision"; `_config/decision-log.md`, D17, D20; `_config/risks.md`, K15, K18.
- `_config/research-question.md`, "Ground truth for variational methods (standard map)": map and tangent map (Equation 1).
- `schema/run-store.md`, root attributes `precision`, `backend`, `device`, `backend_fallback_reason`; `_config/coding-conventions.md`, "Tests".
- `_config/bibliography/benchmarks.md`: row "Lyapunov exponent, single chaotic orbit" (0.298 at K = pi/2; 1.172 at K = 2 pi), SalaArtuso2015, p. 9 as stated in the row, Eq. 5 (McMillan form, converted in Equations 2); row "Standard map largest Lyapunov exponent, large K", Shevchenko2004. Both used for horizon sizing only (K18).
- Derived (no outside source): Equations 2, 3, 4, 5 and the Float32 lattice heuristic; the pi rounding errors were computed from a 60-digit value of pi.
