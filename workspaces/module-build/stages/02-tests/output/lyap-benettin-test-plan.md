# Test plan: lyap-benettin

Spec: `workspaces/module-build/stages/01-spec/output/lyap-benettin-spec.md` (D30) with its Amendment D35 (the bound B_S of S3 with the n^2 term, the state-magnitude factor of AC7 and AC11, the amplification sum Psi_n; the stage 02 relaxations of this plan are withdrawn, see the table `scratchpad/m2/batch-b-bound-table.md`). API contract: M2 Batch B `api-batch-b.md`, final (D34); names tested: `Benettin`, `safe_norm`, `effective_eps`, `reducer_provenance`, `NumericalHealthError.column`, and the Batch A names `propagate`, `propagate_orbit`, `ResidualMonitor`, `initial_deviations`, `StandardMap`, `WrappedStandardMap`, `LinearSymplecticMap`. Tests: `test/lyap-benettin/runtests.jl` including `lb_helpers.jl` (which includes `lb_psi.jl`, the one definition of Psi_n, also included by the householder and two-trajectory suites), `lb_ac01_04.jl`, `lb_ac05_07.jl`, `lb_ac08_11.jl`, `lb_ac17_19.jl` (helpers and test maps prefixed `lb_`, nothing in `src/`).

Conventions. Element types Float64 and Double64 (Float32 dropped, D24). Every tolerance is in units of `effective_eps(T)` (Float64 eps, Double64 2^-100); B_S as in S3 as amended by Amendment D35 (`lb_BS(eps, n, G_first, K)` = 8 eps (n + G_first) + 16 abs(K) n^2 eps, written once; K = 0 where the tangent is constant, i.e. AC1 to AC4). Seeds: base 20261012 (`numerics-defaults.toml`) plus offsets 21000 to 21999 (`LB_OFF`): 21010 AC1, 21020 AC2, 21030 AC3, 21040 AC4, 21050 + case AC5, 21060 + case AC6, 21070 S7 regular-class screening (shared by AC7, AC11 and AC17), 21080 AC8, 21090 AC9, 21100 AC10, 21110 AC11, 21190 AC19, 21200 API, 21210 prefix invariance. BigFloat references at 256 bits inside `setprecision(BigFloat, 256) do`. Orbit ids are 1..N; the deviation vector of orbit j is `initial_deviations(seed, j, Val(2), Val(1), T)`.

Test-side tools (all in `lb_helpers.jl`): `lb_Probe` is a wrapper reducer that forwards `reducer_init/step/failed/error/final` to `Benettin` and records the Benettin final at chosen handoffs (the stream has no emission channel, Q9) and the norm alpha_j of the incoming deviation vector; `lb_run` runs `propagate_orbit` over a vector of orbits (same body as `propagate`, no kernel compile; the kernel path is exercised in AC9, AC10 and the provenance test); `lb_Fault` (zero or NaN tangent for flagged orbits), `lb_Boom`, `lb_Trace` (records the orbit at every handoff).

## Criterion to testset

| AC | Testset | Sample | Tolerance |
|----|---------|--------|-----------|
| 1 | "AC1 Y5 hyperbolic map" | A = [[2,1],[1,1]] via `LinearSymplecticMap`, start at the origin, ids 1..100 (Float64) and 1..50 (Double64), tau in {1,5,10}, n in {20,100,1e4}; ids with abs(c_1) < 0.05 excluded (at most 10 percent of the ids) | (a) vs ln(norm(A^n w0)/norm(w0)) from BigFloat A^n; (b) vs the closed form of E10 (last term kept); both <= B_S(n), G_first = squared Frobenius norm of A^tau; also `n` of the final = N |
| 2 | "AC2 rotation and elliptic point" | (a) rot4, w0 = (1,0), tau in {1,3}, n in {12,120,1200,9996}; (b) K = -2 at (0,0), tau 1 with n in {1,2,3,4,10,11,100,101,1000,1001}, tau 2 with n in {2,4,10,100,1000}, ids 1..100 (Float64) and 1..50 (Double64) | (a) sum = comp = lambda = 0 exactly; (b) vs BigFloat M^n w0; even n |nX| <= B_S; odd n |nX| <= ln((3+sqrt 5)/2) + B_S; G_first 7 (tau 1), 2 (tau 2) |
| 3 | "AC3 K = 0 shear" | seeded dyadic starts, tau in {1,10}, n in {10,100,1e4} with 20 ids (Double64 at n = 1e4: 8) and n = 1e6 with 3 ids (Float64 only); plus w0 = (0,1), n = 1e4 | vs ln(norm(w_n)/norm(w_0)), w_n = (a + n b, b); <= B_S, G_first = 2 + tau^2; nX < ln(n+1); (0,1) case <= 8 n eps (tangent-core AC3) |
| 4 | "AC4 hyperbolic fixed point" | (0,0), K in {1,2,6}, tau in {1,5}, ids 1..100, n in {20,100,1000}, both types | vs the exact closed form w_n = c_1 mu^n v_u + c_2 mu^-n v_s (E12 with the O(mu^-2n) cross term kept, so no exclusion of small c_1 is needed); <= B_S |
| 5 | "AC5 same-orbit tangent propagation" | K = 0.5 (n = 1e4), 2 pi (n = 1e4), 10 (n = 1e3); tau in {1,4}; 20 starts (Double64: 6 at n = 1e4, 10 at n = 1e3, runtime cut); checks at every handoff multiple of 100 | E1 at 256 bits at the realized points for every n (S3 as amended); <= B_S; orbit equals the forward orbit bitwise; G_first > 1e4 counted with `@info` |
| 6 | "AC6 tau invariance" | tau in {2,5,10} against tau = 1, K in {0.5, 2 pi}, n in {100,1000}, 20 starts, both types | <= B_S^tau(n) + B_S^1(n) (S3 as amended) |
| 7 | "AC7 full BigFloat orbit" | K = 0.5, S7 screening (X(1e5) <= 5e-4, replaced starts counted with `@info`), 50 starts, n in {10,100}, both types | <= 32 K M_state n^2 eps rho_n (Amendment D35; M_state of S1, `lb_ref_full`); rho_n > 50 reported and at most 5 excluded |
| 8 | "AC8 norm evaluation" | (3c,4c) for c = 2^600, 2^-600, 1, 2^1000, 2^-1000; 200 seeded random vectors per type (d = 2 and 4, scale 2^e, abs(e) <= 900); zero; NaN and Inf | relative error <= 2 eps against BigFloat; zero gives 0; non-finite input gives a non-finite norm |
| 9 | "AC9 numerical health" | (a) zero tangent for flagged orbits: orbit 4 of 6, chunk 2, tau 3 (`propagate`, Float64) and the per-orbit body (Double64); lowest id first with two flagged; (b) NaN tangent: no masking by the chain; `Benettin` with k = 2 | what "deviation matrix", orbit, step = tau, backend :cpu, `column` 1 and value 0.0 for the zero case (threshold not asserted, Q11), `column` 0, threshold Inf for NaN; nothing stored; ArgumentError message "Benettin needs k = 1, got k = 2" before any work |
| 10 | "AC10 orbit and execution invariance" | (a) 5 starts, tau in {1,4}, chains {Trace}, {Benettin, Trace}, {Benettin, monitor, Trace}; (b) 4096 orbits, K = 2 pi, n = 1e3, tau in {1,4}, `propagate` with chunk 64 and 1000 against the serial map of `propagate_orbit` (Float64); (c) `@allocated`, n = 1e4, chain {Benettin, monitor}, tau in {1,4}, both types | (a) orbit points bitwise equal; (b) x, W, finals bitwise (`isequal`), `provenance.estimator == "benettin"`; (c) 0 bytes |
| 11 | "AC11 wrapped versus lifted state" | (a) K = 0, 50 dyadic starts, tau in {1,3}, n = 300; (b) K = 0.5 regular class, 50 starts, n in {10,100}; Float64 | (a) finals bitwise equal; (b) <= 64 K M_state n^2 eps rho_n (Amendment D35) |
| 12 | stage 04 (V) | precision agreement of the series (Float64 vs Double64 vs BigFloat, 1000 starts) | ensemble statistics, `verify/` |
| 13 to 15 | stage 04 (V) | external tangent-QR oracle (ORACLE-CT1, ORACLE-CT2, strong chaos) | needs the external package |
| 16 | stage 04 (G) | CPU vs GPU, tolerances tangent-core AC12 | needs the GPU |
| 17 | "AC17 regular-orbit decay" | 20 regular-class starts (S7), K = 0.5, tau = 1, n up to 1e6, 21 log-spaced n in [1e4, 1e6] | least-squares slope of ln X against ln n in [-1.1, -0.8] for at least 16 of 20 starts |
| 18 | "AC18 range invariant" | alpha_j of every run of AC5 and AC6 (`lb_Probe` records the norm of the incoming vector) | alpha_j within sigma(K)^+-tau (1 +- 16 eps tau Phi_j), Phi_j = squared Frobenius norm of the window product (Float64, from the stored tangents) |
| 19 | "AC19 symplectic dimension 2" | (a) chain {Benettin, ResidualMonitor(sys, 4, Inf)} on the AC5 cases, tau in {1,4} (Float64, 20 starts), tau = 1 (Double64, 8 starts); (b) LinearSymplecticMap [[2 + delta, 1], [1, 1]] (det 1 + delta), both sides of the det threshold, both types, per-orbit body | (a) no health error, det_max <= 4 eps (2 + K) and > 0; (b) delta in [2 theta, 4 theta): "det" error, value delta, orbit 1, step 1; delta in (theta/4, theta/2]: no error over 100 steps |
| API | "API reducer protocol and provenance" | `reducer_init/step/update/final/failed/error` on one handoff and on a failing handoff, both types; provenance of a `propagate` run | n = 0 final (zeros, `SVector{1,T}` fields); step arithmetic against ln 5, ln 5 + ln 10; failure recorded without throwing in `reducer_step`, W and accumulators unchanged, thrown by `reducer_update` (orbit 0, :cpu); `reducer_error` (column 1); `reducer_provenance(Benettin()) == (estimator = "benettin",)` merged into `propagate` provenance; `lambda == (sum + comp)/n` |
| API | "API prefix invariance" | all-handoff snapshots of a run of 400 steps against separate runs of length j tau, tau in {1,4}, both types | bitwise equal sum and comp (Q9) |

Testset count: 16 (AC1 to AC11, AC17 to AC19, two API). AC12 to AC16 are stage 04. 15 009 assertions when green against the throwaway stub (below).

## Deviations from the spec

1. AC5 (withdrawn). Stage 02 had replaced the reference for n > 1000 by the tangent matrices as stored in T and moved K = 0.5 at n = 1e4 to stage 04 because B_S failed on correct code. Amendment D35 adds the n^2 term to B_S (S3, the single home); the test uses the E1 reference at 256 bits for every n and B_S as amended, and K = 0.5 at n = 1e4 is in the package suite. No bound of this plan differs from the spec as amended.
2. AC7 and AC11 (withdrawn). The test factor M is the spec's M_state (Amendment D35, S1 symbol), not an own form: `lb_ref_full` returns max(1, |K|, largest lifted coordinate).
3. Package-suite sample reductions (runtime, owner question on the Double64 coverage): AC3 at n = 1e6 only for Float64 (Double64 `sin`/`cos` at |x| ~ 1e6 cost 10 us per step: 30 s for three orbits) and 3 ids; AC1 and AC2 Double64: ids 1..50; AC5 Double64: 10 starts at n = 1e3 and 6 at n = 1e4; AC8: 200 random vectors per type; AC10(b) Float64 only; AC11 Float64 only; AC19(a) Double64 at tau = 1 with 8 starts; the S7 screening is done once (first 50 regular starts, shared by AC7, AC11 and AC17, the latter using the first 20). The full-size versions go to the stage 04 list (`verify/`): AC3 Double64 n = 1e6, AC5 Double64 20 starts at every n, AC10(b) Double64 thread invariance.

## Spec ambiguities resolved

- AC1(b) and AC4: the closed forms keep the exponentially small cross and last terms, so they are exact (not asymptotes); small c_1 then needs no exclusion in AC4; AC1 keeps the spec exclusion of abs(c_1) < 0.05 (counted, at most 10).
- Reference normalisation. n X(n) is compared with ln(norm(Y_n w0)/norm(w0)) where w0 is the stream's seeded vector widened exactly (its norm is 1 up to eps(T)); the ratio form removes that rounding from the comparison.
- AC9(b) non-finite entries are tangent-core AC8; only the check that the Benettin chain does not change that verdict (what "deviation matrix", `column` 0) is added. The threshold of the zero-norm case is not asserted (Q11).
- AC10(b) "1 vs N threads": the test never changes thread settings (Batch A rule); it compares `propagate` with the serial `map` over `propagate_orbit`, which is thread-count independent by construction; the run's `Threads.nthreads()` is whatever the runner uses.
- AC18: alpha_j is the norm of the deviation vector the estimator receives (read by `lb_Probe` before the inner step), Phi_j is computed in Float64 from the stored tangents; sigma(K) is evaluated in T.
- AC17: n sampled at 21 log-spaced values (the series has no sink, Q9; `lb_Probe` snapshots Benettin's final at those handoffs of one run). A non-positive X(n) counts as a start outside the interval.
- AC19(a): the monitor is chained with orth_tol = Inf (k = 1 has no orthogonality check anyway); the det bound uses max|DT| <= 1 + K.
- Series: prefix invariance (Q9) is tested directly ("API prefix invariance") because every series test relies on it being the same as separate runs.
- Provenance: `estimator = "benettin"` is asserted in the `propagate` runs of AC10(b) and the API test; the schema v2 keys beyond that are not Benettin's.

## RED log

Re-run after the Amendment D35 edits (new shared include `lb_psi.jl` included): 0 passed, 0 failed, 16 errored in 16 testsets (51 s), the same table; nothing comes from a helper or an include.

Command (repo root, test environment with Chaoslab dev'd at the unmodified repo, 4 threads, whole suite inside one testset): `julia -t 4 --project=<scratchpad>/m2/testenv -e 'using Test, Chaoslab; @testset verbose=true "lyap-benettin" begin include("test/lyap-benettin/runtests.jl") end'`. Log: `scratchpad/m2/red-lyap-benettin.log`. Result: 0 passed, 0 failed, 16 errored, 16 testsets; no parse error and no error inside a helper. Every testset errors on its first missing Chaoslab name: `UndefVarError: effective_eps` (AC1 to AC6, AC8, AC18, API protocol; the Foundation amendment A1b) or `UndefVarError: Benettin` (AC7, AC9, AC10, AC11, AC17, AC19, API prefix invariance). Both are names the contract adds, so the verdict is the same for each; after `effective_eps` exists the next missing name is `Benettin` in every one of them.

| Testset | Error | Fail | First error |
|---|---|---|---|
| AC1 Y5 hyperbolic map | 1 | 0 | UndefVarError `effective_eps` |
| AC2 rotation and elliptic point | 1 | 0 | UndefVarError `effective_eps` |
| AC3 K = 0 shear | 1 | 0 | UndefVarError `effective_eps` |
| AC4 hyperbolic fixed point | 1 | 0 | UndefVarError `effective_eps` |
| AC5 same-orbit tangent propagation | 1 | 0 | UndefVarError `effective_eps` |
| AC6 tau invariance | 1 | 0 | UndefVarError `effective_eps` |
| AC7 full BigFloat orbit | 1 | 0 | UndefVarError `Benettin` |
| AC8 norm evaluation | 1 | 0 | UndefVarError `effective_eps` |
| AC9 numerical health | 1 | 0 | UndefVarError `Benettin` |
| AC10 orbit and execution invariance | 1 | 0 | UndefVarError `Benettin` |
| AC11 wrapped versus lifted state | 1 | 0 | UndefVarError `Benettin` |
| AC17 regular-orbit decay | 1 | 0 | UndefVarError `Benettin` |
| AC18 range invariant | 1 | 0 | UndefVarError `effective_eps` |
| AC19 symplectic dimension 2 | 1 | 0 | UndefVarError `Benettin` |
| API reducer protocol and provenance | 1 | 0 | UndefVarError `effective_eps` |
| API prefix invariance | 1 | 0 | UndefVarError `Benettin` |

## Helper check and runtime

The helpers and every closed form were checked against a throwaway implementation (scratchpad `m2/lb_patch`, never committed): `safe_norm` by max-abs scaling with plain division, `Benettin` exactly as the contract describes it (compensated sum of ln alpha, replacement by the unit vector, code 1 on a zero or non-finite alpha, `reducer_step` non-throwing, `column` 1, provenance), plus the A1 to A3 amendments (Batch A, D34). Result after Amendment D35: all 15 009 assertions pass in 2 min 36 s wall (`scratchpad/m2/green-lyap-benettin-d35.log`; cold start, including all JIT compilation, 4 threads, a machine shared with other test runs). Before the amendment (2 min 7 s wall, 133 s of process CPU) about a third was compile time (a second pass in the same process needed about 30 percent less CPU), mostly one kernel and per-type `propagate_orbit` compilation for each argument-type combination, which is why most testsets use the serial per-orbit body (`lb_run`) and `propagate` is used only in AC9, AC10 and the provenance test. Realized quantities printed with `@info` (stage 04 should record them): worst ratio error/B_S (B_S as amended) in AC5 0.0015, AC6 0.0034 (with the old linear B_S the same runs reached 5.5, see Amendment D35); AC5 takes 31 s of the suite (E1 at 256 bits for every n); 48 of the 164 AC5 runs have G_first > 1e4; the S7 screening replaced 1 candidate (first 50 regular starts).

## Spec items not turned into a package test

- AC12 to AC16: verification environment and local GPU (stage 04), as the spec tags them.
- Items moved to stage 04 for runtime (deviation 3 above).
- The "ensemble statistics" of S5 and the series storage schedule of S8 (no sink in the contract, Q9): only prefix invariance is tested.
