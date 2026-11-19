# Test plan: clv-ginelli

Status: approved at the stage 02 checkpoint, 2026-11-19 (owner: "Approve tests"), with the test-side bounds U3 to U8. RED run of the new suite: 12 pass, 29 error, 41 total; all 29 errors are `UndefVarError` (`ClvGinelli` 28, `clv_initial_coefficients` 1). The 12 passes are the three negative-control testsets, which use no new name. Run against a throwaway scratch implementation outside the repo (session scratch folder `clv-red`, partial: no checkpointing, sink, health errors or ensemble loop): every test of AC1 to AC7, AC9, AC11, AC12 passes except the points that are scratch limits (listed under Calibration); total scratch run about 32 s.

Spec: `workspaces/module-build/stages/01-spec/output/clv-ginelli-spec.md` (approved at the stage 01 checkpoint, D122, open points P1 to P16 as proposed; 14 acceptance criteria, recorded lines R1 to R6). Suite: `test/clv-ginelli/` (`runtests.jl`, `clv_helpers.jl`, `clv_ac01_05.jl`, `clv_ac06_10.jl`, `clv_ac11_14.jl`, 946 lines, 32 testsets).

Conventions copied from continuous-qr: every testset begins by naming the new symbols it uses (`_ = (ClvGinelli, ...)`), so a missing module errors that testset only; no top-level call and no helper calls new API at include time; helpers carry the prefix `clv_`; package sizes sit behind `clv_size(small, full)` with `CLV_FULL_RUN` (default false, set by a later `verify/clv-ginelli.jl`); test-side doubles act as negative controls. The independent reference algorithm of the suite (`clv_ref`, `clv_ref_run`) uses `LinearAlgebra.qr` and triangular solves only, in any element type, on given tangent products; at stage 03 it is the oracle for the realized orbit points. Nothing in `src/` is edited at stage 02.

## API proposal for stage 03

Layer `Chaoslab.Lyapunov`, file `src/Lyapunov/clv_ginelli.jl`. The module owns a per-orbit loop (not a reducer of the tangent stream, P10): the orbit is advanced by `forward(map, z)` and the tangent of one step is `tangent(map, z)`, `tau` steps per handoff, tangent product with later steps on the left. CPU only (P9). All arrays are `Array{T}` of the element type `T` of the start state; handoff index `n` is stored at index `n + 1`.

| Name | Meaning |
|------|---------|
| `ClvGinelli(map; tau = 1, h = 1.0, orthonormalizer = HouseholderQR(), storage = :store, block = 0)` | `tau >= 1` steps per handoff; `h` the step length (time per map step, 1 for maps), used only to convert exponents to per unit time (Dt = tau h); `orthonormalizer` any delivered `Orthonormalizer` (P8); `storage` in `:store` or `:checkpoint`, `block` = h_b >= 1 for `:checkpoint`. `tau < 1`, an unknown `storage`, `:checkpoint` with `block < 1` throw `ArgumentError` |
| `propagate_clv_orbit(clv, x0, n_f, n_s, n_b; k, seed, orbit = 1, W0 = nothing, C0 = nothing, cseed = seed, outputs = (), angle_sets = nothing, sink = nothing, backend = :cpu)` | one orbit. `x0::SVector`; `k` in `1:d`; `W0` the starting d x k frame (default `initial_deviations(seed, orbit, Val(d), Val(k), T)`); `C0` the starting k x k upper-triangular matrix (default `clv_initial_coefficients(cseed, orbit, Val(k), T)`). `outputs` a tuple from `:V` (phase-space CLVs), `:G` (frames), `:R` (factors); `angle_sets = (r1, r2)` two column ranges for the principal angles (P11). `ArgumentError` before any work for `k` outside `1:d`, `n_s < 1`, `n_f < 0`, `n_b < 0`; `NumericalHealthError` for health failures |
| result (`ClvResult`) | fields `C` (k, k, n_s + 1) coefficient matrices C_n, n = 0 .. n_s, unit columns, positive diagonal; `lambda` (k, n_s + 1) local exponents ln Lambda_i(n) / (tau h), last column NaN when `n_b = 0`; `cosines` (k, k, n_s + 1) the matrix \|C_n^T C_n\| (V12 a); `V`, `G` (d, k, n_s + 1) and `R` (k, k, n_s + 1, last slice NaN when `n_b = 0`) when requested, else `nothing`; `angles` (min(m1, m2), n_s + 1) all principal angles between the two sets, ascending, when `angle_sets` is given, else `nothing`; `lambda_fwd` (k) forward exponents of the window n = 0 .. n_s - 1, per unit time; `x0` the state at the window start (after n_f handoffs), `x` the state after all n_f + n_s + n_b handoffs; `bytes` the bytes the run holds (stored R, frames, buffers; peak when `sink` is given); `record` NamedTuple of the run parameters (see AC14e) |
| `sink(n, C_n, lambda_n)` | optional callback, called once per sample for n = n_s down to 0 (decreasing, as the backward recording runs); with a sink `C`, `lambda`, `cosines`, `V`, `G`, `R`, `angles` of the result are `nothing` and only `lambda_fwd`, `x0`, `x`, `bytes`, `record` are returned |
| `propagate_clv(clv, xs, n_f, n_s, n_b; k, seed, chunk = nothing, ntasks = nothing, backend = :cpu, outputs = (), angle_sets = nothing)` | many orbits, orbit id j = index in `xs`; returns an object with `results` (vector of `ClvResult`, member j equal to `propagate_clv_orbit(...; orbit = j)` bitwise) and `provenance` (NamedTuple with `backend`, `backend_fallback_reason`, `chunk`, `ntasks`). A `:gpu` or `:auto` request runs on the CPU with `backend = :cpu` and the reason "reducer is CPU-only" (`FALLBACK_REASONS[6]`) |
| `clv_initial_coefficients(seed, orbit, Val(k), T)` | the starting C: upper triangular, positive diagonal, unit columns (to 8 eps), from a stream keyed by (seed, orbit) and different from the frame stream; column j independent of k (P13) |
| `clv_bytes(Val(d), Val(k), T; n_s, n_b, storage = :store, block = 0, frames = true)` | V17 byte count: `bytes_store = sizeof(T) (N k(k+1)/2 + (n_s + 1) d k) + bytes_work`, `bytes_work = sizeof(T) (d + 2 d k + k(k+1)/2 + 2 k) + 16`, N = n_s + n_b; checkpointed: `sizeof(T) (n_c (d + d k) + h_b (k(k+1)/2 + d k) + k(k+1)/2) + bytes_work`, `n_c = ceil(N / h_b)`; `frames = false` drops the frame term of the store-all count. 288576 and 24528 for the spec examples |

Health errors (`NumericalHealthError(what, value, threshold, step, orbit, backend, column)`): `what` = "state" (non-finite orbit state or start), "deviation matrix" (non-finite starting frame), "orthonormalisation breakdown" (forward factorization: zero or non-finite pivot, passes, conditioning; `column` the failing column; the delivered kind), "backward breakdown" (zero or non-finite column norm after the back-substitution, P14). `step` counts handoffs from 1 over the forward transient and the recorded window in forward order; the backward breakdown reports the handoff `n` at which the solve broke, counted in the same global numbering. `orbit` the orbit id. Nothing is stored or delivered (the sink is never called). Exports: `ClvGinelli`, `ClvResult`, `propagate_clv_orbit`, `propagate_clv`, `clv_initial_coefficients`, `clv_bytes`.

## Criterion to testset

| AC | Testset (file) | Content |
|----|----------------|---------|
| 1 | AC1, AC1 negative controls (`clv_ac01_05.jl`) | fixtures a to d, 4 frames, tau in {1, 4} (flow tau = 20), n_f = n_b = 100, n_s = 50: delta of every CLV, exponents, column 1 = e_1 and v_1 = g_1 bitwise, positive diagonal; Float64 and Double64 (fixture b). Controls (pass now): Gram-Schmidt vectors, forward step, lower triangle |
| 2 | AC2, AC2 control (`clv_ac01_05.jl`) | I1 (orthonormal G, triangular R and C, unit columns), (a) covariance, (b) product of norms, (c) telescoping and the window mean against the delivered `QRSpectrum` on the identical orbit; fixtures a, b, standard map starts, F3 orbit C and chaotic starts. Control (pass now): 1-norm normalisation |
| 3 | AC3a, AC3b, AC3c (`clv_ac01_05.jl`) | (a) slope of ln error against n_b, fixtures a, b, j >= 2 within 20 percent, column 1 at the floor; (b) 2 frames x 2 C seeds x 2 taus at common samples, standard map, coupled maps, F3 chaotic starts (pair 2, 3 as a span); (c) covariance chain against the 256-bit tangent product, m in {1, 10, 50}, Float64 and Double64 |
| 4 | AC4 (`clv_ac01_05.jl`) | order 2 and 4 tables, h in {1/10, 1/20, 1/40}: ratios in [2^p / 2, 2 x 2^p] above the floor; errors reported |
| 5 | AC5, AC5 control (`clv_ac01_05.jl`) | (a) pair identity with ln\|omega\| for a random C and n_b = 1 and after the transients, (b) non-partner omega <= 1e-9, (c) determinant identity; standard map, coupled maps, fixture b, F3 chaotic starts. Control (pass now): wrong Lambda normalisation |
| 6 | AC6a, AC6b, AC6c (`clv_ac06_10.jl`) | (a) cat map pi/2, fixture b arctan(1/3) and pi/4; (b) cosines from C against V = G C within 16 eps; (c) planted tangency delta in {1e-2, 1e-4, 1e-6}, Float64 and Double64, smallest principal angle against the pairwise angle |
| 7 | AC7 (`clv_ac06_10.jl`) | (a) f(z) in span(v_2, v_3), (b) grad H orthogonal to v_1 and v_4, both <= 1e-6; (c) forward exponents pairing for 95 percent of chaotic starts |
| 8 | AC8, AC8b, AC8c, AC8d (`clv_ac06_10.jl`) | (a) bitwise store-all against checkpointed for h_b in {1, 7, 64, 500} (and F3, h_b = 7); (b) counted forward and tangent calls; (c) bytes against `clv_bytes`, 288576 and 24528; (d) streaming sink, order, nothing stored, peak bytes |
| 9 | AC9, AC9c (`clv_ac06_10.jl`) | coupled maps and F3 orbit C, k in {1, 2, 4}, Householder and repeated Gram-Schmidt: bitwise first k columns; k = 0 and k = 5 rejected; starting C properties |
| 10 | AC10a, AC10b, AC10c (`clv_ac06_10.jl`) | (a) 1 against 4 tasks, chunk 64 and 1000, ensemble position, orbit bitwise the integrator and the delivered spectrum; (b) forward exponents bitwise the delivered spectrum (standard map, F3); (c) allocation by differences of two runs, Float64 and Double64 |
| 11 | AC11, AC11b (`clv_ac11_14.jl`) | (a) Float64 and Double64 against the 256-bit reference at the realized orbit points, B_v; (b) BigFloat module run, two starts agree to 1e-10 |
| 12 | AC12a, AC12b (`clv_ac11_14.jl`) | (a) linear fixtures a, b, c with their inverse (index order reversed, with a control); (b) standard map K = 5 in BigFloat(256), reversed orbit from the end point, 1e-12 |
| 13 | none | the static null-space algorithm in 256 bits is the verification environment (stage 04, `verify/clv-ginelli.jl`); the package suite carries the independent reference `clv_ref` as a stand-in for the realized-orbit checks (AC3c, AC11) |
| 14 | AC14a to AC14e (`clv_ac11_14.jl`) | (a) non-finite start, frame, injected NaN state; (b) singular tangent; (c) overflow of the backward solve; (d) argument errors with a call counter proving no work was done; (e) run record keys and the `:gpu` fallback |
| R1 to R6 | none | recorded lines, stage 04; R3 ratios are printed by `@info` lines of AC1, AC2, AC3c, AC5, AC6b, AC11 |

## Tests that pass already

Three testsets use no new name and pass in the RED run (12 tests): "AC1 negative controls" (8), "AC2 negative control" (2), "AC5 negative control" (2). They test the reference algorithm and the doubles, so a wrong helper fails now rather than at stage 04.

## Points found while writing the tests (test-side readings, U)

- U1: the identity frame is not generic for fixture (b): v_4 = S e_3 = (1, 2, 1, 0) lies in span(e_1, e_2, e_3), so the forward 3-frame never leaves it (realized delta 0.745). Fixture (b) uses a fourth seeded frame in place of the identity; the other fixtures keep the identity.
- U2: layout. `C`, `V`, `G`, `R`, `cosines` have n_s + 1 samples, n = 0 .. n_s, stored ascending; `lambda` has n_s + 1 columns (Lambda_i(n) needs R_n); with `n_b = 0` the last column of `lambda` and `R` is NaN. The sink delivers decreasing n (the spec's output order). The spec lists the angles as the smallest principal angle only; the proposal returns all principal angles ascending (AC6a needs both for fixture b).
- U3: n* = ceil(40 / (g_min Dt)) instead of 25 (chaotic prefactor of the convergence: measured 1e5, so 25 gives about 2e-6 on F3, 40 gives below 1e-10 at the starts used); g_min comes from a test-side pilot spectrum (reference algorithm, 300 handoffs), the pair (2, 3) of F3 is skipped; n* is capped at 3000 (package) or 6000 (full). F3 orbit C (sticky, lambda_1 = 0.0127, gap lambda_1 - lambda_2 = 0.007, n* above 3000) cannot meet 1e-10 or 1e-6 at package sizes: criteria that need convergence on F3 (AC3b, AC5b, AC7) use chaotic starts selected by the SALI rule at T = 1000 (henon-heiles criterion 10) from the recorded stream `clv_seed(0, 1)`; orbit C stays in the identities that hold for any start (AC2, AC5a, AC6b, AC8, AC9, AC11).
- U4: B_fix = 1e3 eps + 100 exp(-gamma n) (spec 4 exp(-gamma n)): in Double64 the exponential term is the bound and the realized ratio was 4.9 with prefactor 4, 0.2 with 100 (fixture b, gamma n = 44.6).
- U5: the exponent bound of AC1 adds the rounding of the interval product, 4 tau eps exp((lambda_1 - lambda_d) Dt) / Dt: fixture a at tau = 4 has the interval condition exp(16) and a realized exponent error of 1.8e-13, 3.3 times the bare B_fix / Dt.
- U6: AC3c: ln \|\|M_{m,n} v_i\|\|\| is sensitive to the rounding of v_i along faster directions (growth exp((lambda_1 - lambda_i) m Dt)); the standard map at i = 2 and m = 50 errs by 20 with the bare bound B_cov m. The test adds 256 eps kappa exp((lambda_1 - lambda_i) m Dt), with lambda from the forward exponents; at i = 1 and m = 1 it is 5 B_cov. Realized ratio 0.24 (F3 Float64 alone exceeds the bare bound by 1.27 at m = 1).
- U7: AC5 (a) uses 512 eps (1 + ...) instead of 64 eps: realized 0.95 of 64 eps (omega is formed from the delivered V, absolute error eps, relative eps / sin theta). Realized 0.12 of 512 eps.
- U8: B_v = c_v eps (1 + G_max) with c_v = 1e4 (spec 1e3): realized ratio 1.6 at c_v = 1e3 on F3 orbit C, 0.16 at 1e4 (Float64); the others at or below 2.4e-4.
- U9: `lambda_fwd` is defined as (sum of ln R_ii over the window / (n_s tau)) / h with the compensated sum of the delivered spectrum, so that it equals `QRSpectrum(alg).lambda ./ h` bitwise (the delivered lambda is per step, not per handoff: the first draft divided by tau h and was wrong by tau for F3). AC10b tests this on the standard map and on F3 for Householder and repeated Gram-Schmidt, with n_f = 0 and the frame handed in.
- U10: AC7 (c) at package size: 10 chaotic starts at n_s = 300 time units (10 of 10 pass; 4 of 10 at 150); the spec's 200 starts at 2000 are the full size.
- U11: counting and allocation. AC8b counts calls of `forward` and of `tangent` through a wrapper, one call of each per map step (tau = 1 map, the coupled maps); a module that fuses the two must keep the counts as stated. AC10c measures allocation as the difference of two runs (200 more forward handoffs: equal bytes; 200 more backward-transient handoffs in checkpointed mode: at most one extra checkpoint), not a counter inside the module.
- U12: AC14. `step` convention: global handoff count from 1; the NaN map returns NaN from its 5th forward call, expected `step == 5` for `n_f = 3`. AC14c relies on a tangent diag(1, 1e-310) (subnormal pivot, back-substitution overflows) so the forward QR must accept the subnormal pivot; if stage 03 finds the forward conditioning rule fires first, replace the fixture, not the kind.
- U13: criterion I1 states v_1 = g_1 exactly, i.e. column 1 of C is e_1 as stored values; in Double64 the plain normalisation x / \|x\| of a one-entry column is not exactly one (0.99999...923), so stage 03 must set column 1 explicitly. The tests assert exact equality.
- U14: AC6a and AC6c in Double64 need the angles computed in the element type (an `acos` through Float64 misses the bound by 1e11).
- U15: criterion 9 is encoded as bitwise (the spec's first reading); the 64 eps fallback is a stage 03 decision on realized differences.

## Not encoded

- Criterion 13 and R1 to R6 (stage 04).
- Criterion 14 (d) "tau h incompatible with the integrator": in this API `h` is only the unit conversion, the integrator step is the map's own, so there is nothing to check; no test.
- The sign-skipped control of criterion 2 (a double with negative c_ii): its observable is the sign of Lambda, which the module's formula Lambda = 1 / \|\|R^-1 c\|\| makes positive by construction; the 1-norm control and the wrong-Lambda control (AC5) are kept.
- Float32 output (D24, spec preamble): no criterion tests it.
- The angles-only checkpointed byte count: the spec gives no formula for it; only the store-all `frames = false` count is tested (AC8c).
- Seeded fault injection for `n_b = 0` Lambda: the NaN last column is a convention (U2), tested only by not reading it.

## Calibration (scratch implementation outside the repo, session scratch folder `clv-red`)

The scratch implemented the reference-style algorithm (about 55 lines, Float64 `acos` for the angles, plain sums) and ran the suite. Realized values at the package sizes:

| Test | Bound | Realized |
|------|-------|----------|
| AC1 worst delta / B_fix, exponent / bound | 1 | 0.195 (Double64), 0.039 |
| AC2 (a), (b), (c) / bound | 1 | 0.013, 0.020, 0.020 |
| AC2 window mean against the delivered spectrum / allowance | 1 | 0.989 (identity-like, boundary term) |
| AC3a worst relative slope deviation | 0.2 | 0.0064 |
| AC3b worst distance at n* = 40 / (g Dt) | 1e-10 | 1.6e-12 (n* = 5, 26, 1420, 1066) |
| AC3c ratio error / bound | 1 | 0.24 |
| AC4 errors h = 1/10, 1/20, 1/40 | ratios 2^p within [1/2, 2] | order 2: 3.9e-4, 9.8e-5, 2.4e-5; order 4 pass |
| AC5 (a) / 512 eps, (b) / 1e-9, (c) / 64 eps | 1 | 0.12, 9e-4, 0.12 |
| AC6b / 16 eps | 1 | 0.625 |
| AC7 worst sine and cosine | 1e-6 | 1.0e-13 |
| AC11 error / B_v at c_v = 1e4 | 1 | 0.16 (F3 C, Float64); others below 3e-4 |
| AC12a, AC12b | B_fix, 1e-12 | pass |

Scratch limits (these failed or errored only because the scratch lacks the feature): the exact e_1 column and the Double64 angles in Double64 (U13, U14), bitwise `lambda_fwd` (plain sums), call counts of the checkpointed mode, the allocation tests, the sink, the health errors, the argument errors, the run record. Cost: the scratch ran the whole suite in about 32 s; the folder ceiling is 180 s (D38).

## Seeds

StableRNGs. `clv_seed(ac, sub) = 20261012 + 45000 + 50 ac + sub`; frames through `initial_deviations(clv_seed(ac, sub), 1, Val(d), Val(k), T)`, starting C through `clv_rand_C(StableRNG(clv_seed(ac, sub')), k, T)` (test side) or `clv_initial_coefficients`; random starts through `clv_starts(StableRNG(clv_seed(ac, sub)), d, m, T)` (uniform in [0, 2 pi)^d); F3 chaotic starts from `random_energy_start` with the stream `clv_seed(0, 1)` (AC2 to AC7) and `clv_seed(0, 2)` (AC7 c); orbit C at (q2, p2) = (-0.015, 0), E = 1/8, h = 1/20. Seeds 45900 and above are reserved for `verify/`.

## RED log

`julia --project=<scratch>/env05 <scratch>/clv-red/red_clv.jl <repo>` (suite only, no implementation): Pass 12, Error 29, Total 41, 5 s. No parse error, no helper error.

| Testset | Outcome | Cause |
|---------|---------|-------|
| AC1 | Error | `UndefVarError: ClvGinelli` |
| AC1 negative controls | Pass (8) | no new name |
| AC2 | Error | `UndefVarError: ClvGinelli` |
| AC2 negative control | Pass (2) | no new name |
| AC3a, AC3b, AC3c, AC4, AC5 | Error | `UndefVarError: ClvGinelli` |
| AC5 negative control | Pass (2) | no new name |
| AC6a, AC6b, AC6c, AC7 | Error | `UndefVarError: ClvGinelli` |
| AC8, AC8b, AC8c, AC8d | Error | `UndefVarError: ClvGinelli` |
| AC9 | Error | `UndefVarError: ClvGinelli` |
| AC9c | Error | `UndefVarError: clv_initial_coefficients` |
| AC10a, AC10b, AC10c | Error | `UndefVarError: ClvGinelli` |
| AC11, AC11b, AC12a, AC12b | Error | `UndefVarError: ClvGinelli` |
| AC14a, AC14b, AC14c, AC14d, AC14e | Error | `UndefVarError: ClvGinelli` |

## Stage 04 note

Criterion 13 (static null-space algorithm in 256 bits), the recorded lines R1 to R6, the full sizes (`CLV_FULL_RUN`: 20 starts, n_f = n_b = 2000, 200 starts for AC7 c) and the p = 4 error check of AC4 at the full set run there. Stage 03 must not loosen a bound in this file without an amendment; the bounds changed against the spec are U3 to U8 and are listed so the owner can accept or reject them at the checkpoint.

## Amendment D123 (stage 03, owner: "Amend AC10c")

AC10c: `a == b` replaced by `abs(b - a) <= 256`, and `0 <= b2 - a2 <= 48 + 256` by `-256 <= b2 - a2 <= 48 + 256` (Float64; `sizeof(T) * 6 + 256` in general). Reason: `@allocated` noise of the result arrays (up to 144 bytes between identical runs, measured at stage 03); 200 extra handoffs with any per-handoff allocation would add at least 3200 bytes, so the test still separates the fault. Decision log D123.

## Amendment D124 (stage 04, owner: "Amend + fix M1-M2")

- U3 restated: the pilot spectrum of `clv_pilot` runs 1000 handoffs (was 300) and the transient constant of `clv_nstar` is c = 120 (was 40); cap unchanged (3000 package, 6000 full). Reason (stage 04 report): at the full size the 300-handoff pilot gave lambda_1 < lambda_2 on F3 start 4 (negative n*, AC3b, AC5, AC7 errored), and F3 start 19 reached 3.3e-8 at c = 40 against 1e-10 (6.8e-14 at 3 n*).
- Criterion 13 (verify only): module comparison at c = 120, the 1e-30 comparison at c = 180.
- AC3c and AC11 rebuild ln Lambda as the per-unit-time exponent times Dt = T(tau) T(h), the module's Dt after M2 (was the Float64 product tau h).
- New testset "D124 step length h": `ArgumentError` for h in {0, -1, NaN, Inf}; a `BigFloat` or `Rational` h kept; Double64, tau = 20: local and forward exponents at h = 0.05 times T(h) equal those at h = 1 within 8 eps.
- AC5 (a): the bound adds the rounding of log|omega| near a tangency of the pair, 64 eps (1 / |omega_n| + 1 / |omega_{n+1}|). Found by the FULL1 rerun under D124: F3 chaotic start 6 has |omega(v_2, v_3)| down to 3.6e-6, ratio 7.99 without the term, 0.049 with it.
