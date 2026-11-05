# Test plan: gali

Spec: `workspaces/module-build/stages/01-spec/output/gali-spec.md` (approved at the stage 01 checkpoint, D106: contract G1, P1 to P9 as proposed; Amendment D107 approved at the stage 02 checkpoint, 2026-11-05, settles the deviations below). Shared symbols are cited from the sali spec (A2 to A16, Amendments D67, D68), the coupled-maps spec (items 1, 3, 7 to 10, criterion 13, Amendment D105) and the Benettin spec (L1 S2, S3, S7, Amendment D35). Suite: `test/gali/` (`runtests.jl`, `ga_helpers.jl`, `ga_ac01_04.jl`, `ga_ac05_08.jl`, `ga_ac09_12.jl`, `ga_ac13_16.jl`, `ga_ac17_20.jl`), registered automatically by `test/runtests.jl`.

## API proposal for stage 03

Layer `Chaoslab.Indicators` (`src/Indicators/`), the reducer protocol of `Chaoslab.Tangent` (`reducer_init`, `reducer_step`, `reducer_update`, `reducer_final`, `reducer_failed`, `reducer_error`, `reducer_provenance`), as `SALI`. Names the tests use:

| Name | Signature and meaning |
|------|----------------------|
| `GALI()` | reducer for `propagate` and `propagate_orbit` with `k = K` columns, `2 <= K <= d`. `reducer_init(GALI(), T, Val(d), Val(K))` throws `ArgumentError` for `K < 2` or `K > d`. Isbits state (GPU-compatible); the one replacer of its chain, placed last (D71). |
| `reducer_step(GALI(), st, n, x, W)` | at the handoff, for each column `i = 1..K`: scaled norm `alpha_i`; a zero or non-finite `alpha_i` records the failure (kind `"deviation matrix"`, the column, value, step), leaves sums and `W` unchanged, no throw; else adds `log(alpha_i)` to the compensated sum `S_i` and replaces column `i` by the unit column. Then evaluates `G_k`, `k = 2..K`, from the stored unit columns by the nested route (spec Algorithm 3). |
| `reducer_update`, `reducer_failed`, `reducer_error` | as `SALI`: `reducer_update` throws `reducer_error(r, st, 0, :cpu)`; `reducer_error` is `NumericalHealthError("deviation matrix", value, threshold, step, orbit, backend, column)`. |
| `reducer_final(GALI(), st)` | `(n, gali, L, sum, comp, floor)`: `n` last successful handoff; `gali::SVector{K-1,T}` = `(G_2, ..., G_K)` at `n`; `L::T` = `-(compensated left-to-right total of S_1..S_d)` when `K = d` (spec Algorithm 4; for `d = 2` bitwise the sali `L`), `nothing` when `K < d`; `sum`, `comp` `SVector{K,T}` (`S_i = sum[i] + comp[i]`); `floor::SVector{K-1,T}` = `gali_floor(T, k)` for `k = 2..K`. At `n = 0`: `gali` zero, sums zero, `L` zero (K = d), `floor` set. |
| `reducer_provenance(GALI())` | `(indicator = "gali",)`, merged into the run record next to `tau`, `k`, `precision`. |
| `gali_values(W)` | pure function: `SVector{K-1}` of `G_k`, `k = 2..K`, of the `d x K` matrix `W` (`SMatrix` or `Matrix`): each column divided by its scaled norm, then the nested evaluation (B2 iv). Any element type, BigFloat included. |
| `gali_value(W, k)` | `G_k` of the first `k` columns; `ArgumentError` for `k` outside `2..size(W, 2)`. |
| `gali_floor(T, k)` | `B_meas(T, k) = 2 k effective_eps(T)` (B9). |
| `check_frame(GALI(), W0; orbit = 0)` | frame check (spec Algorithm 1): returns `G_K(0)` of the `K` columns if `>= 10 gali_floor(T, K)`, else throws `NumericalHealthError("initial frame", value, threshold, 0, orbit, :cpu, 0)`. Called by the Experiments layer for a supplied `W0` before the first chunk (as `check_frame(SALI(), ...)`; `propagate` does not call it, Amendment D107e). |

Series are read with a test-side probe reducer that records `reducer_final` at chosen handoffs (as the sali suite).

Two remarks on the proposal as written (no change made, see "Deviations"): the last sentence of the `check_frame` row does not match the code: `propagate` and `propagate_orbit` do not call `check_frame` (the Experiments layer does, `src/Experiments/run_indicator_record.jl` lines 234 and 359, before the first chunk), so the tests call it directly and check the "before any work" order by calling it first, as the sali suite does; and `reducer_final` has the two-argument form only (the workspace form defaults to it).

## Conventions

Element types Float64 and Double64 (D24; Double64 eps = 2^-100 through `effective_eps`), BigFloat(256) for references and for the oracle-grade module runs (T = BigFloat in AC5, AC6(e), AC11, AC13, AC15(b)); Float32 is out.

Test-side constants (`ga_helpers.jl`): `B_meas(T, k) = 2 k eps(T)` (`ga_meas`; `gali_floor` is compared with it); `B_G^k(n) = 4 k eps [Lambda_n^(k) + G_k,ref(n) sum_{i<=k} Psi_n^(i)] + B_meas(T, k)` (`ga_BG`); `B_L(n) = 32 eps sum_{i<=d} Psi_n^(i)` (`ga_BL`, multiplier 32 where the spec proposes 8, see Deviations 8). `Lambda_n^(k)` from the reference at every map step (`ga_Lambda`, ratio clipped at e^700, a scale). `Psi_n^(i)` for d-dimensional tangents (`ga_psis`, the lb_Psi recursion with the spectral norm `opnorm`, Float64 scales, one sweep per checkpoint for all columns).

References. `ga_ref`: the BigFloat(256) product `Y_n W_0` of the one-step tangents (`tangent` of the system in BigFloat at the widened realised orbit points, `ga_tanB`), column-normalised at every step, volumes `G_k` by a twice-orthogonalised modified Gram-Schmidt of the unit matrix (route B2 iv, `ga_vols`), Float64 `ln G_k` at every step (for Lambda), BigFloat `G_k` and column sums `S_i` at the checkpoints; for K = d the log reference `L_ref = ln abs(det W_0) + sum_i ln abs(det M_i) - sum_i S_i` (exact for the product, valid at any depth; the volume of unit columns rounded to 256 bits is meaningless below 2^-256). BigFloat module runs use a 512-bit reference. The four routes of B2 (minors, Gram determinant, singular values by Jacobi, QR) are `ga_route_*`. Checkpoints: multiples of 10 at 10 per decade (`ga_checkpoints`, D38), every multiple of 10 with `GA_FULL`. Exact spectra at fixed points: `ga_exponents` (Laplacian eigenvalues in closed form, `lambda + 1/lambda = abs(2 + g)`), `E_k = (k-1) l_1 - sum_{j=2..k} l_j` (`ga_Ek`). Singular values of long 256-bit products (AC13): `ga_logsv` accumulates the product of the 256-bit tangents at `hp` bits chosen from the Float64 growth, Jacobi eigenvalues of `Y' Y`.

Systems: `GACase` (Y3 `FroeschleMap`, Y4 `StandardMapChain`, fixed-size SVector states, D105), orbits C = (0.5, 0.5, 0, 0) and D = (3, 0.5, 0, 0) at (0.5, 0.1, 1e-3) (cm item 3, item 1 labelling), the chaotic panel (2, 0.5, 0.1) (`GA_Y3B`). Frames: the first K columns of `initial_deviations(seed, j, Val(d), Val(d), T)` (nest bitwise, tested in AC8a). Probe reducer `ga_Probe` wraps the module reducer, records `reducer_final` at chosen handoffs (and whether the first two replaced columns are bitwise equal, `eq12`); `ga_Trace` reads the orbit. Test maps: `ga_Fault` (4D, zero tangent, collapsing column 3, NaN), `ga_Boom`, `ga_Seq` (time-dependent tangent for products of symplectic matrices). Negative-control doubles: `ga_ctl_noscale`, `ga_ctl_gram`, `ga_sali`, plus the partial-sum log form and the k < d identity inside their testsets.

Missing-module behaviour: every module-dependent testset starts with `_ = (GALI, ...)` naming the new symbols (`GALI`, `gali_values`, `gali_value`, `gali_floor`), so without the module it errors on its first statement (UndefVarError); no top-level code and no helper calls a new name at include time; the 5 negative-control testsets contain no new symbol and pass now.

Seeds (base 20261012 + 40000, `ga_seed(ac, sub) = base + 50 ac + sub`, range 40000..40999, StableRNGs only, every draw seeded): testset AC<n> uses `ga_seed(n, sub)` for sub = 0..6 and, as frame seeds, sub up to 7; AC20 uses n = 19 (n = 20 would leave the range; criterion 19 has no testset), the API tests n = 0 with sub 10 to 12. Frames are `ga_frame(seed, j, T, d, K)`.

Sizes: `ga_size(small, full)` with `GA_FULL = isdefined(Main, :GA_FULL_RUN) && Main.GA_FULL_RUN`; `verify/gali.jl` (stage 04) sets `GA_FULL_RUN = true` before including `runtests.jl`. Every reduction is marked "D38" at its use; the list is in "Runtime".

## Criterion to testset

Tags: S package suite, V verification environment, G local GPU. Tolerance column: bound used (B_meas = `ga_meas`).

| AC | Tag | Testset (file) | Sample | Tolerance |
|----|-----|----------------|--------|-----------|
| 1 | S | AC1 (ga_ac01_04); AC1c negative control (noscale); the Gram control is AC6d | 9 (d, K) pairs, 100 sets per pair and type (spec 1000), routes 20 sets per pair; control 180 sets | (a) relative 1e-60; (b) 2 B_meas(T, k) |
| 2 | S | AC2 | k = 2: 56 angles x d in {2, 4, 8}; k = 3: 12 cosines x 3 d; k = 4: 3 x 2; k = d: 15 random sets; per type | 3 B_meas(T, k) |
| 3 | S | AC3 | 30 sets per pair and type (spec: the AC1 sets), scale test Float64 | range 1 + B_meas; 2 B_meas(T, k) |
| 4 | S | AC4; AC4f control | d in {2, 3, 4, 8}, 300 random pairs per d (spec 1e4) + 56 constructed + 18 threshold pairs; per type | 2 B_meas(T, 2) + sali_floor(T), Deviations 3 |
| 5 | S | AC5a (standard map), AC5b (symplectic products), AC5c (Y3, Y4); AC5d control | BigFloat; a: 3 starts x 3 K x 2 tau x 2 frames; b: 6 products per d x 2 frames; c: 5 cases at n in {5, 18, 100}; d: 60 products | B_G^d + 8 eps G Sum Psi (Deviations 5); d: factor 10 for 90 percent |
| 6 | S | AC6; AC6d control | 4 maps x 40 values of m x 3 orders x 2 types; (e) 4 maps x 14 m x 3 k | 2 B_meas(T, k); m_floor; sharper bound 99 percent; (e) BigFloat |
| 7 | S | AC7a (argument errors), AC7b (non-finite), AC7c frame check, AC7c before any work; the run record (d) is in API provenance | K in {2, 3, 4}; 1000 seeded ids Float64, 200 Double64 (spec 1e4); 200 frames against sali | 5 and 20 B_meas; decision window 2 B_meas (Deviations 2) |
| 8 | S | AC8a (nesting), AC8b (sali sums, B4), AC8c (K = d = 2), AC8d (Y4) | a: 7 orbits, n in {10, 100, 1000}; b: 5 orbits x 2 types; c: 4 starts x 3 K x 2 types; d: N in {2, 5} | bitwise (sums, L); 2 B_meas (bitwise reported) |
| 9 | S | AC9; AC9e control | 9 cases (C, D, 2 chaotic Y3, 2 chaotic Y4 N = 2, Y4 N = 3 ring, regular Y4 N = 2), K = 4 and K = d, both types, about 20 checkpoints | B_G^k, B_L, 9(c) form, floor 2 B_meas, control 10 B_L |
| 10 | S | AC10a (tau), AC10b (K = 0 shear and seeded frame) | a: 4 orbits, tau in {1, 2, 5, 10}, n in {10, 100}; b: tau in {1, 10}, n up to 1e5 (tau = 10: 1e6), Double64 to 1e4 | 2 B_G^k, 2 B_L; shear B_G^k with closed-form Lambda, Psi |
| 11 | S | AC11 | 2 Y3 + Y4 systems (N in {2, 3, 5}, open and ring), 2 frames; BigFloat (d <= 6), Double64, Float64 log form | 0.5 (+ B_L at the log form), window n = (100, 200) |
| 12 | S | AC12a (stable fixed points), AC12b (neighbour) | 6 system-frames, n <= 2000; orbit C to n = 2e4 | G in [0.05, 1 + B_meas], slope 0.2; band of AC14 |
| 13 | S | AC13 | orbit D + 2 chaotic Y3 + Y4 N = 2, 3, BigFloat module, n up to 1600 (D) or 800 | offsets 1.0 and 15; orbit D exponents loose |
| 14 | S | AC14 | orbit C both types; 18 pooled other starts (6 per system) | bands of the spec; 50 percent in band at D38 (70 with GA_FULL) |
| 15 | S | AC15a (floors), AC15b (BigFloat), AC15c (log form means), AC15d (representation) | a: 20 starts, n = 600; b: 4 orbits, n = 150; c: 30 starts, n = 500; d: 3 starts, n = 2000 | median ratio in [1.5, 2.5]; B_G^k, B_L; 4 standard errors |
| 16 | S | AC16a (zero tangent, column 3), AC16b (non-finite), AC16c (overflow), AC16d (equal columns), AC16e (absorbing zero) | 6 orbits per fault map; 100 starts n = 400 for (e) | exact error fields; 2 B_meas |
| 17 | S | AC17a (orbit), AC17b (propagate vs serial), AC17c (allocation), AC17d (reader first) | a: 2 tau x 5 starts, 5 chains; b: 512 orbits, n = 200, chunk 64 and 1000, tau in {1, 4}; c: Y3 n = 1e4, Y4 N = 100 n = 300 | bitwise; 0 bytes |
| 18 | S | AC18a (shear, bitwise), AC18b (wrapped vs lifted) | a: 10 dyadic starts x 2 systems x 2 types x 2 tau; b: 7 regular and orbit starts at n in {10, 100}, 6 chaotic at n = 10 | bitwise; 64 K_sys M n^2 eps Sum rho |
| 19 | V | none (stage 04) | ChaosTools `gali`, k = 2 and k = 3, 4 | B_G^k + B_oracle |
| 20 | S | AC20 (package part) | Y4 N = 100 open and ring, K = 4, 4 orbits, n = 100 through `propagate`; N = 10 chain K = 4 and K = 20, n <= 50 | range, nesting; B_G^k, B_L |
| 20 | V | none (stage 04) | time per orbit-step against the L1 pass | ratio reported |
| 21 | G | none (stage 04) | CPU against GPU, B11 | B11 |
| API | S | API init and the final at n = 0; one handoff against hand values; failure in a column; provenance and the run record; gali_floor; gali_values and gali_value | hand matrix columns (3, 4, 0, 0), (0, 5, 0, 0), (0, 0, 2, 0), (0, 0, 0, 7): G_2 = G_3 = G_4 = 3/5, S = ln 5, ln 5, ln 2, ln 7 | 2 B_meas, 4 eps |

Not tested in the package suite (stage 04): 1 against N threads (AC17(b) uses the run's thread count), AC19, the timing part of AC20, AC21. Recorded lines R1 to R5 are not tests; AC4, AC6, AC8, AC9, AC11, AC12, AC13, AC14, AC15 print realised ratios, slopes and counts with `@info`.

## Deviations and spec items not testable as written

1. `check_frame` is not called by `propagate` or `propagate_orbit` (the API section above says it is): the Experiments layer calls it before the first chunk. The tests call `check_frame(GALI(), W0; orbit)` directly, and AC7c "before any work" calls it before the map (`ga_Boom`). Criterion 7(d) (run record): the record can carry `tau`, `k`, `precision`, the indicator name and the floors; "chain" and "frame source" are recorded by Experiments (as sali Amendment D67); not tested here.
2. Criterion 7(c), agreement with the sali frame check "except within a factor 1.0013 of the threshold": the computed values near the threshold 40 eps carry an absolute error of order eps (2.5 percent of the threshold), so the factor is untestable; the test allows disagreement within 2 B_meas(T, 2) of the threshold. Not tested: `gali_values` of a one-column matrix (K < 2 is tested through `reducer_init`).
3. Criterion 4: the SALI input carries its own floor (4 eps), so the tolerance of (a) to (e) is `2 B_meas(T, 2) + sali_floor(T)`. In (d) the relative statement `SALI/g <= 1.0013` is untestable at g = 10 B_meas (g itself is uncertain by 10 percent); tested as `SALI <= 1.0013 g + 2 tolerance`.
4. Criterion 3(e): columns are scaled by 2^+-498 (about 1e+-150) so that the scaling is exact, on unit-scale columns, Float64 only.
5. Criterion 5: the stated bound `2 B_meas(BigFloat, 2) / G_ref` for the logarithm holds for one evaluation, not after the transport of n handoffs: the reference stub exceeds it by a factor 1.69 at d = 2, K = 2 pi, tau = 5, n = 40 (G about 1e-35). The identity is tested as the absolute bound of criterion 9(c), `abs(G_d - abs(det W0) exp(L)) <= B_G^d(n) + 8 eps G_ref sum_i Psi_n^(i)`, with a 512-bit reference. In 5(d) the factors have amplitude 6 (amplitude 1.5 gives a factor 10 in only 37 percent of the products).
6. Criterion 6(d) and AC1(c): the Gram-determinant double fails the bound at m = 9 to 13 in Float64 (value 0 or 1.5e-8); beyond m = 15 its value is 0 and the error 1e-m is below the bound. The control asserts failure in at least half of the sets at m = 9..13 (stub: 150 of 150).
7. Criterion 9(b) reference: `L_ref` is `ln abs(det W0) + sum ln abs(det M_i) - sum S_i` (Conventions), not the logarithm of the 256-bit volume, which saturates near 2^-256 (it gave error / B_L up to 1e23 for a correct log form).
8. `B_L = 8 eps sum Psi` (B10) is exceeded by a faithful implementation: the reference stub reaches error / (8 eps sum Psi) = 1.54 on a strongly chaotic Y4 N = 2 start (0.94 on a Y3 start), because the per-handoff rounding of `ln alpha_i` is a few eps per column. The tests gate with the multiplier 32. Stage 04 reports the realised ratio; the spec constant should become 16 to 32.
9. Criterion 11: the window n_1 = 40, n_2 = 80 is still in the transient for the tight spectra of N = 5 (log form deviation 1.19 and 0.59), where the spec's drafting values exist for N = 3 only; the gated window is (100, 200) (allowed: n_1 >= 40, n_2 >= n_1 + 40); the [40, 80] deviations of the direct orders are counted and printed. BigFloat(256) resolves only about 1e-76, so the high orders of fast hyperbolic points are carried by the log form (K = d); unresolved orders are counted.
10. Criterion 13: "every order is resolved in BigFloat (always)" is false for strongly chaotic starts (G_4 below 1e-76 after a few hundred steps); the orders with a resolved value are tested and k = d goes through the log form. The singular values of the 256-bit product are not computable at 256 bits for condition numbers of 1e500: the tangents are 256-bit, the product and the Jacobi eigenvalue solver run at `hp` bits (`ga_logsv`). D38: 1600 only for orbit D, 400 and 800 for the other starts.
11. Criterion 14: the spec's pass rate of 70 percent is at the edge of the realised rate (stub over 60 seeded starts: about 40 in band, 0.67; the k = 2 band of 0.4 is missed by starts at 0.41 to 0.58), so a sample of 18 would fail it by chance. The package suite uses 50 percent with 6 starts per system (pooled 18); `GA_FULL` uses the spec's 70 percent with 20 per system. Stage 04 reports the full-size rate; the proposal (P6) may need 60 percent or a k = 2 band of 0.6.
12. Criterion 15(a): the ratio of floor handoffs is taken on the chaotic panel (2, 0.5, 0.1), not on orbit D (orbit D does not reach the Double64 floor of k = 3 within 3000 steps); starts reaching both floors are used, at least half required.
13. Criterion 18(b): the bound has no orbit-divergence term, so it is applied where wrapped and lifted orbits stay together: orbit C, orbit D, regular starts at n in {10, 100}, strongly chaotic starts at n = 10 only (at n = 100 the ratio is 2e7).
14. Criteria 17(c) and 20: for d = 200 the frames are Gaussian frames from a seeded StableRNG (the 200 x 200 static QR of `initial_deviations` is out of scope).
15. Criterion 16(c): the 4D overflow map is diag(mu, 1, 1/mu, 1) (symplectic in the (q, p) pairs), the transposition of diag(mu, 1/mu) of sali 13(d). Criterion 16(e): the bitwise equality of two stored columns is read by the probe wrapper from the replaced W (the stream has no other channel); the number of starts reaching it is reported (stub: 41 of 100).
16. Criterion 17(a): the read-only trace reducer precedes the replacer, so the chain rule (replacer last, D71) is not violated; 17(b): the thread count is the run's.
17. Criterion 8: bitwise equality of the sums with the sali sums and of L (d = 2) is asserted as the spec states; the nesting of G_2 and G_3 is gated by the bound and reported bitwise (stub: 63 of 63 bitwise, the SMatrix product being column independent). Stage 03 must keep that.
18. Criterion 12(a) uses 2 seeded frames per system; `G_k >= 0.05` holds for every frame tried; 12(b) gates orbit C only, the Y4 neighbour is reported (it can be a lower-dimensional torus).

## Runtime

Package suite sizes (D38, 3 minute ceiling), reductions marked "D38" in the test files: AC1 routes 20 sets per pair and 100 per pair and type (spec 1000); AC3 30 sets per pair; AC4 300 pairs per d (1e4); AC5 starts 3 (20), products 6 per d, n = 100 for K = 0.5 (1000), n = 18 for N = 3, 5; AC6 4 maps (100), (e) 4 maps (200); AC7c ids 1000 and 200 (1e4); AC8 orbits 5 to 7 (40), starts 4 (20), Y4 n = 100; AC9 nine cases (spec: 20 starts per class and horizons 1e3 and 1e4), n = 1000 (C), 300 (D), 150 (chaotic), checkpoints 10 per decade, Y4 N = 5 only with GA_FULL; AC10 two chaotic starts (20), shear to 1e5 (tau = 1) and 1e6 (tau = 10), Double64 to 1e4; AC11 2 frames (5), N = 10 only with GA_FULL, BigFloat only for d <= 6; AC12 2 frames (5); AC13 2 starts (20), horizons 800 (1600 for orbit D); AC14 6 starts per system (20), pass rate 50 percent; AC15 20 starts and n = 600 (100, 3000), 2 starts (20), 30 starts and n = 500 (100, 1000), n = 2000 (1e4); AC16(e) 100 starts, n = 400 (1000, 1000); AC17 512 orbits, n = 200 (4096, 1000), Y4 N = 100 n = 300 (1e4); AC18 10 starts (50), 4 and 2 chaotic starts; AC20 4 orbits, n = 100 (200, 1000).

Measured with a throwaway stub of the API (below): 131 s wall on 2 threads including compilation, of which 79 s are the allocation test AC17c (the stub allocates in every step and the d = 200 methods compile; with a zero-allocation module this testset is expected to take far less) and about 50 s all the rest. Estimate for the real module: 1 to 2.5 minutes, within the D38 ceiling; the largest items besides AC17c are AC20 (6 s), AC1 and AC3 (9 s together), AC14 (4 s), AC5 (8 s together). The BigFloat module paths (AC5, AC11, AC13, AC15b) depend on the real module compiling for BigFloat; a slow generic QR there is the main risk for the ceiling.

## Dry run against a throwaway stub (not in the repo)

To check the helpers and tolerances, the suite was run against a 110-line stub of the API in the session scratchpad (`scratchpad/ga02/stub.jl`: column normalisation, compensated sums, nested modified Gram-Schmidt, log form, `check_frame` by the same evaluation; nothing of it is in the repository). Result: 2 696 passed, 8 failed, 0 errors, 51 testsets, 131 s. The 8 failures are all AC17c (allocation), expected for an allocating stub. Realised ratios of the stub: AC1 worst error / (2 B_meas) 0.26 (Float64) and 0.010 (Double64); AC9 error / bound 0.07 (G, Float64), 0.002 (G, Double64), 1.54 for the log form against the spec's 8 eps (Deviations 8); AC10a 0.006; AC8 bitwise (63 of 63); AC13 orbit D l_1 = 0.00765 and l_1 - l_2 = 0.00676 at n = 3000; AC15a median floor ratio 1.66 (k = 3) and 1.83 (k = 4); AC16e 41 of 100 starts reach bitwise equal columns; AC6 m_floor 16 to 17 (Float64), 32 to 34 (Double64), sharper bound violated 0 of 160 sets.

## RED log

Command (scratchpad driver `scratchpad/ga02/myred.jl`: `using Test, Chaoslab` and the include of `test/gali/runtests.jl`; environment `scratchpad/cm02/testenv`): `julia --project=<testenv> -t 2 myred.jl <repo>`, 10 s wall (Julia 1.13.2). Log: `scratchpad/ga02/red.out` (not in the repo).

Result for the `gali` suite: 5 Pass, 0 Fail, 46 Error, 51 testsets. The 46 module-dependent testsets each error on their first statement, the line naming the new symbols: `UndefVarError: GALI` in 38, `UndefVarError: gali_values` in 7 (AC1 to AC4, AC6, AC7b, API gali_values and gali_value), `UndefVarError: gali_floor` in 1 (API gali_floor). No MethodError, no parse error, no include-time error, no error inside a helper. The 5 Pass are the negative controls (AC1c noscale double, AC4f SALI as GALI_2, AC5d identity for k < d, AC6d Gram determinant, AC9e partial-sum log form), which do not need the module and pass as controls.

```
Test Summary:  | Pass  Error  Total   Time
red            |    5     46     51  7.8s
  AC1 ... API gali_values and gali_value  |  1 Error each (UndefVarError GALI, gali_values or gali_floor)
  AC1c, AC4f, AC5d, AC6d, AC9e            |  1 Pass each (negative controls)
```

The other suites were not run in this stage (the driver runs the gali suite alone; the main project lacks the test extras).

## Rework D108 (2026-11-05)

Owner approval of spec Amendment D108 at the stage 04 checkpoint. Tests first for the code change (D108j), then the test-side changes (D108k). Nothing else in `test/gali/` changed.

New file `test/gali/ga_d108.jl` (included last from `runtests.jl`), one testset "D108j subnormal column norm is a health failure" with 20 assertions, Float64 and Double64 (BigFloat has no reachable subnormal range: `floatmin(BigFloat)` is 2e-323228497, so the guard is the same code but cannot be exercised): `reducer_step` on a frame with the column `[5e-324, 5e-324, 0]` fails with column 2, kind "deviation matrix", threshold `floatmin(T)` in `reducer_error`, finite value below `floatmin(T)`, `W` unchanged; `reducer_update` throws the same error; `gali_values` returns no finite value above 1 (non-finite for the orders containing that column); `check_frame` rejects the frame ("initial frame", step 0, the orbit id); a column of norm just above `floatmin(T)` (2 `floatmin` times a direction) is accepted, gives the values of the scaled frame within 1e-12 and passes `check_frame`.

RED (against the stage 03 code, standalone, helpers only): 6 pass, 14 fail, 0 error (the 14 are the assertions on the failure code, the threshold, `gali_values` and `check_frame`; the 6 are the accepted-column assertions). GREEN after the code change: 20 of 20.

Test-side changes (D108k):

| Item | Change |
|---|---|
| D108a | `ga_BL(epsT, psi, amax = 0.0) = 32 eps (1 + amax) sum(psi)`; the bundle carries `amax` (largest absolute coordinate of the realised orbit up to each checkpoint); `ga_check_orbit` and AC9(b), AC10a use it (AC11 fixed point: a_max 0) |
| D108b | AC9(d) floor control: computed value at most `2 B_meas + B_G^k` where the reference is below `B_meas` (reference exactly 0 or non-finite bound: neither pass nor fail) |
| D108c | `ga_Lambda` drops terms whose log ratio is not finite; AC9(a) is gated where `G_ref >= GA_GREF_MIN = 1e-70` (skipped handoffs counted in `skipped`), AC9(c) requires `G_d,ref >= max(10 B_meas, 1e-70)`; the log form 9(b) stays gated at every handoff |
| D108d | AC11 window: `(100, 200)` for `d <= 10`, `(400, 800)` for `d > 10` (snapshots at 40, 80, 100, 200, 400, 800 for Y4 N = 10); the other window is counted as before |
| D108e | AC12a: Y4 N = 5 uses eps = 0.2 |
| D108f | AC13: `ga_offsets` returns `valid[n]` (the log sum of the 256-bit singular values within 1e-6 of ln abs(det) of the product); invalid handoffs are counted and reported, not gated; stability bound `abs(c_k(n_2) - c_k(100)) <= 5.0`; `abs(c_k) <= 15` unchanged |
| D108g | AC14: pass rate 0.5 at every size (was 0.7 at the full size) |

Package-size result after the rework: gali suite 2724 of 2724 (2704 before plus the 20 of `ga_d108.jl`), 0 fail. At the package sizes none of the changes alters a verdict (the failures D108 corrects appear at the full sizes of `verify/gali.jl`).

## Rework D109 b to d (2026-11-05)

Review findings 2 to 4 (LOW), spec Amendment D109 b, c, d. Tests first. New file `test/gali/ga_d109.jl` (included last from `runtests.jl`), three testsets, 103 assertions:

| Testset | Content |
|---|---|
| D109b | `gali_values`, `gali_value`, `check_frame` on a 2 by 3 matrix (Float64, Double64, BigFloat; `Matrix` and `SMatrix`) throw `ArgumentError`; one column still rejected, d = K still accepted; the `SMatrix` message is constant (contains no size) |
| D109c | Float64, Double64, BigFloat, d in {4, 8, 200} (BigFloat: 30 instead of 200, same loops), K in {2, 4}: `@inferred gali_values(Matrix)` and `@inferred check_frame`, element type `T`, values bitwise equal to the `SMatrix` path (`===`; `==` for BigFloat), `gali_value(M, k)` equal to the order `k` value |
| D109d | Integer matrices (exact 0/1 entries, general `Int`, `Int32`, `Rational`) give `float(T)` values equal to the Float64 evaluation, on the `Matrix` and `SMatrix` paths; `check_frame` returns a `Float64`, rejects dependent integer columns with `NumericalHealthError` |

RED (against the D108 code, standalone with the helpers, `scratchpad/ga05/red109.out`): 6 pass, 18 fail, 5 error, 29 counted. A testset stops at its first error, so D109c (error at its first `@inferred`) and part of D109d (`effective_eps(Int)` in `check_frame`) are counted short. GREEN: 103 of 103; whole gali suite 2827 of 2827 (2724 + 103).
