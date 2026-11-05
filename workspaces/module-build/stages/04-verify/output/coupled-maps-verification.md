# Verification: coupled-maps

Spec: `../../01-spec/output/coupled-maps-spec.md` (criteria 1 to 21, Amendments D104 and D105). Package tests: `test/coupled-maps/` (630 of 630 pass; whole suite 1 405 539 pass, 0 fail, 1 existing broken, stage 03 rework notes). Oracle script: `verify/coupled-maps.jl`, run from the repo root as `julia --project=verify -t 8 verify/coupled-maps.jl`, then `CM_ONLY=threads CM_APPEND=1 julia --project=verify -t 1 verify/coupled-maps.jl` and `CM_ONLY=compare CM_APPEND=1 ...` (1 against 8 threads). Raw output: `coupled-maps-verify-log.txt` (227 PASS, 0 FAIL, 0 ERROR). GPU: `verify/gpu.jl` with `GPU_ONLY=coupled-maps`, raw output `gpu-coupled-maps-verify-log.txt`. Figure: `coupled-maps-orbits-cd.png`.

Run facts. Julia 1.13.2, 8 threads (1 for the thread comparison), base seed 20261012, offsets 41000 and above (package tests use offsets below 20000). The script includes the package test files with `@testset` and `@test` turned into no-ops, so the measurement functions, assertions and tolerances are the package ones; only sample sizes are the full D38 sizes. Device: NVIDIA GeForce RTX 5060 Laptop GPU, capability 12.0, CUDA.jl 6.4.2. Main run wall time 473 s (timing section 293 s of it).

Sample sizes (full D38): criteria 1, 2: Float64 14 000 points (Y3 five sets plus uniform, Y4 panel of 8 configs N <= 10, 1000 each), Double64 12 000 (N <= 5), N = 100 20 points; 3: 1000 per case; 4, 10: 200 orbits per case (N = 100: 20); 5, 6: 200 orbits per case, 2400 orbits per type; 7: 1000 Float64 and 200 Double64 orbits per system; 9: 200 orbits per case; 13: 10 per case; 14, 15: 200; 16: 200 per case, N = 100 20 orbits; 17: 20 draws per case. Exceptions: none reduced.

## Evidence table

"T" = package test, "V" = `verify/coupled-maps.jl`, "G" = `verify/gpu.jl`. Ratios are error over bound (pass at most 1).

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| 1 | tangent vs ForwardDiff / complex step / plain formula, Float64 N <= 10 | 0.055 / 0.083 / 0.059 | 1 | PASS | T, V |
| 1 | same, N = 100 | 0.039 / 0.078 / 0.039 | 1 | PASS | V |
| 1 | tangent vs ForwardDiff, Double64 (eps 2^-100), 12 000 pts | 0.0014 | 1 | PASS | T, V |
| 2 | M M_inv - I; inverse tangent vs ForwardDiff; vs plain, Float64 | 0.090 / 0.052 / 0.052 | 1 | PASS | T, V |
| 2 | same Double64; Float64 N = 100 | 0.0018 / 0.0014 / 0.0012; 0.050 / 0.039 / 0.019 | 1 | PASS | V |
| 3 | det (256 bits), symplectic, module residuals, Float64 | 0.0046 / 0.015 / 0.0055 | 1 | PASS | T, V |
| 3 | same Double64; block symmetry flag | 1.1e-4 / 4.2e-4 / 1.7e-4; 0 | 1; 0 | PASS | V |
| 4 | round trip over tau(n), 200 orbits per case, Float64 / Double64 / N = 100 | 0.015 / 4.3e-4 / 0.0059 | 1 | PASS | T, V |
| 4 | Chaoslab 256-bit orbit and inverse vs independent plain formulas | 0 | 1 | PASS | T, V |
| 5 | wrapped vs lifted torus distance, lift at 256 bits, Float64 / Double64 | 0.017 / 3.4e-4 | 1 | PASS | T, V |
| 5 | range, integrality, counts vs 256-bit increments (50 teacher-forced steps), 2400 orbits per type; torus_state / lift incl. N = 100 | true; 0 failures | true | PASS | T, V |
| 6 | wrapped inverse torus distance Float64 / Double64; counts exactly equal; overflow guard | 0.032 / 7.5e-4; true; 10 of 10 | 1; true; 10 | PASS | T, V |
| 7 | K = 0 orbit defects (Y3, N = 2, 3, 100 ring), products vs [[I, nI], [0, I]] | 0, 0 | 0 | PASS | T, V |
| 8 | K = 0 exponent law, n up to 1e6, Float64 / Double64; position vector | 0.025 / 9.3e-4; 0 | 1; 0 | PASS | T, V |
| 9 | decoupled limit bitwise: (a) Y3 mu = 0 vs Y1, (b) eps = 0 vs Y1, (c) independence, Float64 and Double64, 200 orbits per case x 1000 steps | 0 defects in all six | 0 | PASS | T, V |
| 10 | symmetries over tau(k): Float64 / Double64 / N = 100 | 0.020 / 3.0e-4 / 0.0049 | 1 | PASS | T, V |
| 11 | Y3 fixed point: origin stays (7 sets, both types), traces, class, exponent (36 pairs), Jordan block | 0 deviations, 0 wrong classes, Jordan 0 | 0; 1e-12 | PASS | T, V |
| 12 | eigenvalues of G(0) vs K + eps l_k (N up to 100, open and ring) | 0.276 | 1 | PASS | T, V |
| 12, D104a | off-block entries / (16 N eps (1 + a_max)); with the base bound 16 eps (1 + a_max) | 0.443; 44.3 | 1; info | PASS | T, V |
| 12, D104b | trace per mode 0.116; exponent 0.134; wrong classes 0; modes inside the parabolic band 10 of 4056 | | 1 | PASS | T, V |
| 13 | pairing and sum, n = 5 / n = n_H (110 orbits) | 1.6e-4, 2.4e-4 / 3.6e-6, 3.5e-6 | 1 | PASS | T, V |
| 14 | total momentum Float64 / Double64 | 0.032 / 7.1e-4 | 1 | PASS | T, V |
| 15 | printed unit-torus map vs chain, 200 orbits per case | 0.010 | 1 | PASS | T, V |
| 15 | finding: determinant of the literal reading of the printed lines at 3 points | -2.78, -6.62, -2.67 | not scored | INFO | V |
| 16 | Float64 vs 256 bits, 2820 orbits; N = 100 only | 4.4e-5; 1.6e-6 | 1 | PASS | T, V |
| 16 | effective Double64 eps (sin, cos), three samples | 6.78e-32 (2^-103.5) | <= 2^-100 | PASS | T, V |
| 16 | Double64 vs 256 bits with the effective eps, 2420 orbits; nominal 2^-104 | 1.8e-6; 2.4e-6 | 1 | PASS | T, V |
| 17 | banded action vs dense product, N = 10, 50, 100, open and ring, 20 draws per case | 0.10 | 1 | PASS | T, V |
| 17 (V) | time ratio N = 100 / N = 50, tangent_action! open / ring | 2.017 / 2.046 | <= 2.6 | PASS | V |
| 17 (V) | inverse_tangent_action! open / ring | 2.006 / 1.925 | <= 2.6 | PASS | V |
| 17 (V) | `tangent(sys, SVector) * SMatrix` (band built per call) open / ring | 2.212 / 2.217 | <= 2.6 | PASS | V |
| 17 (V) | dense M W, quadratic, information | 4.13 / 4.13 | info | INFO | V |
| 18 | Y3 10^5 orbits: threaded loop vs serial; `propagate` x and W vs plain loops; chunk 777; both parameter sets (n_H = 18 and 9); wrapped Y3 10^4 | all bitwise equal | bitwise | PASS | T, V |
| 18, D104c | 1 vs 8 threads, per-orbit hash of x and W: Y3 two sets (forward loop, stream with identity W0, stream with random W0) 10^5 orbits each; wrapped Y3 10^4; Y4 N in {2, 5} lifted and wrapped through the stream 10^4 | 0 orbits differ in all 12 cases | 0 | PASS | V |
| 18 GPU | Y3 (0.5, 0.1, 1e-3), 10^5 orbits, n = 18, GPU vs CPU: max torus distance / tau(n) (max distance 1.5e-12) | 2.1e-7 | 1 | PASS | G |
| 18 GPU | Y3 (2, 0.5, 0.1), 10^5 orbits, n = 9 (max distance 5.3e-11) | 1.7e-5 | 1 | PASS | G |
| 18 GPU | wrapped Y3, both sets, 10^4 orbits (lifted difference at 256 bits) | 1.1e-7 / 1.3e-5 | 1 | PASS | G |
| 18 GPU | Double64 gate, Y3 and Y4 N = 2, 5 | FALLBACK, reason "kernel does not compile for the element type" (probe: `compiles = false`) | compile or reason | FALLBACK (valid result) | G |
| 19 | zero allocation after warm-up (forward, tangent, inverse, wrapped, frame actions), Y3 and Y4 N in {2, 10, 100}, Float64 and Double64 | 0 bytes | 0 | PASS | T |
| 20 (a) | SALI of C, time average over n in [1e4, 1e5] (reported 0.28; single values 0.10 to 0.53) | 0.295 | [0.1, 0.5] | PASS | V |
| 20 (b) | SALI of D first below 1e-12 at n (reported about 4.7e3); SALI at 1e3 / 5e3 | 3947; 1.3e-5 / 1.6e-18 | n in [1e3, 2e4] | PASS | V |
| 20 (c) | L_N of D at n = 1e6 (Skokos2010 Fig. 7 about 8e-3) | 0.0087 | [5e-3, 1.5e-2] | PASS | V |
| 20 (d) | L_N of C at n = 1e6 (reported 1.6e-6 at 1e7) | 1.07e-5 | <= 1e-4 | PASS | V |
| 21 | ChaosTools hand system vs module, Y3 (0.5, 0.1, 1e-3), 100 orbits: orbit torus distance / tau; tangent product / bound | 3.0e-4; 0.027 | 1 | PASS | V |
| 21 | Y3 (2, 0.5, 0.1): same | 0.0012; 0.43 (0.10 with the oracle Jacobian at the module's orbit points) | 1 | PASS | V |
| 21 | Y4 N = 2, K in {1, 3}, eps in {0.2, 0.6102}: same | 0; 0 (bitwise equal orbits and products) | 1 | PASS | V |
| 21 | Y4 N = 3, same four parameter pairs: same | at most 0.0043; at most 0.31 (0 with the oracle Jacobian at the module's points) | 1 | PASS | V |
| D105a | fixed-size vs vector state, bitwise, N in {2, 3, 5, 10, 100}, open and ring, 25 (K, eps) pairs (N = 100: 4), 358 cases per type, Float64 and Double64, lifted and wrapped | 0 mismatches | 0 | PASS | T, V |
| D105b | stream on the CPU backend vs plain serial loop, N in {2, 5}, lifted and wrapped, 2000 orbits, k = 4: orbits with a different final state / frame; fixed-size types kept | 0 / 0 in all four | 0 | PASS | T, V |
| D105c | isbits system values; bytes allocated by forward, inverse and both frame products, N in {2, 10, 100}, Float64 and Double64, lifted and wrapped | 0 | 0 | PASS | T, V |
| D105d | Y4 N = 2 and N = 5 (K = 3, eps = 0.2) and wrapped N = 2, 10^4 orbits, GPU vs CPU through the stream with fixed-size states | after the rework below: 6.1e-6 (N = 2), 9.0e-6 (N = 5), 1.75e-5 (wrapped N = 2); max torus distance 1.0e-11 (first run: fell back to the CPU, "kernel does not compile for the element type") | 1 | PASS (after rework) | G |

Extra oracle items. O1c, O2, O3, O4, O5, O10 are covered by criteria 9, 16, 12, 18 and 13 above. O7: independent re-implementation of Y3 in the labelling of Skokos2003a Eq. (9) with plain `qr`, same initial frame as the module: orbit C at n = 100, 1e3, 1e5, 1e6 and orbit D at n = 100, 1e3 agree to at most 1e-9 (scored, all PASS); 10 seeded random starts at n = 200 agree to at most 1e-9; D at n = 1e5, 1e6 and the same starts at n = 2e4 differ by up to 7.8e-4 and 5.6e-3 (informational: orbits that differ by rounding separate exponentially). O6b: ChaosTools `lyapunovspectrum` against the Chaoslab QR spectrum with the same frame: orbit C at 1e5 and 1e6 agree to 5.4e-15; orbit D differs by 7.8e-4 (1e5) and 3.7e-4 (1e6), loose bound 1e-3 (K18); chain N = 2, K = 3, eps = 0.2, 2e5 steps through the stream: ChaosTools and Chaoslab (0.7502, 0.4607, -0.4607, -0.7502), difference 1.3e-14 (spec drafting check on another orbit: 0.749, 0.461). O9: ForwardDiff of the n-step map vs the product of tangents at n = n_H, 10 002 Y3 points: 0.0081 of the bound 8 n eps prod ||M_k||_F sqrt(d), product against 256 bits 0.61; Y4 N in {2, 5}, 4000 points: 0.0033 and 0.010. The log also records the QR pairing and sum for C and D (pairing defect of the QR exponents 1.9e-5 at 1e5 for C and 3.5e-6 at 1e6 for D, decaying like 1/n; sum at most 1e-16).

## Figure

`coupled-maps-orbits-cd.png`: left, finite-time exponent of C (decays like ln(n)/n) and of D (plateau near 8.7e-3, inside the grey band of criterion 20 c); centre, SALI of C (fluctuates in the blue band) and of D (falls to rounding level at n about 4e3, inside the red range of criterion 20 b); right, orbit D, pairing defects of the QR exponents and sum against n.

## Findings and notes for the owner

1. D105d failed at the first run and passes after the rework (section "Rework after stage 04"). The Y4 fixed-size path did not compile for the GPU kernel. Cause (error reproduced with a scratch script, `src/` not touched): `forward(::StandardMapChain, ::SVector)` in `src/Systems/chain_tangent.jl:79` contains `throw(DimensionMismatch("chain of $(m.N) sites ..."))`; the string interpolation allocates a `String` inside the kernel, which CUDA rejects (InvalidIRError, "unsupported allocation of an object with references"). The same check is in the `inverse` method (line 87) and probably in `tangent`/`inverse_tangent` and the wrapped methods. Y3 does not have it and runs on the device. The runtime fell back to the CPU with the documented reason, so the GPU-vs-CPU comparison for Y4 could not be made. Fix belongs to stage 03 (a constant message, or the check moved to the host side), then rerun `GPU_ONLY=coupled-maps`.
2. n_H of the published Y3 parameters is 18 in the code and tests (sigma = 2.074, 2.074^19 = 1.03e6 > 1e6), the spec says 19. Documentation inconsistency, no effect on results.
3. Criterion 20 (a) to (d) all pass. The text value 5e-2 of Skokos2003a for D is not reproduced (0.0087 at 1e6, consistent with Fig. 7 and the spec drafting check); P10 stands. SALI of D crosses 1e-12 at n = 3947 (reported 4.7e3).
4. Criterion 21: Y4 orbits and products of ChaosTools equal the module bitwise for N = 2 and nearly so for N = 3; the largest product ratio is 0.43 (Y3 at (2, 0.5, 0.1), own orbits), so the spec bound holds without an orbit-divergence term.
5. Criterion 12 base bound without the factor N gives 44.3 at N = 100, as D104a states.
6. The AC17 time measurements are minima of BenchmarkTools samples at the 100 ns timer resolution, on a machine that may have been loaded.
7. Conflicts between documentation and code carried from stages 02 and 03: `_config/architecture.md` names `step`, `invariants` and StaticArrays up to dimension 8, the code uses `forward`, `inverse`, `tangent`, `inverse_tangent`; the `forward` docstring still describes the Y1 state; spec criterion 12 base bound (superseded by D104a).

## Audit

| Check | Result |
|---|---|
| Evidence | Criteria 1 to 21, D104a to D104c, D105a to D105d have numbers |
| Review closed | no CRITICAL or HIGH finding open: D105d fixed and rerun, review findings 2 and 3 fixed, 1, 4, 5, 6 recorded with no code change |

## Independent review

Fresh Sonnet subagent, given the spec and the code (`froeschle_map.jl`, `standard_map_chain.jl`, `wrapped_standard_map_chain.jl`, `chain_tangent.jl`, the coupled-maps parts of `Systems.jl` and `residuals.jl`), not the test plan, notes or this report; read only, small scratch snippets. No CRITICAL or HIGH finding, no math error in Y3 or Y4. Confirmed: Y3 forward, tangent at the pre-step sums, inverse and inverse tangent at the image; the mu = 0 limit equal to Y1 through the shear with K = -nu; Y4 force, G, open and ring bonds; dense tangent and `ChainTangent` entries bitwise equal (including -0.0 off the band), banded action equal to `ChainTangent * W` bitwise; wrapped count recursions and the overflow guard; element types Float64, Float32, Double64, BigFloat, Complex; `@inferred`, isbits and 0 bytes at N = 100.

| # | Finding | Severity | Disposition |
|---|---------|----------|-------------|
| 1 | With eps = 0 the neighbour term is still `eps * s`: a non-finite neighbour (0 * Inf = NaN) poisons a decoupled site, so criterion 9(c) holds for finite perturbations only | LOW | no change: criterion 9(c) is read for finite perturbations (the spec domain makes non-finite input give non-finite output) |
| 2 | Vector forms `forward` and `inverse` of the wrapped chain use `similar(s)`, so a Float64 state with Double64 parameters is narrowed to Float64 (the lifted chain and the fixed-size forms promote) | LOW | fixed: `similar(s, _promoted(w.chain, s))` |
| 3 | The positional constructor `StandardMapChain(N, K, eps, ring)` bypassed the checks of N and of the ring (N = 2 ring double counts the bond) | LOW | fixed: the checks are in the inner constructor |
| 4 | At eps = 0, `K sin + eps s` gives +0.0 where Y1 keeps -0.0 | LOW | no change: the criterion 9 checks pass (0 bitwise defects in this report) |
| 5 | Compile cost of the fixed-size forms grows with N k (3.1 s at N = 100, k = 4) | LOW (risk) | recorded; no change |
| 6 | The GPU compiler may contract multiply-add; the CPU path has none | LOW (note) | no change: the GPU comparison is tolerance based (D105d) |

## Rework after stage 04 (2026-11-05)

Owner decision: "Fix all" (D105d and review findings 2 and 3), tests first. Test `test/coupled-maps/cm_review.jl` (included from `runtests.jl`): RED 7 fail of 13 (4 promotion, 3 constructor), 636 earlier assertions pass. Code (in this stage, not by a subagent): `standard_map_chain.jl` inner constructor with the N and ring checks (finding 3); `wrapped_standard_map_chain.jl` vector `forward` and `inverse` promote (finding 2); the 7 `DimensionMismatch` messages of the fixed-size methods in `chain_tangent.jl` and `wrapped_standard_map_chain.jl` are constant strings (D105d). GREEN: coupled-maps 643 of 643; whole package suite 1 405 552 pass, 0 fail, 1 existing broken (weighted-birkhoff), quality 12 of 12, layering 2 of 2. GPU rerun (`GPU_ONLY=coupled-maps`, log `gpu-coupled-maps-verify-log.txt` overwritten): PASS 7, FALLBACK 3 (the Double64 gate); Y3 rows unchanged. The oracle script `verify/coupled-maps.jl` was not rerun after the rework: the changes touch only error paths, promotion of the wrapped vector forms and the positional constructor, and the package suite covers them.

