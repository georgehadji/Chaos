# Verification: standard-map

Spec: `../../01-spec/output/standard-map-spec.md` (criteria 1 to 19 and Amendment D28, criteria W1 to W8). Package tests: `test/standard-map/` (`sm_*.jl`). Oracle script: `verify/standard-map.jl`, run as `julia --project=verify -t 4 verify/standard-map.jl` from the repo root; full printed output in `standard-map-verify-log.txt` (one line per criterion: `AC<n> | quantity | value | bound | verdict`).

Run facts. Julia 1.13.1, 4 threads, base seed 20261012 with offsets 40000 and above (the package tests use offsets below 20000, so the samples here are fresh). Library state: working tree of branch `m2-p1-modules` including the uncommitted D28 files (`src/Systems/wrapped_standard_map.jl`) at the time of the run. Timings in the log (and the 17 min wall time of the whole script, 1039 s, of which about 370 s is compilation of the Y5 and Double64 paths and 295 s the first CairoMakie call) were taken on a loaded machine (another agent's Julia process was running); they are not benchmarks. Values are "worst ratio of error to bound" unless stated; a criterion passes at ratio <= 1.

Sample sizes: criteria 1, 2, 4: Float64 500 points per K (18 K slots: 17 listed plus uniform, 9000 points), Double64 100 per K (1800); 3, 5, 6: 200 orbits per listed K (3400); 7 as in the package test (1000 and 200 orbits, 20 tangent-product starts); 14: 1000 orbits per listed K (17000) against 256-bit; effective sin/cos error: 259346 evaluations; O6: 10008 points; 15 to 17: 10 draws per N in {1,2,3,5}; 18: 10^5 orbits per K; 19: 100 orbits; W1, W2, W5, W6: 500 (Float64) or 200 (Double64) per K in {0.5, 2, 5, 10}. The measurement functions are those of the package tests (the test files included up to their first `@testset`), so a definition error in a test helper would not be caught here; the independent pieces are the 256-bit references, the ForwardDiff Jacobian, ChaosTools and the plain-formula 256-bit orbit.

## Criterion 19: ChaosTools convention and conversion

The verify environment holds ChaosTools 3.5.4 and DynamicalSystemsBase 3.19.3. Neither contains a predefined standard map (`PredefinedDynamicalSystems` / `Systems.standardmap()` is not in the manifest; the string appears only in a commented line, `ChaosTools/test/stability/fixedpoints.jl:26`). The standard map that ChaosTools itself defines and uses is `standardmap_rule` in `ChaosTools/test/chaosdetection/01test.jl` lines 3 to 13 (`mod2pi` variant in `docs/src/chaos_detection.md` lines 97 to 106), read from `~/.julia/packages/ChaosTools/0SbsT/`:

```
theta = x[1]; p = x[2]            # state (theta, p)
p += par[1]*sin(theta)            # parameter k = par[1]
theta += p
theta, p reduced to [0, 2 pi)     # after every step (while loops, or mod2pi)
```

Canonical form of this module (spec Equations item 1): p' = p + K sin x, x' = x + p'. Conversion: x = theta, p = p, K = k; ordering (x, p) is the same. The conversion is the identity. The only difference is that ChaosTools wraps both variables every step, so the comparison is on the torus: the Chaoslab lifted orbit against the ChaosTools orbit, torus distance computed at 256 bits (difference reduced modulo 2 pi exactly), which is the exact version of "lifted orbit wrapped at step k". Both ChaosTools variants (the `mod2pi` rule from the docs and the while-loop rule from the tests) were run, driven through `DeterministicIteratedMap` and `step!`; initial points (uniform in [-4 pi, 4 pi]^2, unwrapped) go to both systems unchanged.

## Evidence table

"Where": T = package test (`test/standard-map/`), V = `verify/standard-map.jl`. All values from `standard-map-verify-log.txt`.

| Criterion | Quantity | Value | Bound | Verdict | Where |
|---|---|---|---|---|---|
| 1 | tangent vs ForwardDiff, Float64, 9000 pts | 0 | 1 | PASS | T, V |
| 1 | tangent vs complex step 1e-30, Float64 | 0.1485 | 1 | PASS | T, V |
| 1 | tangent vs ForwardDiff, Double64 (eps 2^-100), 1800 pts | 0 | 1 | PASS | T, V |
| 1 (O6) | ForwardDiff n-step Jacobian vs product of analytic DT, 10008 pts, n = n_H(K) | 0.01926 | 1 | PASS | V |
| 1 (O6) | product of analytic DT vs 256-bit product at the same points | 0.01917 | 1 | PASS | V |
| 1 (O6) | ForwardDiff n-step Jacobian vs 256-bit product | 0.02402 | 1 | PASS | V |
| 2 | det and symplecticity residuals, Float64 | 0.7875 | 1 | PASS (marginal) | T, V |
| 2 | same, Double64 (eps 2^-100) / against nominal 2^-104 | 0.04969 / 0.795 | 1 | PASS | T, V |
| 3 | inverse round trip / tau(n), Float64, Double64 | 0.02946, 0.0006006 | 1 | PASS | T, V |
| 4 | inverse tangent, Float64, Double64 | 0.25, 0.005615 | 1 | PASS | T, V |
| 5 | lifted vs wrapped / tau(k), Float64, Double64 | 0.121, 0.002364 | 1 | PASS | T, V |
| 6 | K to -K symmetry / tau(k), Float64, Double64 | 0.07297, 0.0009569 | 1 | PASS | T, V |
| 7 | orbit defects (Float64, Double64), tangent-product defects | 0 / 0, 0 / 0 | 0 | PASS | T, V |
| 8 | K=0 exponent error / 8 eps, n up to 1e6, Float64, Double64 | 0.02472, 0.0009313 | 1 | PASS | T, V |
| 8 | abs(n X / ln n - 1) at n = 1e4 (both types) | 5.429e-10 | 1e-9 | PASS | T, V |
| 8 | max abs X(n) for w0 = (1,0) | 0 | 0 | PASS | T, V |
| 9 | trace, eigenvalue/angle ratios, wrong classes, Jordan, Float64 | 0, 0.01875, 0, 0 | 1, 1, 0, 1e-12 | PASS | T, V |
| 9 | same, Double64 | 6.1e-37, 7.3e-05, 0, 0 | same | PASS | T, V |
| 10 | steps leaving (0,0), 17 K, 1e4 steps, both types | 0 | 0 | PASS | T, V |
| 11 | trace at (pi,0) and (0,0), classes, both types | 0 and 5.9e-37, 0 | 1 | PASS | T, V |
| 12 | eigenvector ratios 0.618 and -1.618: error; eigen residual; printed 0.62 | 0; 2.2e-16; 0.001966 | 1e-12; 1e-12; 0.005 | PASS | T, V |
| 13 | O1 after 2 steps; P1 after 18; P1 after 9 (torus distance) | 2.816e-04; 2.995e-05; 1.786 | <= 6.283e-04; <= 6.283e-04; > 6.283e-04 | PASS | T, V |
| 14 | Float64 vs 256-bit orbit / tau(n_H), 17000 orbits | 0.0001174 | 1 | PASS | T, V |
| 14 | effective Double64 eps (max rel. error of sin, cos, 259346 evaluations) | 1.937e-31 (2^-102.0) | <= 2^-100 = 7.889e-31 | PASS | T (small), V |
| 14 | Double64 vs 256-bit orbit / tau(n_H) with that eps (and with nominal 2^-104) | 1.117e-05 (4.386e-05) | 1 | PASS | T, V |
| 14 | Chaoslab 256-bit orbit vs independent plain 256-bit orbit | 0 | 1 | PASS | T, V |
| 15 | Y5 symplecticity and determinant, A and A_S, Float64 / Double64 | 0.05658 / 0.0052 | 1 | PASS | T, V |
| 16 | exponents, block-diagonal A, Float64 | 0.03906 | 1 | PASS | T, V |
| 16 | exponents, A_S (bound 16 eps + ln kappa(S)/n) | 0.9764 | 1 | PASS (marginal) | T, V |
| 17 | Y5 inverse, Float64 / Double64 | 0.02309 / 0.001535 | 1 | PASS | T, V |
| 18 | CPU-backend kernel (`propagate`) vs plain loop, 10^5 orbits, K = 1 (n 14) and K = 5 (n 6): x and W bitwise equal (torus distance 0) | bitwise | tau(n) | PASS | V |
| 18 | threaded vs serial loop, bitwise; chunk 777 vs default; vs serial `propagate_orbit` | bitwise | bitwise | PASS | T (map level), V |
| 18 | GPU vs CPU; Double64 GPU compile gate | - | tau(n) | STAGE04-GPU (`verify/gpu.jl`) | |
| 18 | literal "1 vs 8 threads" | only serial vs 4 threads run | bitwise | NOTE | V |
| 19 | ChaosTools (`mod2pi` rule) vs Chaoslab, 100 orbits, k <= n_H, exact torus reduction | 0.01923 (per K: 0.019, 0.0074, 0.0089, 0.0020 for K 0.5, 1, 2, 5) | 1 | PASS | V |
| 19 | same, Chaoslab.wrap in Float64 instead of exact reduction | 0.03235 | 1 | PASS | V |
| 19 | ChaosTools (while-loop rule), exact / Float64 wrap | 0.02509 / 0.01688 | 1 | PASS | V |
| 19 | sanity, n = 100 (unscored) | 3.1 (decorrelated chaotic orbits) | unscored | INFO | V |
| W1 | wrapped state vs lifted orbit / tau(k); range and integrality, Float64, Double64 | 0.04475, 0.0006249 | 1 | PASS | T, V |
| W2 | lift at 256 bits vs lifted orbit / tau(k), Float64, Double64 | 0.04662, 0.0005272 | 1 | PASS | T, V |
| W3 | 1e6 steps, K = 5: x, p in [0, 2 pi), integer counts below 2^53 (final m = 94560675, n = 106) | true | true | PASS | T, V |
| W3 | teacher-forced 256-bit increments, first 50 steps; free run to n_H counts equal; free run torus ratio | true; true; 0.002254 | true; true; 1 | PASS | T, V |
| W4 | K = 0 known answer, m_k = floor((x0 + k p0)/2 pi) at 256 bits, k <= 1000 | true | true | PASS | T, V |
| W5 | wrapped tangent product vs 256-bit; vs lifted product (Float64) | 0.07279; 0.0004099 | 1 | PASS | T, V |
| W5 | same, Double64 | 0.003006; 7.973e-06 | 1 | PASS | T, V |
| W6 | round trip, counts exactly equal, Float64, Double64 | 0.07371, 0.001891 | 1 | PASS | T, V |
| W7 | allocations after warm-up, Float64 | (0,0,0,0,0) | 0 | PASS | T, V |
| W7 | allocations after warm-up, Double64 | (0,0,0,0,0) | 0 | PASS (final run; see Review findings) | T, V |
| W8 | overflow guard NaN, health error "state", two below limit exact | true | true | PASS | T, V |

Count: 19 spec criteria plus W1 to W8 (27 criteria). PASS on every measured line (90 PASS lines in the final run); criterion 18 GPU parts are STAGE04-GPU; Lyapunov-estimator comparisons are not in this slug (DEFERRED lyap-benettin 04 for exponent comparisons with ChaosTools).

## Figures

- `standard-map-agreement.png`: Float64 orbit against the 256-bit orbit from identical data, n = 1 to 40, K = 1 (n_H = 14) and K = 5 (n_H = 6), 1000 seeded orbits each (seed offsets 40701, 40702). Solid red: ensemble maximum of the lifted max-norm error; dashed red: median; orange: ensemble maximum of the torus distance; black: tau(n) with the ensemble minimum over orbits (the conservative bound line); dotted: ensemble maximum of tau(n); dash-dot vertical: n_H(K). Log y axis. The measured error stays below the minimum tau over the whole plotted range at both K. At K = 5 the error grows exponentially to about n = 26, then flattens (the Float64 and 256-bit orbits have decorrelated; the torus distance saturates near 1.1), so beyond that point the curve measures decorrelation, not rounding growth. The worst ratio of error to tau for n <= n_H is 0.017 (K = 1) and 0.0064 (K = 5); the bound is conservative by a factor of about 60 to 160 at the horizon.

## Review findings

Code review (two independent reviewers, stage 04): no CRITICAL or HIGH findings; the MEDIUM findings were fixed in commit 1361767. The D28 amendment (wrapped state with winding counts) was reviewed and committed in 55e5987 before the final run.

Findings on the evidence:
- W7 Double64 (resolved). The first full run measured (160, 160, 0, 0, 16000) bytes. A rerun of the identical script order up to W7 against the committed source gave (0, 0, 0, 0, 0), `code_warntype` of `forward` on the wrapped Double64 map is fully inferred, and the final full run (857 s) gives (0, 0, 0, 0, 0) and 90 PASS lines. The likely cause is that the first run loaded the D28 source before the `_pi` constant replaced `oftype(x, pi)`, which allocated for Double64 (noted by the D28 implementer). The log in this folder is the final run.
- Criterion 16, A_S at 0.976: with n >= 10 the term ln(kappa(S))/n is about 10^14 times 16 eps, so the ratio measures how close the worst draw comes to the analytic inequality smax(S^-1 B S) <= kappa(S) smax(B), not rounding. A value near 1 is expected for a draw whose dominant direction aligns with the extreme directions of S; it cannot exceed 1 except by rounding. Not marginal.
- Criterion 18 compares serial with 4 threads (bitwise) rather than 1 with 8 threads; orbits are independent work items, and tangent-core AC6 covers a 1-thread child process against 4 threads. Accepted.
- Winding counts (D28): m grows like n^2/2 on an accelerator-mode orbit and reaches 2^53 after about 10^8 steps; the guard turns this into a health error. Owner decision pending at the Batch A checkpoint.

## Notes

- **W7 Double64 allocation (resolved, see Review findings).** In the first full 17-minute run `smw_allocs` returned (160, 160, 0, 0, 16000) bytes for forward, inverse, tangent, inverse_tangent and `propagate_orbit` (100 steps) with the wrapped map in Double64, i.e. 160 bytes per forward and per inverse step. The same function in a fresh process returns (0, 0, 0, 0, 0): alone (`SM_ONLY=W`), with the AC1 to AC13 sections before it, with O6, AC18 and AC19 before it, with AC14 before it, with AC15 to AC17 before it (the full sequence was not rerun a second time), and in a separate script with ChaosTools, CairoMakie and BenchmarkTools loaded; all give 0 bytes. Float64 was 0 in every run. The original FAIL line is kept in the log unedited. The package test W7 passes in its own process. Treat as: the module allocates nothing in isolation, but the measurement is context dependent, and the cause has not been found. Owner decision: accept the isolated result, or send back to 03.
- **Criterion 2 marginal (Float64): worst ratio 0.7875** of the bound 8 eps (1 + |K|). Double64 against the nominal 2^-104 gives 0.795, so the Double64 symplecticity residual would also meet the nominal eps on this sample (the test plan, "Spec ambiguities" 1, reported 1.4 times the bound at K = 20 on its sample; not reproduced here at 100 per K, sample dependent).
- **Criterion 16, A_S: ratio 0.9764** against 16 eps + ln(kappa(S))/n. Close to 1; the bound for the non-normal variant is derived and loose in the other direction for most draws, but the worst draw uses 98 percent of it. Not a failure.
- **Criterion 14, effective Double64 eps.** The measured maximum relative error of sin and cos in Double64 against BigFloat-256 is 1.937e-31 (2^-102.0, sin at x = -141.3), on 259346 evaluations (119346 orbit arguments, 100000 uniform in [-4 pi, 4 pi], 40000 uniform in [-1e3, 1e3]); 99.9 percent quantiles 2.0e-32 to 7.6e-33. The pass limit is 2^-100 = 7.889e-31. The test plan measured 7.2e-31 on a 2000-point sample; this larger sample does not reproduce that, so the 2^-100 used for all Double64 tolerances in the tests is a factor 4 looser than what DoubleFloats `sin` and `cos` delivered here. The nominal 2^-104 is still not reached (1.9e-31 is 2^-102). Owner confirmation of the eps 2^-100 convention (test plan item 1) stands.
- **Criterion 18, literal "1 and 8 threads".** The script runs with 4 threads; it compares `propagate` (CPU KernelAbstractions backend), the threaded loop and the serial loop, all bitwise equal. A run at 1 and at 8 threads was not made; each orbit is one independent work item, so thread count cannot change an orbit result unless the kernel launch reorders reductions (none here). GPU comparison and the Double64 GPU gate: `verify/gpu.jl`.
- **O6 bound.** The spec states no bound for O6; the script uses the derived `8 n eps ||W0||_F prod ||DT_k||_F` of the tangent-core criterion cited in W5 (||W0||_F = sqrt 2). ForwardDiff and the analytic tangent agree bitwise at one step (criterion 1 ratio exactly 0), because the dual-number arithmetic performs the same floating-point operations as the analytic formula.
- **Criterion 13.** P1 after 9 steps is 1.786 away on the torus (0.284 in unit coordinates), consistent with the 0.28 and 0.26 quoted in the spec.
- **Criterion 19 convention.** No predefined map exists in the environment; the comparison is against the rule ChaosTools uses in its own tests and docs (see the section above). The conversion is the identity; the agreement ratios (at most 0.03) confirm it.
- **Script runtime** 17 min with 4 threads on a loaded machine, above the 10 minute target; most of it is first-call compilation (Y5 Double64 sections 373 s, CairoMakie 295 s, AC14 116 s).
- tau(n) is conservative: every measured orbit ratio is below 0.13 (largest 0.121 for AC5), so these checks would not detect an error below about ten times the rounding level; the exact-answer criteria (7, 8, 9, 10, 11, 12, 13, W4) and the independent 256-bit comparisons carry the sharper evidence.

## GPU

Measured by `verify/gpu.jl` (decisions D29 and D33), raw output `gpu-verify-log.txt`. Machine: GeForce GTX 1050 (sm_61, 3 GB, display GPU, WDDM driver 581.83), CUDA.jl 6.4.2 with runtime 12.9.0, Julia 1.13.1, 4 CPU threads. Windows TDR registry: TdrDelay 8 s, TdrLevel 3 (the default delay is 2 s). Every GPU line below was checked for `provenance.backend == gpu`; an element type or run that ended on the CPU is printed as FALLBACK(reason), which is a result of the fallback path and not a kernel pass. Since D33 a launch is kept under a 1 s estimated budget (gpu_ns_per_step x N x orbits, at least one wave of 1024 orbits); a run beyond it goes to the CPU with the reason "launch exceeds the display watchdog".

| Criterion | Quantity | Value | Bound | Verdict |
|---|---|---|---|---|
| 18 | GPU vs CPU, K = 1, n = n_H = 14, 1e5 seeded orbits in [-4 pi, 4 pi]^2: max over orbits of torus distance / tau(n) (largest distance over the first 2000 orbits 5.587e-13) | 2.854e-05 | 1 | PASS |
| 18 | GPU vs CPU, K = 5, n = n_H = 6, 1e5 orbits (largest distance over the first 2000 orbits 2.011e-11) | 1.878e-04 | 1 | PASS |
| 18 | Double64 kernel gate (compile or reason): an explicit `:gpu` run with `StandardMap{Double64}` ends on the CPU with the reason in provenance | kernel does not compile for the element type | compile or reason | FALLBACK(kernel does not compile for the element type) |

The compared quantity is the final lifted state from `propagate` with `k = 2` and the identity as initial deviation matrix, the bound is the spec's tau(n) = 16 eps M_n n Gamma_n evaluated on the plain CPU orbit of each start, and the distance is the 256-bit torus distance of the spec. The GPU stays within 2e-4 of the bound on both K. The Double64 kernel does not compile for the GPU (cause in the tangent-core GPU section).
