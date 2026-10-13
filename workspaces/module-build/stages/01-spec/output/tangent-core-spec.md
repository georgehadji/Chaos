# Spec: tangent-core

Owner checkpoint 2026-10-13 (decision D24): Float32 (binary32) is dropped from M2. Every Float32 item in this spec (element-type lists, the screening tier, the Float32 GPU gate and Float32 test cases) is inactive until revisited at M4; stage 02 writes no Float32 tests.

Status: draft for the 01-spec checkpoint. Layers: Foundation pieces (backend detection and fallback, numerical health error, compensated summation, provenance fields) plus Tangent (the propagation stream). Notation: d = phase-space dimension, k = number of deviation vectors (1 <= k <= d), tau = handoff interval in map steps, T = element type, eps(T) = machine epsilon of T (spacing of T at 1; 2^-52 for Float64). "Float64", "Double64", "BigFloat" are precision labels as in `schema/run-store.md`.

## Scope

- Catalog rows (`_config/catalogs/lyapunov.md`): intro paragraph (the shared tangent stream that every variational method consumes), L7 (tau as a recorded run parameter, not a module), section D2 (the stream consumes the system's analytic one-step tangent; AD and complex-step validation belong to the systems spec, not here). `_config/architecture.md`: "Compute" and module-map rows Foundation and Tangent.
- Out of scope: any orthonormalizer or exponent estimator (Lyapunov module), diagnostics, flows and their integrators (M5), adaptive tau (L7b), banded tangents (L7c), element-type support beyond what the precision spec adds.
- Prerequisite: the standard-map spec (Y1, Y5) is `planned` in the same milestone, not `verified`. This spec states needs on a system only as a contract (below); its tests use Y5 linear maps and the K = 0 standard map, which depend on nothing else.
- Foundation rule (brief): the systems layer imports nothing from Foundation; this module creates the Foundation pieces it needs.

## Equations

| # | Statement | Source |
|---|-----------|--------|
| E1 | Map of dimension d: x_{n+1} = f(x_n); orbit x_n = f^n(x_0) | Skokos2010 Eq. (4), p.8 |
| E2 | One-step tangent map: w_{n+1} = M_n w_n, M_n = Df(x_n) (evaluated at the point BEFORE the step) | Skokos2010 Eq. (11), p.9 |
| E3 | Product form: w_n = M_{n-1} ... M_0 w_0 =: Y_n w_0; Y_{n+1} = M_n Y_n, Y_0 = I | Skokos2010 Eqs. (12), (13), p.9 |
| E4 | Symplectic maps: M^T J M = J, so every Y_n is symplectic (product of symplectic matrices) and det Y_n = 1 for d = 2 | Skokos2010 Eq. (5), p.8; consequence derived here |
| E5 | k vectors at once: W_{n+1} = M_n W_n, W in R^{d x k}, columns evolve independently because the tangent map is linear in w | Skokos2010 Eq. (53), p.31 (linearity); matrix form derived column by column from E2 |
| E6 | Composition over a block of tau steps: Y_tau(x_{j tau}) = M_{(j+1)tau-1} ... M_{j tau}; block j+1 starts from the replaced matrix | Skokos2010 Eq. (22), p.15 (composition law); Eqs. (57), (58), p.32 (replace-by-normalized vector, local expansion alpha_i) |
| E7 | Replacement algebra (derived): a handoff that replaces W by W' = W S^{-1}, S an invertible k x k matrix, leaves the true deviation matrix recoverable as W' S followed by the remaining steps. For k = 1, S = alpha (norm), the accumulated log of the S factors equals the log of the true norm | derived from E5 and linearity; consumed by Skokos2010 Eq. (59), p.32 (k = 1) and Eq. (79), p.43 (k >= 1) |
| E8 | Standard-map Jacobian (derived, for tests): map p' = p + K sin x, x' = x + p' gives dx'/dx = 1 + K cos x, dx'/dp = 1, dp'/dx = K cos x, dp'/dp = 1, so M = [[1 + K cos x, 1], [K cos x, 1]], det = (1 + K cos x) - K cos x = 1 (matches `_config/research-question.md`, "Ground truth"). For K = 0: M = [[1, 1], [0, 1]] for every x (derived) | derived |
| E9 | Y5 hyperbolic test map A = [[2, 1], [1, 1]] (integer entries, det = 2 - 1 = 1, trace 3): eigenvalues mu_{+-} = (3 +- sqrt 5)/2, exact exponent ln mu_+ = ln((3 + sqrt 5)/2) = 0.9624 (derived); A^n = [[F_{2n+1}, F_{2n}], [F_{2n}, F_{2n-1}]] with F the Fibonacci numbers (derived, induction) | derived |

## Algorithm

**Contract.** Input: a system giving the one-step map and its one-step tangent action (E1, E2); an ensemble of orbits (initial points); k; tau; total steps N (a positive multiple of tau, otherwise rejected before any work); element type T; seed; a chain of reducers. Per orbit, for step n = 0, 1, ..., N-1: apply the tangent at x_n to W, then advance x (E2, E5; tangent uses the OLD point). After every tau steps (steps tau, 2 tau, ...), hand the tuple (orbit id, step index n, x_n, W_n) to the reducer chain in fixed order. A reducer may read x and W and may replace W (E7); later reducers see the replaced matrix. Reducer state is per orbit and lives in the ensemble layout; no reducer sees two orbits. No wrapping of x inside the stream (lifted coordinates, `_config/numerics-policy.md`); wrapping is for output only. The stream applies no min, max or comparison to state values, so NaN and Inf propagate (a non-finite value cannot disappear before the next handoff).

**tau.** A recorded run parameter in units of map steps (L7; `_config/numerics-policy.md`). Skokos2010 fixes no value; it says only to fix "a small time interval tau" (p.31, after Eq. (56)) and lists the renormalization time as an input of the mLCE algorithm (Table 1, input 4, p.34) and of the standard method (Table 2, input 5, p.44).

**Initial deviation vectors.** For orbit id j and run seed s: a d x k matrix of independent standard normal numbers drawn from a stream keyed only by (s, j) with the project's stable seeded generator, orthonormalized with the sign fixed so that the triangular factor has a positive diagonal; the first k columns are then uniform (Haar) on orthonormal k-frames. A random choice gives the largest exponent with probability one (Skokos2010 Sect. 4.2, p.20-21) and a random set of p vectors gives the sum of the p largest (Skokos2010 Eqs. (66), (67), p.39). Generated on the host for every backend and copied to the device, so CPU and GPU start identical. The per-orbit stream never depends on thread, chunk, device or ensemble position. Record the run seed and the rule that maps (s, j) to the per-orbit seed (stored as the orbit attribute `seed`).

**Accumulated logarithms.** Rule: every accumulator that sums logarithms of stretch factors (and any other sum over more than 10^3 terms handed to a reducer by the core's helper) is a running sum plus a compensation term, updated on the branch of the operand with the larger magnitude (Neumaier form); the compensation term is part of the accumulator state; reassociation and fused fast-math are not permitted (`_config/coding-conventions.md`, "Tests"). Why: logs of stretch factors are positive and negative, and n up to 10^9 steps (`_config/architecture.md`, "Scalability"); a plain running sum has relative error that grows with n (measured below: about 1.9e-9 after 10^8 terms in Float64, i.e. 8.5e6 eps). Test: AC9.

**Numerical health.** At every handoff, and again once after the last step, check that x and every entry of W are finite. On failure raise a numerical health error carrying orbit id, the step index of the failing handoff, the backend, and which of {state, deviation matrix} was non-finite. Checking only at handoffs is sufficient because NaN and Inf propagate through the stream's add and multiply. A chunk with an error stores no results. On a GPU the check writes a per-orbit flag that the host reads when the chunk ends, then raises the same error.

**Proposed `renorm_interval_tau = 1`** (map steps; owner may override per study in its preregistration).

Justification, two bounds (Float64 unless stated; derived):
1. Overflow. The norm of a deviation vector grows as exp(lambda n); steps to overflow = ln(Fmax)/lambda with ln(Fmax) = 1024 ln 2 = 709.8 (Float64 range; the Double64 exponent range is to be checked in the precision spec) and 88.7 (Float32). Lambda scale: benchmarks.md, quantity "Standard map largest Lyapunov exponent, large K": lambda ~ ln(K/2) (Chirikov law, asymptotic; Shevchenko2004), cross-checked by quantity "Lyapunov exponent, single chaotic orbit": 1.172 at K = 2 pi (SalaArtuso2015, McMillan form, converted to p' = p + K sin x, x' = x + p' by p = x - y as that row states). Assumed K_max = 2 pi (largest K with a benchmark row; the panel's "large K" is not numerically fixed in `_config/research-question.md`, open question). With lambda <= 1.2: 709.8/1.2 = 591 steps (Float64), 74 steps (Float32). At K = 20, ln(K/2) = 2.30 (illustrative, from the same row): 308 steps (Float64). Overflow is not the binding bound.
2. Alignment (k >= 2). Columns line up like exp(-(lambda_1 - lambda_k) n); the component of the second column independent of the first falls below eps(T) when (lambda_1 - lambda_k) n > ln(1/eps(T)). On the standard map lambda_2 = -lambda_1, so the gap is 2 lambda <= 2.4 at K = 2 pi. Steps to total loss: ln(1/eps)/2.4 = 36.0/2.4 = 15 (Float64), 15.9/2.4 = 6.6 (Float32), 73.7/2.4 = 30 (about 32 digits, Double64); at K = 20 (gap 4.6): 7.8 (Float64), 3.5 (Float32). Skokos2010 p.39 names this small-angle problem for p >= 2.
3. Rule for overrides: tau <= tau_max(T) = floor(ln(1/eps(T)) / (2 (lambda_1 - lambda_k))), which keeps about half the digits (stretch ratio at most eps^-1/2). At K = 2 pi: 7 (Float64), 3 (Float32), 15 (Double64). For k = 1 only overflow applies: tau <= floor(709.8/(2 lambda)) = 295 (Float64, K = 2 pi).
4. tau = 1 is the only value that satisfies the rule for every element type and every K up to 2 pi, is independent of k and of the orthonormalizer (classical Gram-Schmidt loses orthogonality like eps times the condition number of W), and is the baseline of the L7 sweep (tau in {1, 2, 4, 8, ...} up to tau_max). Cost: the maximum handoff rate; studies that need speed raise tau inside tau_max. At K >= 20 in Float32 even tau = 1 is outside the margin (3.5/2 < 2); flagged for the precision spec.

## Invariants

- I1 Order: the tangent used at step n is Df(x_n); orbit and W advance together (E2, E3).
- I2 Handoffs occur exactly at steps that are multiples of tau, with the step index carried in the stream.
- I3 Replacement transparency (E7): accumulated log S_n of the S factors plus the log norm of the final replaced vector equals the log norm of the unreplaced evolution (k = 1).
- I4 Symplecticity of the accumulated product for symplectic systems (E4); the core does not enforce it, the diagnostics module measures it.
- I5 Orbit independence: orbit j's result depends only on (system parameters, initial point, seed, j, k, tau, N, T). Never on thread count, chunk size or position in the ensemble. Across backends it holds only to the tolerance of AC12.
- I6 Finiteness at every handoff (numerical health above).
- I7 Lifted coordinates are never wrapped inside the stream.
- I8 The hot loop allocates nothing after warm-up (AC5).

## Compute

Policy: `_config/architecture.md`, "Compute" (decision D17, GPU first, CPU fallback). Owner decision, M2 kickoff: extended precision is Double64 (DoubleFloats.jl, about 32 digits, isbits); BigFloat only for short CPU oracle runs.

- Data-parallel work: independent orbits (grids, indicator maps, sweeps). Per orbit the work is sequential in n. One KernelAbstractions kernel source serves CPU and GPU; there is no second CPU stepping loop. The CPU backend runs it with worker threads over orbits, dynamic scheduling, and BLAS threads set to 1 when worker threads > 1. CPU-only: BigFloat (not isbits), any reducer needing a dense factorization library call, single long orbits (one orbit does not fill a device).
- Backend resolution (`numerics-defaults.toml` key `backend`): `cpu` always runs on CPU. `auto` and `gpu` both pick GPU when CUDA.jl loads and reports a functional device and the kernel compiles for T; otherwise they run on CPU (an explicit `gpu` request falls back too; the schema says "why auto/gpu ran on CPU"). The resolved value is recorded as `backend` = `gpu` or `cpu` (never `auto`), with `device` (GPU name or `cpu`) and `backend_fallback_reason` (empty when no fallback), `precision`, `threads`, `cpu` (`schema/run-store.md`, root attributes). Allowed reasons, verbatim from architecture.md: "CUDA.jl not loadable", "no functional device", "kernel does not compile for the element type", "device out of memory", "NumericalHealthError raised only on the GPU" (the chunk is rerun on CPU; if the CPU raises too the error propagates and nothing falls back). BigFloat on a GPU request maps to "kernel does not compile for the element type" (owner decision, M2 kickoff). Reasons 1 to 3 hold for the whole run; reasons 4 and 5 apply to the chunk and are recorded per chunk file.
- Proposed additions, need an architecture.md edit and an owner decision: reason "reducer is CPU-only" (a chain containing a CPU-only reducer must resolve to CPU; the five listed reasons do not cover it) and reason "GPU not faster for the element type" (K15: GTX 1050 3 GB, compute capability 6.1, weak Float64 throughput; `auto` picks GPU for T only if the benchmark of AC13 shows GPU time <= CPU time at 10^6 orbits).
- Provenance additions for the store (schema bump, proposal): root attributes for `renorm_interval_tau`, `deviation_vectors_k`, `seed` (run seed) and `chunk_size`; the schema today has only the per-orbit `seed` and the config hash.
- K15 check, read 2026-10-13: the CUDA.jl README support matrix lists Pascal (6.0, 6.1) as supported and Julia 1.10 as CI-tested, Julia 1.13 supported (this machine has 1.13.1); it requires Julia >= 1.10. NVIDIA's CUDA 13.0 release notes drop pre-Turing GPUs in cuFFT and cuSPARSE; `nvidia-smi` here reports driver CUDA 13.0, so CUDA.jl may choose a 13.x runtime. Mitigation: pin a 12.x runtime for this device. Whether a 13.x toolkit compiles kernels for sm_61 was not read; to verify in 04-verify (if it does not compile, reason "no functional device" or "kernel does not compile for the element type" is recorded). CI has no GPU and runs the same kernel on the KernelAbstractions CPU backend.
- Memory per orbit (resident, structure-of-arrays, one array per component across orbits): M(d, k, T) = sizeof(T) (d + d k + 2 k) + 16 bytes, i.e. state d, deviation matrix d k, k log accumulators and k compensation terms, plus an 8-byte step index and an 8-byte orbit id. Reducer state is added as each reducer declares it. Scratch per worker is preallocated per chunk. Worked values: d = 2, k = 2, Float64: 8 (2 + 4 + 4) + 16 = 96 bytes; d = 4, k = 4: 8 (4 + 16 + 8) + 16 = 240 bytes.
- Chunk size from memory at run start, never a constant: chunk = floor(f * free_bytes / (M + reducer bytes per orbit)), free_bytes the free device memory (GPU) or free host memory (CPU) read when the run starts, f = 0.5 (proposed headroom for the driver and output buffers; recorded). Example: free_bytes = 2^30, f = 0.5, d = 2, k = 2: floor(536870912 / 96) = 5592405 orbits.
- CPU vs GPU agreement: never bitwise and never per trajectory (risk K16): chaotic orbits amplify FMA and reduction-order differences. Agreement is tested on accumulated logs divided by n (AC12). That comparison runs locally; its result goes in the module's 04-verify report; CI runs only the CPU backend.

## Oracle

Package test suite (small, self-contained, CI on the CPU backend): analytic Y5 linear maps (exact), the K = 0 standard map (exact), direct matrix product at BigFloat(256) along the realized orbit for short horizons, invariants.

Separate verification environment, driven from module-build stage 04 (owner decision, M2 kickoff): ChaosTools tangent integration (`lyapunov`, `lyapunovspectrum`, catalog D2) on the same system and orbit; ChaosTools cannot be handed the same initial deviation vectors, so it is compared on tau-matched, long-horizon exponent-like quantities, not on W. Settings and tolerances:
- Regular orbit, K = 0.5, a seeded initial point, n = 10^5: |lambda_1(n) Chaoslab - lambda_1(n) ChaosTools| <= 1e-4 (both decay like ln n / n = 1.2e-4, Skokos2010 Eq. (62) row in benchmarks.md, "Finite-time mLCE of a regular orbit"; the start-vector dependence is ln(c)/n with c of order 1).
- Chaotic ensemble, K = 2 pi, 1000 seeded points, n = 10^4: difference of ensemble means <= 4 combined standard errors.
- Full BigFloat orbit (not only along the Float64 orbit) for a regular case, K = 0.5, n = 100: relative error of W_n <= 16 K n^2 eps(T) (derived: an initial error of size eps in p drifts x by n eps, so a tangent entry K cos x errs by K n^2 eps).
- ForwardDiff-vs-analytic tangent along long orbits belongs to the systems spec; TaylorIntegration BigFloat applies to flows (M5).
K18: published lambda(K) values (benchmarks.md) are loose sanity checks only; no acceptance criterion below uses one as ground truth.

## Acceptance criteria

All run on the CPU backend (kernel on CPU) unless marked [local]. Seeds: run seed from `numerics-defaults.toml`; orbit ids 1..N_orb.

1. Y5 exactness. A of E9, W_0 = I, initial point integer-valued, k = 2, tau = 1, n = 1..35 (entries below 2^53 for n <= 38). Pass: W_n equals the integer matrix power bitwise, orbit x_n equals A^n x_0 bitwise. Tolerance 0.
2. Y5 renormalization and exponent. A of E9, k = 1, tau in {1, 5, 10}, n = 10^4, orbit ids 1..100, T = Float64. Pass: (a) |S_n + ln ||w_n^replaced|| - ln ||A^n w_0||_BigFloat(256)| <= 8 n eps(T); (b) |S_n / n - ln((3 + sqrt 5)/2)| <= 2e-3.
3. K = 0 standard map. W_0 = I, x_0 = 0.25, p_0 = 0.5 (dyadic), tau = 1, n in {1, 10, 1000, 10^6}. Pass: W_n = [[1, n], [0, 1]] and x_n = x_0 + n p_0 bitwise (tolerance 0). With k = 1 started at e_2 and replaced every step: |S_n - 0.5 ln(n^2 + 1)| <= 8 n eps(T) at n = 10^4, and S_n / n <= ln(n)/n + 1e-12 (regular-orbit decay O(ln n / n), Skokos2010 Eq. (62) row of benchmarks.md; the margin ln(1 + n^-2)/(2n) = 5e-13 is derived).
4. BigFloat product. Standard map, 20 seeded initial points, T = Float64 and Double64, W_0 random orthonormal, tau = 1: K = 0.5, n = 100, and K = 2 pi, n = 20. Reference: the product of E2 at BigFloat(256) evaluated at the SAME realized orbit points (isolates tangent propagation from orbit sensitivity). Pass: ||W_n - W_n^ref||_F <= 8 n eps(T) (prod_{i<n} ||M_i||_F) ||W_0||_F (derived: standard forward error bound of a matrix product; the relative form ||W_n - W_n^ref|| / ||W_n^ref|| <= 8 n eps is not used because chaotic products align their columns, so ||W_n|| can be far below the product of the norms).
5. Zero allocations. After one warm-up call, n = 10^4 steps, k = 2, a chain of two test reducers (one replacing W by its column-normalized form, one only reading), T = Float64: the allocation measurement inside the per-orbit stepping region reports exactly 0 bytes.
6. Thread and chunk invariance. K = 2 pi, 4096 seeded orbits, n = 10^3, k = 2, tau = 1; 1 thread vs N threads (N = 4, or all available if fewer, never fewer than 2), and chunk sizes 64 and 1000. Pass: final x, W and accumulated logs bitwise identical in all four runs. Tolerance 0.
7. Seeded initialization. Orbit ids 1..10^4, d = 2, k = 1 and d in {2, 4, 8}, k = d. Pass: (a) ||Q^T Q - I||_max <= 8 k eps(T); (b) the matrix for (s, j) is bitwise identical across thread counts, chunk positions and repeated calls, and differs for j != j'; (c) for d = 2, k = 1 the angle atan2 is uniform on [0, 2 pi): Kolmogorov-Smirnov statistic <= 1.63 / sqrt(10^4) = 0.0163 (1% critical value).
8. Numerical health. (a) A test map that returns Inf at step 7, tau = 5: the error is raised at handoff step 10, carries the orbit id, backend and "state", and no result of the chunk is stored. (b) NaN in one entry of W injected at step 3, tau = 5: same, "deviation matrix", step 5. (c) Overflow: A = diag(mu, 1/mu), mu = e^3 (Float64), W_0 = I, k = 1; with tau = 200 no error at n = 200 (growth e^600 is finite); with tau = 300 the error is raised at handoff 300 (e^711 > 709.8). Tolerance: exact step indices.
9. Compensated summation. (a) Sum of 10^8 copies of fl(0.1) in Float64: |computed - exact| / exact <= eps(T), exact = 10^8 fl(0.1) in rational arithmetic (the correctly rounded result passes; measured here 0.25 eps). (b) Negative control: the plain running sum on the same input has error 1.9e-9 relative (8.5e6 eps, measured) and must fail (a), showing the test discriminates. (c) Double64, 10^7 copies of fl(0.1): relative error <= eps(Double64) against rational exact.
10. Backend resolution and provenance. (a) CI, no CUDA: `auto` gives `backend` = cpu, `device` = cpu, reason "CUDA.jl not loadable" or "no functional device". (b) `cpu` request: reason empty. (c) `gpu` request without a device: falls back, reason recorded. (d) BigFloat with `gpu`: cpu, reason "kernel does not compile for the element type". (e) Each of the five reasons can be forced by a test hook and is recorded as the exact string. (f) The provenance record contains backend, device, reason, precision, threads, cpu, seed, tau, k, chunk size; `backend` is never `auto`. Tolerance: exact string equality.
11. Chunk size. For free_bytes = 2^30, f = 0.5: d = 2, k = 2, Float64 gives 5592405; d = 4, k = 4, Float64 gives 2236962 (floor(536870912/240)). With free_bytes halved, chunk halves (floor). Tolerance 0.
12. [local] CPU vs GPU, Float64, report in 04-verify. Compare S/n (accumulated log divided by n), never per-step values or trajectories. (a) K = 0.5, 1024 seeded orbits, n = 10^5, symmetric form of Amendment D33: each backend against a Double64 reference; the counts of orbits above 1e-6 differ by at most 3; where both agree with the reference the CPU-GPU gap is at most 2e-6. (b) Same-orbit short horizon (orbits not yet diverged): K = 2 pi and K = 10, n = 10: |x_cpu - x_gpu| <= 4 n eps G^n with G = sqrt((1 + K)^2 + K^2 + 2) (derived: Frobenius norm of E8, an upper bound of the Jacobian's 2-norm, so initial perturbations grow at most by G per step). (c) Chaotic ensemble: K = 10, 1024 seeded orbits, n = 10^4: sigma_FT = standard deviation of S/n on the CPU ensemble; per orbit |S_cpu/n - S_gpu/n| <= 6 sigma_FT; difference of ensemble means <= 4 combined standard errors. Tolerances (a) and (c) are conservative proposals to be preregistered for the study; the observed maxima go in the report.
13. [local] Performance record for K15: CPU (1 thread and all threads) vs GPU wall time for 10^6 orbits, d = 2, k = 2, N = 10^4, Float64 and Float32, with the CUDA runtime version and device name. Pass: the table is in the 04-verify report; GPU is the default for T only if GPU time <= CPU time (then `auto` may resolve GPU for T).
14. Single kernel. Stage 04 inspection: exactly one per-step update implementation exists; AC1 to AC9 exercise it through the CPU backend. Pass/fail by inspection, recorded in the 04-verify report.

## Sources

Verified bibliography key (all read in `data/literature/lyapunov-algorithms/Skokos2010.txt`, arXiv:0811.0882, page numbers are arXiv PDF pages):
- Skokos2010: Eq. (4) p.8; Eq. (5) p.8; Eqs. (11)-(13) p.9; Eq. (22) p.15; Eq. (53), overflow sentence and finite-time mLCE Eqs. (54)-(56) p.31; Eqs. (57)-(59) p.32; Table 1 p.34; Sect. 4.2 (random initial vector) p.20-21; Eqs. (66), (67) and small-angle remark p.39; Eq. (79) p.43; Table 2 p.44; Eq. (62) via the benchmarks.md row.
- Shevchenko2004 and SalaArtuso2015 only through the benchmarks.md rows quoted above ("Standard map largest Lyapunov exponent, large K"; "Lyapunov exponent, single chaotic orbit"), not read again here. Both are used as the lambda scale for an upper bound, never as a pass condition (K18).
Layer 3 facts: `_config/research-question.md` "Ground truth"; `_config/architecture.md` "Compute", "Layering", "Scalability"; `_config/numerics-policy.md`; `_config/coding-conventions.md` "Tests"; `schema/run-store.md` root attributes; `_config/risks.md` K15, K16, K18.
Non-bibliography facts read 2026-10-13 (for K15 only): CUDA.jl README support matrix (github.com/JuliaGPU/CUDA.jl); CUDA.jl installation overview (cuda.juliagpu.org, documents pinning a runtime version); NVIDIA CUDA Toolkit release notes 13.0, cuFFT and cuSPARSE deprecation entries. Measured here: plain vs compensated sum of 10^8 copies of 0.1 in Float64 (Julia 1.13.1).

## Amendment D28 (2026-10-13)

Decision D28 (wrapped state with winding counts) needs a state vector longer than the phase-space dimension. This amendment generalizes the stream; behaviour for a state of length d stays bit-identical.

### Equations and contract change

- A1 (derived). State length s >= d. The system gives a one-step map on the state (s components) and a one-step tangent action on the phase-space deviation (a d x d matrix, d the system's dimension). The extra s - d components (winding counts) carry no tangent. W is d x k and evolves as W_{n+1} = M_n W_n with M_n the tangent at the OLD state, as in E2.
- A2. The deviation matrix has d rows, with d taken from the system, not from the state length; 1 <= k <= d. The seeded initial vectors are d x k as before (the orbit-keyed stream does not depend on s).
- A3. Numerical health covers all s state components and all d k entries of W at every handoff and after the last step. A non-finite entry of the extra components reports as "state".
- A4. Memory per orbit: M(s, d, k, T) = sizeof(T) (s + d k + 2 k) + 16 bytes; s = d gives the old formula. Worked value: s = 4, d = 2, k = 2, Float64: 8 (4 + 4 + 4) + 16 = 112 bytes. The chunk size formula is unchanged with this M.
- A5. Reducers receive the full state vector (s components) as the point. A reducer written for s = d is unchanged.
- A6. Orbit independence (I5) and thread and chunk invariance (AC6) hold for any s.

### Compute (amendment)

The single stepping body stays the only kernel. Buffers of states hold s components, buffers of deviation matrices d x k. Zero allocation (I8) holds for s > d. Interface: the row count d is read from the system's dimension (every Batch A system defines it); `bytes_per_orbit(Val(d), Val(k), T, s)` takes the state length as an optional fourth argument (default d) and `chunk_size` takes the keyword `s` (default d).

### Acceptance criteria (amendment)

- W9 State longer than the tangent dimension: with a test map of state length 4 and dimension 2 (the last two components are a step counter and a running sum of x; the tangent is the standard-map matrix at x), k in {1, 2}, tau in {1, 3}: (a) the per-orbit propagation returns a state of 4 components equal bitwise to the manual iteration of the map and a 2 x k deviation matrix equal bitwise to the manual tangent product; (b) the chunked propagation over 257 seeded orbits equals the serial map of the per-orbit propagation bitwise (chunk sizes 16 and default); (c) a non-finite value in component 3 or in component 4 raises the numerical health error of kind "state" at the correct handoff; (d) the per-orbit propagation allocates 0 bytes after warm-up; (e) a reducer sees the 4-component state; (f) the memory formula gives 112 bytes for s = 4, d = 2, k = 2 in Float64 and 96 for the default call, and the chunk size follows; (g) k = 3 with d = 2 is rejected before any work; (h) the existing suites (s = d) are unchanged.

## Amendment D29 (2026-10-13)

Decision D29 builds the GPU path in M2: CUDA.jl as a weak dependency, the probe answering per element type, `gpu_faster` taken from a recorded measurement. This amendment fixes the device hook, the per-chunk fallback handling and the record. The kernel, the reducer contract and the CPU results are unchanged.

### Device hook (contract)

- `Foundation` owns `AbstractDevice`, `NoDevice` and the run-wide `DEVICE` (a `Ref`, `NoDevice()` by default). Nothing in `src/` loads CUDA. The hooks are not exported: `probe(device, T)` (the NamedTuple `loadable, functional, device, compiles, gpu_faster`, answered per element type), `free_memory(device)`, `is_oom(exception)`, `is_compile_error(exception)`; the last two default to `false`.
- `resolve_backend(request, T; cpu_only_reducer, probe)`: `probe` is an optional zero-argument function (the test hook of AC10, unchanged); without it the record of `DEVICE[]` for `T` is used. Order of checks unchanged (E10).
- `Tangent` owns `device_backend(device)` (the KernelAbstractions backend of a device; `CPU()` for `NoDevice`) and `workgroup(backend)` (orbits per work group; 16 on the CPU, 256 on CUDA). The CUDA extension `ext/ChaoslabCUDAExt.jl` sets `DEVICE[]` when `CUDA.functional()`, adds the methods and implements `probe`.

### Fallback handling in `propagate`

1. The reducer chain is built before backend resolution; a state that is not `isbits` resolves to the CPU with reason 6 ("reducer is CPU-only").
2. A GPU run allocates the seven chunk buffers on the device backend, fills the initial deviations on the host and copies them (and the slice of `x0`) to the device, and brings status, state, deviation, reducer-state and failing-step buffers back to the host after `synchronize`. Everything after that is the CPU code. The chunk size comes from `free_memory(DEVICE[])`.
3. A device chunk is wrapped in one `try`: `is_oom` reruns the chunk on the CPU and keeps the rest of the run there, reason 4; `is_compile_error` the same with reason 3; a health error (non-OK status) reruns the chunk on the CPU, and if the CPU passes the run continues there with reason 5, otherwise the CPU error is thrown (its `backend` field is `:cpu`). Any other exception is rethrown. Provenance records the final backend and its reason; `device` is `cpu` after a fallback.

### Record

`data/gpu-benchmark.toml` under the package directory (git-ignored, per machine): tables `[<device name>."<precision label>"]` with `gpu_faster`, `cpu_ns_per_step`, `gpu_ns_per_step`, `threads`, `date`, and `reason` when the element type did not run on the device. `gpu_faster` is true when the GPU time per orbit-step is at most the CPU time with all threads on the AC13 ensemble; an absent entry reads as false. `verify/gpu.jl` writes it; the extension only reads it. The record file is the only home of this fact (`_config/architecture.md`, "Compute").

### Acceptance criteria (amendment)

- D29a Device path equals the CPU path: with a test device whose backend is the CPU backend, `propagate` with `backend = :gpu` gives bitwise equal results and provenance `gpu` with the device name.
- D29b Reasons 3, 4, 6 and the health path: a device chunk that throws an out-of-memory or compile error reruns on the CPU (first, middle and last chunk) with the exact reason string; a reducer with a non-isbits state records reason 6; a health error that only the first run of a chunk raises records reason 5; a health error the CPU raises too is thrown with `backend = :cpu`; other exceptions propagate.
- D29c [local] The STAGE04-GPU lines of the four Batch A verify logs, printed by `verify/gpu.jl` (spec items AC12, AC13 above and the matching items of the precision, standard-map and diagnostics specs).

## Amendment D33 (2026-10-13)

Decision D33 (owner, after the D29 GPU run and an independent review) adds the launch budget, the lost-device path and the symmetric AC12a. It changes the Amendment D29 text above as follows.

### Launch budget (contract)

- One launch runs all N steps of every orbit in its chunk, and the Windows display watchdog (TDR) kills a launch that takes too long (8 s on the measured machine, 2 s by default). A launch is estimated as `gpu_ns_per_step x N x orbits`, with `gpu_ns_per_step` from the record (`Foundation.gpu_step_ns`), and never less than for one wave of 1024 orbits (`LAUNCH_WAVE`; the launch time is flat below it: 256 and 1024 orbits x 1e4 steps took 0.38 s and 0.30 s). The budget is 1 s (`LAUNCH_BUDGET_S`).
- The device chunk is the smaller of the memory-based or requested chunk and the largest chunk within the budget. If that chunk is below one wave, N alone exceeds the budget: the run uses the CPU and records the new eighth fallback reason "launch exceeds the display watchdog" (appended last to `FALLBACK_REASONS`; architecture "Compute" lists it). The seven reasons of the D29 text above plus this one make eight.
- No record (no `gpu_ns_per_step` for the device and element type, or an unreadable record): the conservative default 100 ns per orbit-step (`DEFAULT_GPU_NS`, about 5 times the value measured on a GTX 1050) is used. With it, one wave of 1024 orbits runs up to N = 9765 steps.
- Step-chunked relaunch (carrying the state between launches) is a named run before M2 exit; until then a long N runs on the CPU.

### Lost device and setup failures

- `Foundation.is_device_lost(exception)` (default false; the CUDA extension answers true for the CuError codes after which the context is unusable, among them the launch timeout 702). A lost device sends the rest of the run to the CPU with the reason "launch exceeds the display watchdog" (the nearest existing string: a launch timeout is exactly that, and other context losses share it until the owner names a separate string, proposed: "device lost") and resets `DEVICE[]` to `NoDevice()`, so later calls do not touch the dead context.
- The first `free_memory` call and the buffer allocation on the device sit inside the same classification as a chunk: an out-of-memory error there is reason 4, a compile error reason 3, a lost device as above.

### Provenance and record

- After a mid-run fallback the CPU chunk size is recomputed from host memory (an explicit `chunk` is kept). Provenance `chunk_size` is the chunk size of the backend that finished the run; `gpu_chunks` counts the chunks that ran on the device first and `gpu_chunk_size` is the device chunk size (0 if the device was never used).
- The record is read inside `try`; any read, parse or shape error gives `gpu_faster = false` and no ns value. `verify/gpu.jl` writes it through a temporary file and a move.

### AC12a, symmetric form

Each backend is compared with a Double64 CPU reference on all 1024 seeded orbits (K = 0.5, tau = 10): the numbers of orbits whose S/n differs from the reference by more than 1e-6 differ between CPU and GPU by at most 3, and on the orbits where both backends agree with the reference within 1e-6 the CPU-GPU gap is at most 2e-6. This replaces the CPU-selected class of the D29 run (which failed at 1.162e-6 on one orbit although 20 orbits were above 1e-6 for each backend). At n = 10^5 one launch of 1024 orbits is estimated at 2.2 s, above the budget, so the run goes to the CPU with the new reason; `verify/gpu.jl` prints that FALLBACK line and evaluates the same criterion at n = 4 x 10^4, the largest n whose one-wave launch fits the budget.

### Acceptance criteria (amendment)

- D33a Budget: with a recorded gpu_ns_per_step the device chunk is capped (exact chunk size); with N beyond one wave the provenance reason is the eighth string, the results equal the CPU run, `gpu_chunks` is 0; without a record the default applies.
- D33b Lost device: a device chunk that throws a lost-device error reruns on the CPU with the eighth reason, `DEVICE[]` is `NoDevice()` afterwards.
- D33c Setup: an out-of-memory error from `free_memory` or from the buffer allocation gives reason 4 and a complete CPU run.
- D33d Provenance: `gpu_chunks`, `gpu_chunk_size`, and `chunk_size` recomputed from host memory after the switch.
