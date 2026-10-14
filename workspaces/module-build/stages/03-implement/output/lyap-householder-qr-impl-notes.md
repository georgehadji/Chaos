# Implementation notes: lyap-householder-qr

Code, all in `src/Lyapunov/`: `householder.jl` (about 135 lines: `HouseholderQR`, `AugmentedQR`, `GivensQR`, `orthonormalize` methods, `_finish`, provenance of the family), `lapack_qr.jl` (about 85 lines: `LapackQR`, `LapackWork`, `LapackState`, the reducer methods for the library path). Bench: `bench/lyap-householder-qr.jl` (9 entries, not run, BenchmarkTools is not a dependency). No function over 50 lines.

Edits outside the new files, all additive and behaviour-neutral for the Gram-Schmidt family:
- `Lyapunov.jl`: the two includes (`lapack_qr.jl` before `householder.jl`, which needs the `LapackQR` type for the provenance union) and the four exports.
- `qr_spectrum.jl`: `reducer_step(::QRSpectrum, ::QRState, ...)` is now a one-line call of `_qr_step(r, st, n, W, fact)` with `fact = orthonormalize(r.alg, W)`; the body is unchanged. `LapackQR` calls the same `_qr_step` with its own factorisation. The `st.code == 0` early return stays in `reducer_step`, so a failed state is not re-factored.
- `orthonormalize.jl`: `window_Gw` rescales a triangle whose largest entry is above 1e100 or below 1e-100 by that entry before the squared sums (G_w is scale invariant). Unrequested robustness change in Gram-Schmidt-owned code, an owner item (revert if unwanted, no test depends on it: AC11(c) reads status, R, Q and lambda, never `G_max`). Without it `reducer_step` at `W = 1e200 Q` gave `G_max = NaN` silently. Triangles inside the range take the old path, so Gram-Schmidt results are bitwise unchanged.

## Status

Final code, Julia 1.13.1, 2 threads: `test/lyap-householder-qr/` 2035 pass, 0 fail, 0 error (with the four plan cuts below; 2196 pass before them, the 161 fewer tests are the cut parts). AC8d passes now because `propagate_shadow` landed from the parallel lyap-two-trajectory work; before that it was the only error (`UndefVarError`, not this slug).
- `lyap-gram-schmidt` on 1.13 (full-suite run, before the last two small edits): 19555 pass, 0 error; the two `gs_api.jl:66` assertions pass.
- Quality (Aqua and JET) 12/12 on 1.13 and on 1.10 with the final code; layering 2/2. An earlier run found an Aqua "unbound type parameters" failure in the first form of the Givens helpers (`NTuple{d, SVector{m, T}}` signatures); fixed by taking `rows::Tuple`. JET only runs after Aqua passes, so it was re-run on the final code.
- Other suites in the 1.13 full package run: diagnostics 245, lyap-benettin 15009, lyap-svd-stabilized 12915, precision 150, standard-map 774, tangent-core 565 all green; lyap-two-trajectory was still RED then (its own module, then incomplete).
- Julia 1.10.12, per folder (the full `Pkg.test` was not run; standard-map is slow on 1.10): lyap-gram-schmidt 19555, lyap-benettin 15009, tangent-core 565, precision 150, diagnostics 245 all pass; lyap-householder-qr (final code, with the cuts) 2035 pass, 0 error, 15 min 33 s testset time. The raw `ccall` of `LapackQR` through `libblastrampoline` and the `@generated` Givens sweep run on 1.10 in the Householder suite.

## Algorithms

- V-H: reflectors on the d x k matrix (E5, `sign(0) = +1`, overflow-safe `c` and `|v|` through `safe_norm`), skipped when `c = 0` or not finite; the swept diagonal is overwritten with `-sign(a_j) c`; `Q` by applying the stored reflectors in reverse order to the first k columns of `I_d`. A negative pivot flips the row of R and the column of Q; a pivot that is not finite and positive gives `QR_PIVOT` with the pivot in `R[column, column]` (same channel as the Gram-Schmidt variants).
- V-A: the same sweep on `[W, I_d]`; `Q` is the transposed first k rows of the augmented block. V-G: rotations on the rows of `[W, I_d]`, the sweep `j = 1..min(k, d-1)`, `i = d..j+1` emitted by a `@generated` function (straight-line code, constant row and column indices). Type-level `Val` recursion for the same sweep did not infer and allocated; a first version with runtime loops over a masked `SMatrix` compiled for 27 s at d = k = 8 (Float64) and ran 72 us; the current one compiles in 4 s and runs 7 us.
- Counts: V-H and V-A report the reflectors actually applied, V-G the rotations applied (equal r(d,k), g(d,k) of E9 on generic matrices; zero vectors are skipped, as the spec says).
- `LapackQR`: `dgeqrf_` and `dorgqr_` called through `libblastrampoline` with a preallocated work array (`LinearAlgebra.LAPACK.geqrf!` allocates `work` per call), R read from the upper triangle before `dorgqr_` overwrites it, same sign fix and status as the kernels. `reducer_init` throws `ArgumentError` for `T != Float64` and sets BLAS threads to 1 when Julia has more than one thread (global setting, done once). State `LapackState{K}` = the shared `QRState{Float64,K}` plus the workspace, non-isbits, so `_propagate` resolves "reducer is CPU-only" with no new code.
- Provenance of the family: `(estimator = "qr", orthonormalizer = "V-H"|"V-A"|"V-G", iteration_path = "diagonal", implementation = "kernel"|"library", orthogonality_tol)`. `checks_conditioning` stays false and the orthogonality threshold uses g = 1 (the generic defaults), as the contract says.

## Deviations and decisions

- `orthonormalize` for `LapackQR` is the pure form with a fresh workspace per call (allocating); the zero-allocation contract is on `reducer_step` with the reused workspace (measured 0 bytes for d = k in {2x1, 2x2, 4x4, 8x8}).
- V-A and V-G produce Q from the transposed rows of S = Q_full' instead of applying the transposed rotations to I in reverse: same exact result, different rounding order. Tolerance-based tests pass.
- The G_w rescale is the only change of shared behaviour; flagged for the owner because it sits in `lyap-gram-schmidt` code.
- Float32 and BigFloat for the kernels work through the generic code (BigFloat allocates, exempt); not specifically tested here.

## Zero allocation

`reducer_step`, Float64 and Double64, V-H, V-A, V-G, (d,k) in {(2,1), (2,2), (4,4), (8,8)}: 0 bytes; `LapackQR` Float64: 0 bytes (measured with the AC10 pattern after a warm-up call). Allocation checks are fragile in two ways found on the way: a global `L = Chaoslab.Lyapunov` (non-const) and indexing a tuple of `SVector` rows with a computed index both show spurious bytes.

## Performance (Julia 1.13.1, 1 thread, the machine heavily loaded by other agents, timer resolution about 0.5 us, indicative only)

`orthonormalize` Float64: 2x2 about 1 us (V-H), 0.5 us (V-A, V-G), LAPACK 5 us; 4x4 about 2 us (V-H), 1.5 us (V-A, V-G), LAPACK 13 us; 8x8 6.5 to 7 us (kernels), LAPACK 46 us. Double64: 4x4 9 to 19 us, 8x8 116 to 161 us. For d <= 8 the kernels are 3 to 10 times faster than the library call (the call overhead dominates), which is the K15 argument for the kernel path.
Compile time per specialisation (first call): 0.2 to 1 s at d <= 4, 1.8 to 3.8 s (Float64) and 4.7 to 6.9 s (Double64) at 8x8.

## Runtime against the 3 minute ceiling (D35): over, cuts applied

`julia -t 2`, package precompiled, test code compiled in the run, other agents' Julia processes on the machine (2 other processes, 13 percent CPU load at the start of the last run; up to 6 to 8 processes in the first two). Testset time: first run 23 min 02 s (before the V-G rewrite), second 12 min 50 s (all tests, no cuts), third 14 min 09 s (with the cuts, AC8d now running, a Julia 1.10 quality run in parallel). Largest testsets in the second run: AC6 2m12s, AC1 1m55s, AC9 1m18s, AC12 1m15s, AC7 1m10s, AC10 54s, AC11 50s, AC4 49s.

The plan's four cuts, in its order, were applied per the brief (the only edits under `test/lyap-householder-qr/`; the stage 02 diff of that folder is exactly these four lines): AC9 chunk 64 dropped (`hh_runtime.jl`, kernel and library), AC7 tau 4 dropped (`hh_panel.jl:52`), AC8b K = 10 dropped (`hh_cross.jl`), AC1 `HH_NS` 200 to 100 (`hh_factor.jl`, AC3 follows). They did not bring the folder under 3 minutes and the time did not fall measurably, because the cost is compile time, not run time: a second pass of AC1 in the same process took 17 s against 4 min 10 s for the first (before the V-G rewrite). What remains is the compile of every (T, d, k, variant) specialisation plus the BigFloat reference code, in AC1, AC3, AC6, AC10, AC11, AC12. Further reduction is a stage 02 decision (fewer (T, d, k, variant) combinations, Double64 8x8 kernels cost 5 to 7 s each to compile, or splitting the folder).

## Known limits

- Kernel Double64 8x8 is 120 to 160 us per factorisation; fine for the CPU path, to be measured on the device in stage 04.
- GPU compilation of the kernels (generated function, `ntuple` with closures, `Base.setindex` free) is untested here; per-chunk fallback with reason 3 applies.
- `LapackQR` states hold an array each; `propagate` over many orbits allocates one workspace per orbit in `reducer_init` (expected for the library path).
- Coverage not measured in this run.
